#!/usr/bin/env python3
"""Steady-state critical cycle of the EMITTED tick (after the emitter's depth_pass), by phase.
Records the real tick block with phase spans, chains NT applications with hand-off marks,
emits through write_circuit (depth_pass on), expands, schedules with the board model, and
buckets the critical path by (tick, phase). Phases are matched to reordered ops by gate identity."""
import sys, os, io, json
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import two_bit, point_add, wide, arith, emit as qemit
from qec.arith import Mod
from qec.builder import Circuit

SECP = (1 << 256) - (1 << 32) - 977
MODS = [two_bit, point_add, wide, arith]
NT = 8
for a in sys.argv[1:]:
    k, v = a.split("=", 1)
    if k == "NT": NT = int(v); continue
    m = next(m for m in MODS if hasattr(m, k)); setattr(m, k, eval(v))

spans = []
def wrap_fn(mod, name, label):
    fn = getattr(mod, name)
    def g(C, *a, **k):
        s = len(C.prims); r = fn(C, *a, **k); spans.append((s, len(C.prims), label)); return r
    setattr(mod, name, g)
def wrap_init(cls, label):
    fn = cls.__dict__["__init__"]
    def g(self, C, *a, **k):
        s = len(C.prims); r = fn(self, C, *a, **k); spans.append((s, len(C.prims), label)); return r
    cls.__init__ = g
def wrap_m(cls, name, label):
    fn = cls.__dict__[name]
    def g(self, *a, **k):
        C = self.C; s = len(C.prims); r = fn(self, *a, **k); spans.append((s, len(C.prims), label)); return r
    setattr(cls, name, g)
wrap_fn(two_bit, "add_or_shift_la", "add")
for nm in ("gt_tree", "gt_sweep"):
    if hasattr(wide, nm): wrap_fn(wide, nm, "compare")
wrap_init(wide.PrefixTree, "prefix_build"); wrap_m(wide.PrefixTree, "write", "prefix_write")
wrap_init(wide.LowAdd, "lowadd_build"); wrap_m(wide.LowAdd, "undo", "lowadd_undo"); wrap_m(wide.LowAdd, "write", "lowadd_write")

C = Circuit(hierarchical=True)
mod = Mod(SECP); M = point_add.DialogMultiplier(C, mod); n = mod.n
assert len(C.prims) == 0
blk = M._bezout_tick(False, False, False, True, True)
tick_spans = list(spans)
(c,) = C.alloc(1); R = C.alloc(n); S = C.alloc(n)
cells = [C.alloc(2) for _ in range(NT + 2)]
hand = C.alloc(1)[0]
for i in range(NT, 0, -1):
    pout = C.alloc(1)
    C.apply(blk, [cells[i + 1], cells[i], cells[i - 1][:1], R, S, [hand], pout])
    C.free([hand]); hand = pout[0]
buf = io.StringIO(); qemit.write_circuit(buf, C, {"control": [c], "x": R, "y": S})
circ = json.loads(buf.getvalue())
assert len(circ["definitions"]) == 1, len(circ["definitions"])
body = circ["definitions"][0]["operations"]

def phase_of(pi):
    best = None
    for s, e, nm in tick_spans:
        if s <= pi < e and (best is None or e - s < best[0]): best = (e - s, nm)
    return best[1] if best else "glue"
local = {w: k for k, w in enumerate(sorted(blk.wires))}
pool = defaultdict(list)
ops0 = blk.ops
for i in range(0, len(ops0), 4):
    k, a, b, c_ = ops0[i:i + 4]
    if k == 0: key = ("X", local[a])
    elif k == 1: key = ("CX", local[a], local[b])
    elif k == 2: key = ("CCX", frozenset((local[a], local[b])), local[c_])
    else: key = ("SWAP", frozenset((local[a], local[b])))
    pool[key].append(phase_of(i))
labels = []
for r in body:
    o = r["op"]
    if o == "X": key = ("X", r["target"])
    elif o == "CX": key = ("CX", r["control"], r["target"])
    elif o == "CCX": key = ("CCX", frozenset((r["control1"], r["control2"])), r["target"])
    else: key = ("SWAP", frozenset((r["wire1"], r["wire2"])))
    labels.append(pool[key].pop(0) if pool[key] else "?")

ready = defaultdict(int); who = {}; rec = []; fin = []; lab = []
for ci, call in enumerate([r for r in circ["operations"] if "call" in r]):
    wm = call["wires"]
    for j, r in enumerate(body):
        o = r["op"]
        if o == "X":
            a = wm[r["target"]]; ia = ready[a]; ready[a] = ia + 1; eff = {a: ia}; f = ia + 1; ws = (a,)
        elif o in ("CX", "SWAP"):
            a, b = (wm[r["control"]], wm[r["target"]]) if o == "CX" else (wm[r["wire1"]], wm[r["wire2"]])
            ia, ib = ready[a], ready[b]; x = max(ia, ib) + (1 if o == "CX" else 3)
            ready[a] = ready[b] = x; eff = {a: ia, b: ib}; f = x; ws = (a, b)
        else:
            a, b, t = wm[r["control1"]], wm[r["control2"]], wm[r["target"]]
            ra, rb, rt = ready[a], ready[b], ready[t]
            eff = {a: ra + 7, b: rb + 9, t: rt + 10}
            rt += 1; x = max(rb, rt) + 1; rb = rt = x; rt += 1
            x = max(ra, rt) + 1; ra = rt = x; rt += 1
            x = max(rb, rt) + 1; rb = rt = x; rt += 1
            x = max(ra, rt) + 1; ra = rt = x
            rb += 1; rt += 2
            x = max(ra, rb) + 1; ra = rb = x; ra += 1; rb += 1
            x = max(ra, rb) + 1; ra = rb = x
            ready[a], ready[b], ready[t] = ra, rb, rt; f = max(ra, rb, rt); ws = (a, b, t)
        bind = max(eff, key=eff.get)
        idx = len(rec); rec.append(who.get(bind)); fin.append(f); lab.append((ci, labels[j]))
        for w in ws: who[w] = idx
total = max(ready.values())
nb = len(body)
ends = [max(fin[k * nb:(k + 1) * nb]) for k in range(NT)]
print(f"emitted tick: {nb} ops; chain total {total}; periods {[ends[i]-ends[i-1] for i in range(1,NT)]}")
i = who[max(ready, key=ready.get)]; path = []
while i is not None: path.append(i); i = rec[i]
path.reverse()
seq = []; prevf = 0
for idx in path:
    key = lab[idx]; d = fin[idx] - prevf; prevf = fin[idx]
    if seq and seq[-1][0] == key: seq[-1][1] += d; seq[-1][2] += 1
    else: seq.append([key, d, 1])
for (tk, ph), d, cnt in seq:
    if NT - 3 <= tk < NT - 1: print(f"  tick{tk} {ph:14s} {d:5d} ticks {cnt:3d} ops")

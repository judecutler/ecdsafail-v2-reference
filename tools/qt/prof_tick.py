#!/usr/bin/env python3
"""Critical-path profile of one standalone Bezout tick: ticks spent along the binding
path, bucketed by the innermost generator phase that emitted each op.
Phases are spans of prim indices recorded by wrapping generator entry points.
Usage: python3 -I .audit/prof_tick.py [KEY=VAL ...]   (two_bit module constants)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import two_bit, wide, arith
from qec.builder import Circuit

SECP = (1 << 256) - (1 << 32) - 977
for arg in sys.argv[1:]:
    if "=" in arg:
        k, v = arg.split("=", 1)
        setattr(two_bit, k, eval(v))

spans = []
CUR = [None]

def wrap_fn(mod, name, label, cidx=0):
    fn = getattr(mod, name)
    def g(*a, **k):
        C = a[cidx]
        C = getattr(C, "C", C)
        s = len(C.prims)
        r = fn(*a, **k)
        spans.append((s, len(C.prims), label))
        return r
    setattr(mod, name, g)

def wrap_method(cls, name, label):
    fn = getattr(cls, name)
    def g(self, *a, **k):
        C = self.C
        s = len(C.prims)
        r = fn(self, *a, **k)
        spans.append((s, len(C.prims), label))
        return r
    setattr(cls, name, g)

wrap_fn(two_bit, "add_or_shift_la", "add")
for nm in ("gt_tree", "gt_sweep"):
    if hasattr(wide, nm):
        wrap_fn(wide, nm, "compare")
for cls in ("PrefixTree", "PrefixLent", "PrefixSplit"):
    if hasattr(wide, cls):
        c = getattr(wide, cls)
        for m, lab in (("__init__", "prefix_build"), ("write", "prefix_write")):
            if m in c.__dict__:
                # __init__ takes C as first arg after self
                if m == "__init__":
                    fn = c.__dict__[m]
                    def mk(fn, lab):
                        def g(self, C, *a, **k):
                            s = len(C.prims); r = fn(self, C, *a, **k); spans.append((s, len(C.prims), lab)); return r
                        return g
                    setattr(c, m, mk(fn, lab))
                else:
                    wrap_method(c, m, lab)
la = wide.LowAdd
fn0 = la.__dict__["__init__"]
def la_init(self, C, *a, **k):
    s = len(C.prims); r = fn0(self, C, *a, **k); spans.append((s, len(C.prims), "lowadd_build")); return r
la.__init__ = la_init
wrap_method(la, "undo", "lowadd_undo")
wrap_method(la, "write", "lowadd_write")

mod = two_bit.TwoBitMod(SECP)
C = Circuit(hierarchical=False)
n = mod.n
prev = C.alloc(2); g = C.alloc(2); nxt = C.alloc(1); R = C.alloc(n); S = C.alloc(n)
mod.tick(C, prev, g, nxt[0], R, S, pebble=False)

prims = C.prims
N = len(prims) // 4
ready = {}; who = {}; rec = [None] * N; fin = [0] * N
for idx in range(N):
    k, a, b, c = prims[4*idx:4*idx+4]
    if k == 0:
        ia = ready.get(a, 0); ready[a] = ia + 1; eff = {a: ia}; f = ia + 1
    elif k in (1, 3):
        ia, ib = ready.get(a, 0), ready.get(b, 0); x = max(ia, ib) + (1 if k == 1 else 3)
        ready[a] = ready[b] = x; eff = {a: ia, b: ib}; f = x
    else:
        ra, rb, rt = ready.get(a, 0), ready.get(b, 0), ready.get(c, 0)
        e = {a: ra + 7, b: rb + 9, c: rt + 10}
        rt += 1; x = max(rb, rt) + 1; rb = rt = x; rt += 1
        x = max(ra, rt) + 1; ra = rt = x; rt += 1
        x = max(rb, rt) + 1; rb = rt = x; rt += 1
        x = max(ra, rt) + 1; ra = rt = x
        rb += 1; rt += 2
        x = max(ra, rb) + 1; ra = rb = x; ra += 1; rb += 1
        x = max(ra, rb) + 1; ra = rb = x
        ready[a], ready[b], ready[c] = ra, rb, rt; eff = e; f = max(ra, rb, rt)
    bind = max(eff, key=eff.get)
    rec[idx] = who.get(bind); fin[idx] = f
    for w in (a, b, c):
        who[w] = idx
total = max(ready.values())
i = who[max(ready, key=ready.get)]
path = []
while i is not None:
    path.append(i); i = rec[i]
path.reverse()

def phase(idx):
    pi = idx * 4; best = None
    for s, e, nm in spans:
        if s <= pi < e and (best is None or e - s < best[0]):
            best = (e - s, nm)
    return best[1] if best else "glue"

from collections import OrderedDict
seq = []  # contiguous phase runs along the path with tick deltas
prevf = 0
for idx in path:
    ph = phase(idx); d = fin[idx] - prevf; prevf = fin[idx]
    if seq and seq[-1][0] == ph:
        seq[-1][1] += d; seq[-1][2] += 1
    else:
        seq.append([ph, d, 1])
ccx = sum(1 for j in range(N) if prims[4*j] == 2)
print(f"tick ticks={total} ccx={ccx} wires={C.num_wires} path_ops={len(path)}")
agg = OrderedDict()
for ph, d, c in seq:
    agg.setdefault(ph, [0, 0]); agg[ph][0] += d; agg[ph][1] += c
for ph, (d, c) in sorted(agg.items(), key=lambda x: -x[1][0]):
    print(f"  {ph:14s} {d:5d} ticks  {c:4d} ops")
print("order:", " > ".join(f"{ph}:{d}" for ph, d, c in seq))

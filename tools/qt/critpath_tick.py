#!/usr/bin/env python3
"""Critical path of one standalone Bezout tick, by phase.
Backtracks the board schedule from the last-finishing wire and buckets the path
ops by which generator phase emitted them (recorded via prim-index ranges)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import two_bit
from qec.builder import Circuit

SECP = (1 << 256) - (1 << 32) - 977
for arg in sys.argv[1:]:
    if "=" in arg:
        k, v = arg.split("=", 1)
        setattr(two_bit, k, eval(v))

# Phase markers: wrap key methods to record (start_prim, end_prim, name).
spans = []
mod = two_bit.TwoBitMod(SECP)
orig = {}
def wrap(name, fn):
    def g(C, *a, **k):
        s = len(C.prims)
        r = fn(C, *a, **k)
        spans.append((s, len(C.prims), name))
        return r
    return g
for nm in ("_add_or_shift", "_chain_up", "_chain_down", "_begin", "_end", "_cdbl"):
    if hasattr(mod, nm):
        orig[nm] = getattr(mod.__class__, nm)

C = Circuit(hierarchical=False)
n = mod.n
prev = C.alloc(2); g = C.alloc(2); nxt = C.alloc(1); R = C.alloc(n); S = C.alloc(n)
# instrument add
import types
real_add = mod._add_or_shift
def add_hook(Cc, J, a, word):
    s = len(Cc.prims); real_add(Cc, J, a, word); spans.append((s, len(Cc.prims), "add_or_shift"))
mod._add_or_shift = add_hook
s0 = len(C.prims)
mod.tick(C, prev, g, nxt[0], R, S, pebble=False)
spans.append((s0, len(C.prims), "TICK_ALL"))

# schedule with op-index tracking
prims = C.prims
ready = {}; who = {}; rec = []
for i in range(0, len(prims), 4):
    k, a, b, c = prims[i:i+4]
    idx = i // 4
    if k == 0:
        ins = {a: ready.get(a, 0)}; ready[a] = ins[a] + 1; eff = {a: ins[a]}
    elif k == 1:
        ia, ib = ready.get(a, 0), ready.get(b, 0); x = max(ia, ib)+1; ready[a]=ready[b]=x; eff={a:ia,b:ib}
    elif k == 3:
        ia, ib = ready.get(a, 0), ready.get(b, 0); x = max(ia, ib)+3; ready[a]=ready[b]=x; eff={a:ia,b:ib}
    else:
        ia, ib, ic = ready.get(a,0), ready.get(b,0), ready.get(c,0)
        rt=ic+1; x=max(ib,rt)+1; ib2=rt2=x; rt2+=1; x=max(ia,rt2)+1; ia2=rt3=x; rt3+=1
        x=max(ib2,rt3)+1; ib3=rt4=x; rt4+=1; x=max(ia2,rt4)+1; ia3=rt5=x; ib4=ib3+1; rt6=rt5+2
        x=max(ia3,ib4)+1; ia4=ib5=x; ia5=ia4+1; ib6=ib5+1; x=max(ia5,ib6)+1
        ready[a]=ready[b]=x; ready[c]=rt6; eff={a:ia+7,b:ib+9,c:ic+10}
    bind = max(eff, key=eff.get)
    rec.append((k, bind, who.get(bind)))
    for w in (a,b,c): who[w]=idx
total = max(ready.values())
last = who[max(ready, key=ready.get)]
path = []
i = last
while i is not None:
    k, bind, prev_i = rec[i]
    path.append(i)
    i = prev_i
path.reverse()
# bucket path ops by span (innermost = last matching)
from collections import Counter
buck = Counter()
for idx in path:
    pi = idx*4
    name = "other"
    for s,e,nm in spans:
        if s <= pi < e and nm != "TICK_ALL":
            name = nm  # last (innermost) wins since spans appended inner-first? add order: add_hook appended after real_add returns, chain spans inside. Use smallest span.
    # choose smallest containing span
    best=None
    for s,e,nm in spans:
        if s<=pi<e:
            if best is None or (e-s)<best[0]: best=(e-s,nm)
    buck[best[1] if best else "other"] += 1
print(f"tick ticks={total} path_ops={len(path)}")
for nm,c in buck.most_common():
    print(f"  {nm}: {c} path ops")

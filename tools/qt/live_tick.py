#!/usr/bin/env python3
"""Allocated-wire count over the standalone tick: after every gate, sample
num_wires - len(free_stack). Reports the max per phase (phase spans from prof_tick's wrappers)."""
import sys, os, io, contextlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import builder
samples = []
for nm in ("x", "cx", "ccx"):
    def mk(fn):
        def g(self, *a):
            fn(self, *a); samples.append(self.num_wires - len(self.free_stack))
        return g
    setattr(builder.Circuit, nm, mk(getattr(builder.Circuit, nm)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    import prof_tick as P
N, spans = P.N, P.spans
base = 2 * P.n + 5
live = [s - base for s in samples[-N:]] if len(samples) >= N else [s - base for s in samples]
print(f"gates {N} samples {len(samples)} num_wires {P.C.num_wires} scratch_total {P.C.num_wires - base} peak_alloc {max(live)}")
tops = sorted([(s // 4, e // 4, nm) for s, e, nm in spans], key=lambda x: (x[0], -x[1]))
shown, cur_end = [], 0
for s, e, nm in tops:
    if s >= cur_end and e > s:
        shown.append((s, e, nm)); cur_end = e
pos = 0
for s, e, nm in shown + [(N, N, "end")]:
    if s > pos:
        print(f"  glue           ops {pos:6d}-{s:6d} max {max(live[pos:s])} min {min(live[pos:s])}")
    if e > s:
        print(f"  {nm:14s} ops {s:6d}-{e:6d} max {max(live[s:e])} min {min(live[s:e])}")
    pos = e

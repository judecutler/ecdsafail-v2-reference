#!/usr/bin/env python3
"""Frozen hillclimb harness: one command emits the metric.

Usage: python3 -I .audit/measure.py [--toy] [--label NAME]

Emits the 256-bit hierarchical circuit at A = G from the working tree's
generator, schedules it with the board-exact tick model (gate_pass.schedule),
and prints: wires, ticks, qubit_time, CCX, T, and the per-top-level-call tick
attribution. With --toy, first runs the fast correctness subset (exhaustive
gate-level check on 3 toy curves plus the model check) and refuses to print a
metric if any lane fails. Deterministic: the generator and the scheduler are
both exact, so one run is the measurement (no sampling, no noise).
Work counter: prints total ops and lanes checked, and error count.
"""
import json, subprocess, sys, time, io, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import emit as emit_mod
from qec.gate_pass import schedule

label = "run"
if "--label" in sys.argv:
    label = sys.argv[sys.argv.index("--label") + 1]

t0 = time.time()
if "--toy" in sys.argv:
    r = subprocess.run([sys.executable, "-I", "scripts/check.py", "--toy-bits", "4", "7", "9",
                        "--toy-per-size", "2"], capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    out = r.stdout + r.stderr
    sys.stdout.write(out[-2000:])
    if r.returncode != 0:
        print(f"METRIC REFUSED: toy correctness subset failed (exit {r.returncode})")
        sys.exit(1)

# Emit hierarchical at G via the repo's own emitter, in-process.
buf = io.StringIO()
import contextlib
argv_save = sys.argv
sys.argv = ["emit.py", "--hier", "--out", "/tmp/qtclimb.json"]
spec = {}
emit_script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "emit.py")
with open(emit_script) as f:
    code = f.read()
with contextlib.redirect_stdout(buf):
    exec(compile(code, emit_script, "exec"), {"__name__": "__main__", "__file__": emit_script})
sys.argv = argv_save

c = json.load(open("/tmp/qtclimb.json"))
c = c.get("circuit", c)
W = c["num_wires"]
# CCX count expanded
bodies = [d["operations"] for d in c.get("definitions", [])]
from functools import lru_cache
@lru_cache(maxsize=None)
def counts(b):
    ops = bodies[b] if b is not None else None
    raise RuntimeError
def count_ops(ops):
    ccx = tot = 0
    for r in ops:
        if "call" in r:
            s_ccx, s_tot = DEF[r["call"]]
            ccx += s_ccx; tot += s_tot
        else:
            tot += 1
            if r["op"] in ("CCX", "CCZ"):
                ccx += 1
    return ccx, tot
DEF = []
for d in c.get("definitions", []):
    DEF.append(count_ops(d["operations"]))
ccx, tot = count_ops(c["operations"])
ticks = schedule(c)
qt = W * ticks
T = 7 * ccx
print(f"[{label}] wires={W} ticks={ticks} qubit_time={qt} CCX={ccx} T={T} ops={tot} "
      f"T_ok={T <= 50465933} W_ok={W <= 4096} elapsed={time.time()-t0:.0f}s")

# Per-top-level-call attribution: schedule prefixes.
ready_hwm = 0
prefix = []
ops_top = c["operations"]
for i, r in enumerate(ops_top):
    if "call" in r:
        prefix.append((i, r["call"], r.get("inverse", False)))
sub = dict(c)
marks = []
# cumulative schedule after each top-level op index of interest (cheap: every call boundary)
idxs = [i for i, _, _ in prefix]
cuts = sorted(set(idxs + [len(ops_top) - 1]))
for cut in cuts:
    sub_ops = ops_top[: cut + 1]
    sub = {"num_wires": W, "definitions": c.get("definitions", []), "operations": sub_ops}
    marks.append((cut, schedule(sub)))
prev = 0
names = {i: f"call(def {d}{', inv' if inv else ''})" for i, d, inv in prefix}
for cut, tk in marks:
    print(f"  top[:{cut+1}] {names.get(cut, 'glue')}: cum_ticks={tk} (+{tk-prev})")
    prev = tk

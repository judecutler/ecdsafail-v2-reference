#!/usr/bin/env python3
"""Coordinate descent over Bezout-tick knobs against a qubit-time proxy:
    qt ~ (W0 + dw) * (T0 + 768 * (ticks - t0))
W0, T0: circuit wires and raw ticks at the base setting; 768 tick calls per circuit.
Tie-break: fewer CCX. Deterministic. Prints every accepted move."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep
from bpass import bpass as run

W0, T0 = 1556, 1695337
base = run({})
t0, c0, w0 = base

def proxy(r):
    return (W0 + r[2] - w0) * (T0 + 768 * (r[0] - t0))

def tree_vals(cur, lo_first):
    a, b, cap, fan, leaf = cur
    out = set()
    for da in (-1, -0.5, -0.2, 0.2, 0.5, 1):
        out.add((round(a + da, 2), b, cap, fan, leaf))
    for db in (-1, -0.5, -0.2, 0.2, 0.5, 1):
        out.add((a, round(b + db, 2), cap, fan, leaf))
    for c in (16, 24, 32, 48, 64, 96, 128):
        out.add((a, b, c, fan, leaf))
    for l in (1, 2, 3):
        out.add((a, b, cap, fan, l))
    return [v for v in out if v[0] > 1 and v[1] > 1]

def cands(k, cur):
    if k in ("WIDE_ADD_TREE_F", "WIDE_ADD_TREE_B"):
        return tree_vals(cur, True)
    if k == "WIDE_CMP_SWEEP":
        a, b, cap, fan, leaf = cur
        return [(x, b, c, fan, l) for x in (2, 3) for c in (64, 128, 256) for l in (1, 2, 3)]
    if k == "WIDE_ADD_COPIES":
        return [32, 40, 48, 56, 64, 72, 80, 96, 128]
    if k == "SWAP_COPIES":
        return [128, 192, 256, 320, 384, 512]
    if k == "P_FAN":
        return [16, 24, 32, 48, 64]
    if k == "IDLE_FRESH":
        return [0, 2, 4, 6, 8, 12, 16, 24]
    if k in ("WIDE_P", "WIDE_M", "WIDE_LOW"):
        return [2, 3, 4, 5, 6]
    if k == "WIDE_CMP":
        return [(3, 4, 6), (2, 3, 4, 5), (3, 4, 5), (2, 4, 6), (3, 5, 7), (4, 6), (3, 4, 6, 8)]
    if k == "IDLE_SLACK":
        return [0, 1, 2, 4, 8]
    if k == "LATE_MIN_COPIES":
        return [2, 4, 6, 8, 12]
    if k in ("K_FAST","ONE_COPY","SWAP_SPREAD","LA_TOP_TEMPS","WIDE_ADD_ROT_TEMPS","CODEC_ASIDE","HAND_OFF_MARK"):
        return [True, False]
    if k == "ENCODE_SETS":
        return [1,2,3,4,6]
    return []

KNOBS = ["WIDE_ADD_COPIES", "K_FAST", "ONE_COPY", "SWAP_SPREAD", "LA_TOP_TEMPS", "WIDE_ADD_ROT_TEMPS", "CODEC_ASIDE", "HAND_OFF_MARK", "ENCODE_SETS", "WIDE_ADD_TREE_F", "WIDE_ADD_TREE_B", "WIDE_CMP_SWEEP", "WIDE_M", "WIDE_P",
         "WIDE_LOW", "IDLE_FRESH", "SWAP_COPIES", "P_FAN", "IDLE_SLACK", "LATE_MIN_COPIES", "WIDE_CMP"]
state = {}
best = base
bestp = proxy(base)
print("BASE", base, bestp, flush=True)
improved = True
rounds = 0
while improved and rounds < 6:
    improved = False
    rounds += 1
    for k in KNOBS:
        import bpass as _b
        mm = [m for m in _b.MODS if hasattr(m, k)]
        if not mm:
            continue
        cur = state.get(k, getattr(mm[0], k))
        for v in cands(k, cur):
            if v == cur:
                continue
            trial = dict(state); trial[k] = v
            r = run(trial)
            if not isinstance(r[0], int):
                continue
            p = proxy(r)
            if p < bestp or (p == bestp and r[1] < best[1]):
                best, bestp, state = r, p, trial
                cur = v
                improved = True
                print(f"r{rounds} {k}={v}: {r} proxy={p} ({100*(p/proxy(base)-1):+.2f}%)", flush=True)
print("FINAL", state, best, bestp, f"{100*(bestp/proxy(base)-1):+.2f}%")

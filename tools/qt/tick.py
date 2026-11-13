#!/usr/bin/env python3
"""Fast standalone Bezout-tick harness. Builds ONE merged two-bit tick at n=256
with a mid-pass width and schedules it under the board model. Sub-second, so it
is the hillclimb inner loop for the Bezout step (60% of all ticks).

Knob overrides: pass KEY=VALUE pairs; they set qec.two_bit module constants
before the tick is built (eval'd). Example:
    python3 -I .audit/tick.py WIDE_ADD=8 P_FAN=10

Prints: ticks, ccx, wires-touched for a representative interior tick
(prev set, nxt set, pebble off). Correctness is NOT checked here; use
.audit/measure.py --toy for the gated full metric before keeping a change.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import two_bit
from qec.builder import Circuit
from qec.gate_pass import schedule

SECP = (1 << 256) - (1 << 32) - 977

for arg in sys.argv[1:]:
    if "=" in arg:
        k, v = arg.split("=", 1)
        setattr(two_bit, k, eval(v))

def build_tick(width):
    """One interior merged tick at accumulator width = n, Euclid width `width`
    unused here (merged tick is full-width); returns a flat circuit dict."""
    C = Circuit(hierarchical=False)
    mod = two_bit.TwoBitMod(SECP)
    n = mod.n
    # Wire layout: prev cell (2), g cell (2), nxt g0 (1), R (n), S (n).
    prev = C.alloc(2)
    g = C.alloc(2)
    nxt = C.alloc(1)
    R = C.alloc(n)
    S = C.alloc(n)
    base = C.num_wires
    mod.tick(C, prev, g, nxt[0], R, S, pebble=False, pin=None, pout=None)
    return C

def main():
    C = build_tick(256)
    circ = {"num_wires": C.num_wires,
            "operations": _ops(C.prims),
            "definitions": []}
    ticks = schedule(circ)
    ccx = sum(1 for i in range(0, len(C.prims), 4) if C.prims[i] == 2)
    print(f"bezout_tick: ticks={ticks} ccx={ccx} wires={C.num_wires}")

def _ops(prims):
    names = {0: "X", 1: "CX", 2: "CCX", 3: "SWAP"}
    out = []
    for i in range(0, len(prims), 4):
        k, a, b, c = prims[i:i+4]
        if k == 0:
            out.append({"op": "X", "target": a})
        elif k == 1:
            out.append({"op": "CX", "control": a, "target": b})
        elif k == 2:
            out.append({"op": "CCX", "control1": a, "control2": b, "target": c})
        else:
            out.append({"op": "SWAP", "wire1": a, "wire2": b})
    return out

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Faithful Bezout inner loop. Runs the generator's real _bezout_pass (codecs, edge ticks,
hand-off marks) on a fabricated coded tape, shortened to ITERS ticks, built hierarchically,
emitted through the real emitter WITH its depth_pass, then scheduled with the board model.
Reports (ticks per tick, CCX per tick, pass wires). Validated against per-callee
attribution of the full circuit: base 1,198/tick, C1 config 1,249/tick."""
import sys, os, io, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import two_bit, point_add, wide, arith, emit as qemit
from qec.arith import Mod
from qec.builder import Circuit
from qec.gate_pass import schedule

SECP = (1 << 256) - (1 << 32) - 977
MODS = [two_bit, point_add, wide, arith]


def build(over=None, iters=70):
    over = over or {}
    saved = []
    try:
        for k, v in over.items():
            m = next(m for m in MODS if hasattr(m, k))
            saved.append((m, k, getattr(m, k))); setattr(m, k, v)
        C = Circuit(hierarchical=True)
        mod = Mod(SECP)
        M = point_add.DialogMultiplier(C, mod)
        M.iters = iters
        n = mod.n
        (c,) = C.alloc(1)
        Yr, S = C.alloc(n), C.alloc(n)
        nfree = sum(len(M._free_slots(j, M.swept)) for j in range(M.blocks))
        G = M.with_free_slots(C.alloc(2 * iters - nfree))
        M._bezout_pass(G, Yr, S)
        buf = io.StringIO()
        qemit.write_circuit(buf, C, {"control": [c], "x": Yr, "y": S})
        return json.loads(buf.getvalue())
    finally:
        for m, k, v in saved:
            setattr(m, k, v)


def count_ccx(circ):
    memo = {}
    defs = circ["definitions"]
    def cnt(ops):
        t = 0
        for r in ops:
            if "call" in r:
                d = r["call"]
                if d not in memo: memo[d] = cnt(defs[d]["operations"])
                t += memo[d]
            elif r["op"] == "CCX": t += 1
        return t
    return cnt(circ["operations"])


def bpass(over=None, iters=70):
    try:
        circ = build(over, iters)
        return round(schedule(circ) / iters, 1), round(count_ccx(circ) / iters), circ["num_wires"]
    except Exception as e:
        return ("ERR", type(e).__name__, str(e)[:80])


if __name__ == "__main__":
    over = {}; it = 70
    for a in sys.argv[1:]:
        k, v = a.split("=", 1)
        if k == "ITERS": it = int(v)
        else: over[k] = eval(v)
    print(over, bpass(over, it))

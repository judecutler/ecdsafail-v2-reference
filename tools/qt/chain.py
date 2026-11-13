#!/usr/bin/env python3
"""Steady-state Bezout tick period: apply the recorded tick block N times in a row the way
_bezout_pass does (hand-off mark, cells from a tape), on one flat circuit, and report
(T(N) - T(N0)) / (N - N0) plus circuit wires. This is the in-context number; the
standalone tick (all wires ready at 0) is not."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import two_bit, point_add
from qec.arith import Mod
from qec.builder import Circuit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep import sched

SECP = (1 << 256) - (1 << 32) - 977


def chain(over=None, n_ticks=12, n0=4, mods=None):
    over = over or {}
    saved = {}
    try:
        for k, v in over.items():
            m = point_add if hasattr(point_add, k) and not hasattr(two_bit, k) else two_bit
            saved[(m, k)] = getattr(m, k); setattr(m, k, v)
        C = Circuit(hierarchical=False)
        mod = Mod(SECP)
        M = point_add.DialogMultiplier(C, mod)
        n = mod.n
        R, S = C.alloc(n), C.alloc(n)
        cells = [C.alloc(2) for _ in range(n_ticks + 2)]
        blk = M._bezout_tick(False, False, False, True, True)
        hand = C.alloc(1)[0]
        marks = {}
        for j, i in enumerate(range(n_ticks, 0, -1)):
            pout = C.alloc(1)
            C.apply(blk, [cells[i + 1], cells[i], cells[i - 1][:1], R, S, [hand], pout])
            C.free([hand]); hand = pout[0]
            marks[j + 1] = len(C.prims)
        tN = sched(C.prims[:marks[n_ticks]]); t0 = sched(C.prims[:marks[n0]])
        ccx = sum(1 for i in range(marks[n0], marks[n0 + 1], 4) if C.prims[i] == 2)
        return round((tN - t0) / (n_ticks - n0), 1), ccx, C.num_wires
    except Exception as e:
        return ("ERR", type(e).__name__, str(e)[:80])
    finally:
        for (m, k), v in saved.items():
            setattr(m, k, v)


if __name__ == "__main__":
    over = {}
    for a in sys.argv[1:]:
        k, v = a.split("=", 1); over[k] = eval(v)
    print(over, chain(over))

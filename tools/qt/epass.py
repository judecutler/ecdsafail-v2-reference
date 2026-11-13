#!/usr/bin/env python3
"""Faithful Euclid inner loop: the generator's real _euclid_pass (all 384 steps, real width
schedule, tape sweep and encodes), hierarchical, emitted with depth_pass, board-scheduled.
Reports (pass ticks, CCX, wires)."""
import sys, os, io, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qec import two_bit, point_add, wide, arith, emit as qemit
from qec.arith import Mod, load_const
from qec.builder import Circuit
from qec.gate_pass import schedule
from bpass import count_ccx, SECP, MODS


def epass(over=None, inverse=False):
    over = over or {}
    saved = []
    try:
        for k, v in over.items():
            m = next(m for m in MODS if hasattr(m, k))
            saved.append((m, k, getattr(m, k))); setattr(m, k, v)
        C = Circuit(hierarchical=True)
        mod = Mod(SECP)
        M = point_add.DialogMultiplier(C, mod)
        n = mod.n
        (c,) = C.alloc(1)
        Xr = C.alloc(n); U = C.alloc(n)
        pool = list(Xr) + C.alloc(M.tape_external)
        x0, u0 = list(Xr), list(U)
        blk = C.record(M._euclid_pass, [U, pool])
        C.apply(blk, [U, pool], inverse=False)
        buf = io.StringIO()
        qemit.write_circuit(buf, C, {"control": [c], "x": x0, "y": u0})
        circ = json.loads(buf.getvalue())
        return schedule(circ), count_ccx(circ), circ["num_wires"]
    except Exception as e:
        import traceback; traceback.print_exc()
        return ("ERR", type(e).__name__, str(e)[:80])
    finally:
        for m, k, v in saved:
            setattr(m, k, v)


if __name__ == "__main__":
    over = {}
    for a in sys.argv[1:]:
        k, v = a.split("=", 1); over[k] = eval(v)
    print(over, epass(over))

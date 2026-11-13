#!/usr/bin/env python3
"""One-at-a-time knob sweep on the standalone Bezout tick. In-process.
Usage: python3 -I .audit/sweep.py            (built-in grid)
       python3 -I .audit/sweep.py 'K=V' ...  (one combined setting)"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qec import two_bit
from qec.builder import Circuit

SECP = (1 << 256) - (1 << 32) - 977


def sched(prims):
    ready = {}
    g = ready.get
    for i in range(0, len(prims), 4):
        k, a, b, c = prims[i:i + 4]
        if k == 0:
            ready[a] = g(a, 0) + 1
        elif k == 1:
            ready[a] = ready[b] = max(g(a, 0), g(b, 0)) + 1
        elif k == 3:
            ready[a] = ready[b] = max(g(a, 0), g(b, 0)) + 3
        else:
            ra, rb, rt = g(a, 0), g(b, 0), g(c, 0)
            rt += 1; x = max(rb, rt) + 1; rb = rt = x; rt += 1
            x = max(ra, rt) + 1; ra = rt = x; rt += 1
            x = max(rb, rt) + 1; rb = rt = x; rt += 1
            x = max(ra, rt) + 1; ra = rt = x
            rb += 1; rt += 2
            x = max(ra, rb) + 1; ra = rb = x; ra += 1; rb += 1
            x = max(ra, rb) + 1; ra = rb = x
            ready[a], ready[b], ready[c] = ra, rb, rt
    return max(ready.values())


def run(over):
    saved = {k: getattr(two_bit, k) for k in over}
    try:
        for k, v in over.items():
            setattr(two_bit, k, v)
        mod = two_bit.TwoBitMod(SECP)
        C = Circuit(hierarchical=False)
        n = mod.n
        prev = C.alloc(2); g = C.alloc(2); nxt = C.alloc(1); R = C.alloc(n); S = C.alloc(n)
        mod.tick(C, prev, g, nxt[0], R, S, pebble=False)
        ccx = sum(1 for i in range(0, len(C.prims), 4) if C.prims[i] == 2)
        return sched(C.prims), ccx, C.num_wires
    except Exception as e:
        return ("ERR", type(e).__name__, str(e)[:60])
    finally:
        for k, v in saved.items():
            setattr(two_bit, k, v)


def main():
    if len(sys.argv) > 1:
        over = {}
        for a in sys.argv[1:]:
            k, v = a.split("=", 1)
            over[k] = eval(v)
        print(over, run(over))
        return

    base = run({})
    print("BASE", base, flush=True)
    GRID = {
        "WIDE_ADD_COPIES": [32, 48, 96, 128],
        "SWAP_COPIES": [64, 128, 192, 384],
        "P_FAN": [16, 64],
        "IDLE_FRESH": [0, 2, 8, 16],
        "WIDE_ADD_TREE_F": [(1.8, 4, 64, 2, 2), (2.5, 4, 64, 2, 2), (2.0, 3, 64, 2, 2), (2.0, 6, 64, 2, 2),
                            (2.0, 4, 32, 2, 2), (2.0, 4, 128, 2, 2), (2.0, 4, 64, 3, 2), (2.0, 4, 64, 2, 3), (2.0, 4, 64, 4, 2)],
        "WIDE_ADD_TREE_B": [(2, 1.8, 64, 2, 2), (4, 1.8, 64, 2, 2), (3, 1.5, 64, 2, 2), (3, 2.2, 64, 2, 2),
                            (3, 1.8, 32, 2, 2), (3, 1.8, 128, 2, 2), (3, 1.8, 64, 3, 2), (3, 1.8, 64, 2, 3)],
        "WIDE_CMP_SWEEP": [(2, None, 64, 2, 2), (2, None, 256, 2, 2), (3, None, 128, 2, 2), (2, None, 128, 3, 2),
                           (2, None, 128, 2, 3), (2, None, 128, 4, 2), None],
        "WIDE_P": [3, 5, 6],
        "WIDE_M": [2, 4],
        "WIDE_LOW": [2, 4],
        "K_FAST": [False],
        "ONE_COPY": [False],
        "SWAP_SPREAD": [True],
        "EARLY_END": [False],
        "WIDE_TAIL": [False],
        "LA_TOP_TEMPS": [True],
        "WIDE_ADD_ROT_TEMPS": [False],
    }
    for k, vals in GRID.items():
        if not hasattr(two_bit, k):
            print("missing", k); continue
        for v in vals:
            r = run({k: v})
            tag = ""
            if isinstance(r[0], int):
                tag = f"dt={r[0]-base[0]:+d} dccx={r[1]-base[1]:+d} dw={r[2]-base[2]:+d}"
            print(f"{k}={v}: {r} {tag}", flush=True)

if __name__ == '__main__':
    main()

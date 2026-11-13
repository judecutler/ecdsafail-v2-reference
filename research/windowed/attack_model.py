"""Whole-attack cost of Shor ECDLP on secp256k1 by window size w, from measured parts.

Kernel numbers are measured circuits on the board (CCX, wires, ticks under the board's tick
model). Lookup numbers are measured here from qrom.lookup at w <= 12 and scaled per entry.
Per windowed addition: lookup, coherent-addend addition, unlookup. The addend word is
WORD bits: A and 2A (for the doubling exception) plus a valid bit.
adds(w) = ceil(2n / w); the literature's boundary tricks give 2n/w - 4 (28 at w = 16).
"""
import math, random, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qec.builder import Circuit
from qrom import lookup

N_BITS, WORD = 256, 1025


def sched(prims):
    ready = {}
    g = ready.get
    for i in range(0, len(prims), 4):
        k, a, b, c = prims[i:i + 4]
        if k == 0:
            ready[a] = g(a, 0) + 1
        elif k == 1:
            ready[a] = ready[b] = max(g(a, 0), g(b, 0)) + 1
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


def measure_lookup(w, fan):
    rng = random.Random(7)
    C = Circuit()
    addr, out = C.alloc(w), C.alloc(WORD)
    lookup(C, addr, [rng.getrandbits(WORD) for _ in range(1 << w)], out, None, fan)
    n = 1 << w
    return C.ops[0::4].count(2) / n, len(C.ops) / 4 / n, sched(C.ops) / n, C.num_wires - w - WORD

print("lookup per entry at WORD=1025 (measured at w=10): fan, CCX, ops, ticks, scratch")
PER = {}
for fan in (1, 16, 64, 256):
    PER[fan] = measure_lookup(10, fan)
    print(f"  fan {fan:4d}: {PER[fan][0]:.2f} CCX  {PER[fan][1]:.0f} ops  {PER[fan][2]:.1f} ticks  {PER[fan][3]} scratch")

KERNELS = {  # name: (CCX, wires, ticks) of one controlled fixed-addend addition, measured
    "lowt-climb (T leader)": (2_961_008, 4096, 30_279_444),
    "slack-order (depth leader)": (7_194_844, 1556, 1_529_222),
}
FAN = 64
ccx_e, ops_e, ticks_e, scr = PER[FAN]
for name, (pa_ccx, pa_w, pa_t) in KERNELS.items():
    print(f"\n{name}: binary method = 512 additions: T {512*7*pa_ccx/1e9:.2f}e9, ticks {512*pa_t/1e9:.2f}e9, "
          f"qubit-time {512*pa_t*pa_w/1e12:.1f}e12")
    print("   w  adds      T(e9)  lookup%T   ticks(e9) lookup%t  wires  qubit-time(e12)  ops/add(M)")
    best = {}
    for w in range(2, 21):
        adds = math.ceil(2 * N_BITS / w)
        n = 1 << w
        lk_ccx, lk_t, lk_ops = 2 * ccx_e * n, 2 * ticks_e * n, 2 * ops_e * n
        T = 7 * adds * (pa_ccx + lk_ccx)
        t = adds * (pa_t + lk_t)
        wires = pa_w + WORD + w + max(0, w - 1 + FAN - 1 - 0)  # addend word + address + lookup scratch
        qt = t * wires
        row = (w, adds, T / 1e9, 100 * lk_ccx / (pa_ccx + lk_ccx), t / 1e9, 100 * lk_t / (pa_t + lk_t), wires, qt / 1e12, lk_ops / 1e6)
        for key, val in (("T", T), ("ticks", t), ("qt", qt)):
            if key not in best or val < best[key][0]:
                best[key] = (val, row)
        if w in (4, 6, 8, 10, 12, 14, 16, 18, 20):
            print("  %2d  %4d  %9.3f  %7.1f  %10.3f  %7.1f  %5d  %12.2f  %9.1f" % row)
    for key in ("T", "ticks", "qt"):
        r = best[key][1]
        base = {"T": 512 * 7 * pa_ccx / 1e9, "ticks": 512 * pa_t / 1e9, "qt": 512 * pa_t * pa_w / 1e12}[key]
        val = {"T": r[2], "ticks": r[4], "qt": r[7]}[key]
        print(f"  best {key}: w={r[0]} ({r[1]} adds) -> {val:.3f} vs binary {base:.3f}: {base/val:.1f}x")

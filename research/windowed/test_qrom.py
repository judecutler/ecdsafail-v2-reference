"""Exhaustive check of qrom.lookup: every address, random tables, random initial out, with and
without ctrl, with and without fan-out copies; out ^= table[addr], everything else restored,
scratch clean; and lookup twice is the identity. Prints gate counts."""
import random, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qec.builder import Circuit
from qec import sim
from qrom import lookup

rng = random.Random(2026)
fails = cases = 0
for w in range(0, 7):
    for width in (1, 5, 13):
        for ctrl_on in (False, True):
            for fan in (1, 4):
                for n_entries in {1 << w, max(1, (1 << w) - 3)}:
                    table = [rng.getrandbits(width) for _ in range(n_entries)]
                    C = Circuit()
                    addr, out = C.alloc(w), C.alloc(width)
                    ctrl = C.alloc(1)[0] if ctrl_on else None
                    lookup(C, addr, table, out, ctrl, fan)
                    once = len(C.ops)
                    lookup(C, addr, table, out, ctrl, fan)
                    lanes = (1 << w) * (2 if ctrl_on else 1)
                    A = [j % (1 << w) for j in range(lanes)]
                    Cv = [j >> w for j in range(lanes)]
                    O = [rng.getrandbits(width) for _ in range(lanes)]
                    init = {**sim.pack(A, addr), **sim.pack(O, out)}
                    if ctrl_on:
                        init.update(sim.pack(Cv, [ctrl]))
                    full = table + [0] * ((1 << w) - len(table))
                    for stop, twice in ((once, False), (len(C.ops), True)):
                        s = sim.run(C.ops[:stop], C.num_wires, init, lanes)
                        got = sim.unpack(s, out, lanes)
                        want = [o if twice else o ^ (full[a] if (not ctrl_on or c) else 0) for a, c, o in zip(A, Cv, O)]
                        ok = got == want and sim.unpack(s, addr, lanes) == A
                        if ctrl_on:
                            ok = ok and sim.unpack(s, [ctrl], lanes) == Cv
                        live = set(addr) | set(out) | ({ctrl} if ctrl_on else set())
                        ok = ok and all(s[x] == 0 for x in range(C.num_wires) if x not in live)
                        cases += 1
                        fails += not ok
print(f"qrom exhaustive: {cases} cases, {fails} failures")
for w in (4, 8, 12, 16):
    C = Circuit()
    addr, out = C.alloc(w), C.alloc(8)
    table = [rng.getrandbits(8) | 1 for _ in range(1 << w)]
    lookup(C, addr, table, out)
    ccx = C.ops[0::4].count(2)
    print(f"w={w}: {ccx} CCX for {1 << w} entries ({ccx / (1 << w):.3f} per entry), scratch {C.num_wires - w - 8}")
sys.exit(1 if fails else 0)

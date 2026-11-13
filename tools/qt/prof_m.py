#!/usr/bin/env python3
"""Per-callee tick attribution inside one definition of an emitted hierarchical circuit.
Schedules the definition's body alone (all wires ready at 0), forward or inverted, with the
board model, and after each depth-1 record logs the critical time. Increments are summed per
callee definition. usage: prof_m.py circuit.json [def_id] [inv]"""
import json, sys
from collections import defaultdict

c = json.load(open(sys.argv[1])); c = c.get("circuit", c)
defs = c["definitions"]
bodies = [d["operations"] for d in defs]
# default: the definition with the most expanded ticks among top-level calls = most-called biggest; pick by arg or largest arity called at top
top_calls = [r["call"] for r in c["operations"] if "call" in r]
did = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != "auto" else max(top_calls, key=lambda d: defs[d]["arity"])
inv = len(sys.argv) > 3 and sys.argv[3] == "inv"
ready = [0] * defs[did]["arity"]

def run(ops, wmap, inverse):
    for idx in (range(len(ops) - 1, -1, -1) if inverse else range(len(ops))):
        r = ops[idx]
        if "call" in r:
            sub = [wmap[w] for w in r["wires"]] if wmap is not None else r["wires"]
            run(bodies[r["call"]], sub, inverse ^ bool(r.get("inverse", False)))
            continue
        o = r["op"]
        m = (lambda w: wmap[w]) if wmap is not None else (lambda w: w)
        if o == "X":
            ready[m(r["target"])] += 1
        elif o == "CX":
            a, t = m(r["control"]), m(r["target"]); ready[a] = ready[t] = max(ready[a], ready[t]) + 1
        elif o == "SWAP":
            a, b = m(r["wire1"]), m(r["wire2"]); ready[a] = ready[b] = max(ready[a], ready[b]) + 3
        else:
            a, b, t = m(r["control1"]), m(r["control2"]), m(r["target"])
            ra, rb, rt = ready[a], ready[b], ready[t]
            rt += 1; x = max(rb, rt) + 1; rb = rt = x; rt += 1
            x = max(ra, rt) + 1; ra = rt = x; rt += 1
            x = max(rb, rt) + 1; rb = rt = x; rt += 1
            x = max(ra, rt) + 1; ra = rt = x
            rb += 1; rt += 2
            x = max(ra, rb) + 1; ra = rb = x; ra += 1; rb += 1
            x = max(ra, rb) + 1; ra = rb = x
            ready[a], ready[b], ready[t] = ra, rb, rt

ops = bodies[did]
agg = defaultdict(lambda: [0, 0]); prev = 0; seq = []
order = range(len(ops) - 1, -1, -1) if inv else range(len(ops))
for idx in order:
    r = ops[idx]
    if "call" in r:
        run(bodies[r["call"]], r["wires"], inv ^ bool(r.get("inverse", False)))
        key = (r["call"], bool(r.get("inverse", False)) ^ inv)
    else:
        run([r], None, False); key = ("prim", False)
    now = max(ready)
    agg[key][0] += now - prev; agg[key][1] += 1; prev = now
print(f"def {did} arity {defs[did]['arity']} inv={inv} total={prev}")
rows = sorted(agg.items(), key=lambda kv: -kv[1][0])
for (k, i), (t, n) in rows[:14]:
    ar = defs[k]["arity"] if k != "prim" else "-"
    print(f"  callee {k}{' inv' if i else ''}: arity {ar} calls {n} ticks {t} per {t/n:.0f}")

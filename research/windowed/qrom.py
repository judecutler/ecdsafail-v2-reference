"""Exact, measurement-free table lookup (QROM) by unary iteration, gates X/CX/CCX only.

    lookup(C, address, table, out, ctrl=None, fan=1)

out ^= table[address] (under ctrl when given). address is little-endian. The walk visits the
address tree depth first with one AND wire a level: a node's wire is made for its left child
with one CCX, turned into the right child's with one CX from its parent, and cleared with one
CCX, so an N-entry table costs 2(N - 1) CCX uncontrolled (2N - 2 + 2 controlled, less for a
ragged table), and w - 1 clean wires (w with ctrl). The map is an XOR, so the same call undoes it:
no measurement is needed to unlook.

fan > 1 copies each leaf's control onto fan - 1 clean wires before writing the data bits, so the
CX fan-out of a wide word does not queue on one wire under a per-wire depth model.
"""


def lookup(C, address, table, out, ctrl=None, fan=1):
    w = len(address)
    assert len(table) <= 1 << w
    table = list(table) + [0] * ((1 << w) - len(table))
    copies = C.alloc(fan - 1) if fan > 1 else []

    def write(leaf, value):
        bits = [o for i, o in enumerate(out) if value >> i & 1]
        if not bits:
            return
        if not copies:
            for o in bits:
                C.cx(leaf, o)
            return
        src = [leaf] + copies
        tree = [(src[j - (1 << (j.bit_length() - 1))], src[j]) for j in range(1, len(src))]
        for a, b in tree:
            C.cx(a, b)
        for k, o in enumerate(bits):
            C.cx(src[k * len(src) // len(bits)], o)
        for a, b in reversed(tree):
            C.cx(a, b)

    def empty(base, size):
        return not any(table[base:base + size])

    def walk(node, level, base):
        """node is 1 exactly on the lanes whose address bits above `level` name this subtree."""
        if level < 0:
            write(node, table[base])
            return
        half = 1 << level
        b = address[level]
        left, right = not empty(base, half), not empty(base + half, half)
        if not (left or right):
            return
        (t,) = C.alloc(1)
        if left:
            C.x(b)
            C.ccx(node, b, t)
            C.x(b)
            walk(t, level - 1, base)
            if right:
                C.cx(node, t)
                walk(t, level - 1, base + half)
                C.ccx(node, b, t)
            else:
                C.x(b)
                C.ccx(node, b, t)
                C.x(b)
        else:
            C.ccx(node, b, t)
            walk(t, level - 1, base + half)
            C.ccx(node, b, t)
        C.free([t])

    if w == 0:
        if ctrl is None:
            for i, o in enumerate(out):
                if table[0] >> i & 1:
                    C.x(o)
        else:
            write(ctrl, table[0])
    elif ctrl is not None:
        walk(ctrl, w - 1, 0)
    else:
        # The top bit is its own node: no AND for the first level.
        top = address[w - 1]
        half = 1 << (w - 1)
        if not empty(0, half):
            C.x(top)
            walk(top, w - 2, 0)
            C.x(top)
        if not empty(half, half):
            walk(top, w - 2, half)
    C.free(copies)

"""Exact reversible integer and modular arithmetic. Registers are wire lists,
least significant bit first. Every routine returns its ancillas clean on every
input inside its stated domain."""


def _maj(C, c, b, a):
    C.cx(a, b)
    C.cx(a, c)
    C.ccx(c, b, a)


def _uma(C, c, b, a):
    C.ccx(c, b, a)
    C.cx(a, c)
    C.cx(c, b)


def _cuma(C, ctrl, c, b, a):
    C.ccx(c, b, a)
    C.ccx(ctrl, c, b)
    C.cx(a, c)
    # Left ungated because it restores b when ctrl is 0, so _cuma must not become a gated _uma.
    C.cx(a, b)


def _cx(C, ctrl, src, dst):
    if ctrl is None:
        C.cx(src, dst)
    else:
        C.ccx(ctrl, src, dst)


def add(C, a, b, cout=None, ctrl=None):
    """b += ctrl*a mod 2^w (CDKM). cout ^= ctrl*carry out."""
    w = len(a)
    assert len(b) == w
    assert ctrl is None or ctrl not in (*a, *b, cout)
    if w == 0:
        return
    (cin,) = C.alloc(1)
    carries = [cin] + a[:-1]
    for i in range(w):
        _maj(C, carries[i], b[i], a[i])
    if cout is not None:
        _cx(C, ctrl, a[-1], cout)
    for i in reversed(range(w)):
        if ctrl is None:
            _uma(C, carries[i], b[i], a[i])
        else:
            _cuma(C, ctrl, carries[i], b[i], a[i])
    C.free([cin])


def sub(C, a, b, borrow=None, ctrl=None):
    """b -= ctrl*a mod 2^w. borrow ^= ctrl*[a > b]."""
    for w in b:
        C.x(w)
    add(C, a, b, borrow, ctrl)
    for w in b:
        C.x(w)


def gt_into(C, a, b, t, ctrl=None):
    """t ^= ctrl*[a > b], a and b restored. Carry out of (~b + a) is [a > b]."""
    assert ctrl is None or ctrl not in (*a, *b, t)
    w = len(a)
    if w == 0:
        return
    (cin,) = C.alloc(1)
    for x in b:
        C.x(x)
    carries = [cin] + a[:-1]
    for i in range(w):
        _maj(C, carries[i], b[i], a[i])
    _cx(C, ctrl, a[-1], t)
    for i in reversed(range(w)):
        C.ccx(carries[i], b[i], a[i])
        C.cx(a[i], carries[i])
        C.cx(a[i], b[i])
    for x in b:
        C.x(x)
    C.free([cin])


def load_const(C, reg, value, ctrl=None):
    for i, w in enumerate(reg):
        if value >> i & 1:
            if ctrl is None:
                C.x(w)
            else:
                C.cx(ctrl, w)


def _carry_into(C, b, bit, c, t, ctrl=None):
    """t ^= MAJ(b, ctrl*bit, c), where c is a wire or None for 0. With ctrl, c
    must imply ctrl, which holds for every carry of reg + ctrl*value."""
    if not bit:
        if c is not None:
            C.ccx(b, c, t)
    elif c is None:
        _cx(C, ctrl, b, t)
    elif ctrl is None:
        C.ccx(b, c, t)
        C.cx(b, t)
        C.cx(c, t)
    else:
        # MAJ(b, ctrl, c) = c ^ b*(ctrl ^ c) when c implies ctrl.
        C.cx(c, ctrl)
        C.ccx(b, ctrl, t)
        C.cx(c, ctrl)
        C.cx(c, t)


def _carries(C, reg, value, ctrl=None):
    """cs[i] is the carry into reg[i] of reg + ctrl*value. Each entry is None
    when the carry is always 0, reg[i - 1] itself when it equals that bit, or a
    fresh wire. Release the fresh wires with _uncarry."""
    cs = [None]
    for i in range(len(reg) - 1):
        bit, c = value >> i & 1, cs[i]
        if c is None and not bit:
            cs.append(None)
        elif c is None and ctrl is None:
            cs.append(reg[i])
        else:
            (t,) = C.alloc(1)
            _carry_into(C, reg[i], bit, c, t, ctrl)
            cs.append(t)
    return cs


def _uncarry(C, reg, value, cs, i, ctrl=None):
    """Clear and free cs[i]. reg[i - 1] and cs[i - 1] must still hold what they held in _carries."""
    if i and cs[i] is not None and cs[i] != reg[i - 1]:
        _carry_into(C, reg[i - 1], value >> (i - 1) & 1, cs[i - 1], cs[i], ctrl)
        C.free([cs[i]])


def add_const(C, reg, value, ctrl=None, cout=None):
    """reg += ctrl*value mod 2^w, cout ^= ctrl*carry out. The carry chain reads
    the constant's bits, so 0 bits below the first 1 cost nothing."""
    w = len(reg)
    value %= 1 << w
    assert ctrl is None or ctrl not in (*reg, cout)
    cs = _carries(C, reg, value, ctrl)
    if cout is not None:
        _carry_into(C, reg[-1], value >> (w - 1) & 1, cs[-1], cout, ctrl)
    for i in reversed(range(w)):
        if value >> i & 1:
            if ctrl is None:
                C.x(reg[i])
            else:
                C.cx(ctrl, reg[i])
        if cs[i] is not None:
            C.cx(cs[i], reg[i])
        _uncarry(C, reg, value, cs, i, ctrl)


def gt_const_into(C, value, b, t, ctrl=None):
    """t ^= ctrl*[value > b], b restored. Carry out of (~b + value) is [value > b]."""
    assert ctrl is None or ctrl not in (*b, t)
    assert t not in b
    assert value < 1 << len(b)
    for x in b:
        C.x(x)
    cs = _carries(C, b, value, ctrl)
    _carry_into(C, b[-1], value >> (len(b) - 1) & 1, cs[-1], t, ctrl)
    for i in reversed(range(len(b))):
        _uncarry(C, b, value, cs, i, ctrl)
    for x in b:
        C.x(x)


def _reduce(C, r, value):
    """r holds s on m wires, value has bit m - 1 set, and s's top bit and the
    carry into it are never both 1. Leaves (s + flag*value) mod 2^(m-1) in
    r[:-1] and flag = [s + value >= 2^m] in r[-1]. One carry chain serves both
    the comparison and the conditional add."""
    assert value >> (len(r) - 1) & 1
    cs = _carries(C, r, value)
    flag = r[-1]
    if cs[-1] is not None:
        C.cx(cs[-1], flag)
    for i in reversed(range(len(r) - 1)):
        _uncarry(C, r, value, cs, i + 1)
        if cs[i] is not None:
            C.ccx(flag, cs[i], r[i])
        if value >> i & 1:
            C.cx(flag, r[i])


def cswap(C, ctrl, a, b):
    for ai, bi in zip(a, b):
        C.cx(bi, ai)
        C.ccx(ctrl, ai, bi)
        C.cx(bi, ai)


class Mod:
    """Arithmetic modulo an odd q < 2^n on canonical values in [0, q)."""

    def __init__(self, q):
        self.q = q
        self.n = q.bit_length()

    def _wrap(self):
        """2^(n+1) - q, so s + _wrap() carries out of n + 1 bits exactly when s >= q."""
        return (1 << (self.n + 1)) - self.q

    def modadd(self, C, a, b, ctrl=None):
        """b = (b + ctrl*a) mod q, for b < q, and a < q when ctrl is 1."""
        (top,) = C.alloc(1)
        add(C, a, b, top, ctrl)
        # s = b + ctrl*a < 2q, so s's top bit and the carry into it are never both 1.
        _reduce(C, b + [top], self._wrap())
        # top = [q was subtracted] = ctrl*[a > b]; ctrl = 0 means s = b < q and top = 0.
        gt_into(C, a, b, top, ctrl)
        C.free([top])

    def moddbl(self, C, s):
        """s = 2s mod q. 2s is a relabel with a zero low bit, so the reduction
        runs on s against _wrap() >> 1, and the low bit of the result is the flag."""
        _reduce(C, s, self._wrap() >> 1)
        s[:] = s[-1:] + s[:-1]

    def dbladd(self, C, R, S, ctrl):
        """S = (2S + ctrl*R) mod q, for S < q and R < q.

        T = 2S + ctrl*R < 3q is built as an (n+2)-bit integer and reduced by
        kappa*q, kappa = floor(T/q) in {0, 1, 2}, in one pass. With
        c = 2^n - q and K = bitlen(2c), T >= q and T >= 2q are the carries out
        of T + c and T + 2c. Below bit K those carries are two short chains;
        above it both need the same run of ones in T, held once in Q. The
        flags k1 = [kappa >= 1], k2 = [kappa = 2] are cleared by parity and
        one controlled comparison against R."""
        n, q = self.n, self.q
        c = (1 << n) - q
        K = (2 * c).bit_length()
        m = n - K
        z, top = C.alloc(2)
        T = [z] + S + [top]
        L, h0, h1 = T[:n], T[n], T[n + 1]
        (zero,) = C.alloc(1)
        add(C, R + [zero], T[:n + 1], cout=top, ctrl=ctrl)
        C.free([zero])

        A = _carries(C, T[:K + 1], c)
        B = _carries(C, T[:K + 1], 2 * c)
        Q = [None]
        for j in range(m):
            if Q[j] is None:
                Q.append(L[K + j])
            else:
                (w,) = C.alloc(1)
                C.ccx(Q[j], L[K + j], w)
                Q.append(w)

        def and_into(a, b, t):
            if b is None:
                C.cx(a, t)
            else:
                C.ccx(a, b, t)

        k1, k2, g, nz = C.alloc(4)
        # k1 = [T >= q] = h0 | h1 | (carry of L + c); k2 = [T >= 2q] = h1 | h0 & (carry of L + 2c).
        C.cx(h0, k1)
        C.cx(h1, k1)
        C.x(h0)
        C.x(h1)
        C.ccx(h0, h1, nz)
        and_into(A[K], Q[m], g)
        C.ccx(nz, g, k1)
        and_into(A[K], Q[m], g)
        C.ccx(h0, h1, nz)
        C.x(h0)
        C.x(h1)
        C.cx(h1, k2)
        and_into(B[K], Q[m], g)
        C.ccx(h0, g, k2)
        and_into(B[K], Q[m], g)
        C.free([nz])

        # k1 becomes [kappa = 1]; add kappa*c below bit K, carry cin into the run of ones above.
        C.cx(k2, k1)
        cin = g
        C.ccx(k1, A[K], cin)
        C.ccx(k2, B[K], cin)
        (u,) = C.alloc(1)
        and_into(cin, Q[m], u)
        C.ccx(h0, u, h1)
        C.cx(u, h0)
        and_into(cin, Q[m], u)
        C.free([u])
        # (T + kappa*c) >> n = kappa, so the top two bits now spell kappa.
        C.cx(k1, h0)
        C.cx(k2, h1)
        for j in reversed(range(m)):
            if Q[j] is None:
                C.cx(cin, L[K + j])
            else:
                C.ccx(cin, Q[j], L[K + j])
        # Q[j+1] = Q[j] & old L[K+j], and old L[K+j] = new L[K+j] ^ cin wherever Q[j] = 1.
        for w in L[K:]:
            C.cx(cin, w)
        for j in reversed(range(1, m)):
            C.ccx(Q[j], L[K + j], Q[j + 1])
            C.free([Q[j + 1]])
        for w in L[K:]:
            C.cx(cin, w)
        C.ccx(k1, A[K], cin)
        C.ccx(k2, B[K], cin)
        for i in reversed(range(K)):
            _uncarry(C, T, c, A, i + 1)
            _uncarry(C, T, 2 * c, B, i + 1)
            if A[i] is not None:
                C.ccx(k1, A[i], L[i])
            if B[i] is not None:
                C.ccx(k2, B[i], L[i])
            if c >> i & 1:
                C.cx(k1, L[i])
            if (2 * c) >> i & 1:
                C.cx(k2, L[i])
        C.cx(k2, k1)

        # kappa is odd iff r_0 != ctrl*R_0. kappa = 2 iff k1 and kappa even, and
        # then also iff ctrl and R > r, since r = 2S + R - 2q < R exactly then.
        p = cin
        C.cx(L[0], p)
        C.ccx(ctrl, R[0], p)
        C.x(p)
        C.ccx(k1, p, k2)
        C.x(p)
        C.cx(p, k1)
        (mk,) = C.alloc(1)
        C.x(p)
        C.ccx(ctrl, p, mk)
        gt_into(C, R, L, k1, ctrl=mk)
        C.ccx(ctrl, p, mk)
        C.x(p)
        C.ccx(ctrl, R[0], p)
        C.cx(L[0], p)
        C.free([mk, p, k1, k2, h0, h1])
        S[:] = L

    def cmodadd_const(self, C, ctrl, value, b):
        """b = (b + ctrl*value) mod q."""
        value %= self.q
        if value == 0:
            return
        (top,) = C.alloc(1)
        add_const(C, b, value, ctrl, top)
        _reduce(C, b + [top], self._wrap())
        gt_const_into(C, value, b, top, ctrl)
        C.free([top])

    def cneg(self, C, ctrl, b):
        """b = q - b when ctrl, valid for b in [1, q). A bijection on all of
        {0,1}^n either way, since it is ~(b + 2^n - q - 1) mod 2^n."""
        add_const(C, b, (1 << self.n) - self.q - 1, ctrl=ctrl)
        for w in b:
            C.cx(ctrl, w)

"""Exact, complete, controlled fixed-addend point addition:
|c>|R>|0> -> |c>|R + c*A>|0> for every R in E(F_q) encoded as little-endian
affine (x, y) with infinity = (0, 0)."""
from .arith import Mod, cswap, gt_into, load_const, sub
from .builder import Circuit
from .curve import Curve
from .model import iteration_bound, width_schedule


class DialogMultiplier:
    """M(x, y) = (x, x*y mod q) for x != 0, M(0, y) = (0, y); see model.py."""

    def __init__(self, C, mod):
        self.C, self.mod = C, mod
        q, n = mod.q, mod.n
        self.iters = iteration_bound(q)
        self.widths = width_schedule(q, self.iters)
        self.zbit = next(j for j in range(1, n) if q >> j & 1)
        self._euclid_steps = {}
        self._bezout_step = None
        self.block = None

    def _euclid_step(self, w):
        if w not in self._euclid_steps:
            C = self.C

            def step(U, V, g):
                g0, g1 = g
                C.cx(V[0], g0)
                gt_into(C, U[:w], V[:w], g1, ctrl=g0)
                cswap(C, g1, U[:w], V[:w])
                sub(C, U[:w], V[:w], ctrl=g0)
                V[:] = V[1:] + V[:1]

            n = self.mod.n
            self._euclid_steps[w] = C.record_shape(step, [n, n, 2])
        return self._euclid_steps[w]

    def _bezout(self):
        if self._bezout_step is None:
            C, mod, n = self.C, self.mod, self.mod.n

            def step(g, R, S):
                g0, g1 = g
                mod.moddbl(C, S)
                mod.modadd(C, R, S, ctrl=g0)
                cswap(C, g1, R, S)

            self._bezout_step = C.record_shape(step, [2, n, n])
        return self._bezout_step

    def _body(self, Xr, Yr):
        C, q, n, T = self.C, self.mod.q, self.mod.n, self.iters
        U = C.alloc(n)
        load_const(C, U, q)
        G = C.alloc(2 * T)

        def euclid(U, V, G):
            for i in range(T):
                C.apply(self._euclid_step(self.widths[i]), [U, V, G[2 * i:2 * i + 2]])

        E = C.run(euclid, [U, Xr, G])
        S = C.alloc(n)
        bz = self._bezout()
        for i in reversed(range(T)):
            C.apply(bz, [G[2 * i:2 * i + 2], Yr, S])
        cswap(C, U[self.zbit], Yr, S)
        C.apply(E, [U, Xr, G], inverse=True)
        load_const(C, U, q)
        C.free(U + G + Yr)
        Yr[:] = S

    def apply(self, Xr, Yr, inverse=False):
        if self.block is None:
            self.block = self.C.record(self._body, [Xr, Yr])
        self.C.apply(self.block, [Xr, Yr], inverse)


def eq_const_into(C, ctrl, bits, value, target):
    """target ^= ctrl AND [bits == value]."""
    flips = [w for i, w in enumerate(bits) if not value >> i & 1]
    for w in flips:
        C.x(w)
    chain = []
    acc = ctrl
    for w in bits[:-1]:
        (a,) = C.alloc(1)
        C.ccx(acc, w, a)
        chain.append((acc, w, a))
        acc = a
    C.ccx(acc, bits[-1], target)
    for acc, w, a in reversed(chain):
        C.ccx(acc, w, a)
    C.free([a for _, _, a in chain])
    for w in flips:
        C.x(w)


def encode(P, n):
    x, y = (0, 0) if P is None else P
    return x | y << n


class ControlledPointAdd:
    def __init__(self, curve: Curve, A):
        assert A is not None and curve.on_curve(A)
        self.curve, self.A = curve, A
        self.mod = Mod(curve.q)
        self.n = self.mod.n
        neg = curve.neg
        twoA = curve.add(A, A)
        # (exceptional input, its sum with A); the generic kernel divides by zero on these.
        self.exceptions = [(None, A), (A, twoA), (neg(A), None), (neg(twoA), neg(A))]
        enc = [encode(e, self.n) for e, _ in self.exceptions]
        assert len(set(enc)) == 4, "addend order must exceed 3"

    def _mulsub(self, gen, Yr, Xr):
        """X = X - gen * Y^2 mod q, Y untouched."""
        C, mod, n = self.C, self.mod, self.n
        T = C.alloc(n)
        for y, t in zip(Yr, T):
            C.cx(y, t)
        if self._modadd_blk is None:
            self._modadd_blk = C.record_shape(lambda k, a, b: mod.modadd(C, a, b, ctrl=k[0]), [1, n, n])
            self._dbl = C.record_shape(lambda s: mod.moddbl(C, s), [n])
        for i in range(n):
            (k,) = C.alloc(1)
            C.ccx(gen, Yr[i], k)
            C.apply(self._modadd_blk, [[k], T, Xr], inverse=True)
            C.ccx(gen, Yr[i], k)
            C.free([k])
            C.apply(self._dbl, [T])
        for _ in range(n):
            C.apply(self._dbl, [T], inverse=True)
        for y, t in zip(Yr, T):
            C.cx(y, t)
        C.free(T)

    def _kernel(self, gen, Xr, Yr):
        C, mod, q = self.C, self.mod, self.curve.q
        a, b = self.A
        mod.cmodadd_const(C, gen, q - a, Xr)
        mod.cmodadd_const(C, gen, q - b, Yr)
        self.M.apply(Xr, Yr, inverse=True)
        mod.cmodadd_const(C, gen, 3 * a, Xr)
        self._mulsub(gen, Yr, Xr)
        self.M.apply(Xr, Yr)
        mod.cneg(C, gen, Xr)
        mod.cmodadd_const(C, gen, q - b, Yr)
        mod.cmodadd_const(C, gen, a, Xr)

    def build(self, hierarchical=False):
        self.C = C = Circuit(hierarchical)
        n = self.n
        (c,) = C.alloc(1)
        Xr, Yr = C.alloc(n), C.alloc(n)
        declared = Xr + Yr
        self.M = DialogMultiplier(C, self.mod)
        self._modadd_blk = None

        flags = C.alloc(4)
        for f, (e, _) in zip(flags, self.exceptions):
            eq_const_into(C, c, Xr + Yr, encode(e, n), f)
        e, gen = C.alloc(2)
        for f in flags:
            C.cx(f, e)
        C.x(e)
        C.ccx(c, e, gen)
        for f, (e_in, e_out) in zip(flags, self.exceptions):
            load_const(C, Xr + Yr, encode(e_in, n) ^ encode(e_out, n), f)

        self._kernel(gen, Xr, Yr)

        C.ccx(c, e, gen)
        C.x(e)
        for f in flags:
            C.cx(f, e)
        C.free([e, gen])
        for f, (_, e_out) in zip(flags, self.exceptions):
            eq_const_into(C, c, Xr + Yr, encode(e_out, n), f)
        C.free(flags)

        _restore_layout(C, Xr + Yr, declared)
        self.layout = {"control": [c], "x": declared[:n], "y": declared[n:]}
        return C


def _restore_layout(C, cur, declared):
    """SWAP logical bit i from wire cur[i] onto declared[i]. A displaced wire
    holds either another live bit, which is tracked, or a zero ancilla."""
    cur = list(cur)
    where = {w: i for i, w in enumerate(cur)}
    for i, d in enumerate(declared):
        s = cur[i]
        if s == d:
            continue
        C.swap(s, d)
        j = where.pop(d, None)
        del where[s]
        if j is not None:
            cur[j] = s
            where[s] = j
        cur[i] = d
        where[d] = i

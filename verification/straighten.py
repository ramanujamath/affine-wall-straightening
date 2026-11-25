"""
Reference implementation of finite affine-wall straightening for
rational Hurwitz--Lerch multiple polylogarithms (two-chamber form).

    H_s(z; t) = sum_{n_1 > ... > n_r > 0} prod_j z_j^{n_j} / (n_j - t_j)^{s_j}

The algorithm STRAIGHTEN acts on admissible affine-chain sums

    A = sum_{x_k > L_k, x_k - x_{k+1} > D_k} prod_k w_k^{x_k} R_k(x_k)

and returns a formal K-linear combination of symbols Li_h(beta_C * omega^l).
A symbol is encoded as (h, args) with h a tuple of exponents and args a tuple
of (block interval (lo, hi), l mod q); the argument of (I, l) is
prod_{j in I} alpha_j * omega^l.

Two coefficient rings are provided:
  * Numeric: complex numbers (colours with 0 < |z_j| < 1).
  * Cyclotomic: exact elements of Q[mu_M] (group ring), colours roots of unity.
"""
from fractions import Fraction as Fr
from collections import defaultdict
import itertools, math, cmath

# ----------------------------------------------------------------- rings
class NumRing:
    def __init__(self, z, q, alpha=None):
        self.q = q
        self.z = [complex(x) for x in z]
        self.alpha = [complex(a) for a in alpha] if alpha else [x ** (1.0 / q) for x in self.z]
        self.omega = cmath.exp(2j * math.pi / q)
    zero = staticmethod(lambda: 0j)
    one = staticmethod(lambda: 1 + 0j)
    add = staticmethod(lambda a, b: a + b)
    mul = staticmethod(lambda a, b: a * b)
    scal = staticmethod(lambda a, f: a * float(f))
    def zpow(self, I, n):
        v = 1 + 0j
        for j in range(I[0], I[1] + 1):
            v *= self.z[j] ** n
        return v
    def apow(self, I, n):
        v = 1 + 0j
        for j in range(I[0], I[1] + 1):
            v *= self.alpha[j] ** n
        return v
    def opow(self, n):
        return self.omega ** (n % self.q)
    iszero = staticmethod(lambda a: abs(a) < 1e-11)
    def arg(self, I, l):
        return self.apow(I, 1) * self.opow(l)


class CycRing:
    """Exact Q[mu_M], M = q*N.  xi_j = zeta_N^{a_j}, alpha_j = zeta_M^{a_j}."""
    def __init__(self, a, N, q):
        self.q, self.N, self.M = q, N, q * N
        self.a = list(a)
    def zero(self):
        return {}
    def one(self):
        return {0: Fr(1)}
    def add(self, x, y):
        r = dict(x)
        for k, v in y.items():
            r[k] = r.get(k, 0) + v
            if r[k] == 0:
                del r[k]
        return r
    def mul(self, x, y):
        r = {}
        for k1, v1 in x.items():
            for k2, v2 in y.items():
                k = (k1 + k2) % self.M
                r[k] = r.get(k, 0) + v1 * v2
        return {k: v for k, v in r.items() if v != 0}
    def scal(self, x, f):
        return {k: v * f for k, v in x.items() if v * f != 0}
    def _e(self, e):
        return {e % self.M: Fr(1)}
    def zpow(self, I, n):
        return self._e(self.q * n * sum(self.a[I[0]:I[1] + 1]))
    def apow(self, I, n):
        return self._e(n * sum(self.a[I[0]:I[1] + 1]))
    def opow(self, n):
        return self._e(n * self.N)
    def iszero(self, x):
        return len(x) == 0
    def arg_exp(self, I, l):
        return (sum(self.a[I[0]:I[1] + 1]) + l * self.N) % self.M

# ------------------------------------------------------- partial fractions
def partial_fractions(factors):
    """prod (x - p)^(-m)  ->  {(p, h): C} with sum C (x-p)^(-h)."""
    order = defaultdict(int)
    for p, m in factors:
        order[Fr(p)] += m
    poles = list(order)
    out = {}
    for a in poles:
        Ma = order[a]
        ser = [Fr(0)] * Ma
        ser[0] = Fr(1)
        for b in poles:
            if b == a:
                continue
            u, Mb = a - b, order[b]
            fac = [Fr((-1) ** n * math.comb(Mb + n - 1, n)) / u ** (Mb + n) for n in range(Ma)]
            ser = [sum(ser[i] * fac[n - i] for i in range(n + 1)) for n in range(Ma)]
        for k in range(1, Ma + 1):
            c = ser[Ma - k]
            if c != 0:
                out[(a, k)] = c
    return out

def floor(x):
    return math.floor(x)

def defect(L, T):
    """1_{x > L} = 1_{x > T} + sum_a sign * 1_{x = a}  (L, T integers)."""
    if L < T:
        return [(a, 1) for a in range(L + 1, T + 1)]
    if L > T:
        return [(a, -1) for a in range(T + 1, L + 1)]
    return []

# --------------------------------------------------------------- algorithm
class Chain:
    def __init__(self, blocks, letters, L, D, coeff):
        self.blocks = blocks      # list of intervals (lo, hi)
        self.letters = letters    # list of lists of (pole, exponent)
        self.L = L                # lower bounds  x_k > L_k
        self.D = D                # gaps  x_k - x_{k+1} > D_k
        self.coeff = coeff

STATS = defaultdict(int)

def straighten(R, ch, out=None):
    if out is None:
        out = defaultdict(R.zero)
    d = len(ch.blocks)
    if d == 0:
        key = ((), ())
        out[key] = R.add(out[key], ch.coeff)
        return out
    q = R.q
    # admissibility (P): no integer pole exceeds the lower bound
    for k in range(d):
        for p, _ in ch.letters[k]:
            if p.denominator == 1 and p > ch.L[k]:
                raise ValueError("inadmissible chain: pole %s > L=%s" % (p, ch.L[k]))
    pf = [list(partial_fractions(lt).items()) for lt in ch.letters]
    for combo in itertools.product(*pf):
        poles = [c[0][0] for c in combo]
        hs = [c[0][1] for c in combo]
        C = Fr(1)
        for c in combo:
            C *= c[1]
        base_coeff = R.scal(ch.coeff, C)
        Lt = [floor(p) for p in poles]
        Et = [floor(poles[k] - poles[k + 1]) for k in range(d - 1)]
        copt = [[None] + defect(ch.L[k], Lt[k]) for k in range(d)]
        eopt = [[None] + defect(ch.D[k], Et[k]) for k in range(d - 1)]
        for cs in itertools.product(*copt):
            for es in itertools.product(*eopt):
                if all(c is None for c in cs) and all(e is None for e in es):
                    fourier(R, ch, poles, hs, base_coeff, out)
                else:
                    stratum(R, ch, poles, hs, base_coeff, Lt, Et, cs, es, out)
    return out

def fourier(R, ch, poles, hs, coeff, out):
    q, d = R.q, len(poles)
    STATS['fourier'] += 1
    cs, kap = [], []
    pre = coeff
    for k, p in enumerate(poles):
        assert (p * q).denominator == 1, "pole not in (1/q)Z"
        c = int(p * q) % q
        kappa = Fr(c, q) - p
        assert kappa.denominator == 1
        cs.append(c)
        kap.append(int(kappa))
        pre = R.mul(pre, R.zpow(ch.blocks[k], -int(kappa)))
        pre = R.mul(pre, R.apow(ch.blocks[k], c))
    pre = R.scal(pre, Fr(q) ** (sum(hs) - d))
    for ls in itertools.product(range(q), repeat=d):
        e = sum(l * c for l, c in zip(ls, cs))
        key = (tuple(hs), tuple((ch.blocks[k], ls[k]) for k in range(d)))
        out[key] = R.add(out[key], R.mul(pre, R.opow(e)))

def stratum(R, ch, poles, hs, coeff, Lt, Et, cs, es, out):
    d = len(poles)
    sign = 1
    for c in cs:
        if c is not None:
            sign *= c[1]
    for e in es:
        if e is not None:
            sign *= e[1]
    # components of the contracted path
    comps, cur = [], [0]
    for k in range(d - 1):
        if es[k] is not None:
            cur.append(k + 1)
        else:
            comps.append(cur)
            cur = [k + 1]
    comps.append(cur)
    info = []
    for comp in comps:
        off = {comp[0]: 0}
        for k in comp[:-1]:
            off[k + 1] = off[k] - es[k][0]
        m0 = min(off.values())
        lam = {k: off[k] - m0 for k in comp}
        anch = [k for k in comp if cs[k] is not None]
        if anch:
            X = None
            for k in anch:
                v = cs[k][0] - lam[k]
                if X is None:
                    X = v
                elif X != v:
                    return  # incompatible
            for k in comp:
                if cs[k] is None and not (X + lam[k] > Lt[k]):
                    return
            info.append(dict(comp=comp, lam=lam, X=X, anch=True))
        else:
            lb = max(Lt[k] - lam[k] for k in comp)
            info.append(dict(comp=comp, lam=lam, lb=lb, anch=False))
    A = next((i for i, c in enumerate(info) if c['anch']), None)

    def val(k, x):  # summand of variable k at integer value x
        v = R.zpow(ch.blocks[k], x)
        den = Fr(x) - poles[k]
        assert den != 0, "pole hit"
        return R.scal(v, 1 / den ** hs[k])

    if A is None:
        right = R.one()
        nleft = len(info)
    else:
        # enumerate the bounded part to the right of (and including) comp A
        total = R.zero()
        def rec(i, prev_last_x, acc):
            nonlocal total
            if i == len(info):
                total = R.add(total, acc)
                return
            c = info[i]
            first = c['comp'][0]
            ub = None
            if prev_last_x is not None:
                # x_last(prev) - x_first > Et  =>  x_first <= x_last(prev) - Et - 1
                ub = prev_last_x - Et[first - 1] - 1 - c['lam'][first]
            if c['anch']:
                Xs = [c['X']] if (ub is None or c['X'] <= ub) else []
            else:
                Xs = range(c['lb'] + 1, ub + 1)
            for X in Xs:
                a = acc
                for k in c['comp']:
                    a = R.mul(a, val(k, X + c['lam'][k]))
                rec(i + 1, X + c['lam'][c['comp'][-1]], a)
        rec(A, None, R.one())
        right = total
        nleft = A
        if R.iszero(right):
            return
    # left chain (free components 0..nleft-1)
    blocks, letters, L, D = [], [], [], []
    pre = R.scal(R.mul(coeff, right), Fr(sign))
    for i in range(nleft):
        c = info[i]
        comp, lam = c['comp'], c['lam']
        blocks.append((ch.blocks[comp[0]][0], ch.blocks[comp[-1]][1]))
        letters.append([(poles[k] - lam[k], hs[k]) for k in comp])
        for k in comp:
            pre = R.mul(pre, R.zpow(ch.blocks[k], lam[k]))
        lb = c['lb']
        if i == nleft - 1 and A is not None:
            a = info[A]
            xfirstA = a['X'] + a['lam'][a['comp'][0]]
            last = comp[-1]
            lb = max(lb, xfirstA + Et[last] - lam[last])
        L.append(lb)
        if i > 0:
            pc = info[i - 1]
            plast = pc['comp'][-1]
            D.append(Et[plast] - pc['lam'][plast] + lam[comp[0]])
    straighten(R, Chain(blocks, letters, L, D, pre), out)

def straighten_H(R, s, t):
    r = len(s)
    ch = Chain([(j, j) for j in range(r)], [[(Fr(t[j]), s[j])] for j in range(r)],
               [0] * r, [0] * (r - 1), R.one())
    out = straighten(R, ch)
    return {k: v for k, v in out.items() if not R.iszero(v)}

# --------------------------------------------------------- numerics
import numpy as np

def li_num(h, ys, N=4000):
    n = np.arange(1, N + 1, dtype=float)
    inner = np.ones(N, dtype=complex)
    for hk, y in reversed(list(zip(h, ys))):
        term = np.power(complex(y), n) / n ** hk * inner
        inner = np.concatenate([[0], np.cumsum(term)[:-1]])
    return np.sum(term)

def H_num(z, s, t, N=4000):
    n = np.arange(1, N + 1, dtype=float)
    inner = np.ones(N, dtype=complex)
    for zk, sk, tk in reversed(list(zip(z, s, t))):
        term = np.power(complex(zk), n) / (n - float(tk)) ** sk * inner
        inner = np.concatenate([[0], np.cumsum(term)[:-1]])
    return np.sum(term)

def eval_out(R, out, N=4000):
    tot = 0j
    for (h, args), c in out.items():
        if len(h) == 0:
            tot += c
        else:
            tot += c * li_num(h, [R.arg(I, l) for I, l in args], N)
    return tot

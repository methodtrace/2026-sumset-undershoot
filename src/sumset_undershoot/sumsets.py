"""Sumsets h A_0, the intersection form I_h, residue intervals, the spectrum, similar
ordering, the common diagram of four-element sets, corners of the classes, failing levels,
the correction Psi_h and the hinge formulas.

A_0 = {0 = a_0 < a_1 < ... < a_k = d} with gcd 1. S_1 = <a_1, ..., a_k>,
S_2 = <d - a_0, ..., d - a_{k-1}>. h A_0 is the set of sums of at most h positive parts.
"""
from math import gcd
from functools import reduce

from .semigroup import apery, conductor, membership, run_functional, mass_definition


# ---------------------------------------------------------------- basic data

def normalize(parts):
    """(g, A_0) for a physical menu of positive parts: A_0 = {0} u parts/g."""
    ps = sorted(set(int(x) for x in parts if int(x) > 0))
    g = reduce(gcd, ps)
    return g, [0] + [x // g for x in ps]


def S1(A0):
    return sorted(set(a for a in A0 if a > 0))


def S2(A0):
    d = max(A0)
    return sorted(set(d - a for a in A0 if a != d))


class SetData:
    """Cached data of A_0: conductors, masses, membership of S_1, S_2 and min part counts."""

    def __init__(self, A0, upto=None):
        self.A0 = sorted(A0)
        self.d = max(A0)
        self.g1 = S1(A0)
        self.g2 = S2(A0)
        self.c1 = conductor(self.g1)
        self.c2 = conductor(self.g2)
        self.M1 = mass_definition(self.g1)
        self.M2 = mass_definition(self.g2)
        self.upto = 0
        self.mem1 = self.mem2 = None
        self.ensure(upto or (self.c1 + self.c2 + 2 * self.d))

    def ensure(self, upto):
        if upto > self.upto:
            self.upto = upto
            self.mem1 = membership(self.g1, upto)
            self.mem2 = membership(self.g2, upto)
            self.ml = min_parts(self.A0, upto)

    def in_S1(self, n):
        self.ensure(n)
        return n >= 0 and self.mem1[n] == 1

    def in_S2(self, n):
        self.ensure(n)
        return n >= 0 and self.mem2[n] == 1

    def in_sumset(self, n, h):
        """n in h A_0 iff n is a sum of at most h positive parts."""
        self.ensure(n)
        return 0 <= n and self.ml[n] <= h

    def in_I(self, n, h):
        return 0 <= n <= h * self.d and self.in_S1(n) and self.in_S2(h * self.d - n)

    def sumset(self, h):
        self.ensure(h * self.d)
        return [n for n in range(h * self.d + 1) if self.ml[n] <= h]

    def I_set(self, h):
        self.ensure(h * self.d)
        return [n for n in range(h * self.d + 1) if self.in_I(n, h)]


INF = 1 << 60


def min_parts(A0, upto):
    """ml[n] = least number of positive parts of A_0 with sum n (INF if none)."""
    gens = sorted(set(a for a in A0 if a > 0), reverse=True)
    ml = [INF] * (upto + 1)
    ml[0] = 0
    for n in range(1, upto + 1):
        best = INF
        for a in gens:
            if a <= n and ml[n - a] + 1 < best:
                best = ml[n - a] + 1
        ml[n] = best
    return ml


def undershoots(elements_sorted, N):
    """[D_E(T) for T = 0..N] for a set E containing 0."""
    out = []
    es = set(elements_sorted)
    last = 0
    for T in range(N + 1):
        if T in es:
            last = T
        out.append(T - last)
    return out


def D_E(es, n):
    return n - max(x for x in es if x <= n)


# ---------------------------------------------------------------- residue intervals

def moduli(A0):
    sd = SetData(A0)
    return [m for m in range(1, sd.d + 1) if sd.d % m == 0 and sd.in_S1(m) and sd.in_S2(m)]


def uv(A0, m):
    """u_s = (s + w1_{-s}) / m and v_s = (w2_s - s) / m for 0 <= s < m."""
    w1 = apery(S1(A0), m)
    w2 = apery(S2(A0), m)
    u, v = [], []
    for s in range(m):
        a = s + w1[(-s) % m]
        b = w2[s] - s
        assert a % m == 0 and b % m == 0
        u.append(a // m)
        v.append(b // m)
    return u, v


def B_set(u, k):
    return [s for s in range(1, len(u)) if u[s] > k]


def R_set(v, j):
    return [s for s in range(1, len(v)) if v[s] > j]


def holes(u, v, N, k):
    return sorted(set(B_set(u, k)) | set(R_set(v, N - k)))


def counting_undershoot(u, v, m, N, k, s):
    """Theorem prop:counting(a): #{t in [s, m-1] : [s, t] in H_k}."""
    H = set(holes(u, v, N, k))
    t = s
    while t <= m - 1 and t in H:
        t += 1
    return t - s


def spectrum(u, v):
    U, V = max(u), max(v)
    beta = [run_functional(B_set(u, k)) for k in range(1, max(U, 1) + 1)]
    gamma = [run_functional(R_set(v, j)) for j in range(0, V + 1)]
    lam = sorted([x for x in beta + gamma if x > 0], reverse=True)
    return beta, gamma, lam


def similarly_ordered(u, v):
    m = len(u)
    for s in range(m):
        for t in range(m):
            if u[s] > u[t] and v[s] < v[t]:
                return False
    return True


def inversion_intervals(u, v):
    out = set()
    m = len(u)
    for s in range(m):
        for t in range(m):
            if u[s] > u[t] and v[s] < v[t]:
                out.add((u[t] + v[s], u[s] + v[t] - 2))
    return sorted(out)


# ---------------------------------------------------------------- four-element diagram

def diagram(A0):
    """Lemma lem:diagram data for a four-element set: list of (r, t_r, x_r, y_r) with
    x_r, y_r rationals; integral and nonnegative iff A_0 is CM."""
    from fractions import Fraction as F
    assert len(A0) == 4
    _, a1, a2, d = A0
    f = apery(S1(A0), d)
    g = apery(S2(A0), d)
    out = []
    for r in range(d):
        t = F(f[r] + g[(-r) % d], d)
        y = (f[r] - a1 * t) / (a2 - a1)
        x = t - y
        out.append((r, t, x, y))
    return out


# ---------------------------------------------------------------- corners and failing levels

def class_corners(sd, s, kmax):
    """Minimal points (kappa, lambda) of C_s = {(k, l) : k d - s in (k + l) A_0}, k <= kmax."""
    d = sd.d
    sd.ensure(kmax * d)
    pts = []
    best_l = INF
    for k in range(1, kmax + 1):
        n = k * d - s
        if n < 0:
            continue
        mlen = sd.ml[n]
        if mlen >= INF:
            continue
        l = max(0, mlen - k)
        if l < best_l:
            pts.append((k, l))
            best_l = l
    return pts


def corner_kmax(sd):
    """Beyond this k every class's least lambda is stable (each minimal representation of a
    large target uses a part d), so all corners have kappa <= kmax."""
    gens = sd.g1
    second = gens[-2] if len(gens) > 1 else 0
    return ((sd.d - 1) * second + sd.d - 1) // sd.d + 2


def all_corners(sd):
    km = corner_kmax(sd)
    return {s: class_corners(sd, s, km) for s in range(1, sd.d)}


def failing_levels_by_class(sd, corners):
    """Proposition prop:failing-interval: failing levels of class s as an interval (or None)."""
    out = {}
    for s in range(1, sd.d):
        C = corners[s]
        if len(C) < 2:
            out[s] = None
            continue
        lo = C[0][0] + C[-1][1]
        hi = max(C[j + 1][0] + C[j][1] - 2 for j in range(len(C) - 1))
        out[s] = (lo, hi)
    return out


# ---------------------------------------------------------------- windows

def o2(sd, y):
    j = 0
    while not sd.in_S2(y + j):
        j += 1
    return j


def u1(sd, n):
    return n - max(x for x in range(n + 1) if sd.in_S1(x))


def psi(sd, h, W):
    """Proposition prop:correction: Psi_h(W) = sum_{z in F_h} #(W cap [z, z+)) (z - z-)."""
    hd = h * sd.d
    I = sd.I_set(h)
    Iset = set(I)
    A = set(sd.sumset(h))
    F_h = sorted(Iset - A)
    Wset = set(W)
    total = 0
    for z in F_h:
        zp = min(x for x in I if x > z)
        zm = max(x for x in A if x < z)
        total += sum(1 for n in range(z, zp) if n in Wset) * (z - zm)
    assert all(0 <= n <= hd for n in W)
    return total


def minimal_points(points):
    pts = sorted(set(points))
    out = []
    best = INF
    for (k, l) in pts:
        if l < best:
            out.append((k, l))
            best = l
    return out


def hinge_count(C, h, l0, q0):
    """Lemma lem:hinge: H_C(h; l0, q0) for minimal points listed with kappa increasing."""
    def pos(x):
        return x if x > 0 else 0
    total = sum(pos(min(q0, h - lam) - max(l0, kap) + 1) for (kap, lam) in C)
    total -= sum(pos(min(q0, h - C[j - 1][1]) - max(l0, C[j][0]) + 1) for j in range(1, len(C)))
    return total


def hinge_window(sd, corners, h, L, U, kind="B"):
    """Theorem thm:hinge: (window total, Psi_h([L, U])) from corners."""
    d = sd.d
    u, v = uv(sd.A0, d)
    total = 0
    psi_total = 0
    for s in range(1, d):
        ls = max(1, -((-(L + s)) // d))
        qs = min(h, (U + s) // d)
        unionB = []
        unionI = []
        for t in range(s, d):
            unionB.extend(corners[t])
            unionI.append((u[t], v[t]))
            CB = minimal_points(unionB)
            CI = minimal_points(unionI)
            hb = hinge_count(CB, h, ls, qs)
            hi = hinge_count(CI, h, ls, qs)
            total += max(qs - ls + 1, 0) - hb
            psi_total += hi - hb
    return total, psi_total


def hinge_sum(sd, corners, h, L, U, at_level):
    """sum_{s <= t} H_{C^B_{s,t}}(at_level; l_s, q_s) with l_s, q_s computed at cap h."""
    d = sd.d
    out = 0
    for s in range(1, d):
        ls = max(1, -((-(L + s)) // d))
        qs = min(h, (U + s) // d)
        unionB = []
        for t in range(s, d):
            unionB.extend(corners[t])
            out += hinge_count(minimal_points(unionB), at_level, ls, qs)
    return out


def four_corners_formula(A0):
    """Corners of every class of a four-element set from (S_1, S_2, d) alone (not stated in the paper)."""
    from math import gcd as g_
    _, a1, a2, d = A0
    out = {}
    for s in range(1, d):
        pts = []
        for i in range(d // g_(a1, d)):
            for j in range(d // g_(a2, d)):
                tot = i * a1 + j * a2 + s
                if tot % d == 0:
                    kap = tot // d
                    pts.append((kap, i + j - kap))
        out[s] = minimal_points(pts)
    return out

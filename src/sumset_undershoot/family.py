"""The family A_0 = {0, a, b, 2b} (coprime 0 < a < b): spectrum, totals, end regimes,
the three-distance description of phi = phi_{b,a}, quadratic pieces and Euclidean evaluation.
"""
from math import comb

from .twogen import phi_sequence, Phi, mass, P_closed


def phi_family(a, b):
    """phi(j) = phi_{b,a}(j) for j = 1..b, with phi(j) = 0 for j >= b (list indexed from 1)."""
    seq = phi_sequence(b, a)
    return [None] + seq + [0, 0]


def M4(a, b):
    return mass(a, b) + mass(b, 2 * b - a)


def Tot_formula(a, b, h):
    """Equation eq:family-tot; Tot(h) = M4 for h >= b - 1."""
    ph = phi_family(a, b)
    if h == 0:
        return 0
    if h >= b - 1:
        return comb(b, 2) + 2 * sum(ph[j] for j in range(2, b)) if b >= 2 else 0
    return comb(b, 2) + 2 * sum(ph[j] for j in range(2, h + 1)) + ph[h + 1]


def end_regime_one(a, b, h):
    return (comb(b, 2) + 2 * sum(comb(b - l * a, 2) for l in range(1, h))
            + comb(b - h * a, 2) + h * h * comb(a, 2))


def g_chain(n, e):
    Q, R = divmod(n, e)
    return R * comb(Q + 2, 2) + (e - R) * comb(Q + 1, 2)


def end_regime_two_gap(a, b, h):
    abar = pow(a, -1, b)
    e = min(abar, b - abar)
    t = b - 1 - h
    return g_chain(t, e) + 2 * sum(g_chain(n, e) for n in range(1, t)), e


def convergent_data(a, b):
    """Partial quotients c_1..c_K of a/b, continuants q_{-1..K}, eps_{-1..K}."""
    cs = []
    x, y = b, a
    while y:
        c, rem = divmod(x, y)
        cs.append(c)
        x, y = y, rem
    K = len(cs)
    p = {-1: 1, 0: 0}
    q = {-1: 0, 0: 1}
    for j in range(1, K + 1):
        p[j] = cs[j - 1] * p[j - 1] + p[j - 2]
        q[j] = cs[j - 1] * q[j - 1] + q[j - 2]
    eps = {j: abs(q[j] * a - p[j] * b) for j in range(-1, K + 1)}
    return cs, q, eps


def cf_index(a, b, h):
    """Equation eq:cf-index: (k, m, r) with h - 1 = m q_k + q_{k-1} + r."""
    cs, q, _ = convergent_data(a, b)
    K = len(cs)
    for k in range(0, K):
        for m in range(1, cs[k] + 1):
            base = m * q[k] + q[k - 1]
            if base <= h - 1 < base + q[k]:
                return k, m, h - 1 - base
    raise ValueError("no index")


def three_distance_gaps(a, b, h):
    """Theorem thm:three-distance: {gap length: count} for 2 <= h <= b."""
    _, q, eps = convergent_data(a, b)
    k, m, r = cf_index(a, b, h)
    out = {}
    for length, count in ((eps[k], h - q[k]), (eps[k - 1] - m * eps[k], r + 1),
                          (eps[k - 1] - (m - 1) * eps[k], q[k] - r - 1)):
        if count:
            out[length] = out.get(length, 0) + count
    return out


def family_euclid_mass(a, b):
    """M4(a, b) by Euclid's algorithm along the continued fraction of a/b (not stated in the paper)."""
    cs, q, eps = convergent_data(a, b)
    K = len(cs)
    total = sum(P_closed(eps[i - 1], eps[i], cs[i - 1]) * q[i - 1] ** 2 for i in range(1, K + 1))
    return 2 * total - comb(b, 2)


def family_step(a, r, q_):
    """M4(a, qa + r) as an exact rational, a cubic in q (not stated in the paper)."""
    from fractions import Fraction as F
    base = M4(a % r, r) if r >= 2 else 0
    return (base + F(a * a, 3) * q_ ** 3 + (2 * Phi(a, r) + F(a * (2 * r - 1), 2)) * q_ ** 2
            + (F(a * a, 6) + r * r - r + 4 * mass(r, a)) * q_)

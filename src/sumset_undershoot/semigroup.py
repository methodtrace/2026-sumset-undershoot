"""Numerical semigroups: membership, Apery sets, conductor, gaps and the mass.

The mass of S is M(S) = sum over gaps x of (x - max{s in S : s < x}); a maximal run of L
consecutive gaps contributes L(L+1)/2, so M(S) = f(N minus S) with f the run functional.
"""
from math import gcd
from functools import reduce


def gcd_list(values):
    return reduce(gcd, values, 0)


def apery(gens, m):
    """Apery set (w_r) of <gens> with respect to m: w_r = least element congruent to r.

    Computed by a label-correcting relaxation over residues (Bellman-Ford style, no heap),
    so that it is independent of the Dijkstra routine in the intersection engine.
    """
    INF = None
    w = [INF] * m
    w[0] = 0
    changed = True
    while changed:
        changed = False
        for r in range(m):
            if w[r] is None:
                continue
            for g in gens:
                r2 = (r + g) % m
                c = w[r] + g
                if w[r2] is None or c < w[r2]:
                    w[r2] = c
                    changed = True
    assert all(x is not None for x in w), "generators have a common factor"
    return w


def membership(gens, upto):
    """bytearray mem with mem[n] = 1 iff n in <gens>, for 0 <= n <= upto (direct closure)."""
    mem = bytearray(upto + 1)
    mem[0] = 1
    g = sorted(set(x for x in gens if x > 0))
    for n in range(1, upto + 1):
        for a in g:
            if a > n:
                break
            if mem[n - a]:
                mem[n] = 1
                break
    return mem


def conductor(gens):
    """c(S) = Frobenius number + 1 (0 for S = N)."""
    gens = [x for x in gens if x > 0]
    m = min(gens)
    w = apery(gens, m)
    return max(w) - m + 1


def gaps(gens):
    c = conductor(gens)
    mem = membership(gens, max(c, 0))
    return [n for n in range(c) if not mem[n]]


def run_functional(positions):
    """f(X): number of pairs s <= t with [s, t] contained in X, for a finite set of integers."""
    xs = sorted(set(positions))
    total = 0
    run = 0
    prev = None
    for x in xs:
        run = run + 1 if prev is not None and x == prev + 1 else 1
        total += run
        prev = x
    return total


def mass_definition(gens):
    """M(S) straight from the definition: distance from each gap down to S."""
    c = conductor(gens)
    mem = membership(gens, max(c, 0))
    total = 0
    last = 0
    for n in range(c):
        if mem[n]:
            last = n
        else:
            total += n - last
    return total


def mass_apery(gens, m):
    """M(S) from one Apery set (Theorem thm:apery-mass):
    M(S) = sum over 1 <= s <= t <= m-1 of (min_{s<=rho<=t} u_rho - 1),
    u_t = (t + w_{-t mod m}) / m."""
    w = apery(gens, m)
    u = [None] + [(t + w[(-t) % m]) // m for t in range(1, m)]
    for t in range(1, m):
        assert (t + w[(-t) % m]) % m == 0 and u[t] >= 1
    total = 0
    for s in range(1, m):
        cur = None
        for t in range(s, m):
            cur = u[t] if cur is None or u[t] < cur else cur
            total += cur - 1
    return total


def two_generator_membership(p, q, upto):
    """Membership array of <p, q> on [0, upto] built as the union of the progressions
    j*q + p*N, 0 <= j < p (every element i*p + j*q reduces to one with j < p)."""
    mem = bytearray(upto + 1)
    for j in range(p):
        start = j * q
        if start > upto:
            break
        count = (upto - start) // p + 1
        mem[start::p] = b"\x01" * count
    return mem


def two_generator_runs(p, q):
    """Multiset {length: count} of maximal runs of gaps of <p, q>, from the definition."""
    if p == 1 or q == 1:
        return {}
    upto = (p - 1) * (q - 1)
    mem = two_generator_membership(p, q, upto)
    runs = {}
    for chunk in bytes(mem).split(b"\x01"):
        L = len(chunk)
        if L:
            runs[L] = runs.get(L, 0) + 1
    return runs


def mass_from_runs(runs):
    return sum(runs[L] * L * (L + 1) // 2 for L in runs)


def two_generator_mass(p, q):
    """M(p, q) = M(<p, q>) by the definition (via the run multiset)."""
    return mass_from_runs(two_generator_runs(p, q))

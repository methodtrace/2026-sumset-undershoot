"""Engine C (two numerical semigroups).

Gap sets of numerical semigroups via the Apery set with respect to the smallest generator
(shortest paths on residues), then
    I_h = {n in [0, hd] : n not a gap of S_1 = <A>, hd - n not a gap of S_2 = <d - A_0>}.
For Cohen-Macaulay sets I_h equals h A_0 at every level.
"""
import heapq


def apery_numerical(gens):
    gens = sorted(set(gens))
    m = gens[0]
    w = [None] * m
    w[0] = 0
    pq = [(0, 0)]
    while pq:
        c, r = heapq.heappop(pq)
        if c != w[r]:
            continue
        for g in gens[1:]:
            r2 = (r + g) % m
            c2 = c + g
            if w[r2] is None or c2 < w[r2]:
                w[r2] = c2
                heapq.heappush(pq, (c2, r2))
    if any(x is None for x in w):
        raise ValueError("gcd of generators is not 1")
    return w, m


def gaps(gens):
    if 1 in gens:
        return []
    w, m = apery_numerical(gens)
    out = []
    for x in w:
        n = x - m
        while n > 0:
            out.append(n)
            n -= m
    return sorted(out)


def frobenius(gens):
    w, m = apery_numerical(gens)
    return max(w) - m


def reflected_generators(A0):
    d = max(A0)
    return sorted(set(d - a for a in A0 if a != d))


def sumset_bits_intersection(A0, h, G_bottom=None, G_top=None):
    d = max(A0)
    pos = [a for a in A0 if a > 0]
    Gb = gaps(pos) if G_bottom is None else G_bottom
    Gt = gaps(reflected_generators(A0)) if G_top is None else G_top
    N = h * d
    bits = (1 << (N + 1)) - 1
    for g in Gb:
        if g <= N:
            bits &= ~(1 << g)
    for g in Gt:
        if g <= N:
            bits &= ~(1 << (N - g))
    return bits


def undershoot_mass(gap_list):
    """Sum over the gaps of the distance down to the nearest semigroup element."""
    gs = set(gap_list)
    total = 0
    for g in gap_list:
        k = 1
        while (g - k) in gs:
            k += 1
        total += k
    return total

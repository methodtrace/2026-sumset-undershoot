#!/usr/bin/env python3
"""Sumsets: residue intervals and exact counting, the spectrum and monotone totals, strict
levels, the threshold, and similarly ordered four-element sets.

Paper results: thm:cm (as the Cohen-Macaulay test, three engines), lem:residue,
prop:counting (a)-(c), lem:spectrum-mass, lem:merge, thm:spectrum (1)-(3),
prop:strict-levels, ex:spectrum (a) and (b) {0,2,8,11}, eq:reg-cm,
prop:threshold and the {0,1,2,5} and {0,3,4,6} remarks, lem:diagram, thm:similar,
prop:shape, prop:families (a)-(d), and the three examples of the introduction
({0,5,7,14}, {0,2,8,11}, {0,4,9,13}).
"""
import os
import sys
from itertools import combinations
from math import comb, gcd
from functools import reduce

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))

from sumset_undershoot.report import Checker, standard_parser
from sumset_undershoot.engines import direct as engine_a, homogenized as engine_b, intersection as engine_c
from sumset_undershoot.semigroup import run_functional
from sumset_undershoot.sumsets import (SetData, moduli, uv, holes, counting_undershoot, spectrum,
                                    similarly_ordered, inversion_intervals, diagram, min_parts, INF)


def sets_with(k, dmax, dmin=2):
    """All A_0 = {0 < a_1 < ... < a_{k-1} = d} with gcd 1 and dmin <= d <= dmax, in order."""
    for d in range(dmin, dmax + 1):
        for mid in combinations(range(1, d), k - 2):
            A0 = [0] + list(mid) + [d]
            if reduce(gcd, A0) == 1:
                yield A0


class Totals:
    """Tot(h) and D_h from engine A; I_h from engine C."""

    def __init__(self, A0, H):
        self.A0 = A0
        self.d = max(A0)
        self.series = engine_a.sumset_series(A0, H)
        self.gb = engine_c.gaps([a for a in A0 if a > 0])
        self.gt = engine_c.gaps(engine_c.reflected_generators(A0))

    def D(self, h):
        return engine_a.undershoot_scan(self.series[h], 0, h * self.d)

    def Tot(self, h):
        return sum(self.D(h))

    def I_bits(self, h):
        return engine_c.sumset_bits_intersection(self.A0, h, self.gb, self.gt)

    def D_I(self, h):
        return engine_a.undershoot_scan(self.I_bits(h), 0, h * self.d)


def is_cm_three_ways(ck, A0, H):
    ap = engine_b.apery(A0)
    tb = engine_b.is_cm(A0, ap)
    sa = [engine_a.sumset_bits(A0, h) for h in range(H + 1)]
    ta = all(sa[h] == engine_c.sumset_bits_intersection(A0, h) for h in range(1, H + 1))
    if tb:
        ck.check("thm:cm engine B CM implies h A_0 = I_h (engines A, C) for h <= bound", ta, A0)
    # engine B rebuilds h A_0 and the undershoot from the homogenized Apery set alone
    sb = engine_b.sumset_series_from_apery(A0, H, ap)
    ck.check("engine B h A_0 from the homogenized Apery set equals engine A, h <= bound", sb == sa, A0)
    d = max(A0)
    for h in range(1, H + 1):
        scan = engine_a.undershoot_scan(sa[h], 0, h * d)
        ck.check("engine B undershoot formula equals the engine A scan on [0, hd]",
                 all(engine_b.undershoot_formula(A0, h, T, ap) == scan[T] for T in range(h * d + 1)), (A0, h))
    return tb


def sigma_by_definition(sd, Tt, Hs):
    """Least h0 such that for every h0 <= h <= Hs the regularity decomposition holds with a
    nonempty central interval (h A_0 from engine A)."""
    d, c1, c2 = sd.d, sd.c1, sd.c2
    low = {n for n in range(0, max(c1 - 1, 0)) if sd.in_S1(n)}
    top = {n for n in range(0, max(c2 - 1, 0)) if sd.in_S2(n)}

    def holds(h):
        if h * d < c1 + c2:
            return False
        target = set(low) | set(range(c1, h * d - c2 + 1)) | {h * d - x for x in top}
        return set(engine_a.elements(Tt.series[h])) == target
    ok = [holds(h) for h in range(Hs + 1)]
    assert ok[Hs]
    h0 = Hs
    while h0 > 0 and ok[h0 - 1]:
        h0 -= 1
    return h0


def check_cm_set(ck, A0, H):
    sd = SetData(A0)
    d = sd.d
    ud, vd = uv(A0, d)
    HB = max(H, -((-(sd.c1 + sd.c2)) // d) + 2, max(ud) + max(vd) + 1)
    Tt = Totals(A0, HB)
    for m in moduli(A0):
        u, v = uv(A0, m)
        U, V = max(u), max(v)
        ck.check("u_0 = v_0 = 0, u_s >= 1, v_s >= 0", u[0] == 0 and v[0] == 0 and min(u[1:] or [1]) >= 1 and min(v) >= 0, (A0, m))
        beta, gamma, lam = spectrum(u, v)
        ck.equal("lem:spectrum-mass sum beta = M(S_1)", sum(beta), sd.M1, (A0, m))
        ck.equal("lem:spectrum-mass sum gamma = M(S_2)", sum(gamma), sd.M2, (A0, m))
        so = similarly_ordered(u, v)
        inv = inversion_intervals(u, v)
        for h in range(1, H + 1):
            N = h * d // m
            elems = set(engine_a.elements(Tt.series[h]))
            Ibits = Tt.I_bits(h)
            for n in range(0, h * d + 1):
                k, s = -((-n) // m), (-n) % m
                inside = u[s] <= k <= N - v[s]
                if not ck.equal("lem:residue n in I_h iff u_s <= k <= N - v_s", inside, bool((Ibits >> n) & 1), (A0, m, h, n)):
                    break
                ck.equal("lem:residue for CM sets also h A_0", inside, n in elems, (A0, m, h, n))
            Dh = Tt.D(h)
            tot = 0
            for k in range(1, N + 1):
                Hk = holes(u, v, N, k)
                tot += run_functional(Hk)
                for s in range(1, m):
                    ck.equal("prop:counting (a)", counting_undershoot(u, v, m, N, k, s), Dh[k * m - s], (A0, m, h, k, s))
            Tot = sum(Dh)
            ck.equal("prop:counting (b) Tot(h) = sum f(H_k)", tot, Tot, (A0, m, h))
            for (L, Uw) in {(0, h * d), (min(d - 1, h * d), h * d), (h * d // 3, max(h * d // 3, h * d // 2))}:
                triples = 0
                for k in range(1, N + 1):
                    Hk = set(holes(u, v, N, k))
                    for s in range(1, m):
                        if L <= k * m - s <= Uw:
                            t = s
                            while t <= m - 1 and t in Hk:
                                triples += 1
                                t += 1
                ck.equal("prop:counting (c) window sums", triples, sum(Dh[L:Uw + 1]), (A0, m, h, L, Uw))
            bound = sum(lam[:N])
            # lem:merge
            merged = sum(max(beta[k - 1] if k - 1 < len(beta) else 0,
                             gamma[N - k] if N - k < len(gamma) else 0) for k in range(1, N + 1))
            ck.equal("lem:merge", merged, bound, (A0, m, h))
            ck.check("thm:spectrum (1) Tot(h) >= sum of N largest", Tot >= bound, (A0, m, h))
            strict = any(lo <= N <= hi for (lo, hi) in inv)
            ck.equal("prop:strict-levels", Tot > bound, strict, (A0, m, h))
            if so:
                ck.equal("thm:spectrum (2) equality under similar ordering", Tot, bound, (A0, m, h))
                ck.equal("thm:spectrum (3) Tot = mass iff N >= max(U-1,0) + V",
                         Tot == sd.M1 + sd.M2, N >= max(U - 1, 0) + V, (A0, m, h))
        if so:
            tots = [0] + [Tt.Tot(h) for h in range(1, H + 1)]
            ck.check("thm:spectrum (3) nondecreasing, concave, at most the mass",
                     all(tots[i] <= tots[i + 1] for i in range(H)) and
                     all(tots[i + 1] - tots[i] <= tots[i] - tots[i - 1] for i in range(1, H)) and
                     tots[-1] <= sd.M1 + sd.M2, (A0, m))
        if so and m >= 2:
            s1 = min(s for s in range(m) if u[s] == U)
            s2 = max(s for s in range(m) if v[s] == V)
            ck.check("prop:threshold c_1 + c_2 = (U + V - 2) m + s_2 - s_1 + 2 with 2 <= s_2 - s_1 + 2 <= m",
                     sd.c1 + sd.c2 == (U + V - 2) * m + s2 - s1 + 2 and 2 <= s2 - s1 + 2 <= m, (A0, m))
            sig_m = -((-(sd.c1 + sd.c2)) // d)
            for h in range(1, HB + 1):
                ck.equal("prop:threshold at the modulus m: Tot(h) = mass iff h >= sigma",
                         Tt.Tot(h) == sd.M1 + sd.M2, h >= sig_m, (A0, m, h))
        if m >= 2 and any(u[s] == U and v[s] == V for s in range(m)):
            ud, vd = uv(A0, d)
            reg_d = max(ud[s] + vd[s] for s in range(d))
            ceil_m = -((-(sd.c1 + sd.c2)) // m)
            ceil_d = -((-(sd.c1 + sd.c2)) // d)
            ck.equal("prop:threshold at the modulus m: an index attaining U and V gives U + V - 1 = ceil((c_1+c_2)/m) and reg = sigma + 1",
                     (U + V - 1, reg_d), (ceil_m, ceil_d + 1), (A0, m))
        if m == d:
            ck.equal("thm:spectrum (2) at m = d: equality at every level iff similarly ordered",
                     all(Tt.Tot(h) == sum(lam[:h]) for h in range(1, HB + 1)), so, A0)
            reg = max(u[s] + v[s] for s in range(d))
            ceil_c = -((-(sd.c1 + sd.c2)) // d)
            sig = sigma_by_definition(sd, Tt, HB)
            ck.equal("sigma(A) = ceil((c_1+c_2)/d) for CM sets", sig, ceil_c, A0)
            if d >= 2 and any(u[s] == U and v[s] == V for s in range(d)):
                ck.equal("remark after prop:threshold: an index attaining U and V gives U + V - 1 = sigma and reg = sigma + 1",
                         (U + V - 1, reg), (sig, sig + 1), A0)
            if so and d >= 2:
                ck.equal("eq:reg-cm max(u+v) = U + V under similar ordering", reg, U + V, A0)
                ck.equal("prop:threshold U + V - 1 = ceil((c_1+c_2)/d) = sigma(A)", (U + V - 1, ceil_c), (sig, sig), A0)
                ck.equal("prop:threshold reg = sigma + 1", reg, sig + 1, A0)
                for h in range(1, HB + 1):
                    ck.equal("prop:threshold Tot(h) = mass iff h >= sigma",
                             Tt.Tot(h) == sd.M1 + sd.M2, h >= sig, (A0, h))


def check_diagram(ck, A0, cm):
    _, a1, a2, d = A0
    dg = diagram(A0)
    integral = all(x.denominator == 1 and y.denominator == 1 and x >= 0 and y >= 0 for (_, _, x, y) in dg)
    ck.equal("lem:diagram CM iff all x_r, y_r nonnegative integers", integral, cm, A0)
    if not cm:
        return None
    pts = {(int(x), int(y)) for (_, _, x, y) in dg}
    ck.equal("lem:diagram d distinct points", len(pts), d, A0)
    ck.check("lem:diagram order ideal", all((x - 1, y) in pts for (x, y) in pts if x > 0) and
             all((x, y - 1) in pts for (x, y) in pts if y > 0), A0)
    u, v = uv(A0, d)
    for (r, t, x, y) in dg:
        s = (-r) % d
        val = a1 * x + a2 * y
        ck.equal("eq:uv-diagram u_s = ceil((a_1 x + a_2 y)/d), v_s = x + y - u_s (including r = 0)",
                 (u[s], v[s]), (-((-val) // d), x + y - (-((-val) // d))), (A0, r))
    so = similarly_ordered(u, v)
    layers = {}
    for (x, y) in pts:
        layers.setdefault(x + y, []).append((x, y))
    cond2 = all(len({-((-(a1 * x + a2 * y)) // d) for (x, y) in layers[t]}) == 1 for t in layers)
    cond3 = all(-((-min(a1 * x + a2 * y for (x, y) in layers[t])) // d) ==
                -((-max(a1 * x + a2 * y for (x, y) in layers[t])) // d) for t in layers)
    ck.equal("thm:similar (1) iff (2)", cond2, so, A0)
    ck.equal("thm:similar (2) iff (3)", cond3, so, A0)
    # prop:shape
    L = max(x for (x, _) in pts) + 1
    Hh = max(y for (_, y) in pts) + 1
    if len(pts) == L * Hh:
        ck.check("prop:shape rectangle", True)
    else:
        W = min(x for x in range(L + 1) if (x, Hh - 1) not in pts)
        Zz = max(y for (x, y) in pts if x == L - 1) + 1
        shape = {(x, y) for x in range(W) for y in range(Hh)} | {(x, y) for x in range(W, L) for y in range(Zz)}
        ck.equal("prop:shape L-shape", pts, shape, A0)
        ck.equal("prop:shape a_1 W + a_2 Z = 0 mod d", (a1 * W + a2 * Zz) % d, 0, A0)
        if L - W >= 2 and Hh - Zz >= 2:
            ck.check("prop:shape both notch dimensions >= 2 forces not similarly ordered", not so, A0)
    return so, pts


def main():
    args = standard_parser(__doc__).parse_args()
    if args.quick:
        D4, D5, H, D4CENSUS, D5CENSUS, DFAM = 9, 7, 5, 14, 9, 14
    else:
        D4, D5, H, D4CENSUS, D5CENSUS, DFAM = 20, 12, 8, 30, 13, 40
    ck = Checker("sumsets", "Residue intervals, exact counting, spectrum, strict levels, threshold, similar ordering",
                 {"four_element_d_max": D4, "five_element_d_max": D5, "levels_h_max": H,
                  "census_four_element_d_max": D4CENSUS, "census_five_element_d_max": D5CENSUS,
                  "families_d_max": DFAM})

    for k, dmax in ((3, D4), (4, D4), (5, D5)):
        for A0 in sets_with(k, dmax):
            cm = is_cm_three_ways(ck, A0, H)
            if cm:
                check_cm_set(ck, A0, H)
            else:
                # counting holds for I_h itself at every level
                sd = SetData(A0)
                d = sd.d
                Tt = Totals(A0, H)
                u, v = uv(A0, d)
                for h in range(1, H + 1):
                    DI = Tt.D_I(h)
                    for kk in range(1, h + 1):
                        for s in range(1, d):
                            ck.equal("prop:counting (a) for I_h at every level (non-CM sets)",
                                     counting_undershoot(u, v, d, h, kk, s), DI[kk * d - s], (A0, h, kk, s))
            if k == 4:
                check_diagram(ck, A0, cm)

    # Census of four- and five-element sets (statistics recorded for the computation section, informational).
    census = {}
    for k, dmax in ((4, D4CENSUS), (5, D5CENSUS)):
        for A0 in sets_with(k, dmax):
            d = A0[-1]
            if not engine_b.is_cm(A0):
                continue
            u, v = uv(A0, d)
            so = similarly_ordered(u, v)
            for bound in sorted({dmax, 22, 13, 16, 30}):
                if d <= bound <= dmax:
                    key = "%d-element CM sets with d <= %d" % (k, bound)
                    c = census.setdefault(key, [0, 0])
                    c[0] += 1
                    c[1] += 1 if so else 0
    ck.record("census [CM count, similarly ordered at m = d]", census)

    # ex:spectrum (a)
    A0 = [0, 2, 5, 10]
    u, v = uv(A0, 10)
    _, _, lam = spectrum(u, v)
    Tt = Totals(A0, 4)
    ck.stated("ex:spectrum (a) {0,2,5,10} m = 10 positive spectrum", lam, [14, 6, 2, 2])
    ck.stated("ex:spectrum (a) bound at h = 3", sum(lam[:3]), 22)
    ck.stated("ex:spectrum (a) Tot(3)", Tt.Tot(3), 23)
    ck.record("ex:spectrum (a) u at m = 10", u)
    ck.record("ex:spectrum (a) v at m = 10", v)
    u5, v5 = uv(A0, 5)
    ck.stated("ex:spectrum (a) similarly ordered at m = 5", similarly_ordered(u5, v5), True)
    ck.stated("ex:spectrum (a) similarly ordered at m = 10", similarly_ordered(u, v), False)
    dg = {(int(x), int(y)) for (_, _, x, y) in diagram(A0)}
    ck.stated("{0,2,5,10} diagram is [0,4] x [0,1]", sorted(dg), sorted((x, y) for x in range(5) for y in range(2)))
    ck.stated("{0,2,5,10} layer 4 values", sorted(2 * x + 5 * y for (x, y) in dg if x + y == 4), [8, 11])
    # ex:spectrum (b), the example {0,2,8,11}
    A0 = [0, 2, 8, 11]
    sd = SetData(A0)
    ck.stated("{0,2,8,11} is CM (engine B) and h A_0 = I_h for h <= 15",
              engine_b.is_cm(A0) and all(engine_a.sumset_bits(A0, h) == engine_c.sumset_bits_intersection(A0, h)
                                         for h in range(1, 16)), True)
    ck.stated("{0,2,8,11} M(S_1) + M(S_2)", sd.M1 + sd.M2, 18)
    Tt = Totals(A0, 12)
    ck.stated("{0,2,8,11} Tot(1..12)", [Tt.Tot(h) for h in range(1, 13)], [19, 19] + [18] * 10)
    u, v = uv(A0, 11)
    ck.stated("{0,2,8,11} (u_6, v_6) and (u_1, v_1)", [(u[6], v[6]), (u[1], v[1])], [(2, 0), (1, 1)])
    ck.stated("{0,2,8,11} interval of the inversion (s,t) = (6,1)", (u[1] + v[6], u[6] + v[1] - 2), (1, 1))
    ck.record("{0,2,8,11} all inversion intervals (strict levels N = 1, 2)", inversion_intervals(u, v))
    invs = [(s_, t_) for s_ in range(11) for t_ in range(11) if u[s_] > u[t_] and v[s_] < v[t_]]
    ck.stated("{0,2,8,11} number of inversions", len(invs), 5)
    ck.stated("{0,2,8,11} union of the inversion intervals",
              sorted({N for (s_, t_) in invs for N in range(u[t_] + v[s_], u[s_] + v[t_] - 1)}), [1, 2])
    ck.stated("{0,2,8,11} moduli", moduli(A0), [11])
    ck.stated("{0,2,8,11} positive spectrum at m = 11", spectrum(u, v)[2], [10, 5, 3])
    # opening of Section 3: {0,2,3,5}
    A0 = [0, 2, 3, 5]
    u, v = uv(A0, 5)
    ck.stated("{0,2,3,5} is CM, its only modulus is 5", (engine_b.is_cm(A0), moduli(A0)), (True, [5]))
    ck.stated("{0,2,3,5} Tot(h) = 2 for h = 1..8", [Totals(A0, 8).Tot(h) for h in range(1, 9)], [2] * 8)
    ck.stated("{0,2,3,5} pairs (1,1) and (2,0) occur", ((1, 1) in zip(u, v), (2, 0) in zip(u, v)), (True, True))
    # opening of the introduction: {0,5,7,14} rises to its mass, {0,4,9,13} climbs past it
    A0 = [0, 5, 7, 14]
    sd = SetData(A0)
    u, v = uv(A0, 7)
    ck.stated("{0,5,7,14} is CM and similarly ordered at m = 7",
              (engine_b.is_cm(A0), similarly_ordered(u, v)), (True, True))
    ck.stated("{0,5,7,14} M(S_1) + M(S_2)", sd.M1 + sd.M2, 65)
    u14, v14 = uv(A0, 14)
    ck.stated("{0,5,7,14} similarly ordered at m = 14, positive spectrum",
              (similarly_ordered(u14, v14), spectrum(u14, v14)[2]), (True, [32, 16, 8, 5, 3, 1]))
    def circular_phi(b_, a_, j):
        pts = sorted({(i * a_) % b_ for i in range(j)})
        gaps_ = [(pts[(i + 1) % len(pts)] - pts[i]) % b_ or b_ for i in range(len(pts))]
        return sum(comb(g, 2) for g in gaps_)
    ck.stated("{0,5,7,14} phi_{7,5}(1..7) from circular gaps", [circular_phi(7, 5, j) for j in range(1, 8)],
              [21, 11, 5, 3, 2, 1, 0])
    ck.stated("{0,5,7,14} Tot(1..10)", [Totals(A0, 10).Tot(h) for h in range(1, 11)],
              [32, 48, 56, 61, 64] + [65] * 5)
    A0 = [0, 4, 9, 13]
    sd = SetData(A0)
    Tt = Totals(A0, 12)
    failing = [h for h in range(1, 13) if engine_a.sumset_bits(A0, h) != engine_c.sumset_bits_intersection(A0, h)]
    ck.stated("{0,4,9,13} is not CM; failing levels for h <= 12", (engine_b.is_cm(A0), failing), (False, [4, 5, 6, 7]))
    ck.stated("{0,4,9,13} M(S_1) + M(S_2)", sd.M1 + sd.M2, 40)
    ck.stated("{0,4,9,13} Tot(1..12)", [Tt.Tot(h) for h in range(1, 13)],
              [22, 36, 45, 48, 49, 46, 44] + [40] * 5)
    ck.record("{0,4,9,13} conductors (c_1, c_2)", (sd.c1, sd.c2))
    # fig:diagrams: the Apery diagrams of {0,5,7,14} and {0,2,8,11} and the layers that cross a multiple of d
    for A0, pts_expected, bad_layers in (
            ([0, 5, 7, 14], sorted((x, y) for x in range(7) for y in range(2)), []),
            ([0, 2, 8, 11], sorted([(x, y) for x in range(3) for y in range(3)] + [(3, 0), (3, 1)]), [2, 3])):
        _, a1_, a2_, d_ = A0
        pts = sorted((int(x), int(y)) for (_, _, x, y) in diagram(A0))
        crossing = sorted({x + y for (x, y) in pts
                           if len({-((-(a1_ * x2 + a2_ * y2)) // d_) for (x2, y2) in pts if x2 + y2 == x + y}) > 1})
        ck.stated("fig:diagrams %s diagram points" % A0, pts, pts_expected)
        ck.stated("fig:diagrams %s layers whose values cross a multiple of d" % A0, crossing, bad_layers)
    # remarks after prop:threshold
    u, v = uv([0, 1, 2, 5], 5)
    ck.stated("{0,1,2,5} pairs at s = 1, 2", [(u[1], v[1]), (u[2], v[2])], [(1, 1), (1, 1)])
    ck.stated("{0,1,2,5} Apery points (u_s d - s, v_s d + s) for s = 1, 2",
              [(u[s] * 5 - s, v[s] * 5 + s) for s in (1, 2)], [(4, 6), (3, 7)])
    A0 = [0, 3, 4, 6]
    sd = SetData(A0)
    u, v = uv(A0, 6)
    ck.stated("{0,3,4,6} is CM", engine_b.is_cm(A0), True)
    ck.stated("{0,3,4,6} similarly ordered", similarly_ordered(u, v), True)
    ck.stated("{0,3,4,6} (c_1, c_2)", (sd.c1, sd.c2), (6, 2))
    ck.record("{0,3,4,6} u", u)
    ck.record("{0,3,4,6} v", v)
    low = {n for n in range(0, sd.c1 - 1) if sd.in_S1(n)}
    top = {n for n in range(0, sd.c2 - 1) if sd.in_S2(n)}
    dec1 = low | set(range(sd.c1, 6 - sd.c2 + 1)) | {6 - x for x in top}
    ck.stated("{0,3,4,6} decomposition with empty central interval holds at h = 1",
              (dec1 == set(engine_a.elements(engine_a.sumset_bits(A0, 1))), sd.c1 > 6 - sd.c2), (True, True))
    ck.stated("{0,3,4,6} ceil((c_1+c_2)/d)", -((-(sd.c1 + sd.c2)) // 6), 2)

    # prop:families
    for d in range(3, DFAM + 1):
        for b in range(2, d):
            A0 = [0, 1, b, d]
            q_, r_ = divmod(d, b)
            cm = engine_b.is_cm(A0)
            ck.equal("prop:families (a) {0,1,b,d} CM iff r = 0 or q + r >= b", cm, r_ == 0 or q_ + r_ >= b, A0)
            if cm:
                ck.check("prop:families (a) similarly ordered", similarly_ordered(*uv(A0, d)), A0)
                ref = sorted(d - x for x in A0)
                ck.equal("prop:families (b) reflection preserves CM and similar ordering",
                         (engine_b.is_cm(ref), similarly_ordered(*uv(ref, d))), (True, True), A0)
        for a in range(1, d):
            if 0 < a < d - a and gcd(a, d) == 1:
                A0 = [0, a, d - a, d]
                if engine_b.is_cm(A0):
                    ck.equal("prop:families (c) {0,a,d-a,d} similarly ordered iff d = 3",
                             similarly_ordered(*uv(A0, d)), d == 3, A0)
            if a + 1 < d:
                A0 = [0, a, a + 1, d]
                ck.check("prop:families (d) {0,a,a+1,d} is CM", engine_b.is_cm(A0), A0)
                if d <= 30:
                    res = check_diagram(ck, A0, True)
                    so, pts = res
                    L = max(x for (x, _) in pts) + 1
                    Hh = max(y for (_, y) in pts) + 1
                    if len(pts) == L * Hh:
                        expect = True
                    else:
                        W = min(x for x in range(L + 1) if (x, Hh - 1) not in pts)
                        Zz = max(y for (x, y) in pts if x == L - 1) + 1
                        expect = (L - W == 1) or (Hh - Zz == 1)
                    ck.equal("prop:families (d) similarly ordered iff rectangle or a notch of width 1", so, expect, A0)
        for A0 in sets_with(4, d, d):
            if not engine_b.is_cm(A0):
                continue
            ref = sorted(d - x for x in A0)
            so = similarly_ordered(*uv(A0, d))
            ck.equal("prop:families (b) A_0 CM and similarly ordered iff d - A_0 is", (engine_b.is_cm(ref),
                     similarly_ordered(*uv(ref, d))), (True, so), A0)
            if A0[2] == d - 1:
                ck.check("prop:families (b) every CM {0,a,d-1,d} is similarly ordered", so, A0)
    # eq:lambda-test (5): a four-element set is CM iff len(w + kd) - len(w) >= k (checked for k <= 2d)
    for d in range(3, D4 + 1):
        for A0 in sets_with(4, d, d):
            _, a1, a2, _ = A0
            top = (a1 * a2 + 2 * d) * 3 + 3 * d * d
            ln = min_parts([0, a1, a2], top)
            ml = min_parts(A0, top)
            ap1 = {}
            for n in range(top + 1):
                if ml[n] < INF and n % d not in ap1:
                    ap1[n % d] = n
            test = all(ln[w + k * d] - ln[w] >= k for w in ap1.values() for k in range(0, 2 * d + 1)
                       if w + k * d <= top and ln[w + k * d] < INF)
            ck.equal("eq:lambda-test CM iff len(w + kd) - len(w) >= k for k <= 2d", test, engine_b.is_cm(A0), A0)
    u, v = uv([0, 1, 2, 3], 3)
    ck.stated("{0,1,2,3} pairs", list(zip(u, v)), [(0, 0), (1, 0), (1, 0)])
    ck.finish()


if __name__ == "__main__":
    main()

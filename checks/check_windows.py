#!/usr/bin/env python3
"""Corrections and windows: undershoot as overshoot, the correction for failing levels and the
hinge formulas; and the reduction by the gcd (fact (ii) of the preliminaries).

Paper results: prop:overshoot-form, prop:correction, lem:hinge, thm:hinge; fact (ii) of the
preliminaries.
"""
import os
import sys
from itertools import combinations
from math import gcd
from functools import reduce
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))

from sumset_undershoot.report import Checker, standard_parser
from sumset_undershoot.engines import direct as engine_a
from sumset_undershoot.semigroup import apery, conductor
from sumset_undershoot.sumsets import (SetData, normalize, psi, o2, u1, all_corners, hinge_window,
                                    hinge_sum, hinge_count, minimal_points, D_E, min_parts)



def sets_with(k, dmax, dmin=2):
    for d in range(dmin, dmax + 1):
        for mid in combinations(range(1, d), k - 2):
            A0 = [0] + list(mid) + [d]
            if reduce(gcd, A0) == 1:
                yield A0


def physical_undershoot(series, A0, g, h, T):
    n = T // g
    bits = series[h]
    below = bits & ((1 << (n + 1)) - 1)
    return (T % g) + g * (n - (below.bit_length() - 1))



def main():
    args = standard_parser(__doc__).parse_args()
    if args.quick:
        DC, DH, HH = 7, 7, 3
    else:
        DC, DH, HH = 11, 10, 4
    ck = Checker("windows", "Overshoot form, corrections, hinge formulas, gcd reduction",
                 {"overshoot_form_and_correction_d_max": DC, "hinge_d_max": DH, "hinge_h_max": HH})

    # prop:overshoot-form, prop:correction
    for k in (3, 4, 5):
        for A0 in sets_with(k, DC if k < 5 else min(DC, 8)):
            sd = SetData(A0)
            d, c1, c2 = sd.d, sd.c1, sd.c2
            for h in range(0, 6):
                hd = h * d
                I = sd.I_set(h) if h > 0 else [0]
                A = sd.sumset(h)
                DI = [D_E(I, n) for n in range(hd + 1)]
                DA = [D_E(A, n) for n in range(hd + 1)]
                for n in range(0, hd - c2 + 1):
                    ck.equal("prop:overshoot-form (a) D_I(n) = D_{S_1}(n) below hd - c_2", DI[n], u1(sd, n), (A0, h, n))
                Iset = set(I)
                for x in [y for y in I if y >= c1]:
                    for n in range(x, hd + 1):
                        ck.equal("prop:overshoot-form (b) D_I(n) = o_2(hd - n) above every anchor x >= c_1 in I_h",
                                 DI[n], o2(sd, hd - n), (A0, h, x, n))
                if hd >= c1 and (c1 == 0 or hd >= c1 + c2):
                    ck.check("prop:overshoot-form (b) c_1 is an anchor when c_1 = 0 or hd >= c_1 + c_2", c1 in Iset, (A0, h))
                if hd >= c1 > 0:
                    Nh = max(x for x in I if x <= c1 - 1)
                    for n in range(c1, hd + 1):
                        j = o2(sd, hd - n)
                        ck.equal("prop:overshoot-form (c)", DI[n], j if j <= n - c1 else n - Nh, (A0, h, n))
                if hd >= c1 + c2 - 1 and h > 0:
                    ck.check("prop:overshoot-form (d) (a) and (b)/(c) cover [0, hd]", hd - c2 >= c1 - 1, (A0, h))
                if h > 0:
                    windows = [list(range(0, hd + 1)), list(range(c1 if c1 <= hd else 0, hd + 1)),
                               list(range(hd // 3, hd // 2 + 1)), list(range(0, hd + 1, 2))]
                    for W in windows:
                        ck.equal("prop:correction sum_W D_h = sum_W D_I + Psi_h(W)", sum(DA[n] for n in W),
                                 sum(DI[n] for n in W) + psi(sd, h, W), (A0, h, len(W)))
                    full = psi(sd, h, list(range(hd + 1)))
                    ck.equal("prop:correction Psi_h([0,hd]) > 0 iff failing level", full > 0, set(I) != set(A), (A0, h))
    # fact (ii) of the preliminaries: physical undershoot (T mod g) + g D_h(floor(T/g))
    for menu in ([15, 30], [6, 15, 30], [10, 15, 20, 30], [4, 6, 10], [9, 21, 33]):
        g, A0 = normalize(menu)
        for h in range(0, 6):
            sums = {0}
            for _ in range(h):
                sums = sums | {x + y for x in sums for y in menu}
            series = engine_a.sumset_series(A0, h)
            for T in range(0, h * max(menu) + 2 * g + 1):
                direct = T - max(x for x in sums if x <= T)
                ck.equal("fact (ii) physical undershoot = (T mod g) + g D_h(floor(T/g))", direct,
                         physical_undershoot(series, A0, g, h, T), (menu, h, T))
    # lem:hinge directly: H_C counts the integers of [l0, q0] in the union of [kappa_j, h - mu_j]
    for k in (3, 4, 5):
        for A0 in sets_with(k, DH, 3):
            sd = SetData(A0)
            corners = all_corners(sd)
            d = sd.d
            for s in range(1, d):
                union = []
                for t in range(s, d):
                    union.extend(corners[t])
                    C = minimal_points(union)
                    for h in range(1, HH + 1):
                        for l0 in range(0, h + 2):
                            for q0 in range(l0 - 1, h + 2):
                                direct = sum(1 for kk in range(l0, q0 + 1)
                                             if any(kap <= kk <= h - mu for (kap, mu) in C))
                                ck.equal("lem:hinge H_C equals the direct count of the union", hinge_count(C, h, l0, q0),
                                         direct, (A0, s, t, h, l0, q0))
    # thm:hinge
    for k in (3, 4, 5):
        for A0 in sets_with(k, DH, 3):
            sd = SetData(A0)
            d = sd.d
            corners = all_corners(sd)
            for h in range(1, HH + 1):
                hd = h * d
                A = sd.sumset(h)
                DA = [D_E(A, n) for n in range(hd + 1)]
                for (L, U) in ((0, hd), (1, hd - 1), (d - 1, hd), (hd // 3, (2 * hd) // 3), (hd // 2, hd // 2)):
                    if not 0 <= L <= U <= hd:
                        continue
                    tot, ps = hinge_window(sd, corners, h, L, U)
                    ck.equal("thm:hinge window total", tot, sum(DA[L:U + 1]), (A0, h, L, U))
                    ck.equal("thm:hinge Psi_h([L,U])", ps, psi(sd, h, list(range(L, U + 1))), (A0, h, L, U))
                    for c in range(1, d):
                        if c in A0:
                            continue
                        A1 = sorted(A0 + [c])
                        sd1 = SetData(A1)
                        c1n = all_corners(sd1)
                        Dp = [D_E(sd1.sumset(h), n) for n in range(hd + 1)]
                        Dh1 = [D_E(sd.sumset(h + 1), n) for n in range(hd + 1)]
                        lhs = hinge_sum(sd1, c1n, h, L, U, h)
                        rhs = hinge_sum(sd, corners, h, L, U, h + 1)
                        a_new = sum(Dp[L:U + 1])
                        a_cap = sum(Dh1[L:U + 1])
                        ck.equal("thm:hinge comparison: adding c beats raising the cap iff sums compare",
                                 (a_new < a_cap, a_new == a_cap), (lhs > rhs, lhs == rhs), (A0, c, h, L, U))
    ck.finish()


if __name__ == "__main__":
    main()

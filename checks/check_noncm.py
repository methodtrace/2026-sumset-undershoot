#!/usr/bin/env python3
"""Non-Cohen-Macaulay sets: failing levels, the overshoot of the mass, the four-element
consequence of the cited regularity inequality, and the five-element example.

Paper results: prop:failing-interval (and the corner formula for tau), prop:overshoot,
the cited inequality (tau - 1) d >= c_1 + c_2 + d (finite evidence only; the proof is cited),
cor:four-overshoot, ex:five-element.
"""
import os
import sys
from itertools import combinations
from math import gcd
from functools import reduce

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))

from sumset_undershoot.report import Checker, standard_parser
from sumset_undershoot.engines import direct as engine_a, homogenized as engine_b, intersection as engine_c
from sumset_undershoot.sumsets import SetData, uv, all_corners, failing_levels_by_class


def sets_with(k, dmax):
    for d in range(2, dmax + 1):
        for mid in combinations(range(1, d), k - 2):
            A0 = [0] + list(mid) + [d]
            if reduce(gcd, A0) == 1:
                yield A0


def analyse(ck, A0, stats):
    sd = SetData(A0)
    d = sd.d
    corners = all_corners(sd)
    u, v = uv(A0, d)
    for s in range(1, d):
        C = corners[s]
        ck.equal("prop:failing-interval kappa_1 = u_s and mu_t = v_s", (C[0][0], C[-1][1]), (u[s], v[s]), (A0, s))
    intervals = failing_levels_by_class(sd, corners)
    tau_formula = 0
    for s in intervals:
        if intervals[s] is not None:
            tau_formula = max(tau_formula, intervals[s][1] + 1)
    Hmax = max(tau_formula + 2, -((-(sd.c1 + sd.c2)) // d) + 2, 3)
    series = engine_a.sumset_series(A0, Hmax)
    gb = engine_c.gaps(sd.g1)
    gt = engine_c.gaps(sd.g2)
    failing_direct = {s: [] for s in range(1, d)}
    any_fail = []
    mass = sd.M1 + sd.M2
    for h in range(0, Hmax + 1):
        A = series[h]
        I = engine_c.sumset_bits_intersection(A0, h, gb, gt) if h > 0 else 1
        ck.check("h A_0 is contained in I_h", A & ~I == 0, (A0, h))
        Fh = I & ~A
        if Fh:
            any_fail.append(h)
            for n in engine_a.elements(Fh):
                s = (-n) % d
                if not failing_direct[s] or failing_direct[s][-1] != h:
                    failing_direct[s].append(h)
        if h >= 1:
            tot = engine_a.total_undershoot(A, 0, h * d)
            totI = engine_a.total_undershoot(I, 0, h * d)
            if h * d >= sd.c1 + sd.c2:
                ck.equal("prop:overshoot total of I_h = M(S_1) + M(S_2) when hd >= c_1 + c_2", totI, mass, (A0, h))
                ck.check("prop:overshoot Tot(h) >= mass", tot >= mass, (A0, h))
                ck.equal("prop:overshoot equality iff h is not a failing level", tot == mass, Fh == 0, (A0, h))
            ck.equal("fact (iii): Tot of h A_0 exceeds that of I_h iff F_h nonempty", tot > totI, Fh != 0, (A0, h))
    for s in range(1, d):
        got = failing_direct[s]
        iv = intervals[s]
        expect = list(range(iv[0], iv[1] + 1)) if iv is not None else []
        ck.equal("prop:failing-interval levels of class s = corner interval", got, expect, (A0, s))
    tau_direct = (max(any_fail) + 1) if any_fail else 0
    ck.equal("tau - 1 = max over classes of (kappa_{j+1} + mu_j - 2)", tau_direct, tau_formula, A0)
    cm = engine_b.is_cm(A0)
    ck.equal("CM iff no failing level", cm, not any_fail, A0)
    if len(A0) == 4:
        stats["four"] += 1
        if not cm:
            stats["four_noncm"] += 1
            stats["four_noncm_failing_levels"] += len(any_fail)
            ck.check("cited inequality (tau - 1) d >= c_1 + c_2 + d (finite evidence)",
                     (tau_direct - 1) * d >= sd.c1 + sd.c2 + d, A0)
            tot = engine_a.total_undershoot(series[tau_direct - 1], 0, (tau_direct - 1) * d)
            ck.check("cor:four-overshoot Tot(tau - 1) > M(S_1) + M(S_2)", tot > mass, A0)
    return sd, tau_direct, series


def main():
    args = standard_parser(__doc__).parse_args()
    D4, D5 = (12, 8) if args.quick else (23, 14)
    ck = Checker("noncm", "Failing levels, overshoot of the mass, four-element consequence, five-element example",
                 {"four_element_d_max": D4, "five_element_d_max": D5})
    stats = {"four": 0, "four_noncm": 0, "four_noncm_failing_levels": 0}
    for A0 in sets_with(4, D4):
        analyse(ck, A0, stats)
    for A0 in sets_with(5, D5):
        analyse(ck, A0, stats)
    ck.record("four-element sets checked", stats["four"])
    ck.record("non-CM four-element sets", stats["four_noncm"])
    ck.record("failing levels over non-CM four-element sets", stats["four_noncm_failing_levels"])
    # ex:five-element
    A0 = [0, 5, 9, 11, 20]
    sd, tau, series = analyse(ck, A0, {"four": 0, "four_noncm": 0, "four_noncm_failing_levels": 0})
    ck.stated("ex:five-element (c_1, c_2)", (sd.c1, sd.c2), (18, 44))
    ck.stated("ex:five-element M(S_1) + M(S_2)", sd.M1 + sd.M2, 74)
    fails = [h for h in range(0, 9) if engine_c.sumset_bits_intersection(A0, h) != engine_a.sumset_bits(A0, h)]
    ck.stated("ex:five-element failing levels", fails, [3])
    ck.stated("ex:five-element tau", tau, 4)
    tots = [engine_a.total_undershoot(engine_a.sumset_bits(A0, h), 0, 20 * h) for h in range(1, 9)]
    ck.stated("ex:five-element Tot(1..3)", tots[:3], [53, 69, 74])
    ck.stated("ex:five-element Tot(h) = 74 for h = 4..8", tots[3:], [74] * 5)
    ck.stated("ex:five-element 3d < c_1 + c_2 (60 < 62)", (60, sd.c1 + sd.c2), (60, 62))
    ck.stated("ex:five-element is not CM (engine B)", engine_b.is_cm(A0), False)
    ck.finish()


if __name__ == "__main__":
    main()

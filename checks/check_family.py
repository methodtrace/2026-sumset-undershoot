#!/usr/bin/env python3
"""The family A_0 = {0, a, b, 2b}: structure, spectrum, totals, the three-distance
description and the quadratic pieces of Tot (including the example a = 2, b = 9).

Paper results: lem:family-ci, lem:family-residue, thm:family-runs, thm:family-spectrum
(eq:family-tot), cor:family-mass, thm:three-distance, cor:quadratic-pieces.
"""
import os
import sys
from math import comb, gcd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))

from sumset_undershoot.report import Checker, standard_parser
from sumset_undershoot.engines import direct as engine_a, homogenized as engine_b
from sumset_undershoot.semigroup import conductor, two_generator_mass
from sumset_undershoot.sumsets import uv, similarly_ordered, spectrum, normalize
from sumset_undershoot.twogen import spacings, Phi_insertion, phi_definition, Phi_definition
from sumset_undershoot.family import phi_family, M4, Tot_formula, convergent_data, cf_index, three_distance_gaps


def runs_outside(bits, N):
    out = {}
    run = 0
    for n in range(N + 1):
        if (bits >> n) & 1:
            if run:
                out[run] = out.get(run, 0) + 1
            run = 0
        else:
            run += 1
    if run:
        out[run] = out.get(run, 0) + 1
    return out


def level_multiset(b, a, j):
    out = {}
    for g in spacings(b, a, j):
        if g >= 2:
            out[g - 1] = out.get(g - 1, 0) + 1
    return out


def main():
    args = standard_parser(__doc__).parse_args()
    if args.quick:
        BENG, BFORM, BTD = 12, 30, 30
    else:
        BENG, BFORM, BTD = 25, 60, 80
    ck = Checker("family", "The family {0,a,b,2b}: spectrum, totals, three distances, quadratic pieces",
                 {"engine_checked_b_max": BENG, "formula_checked_b_max": BFORM, "three_distance_b_max": BTD})
    for b in range(2, BFORM + 1):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            A0 = [0, a, b, 2 * b]
            d = 2 * b
            ph = phi_family(a, b)
            M4v = M4(a, b)
            ck.equal("cor:family-mass M4 = 2 Phi(b,a) - binom(b,2)", M4v, 2 * Phi_insertion(b, a) - comb(b, 2), (a, b))
            ck.equal("phi and Phi from the definition (circular spacings) equal the library values",
                     ([phi_definition(b, a, j) for j in range(1, b + 1)], Phi_definition(b, a)),
                     ([ph[j] for j in range(1, b + 1)], Phi_insertion(b, a)), (a, b))
            ck.check("phi nonincreasing", all(ph[j] >= ph[j + 1] for j in range(1, b)), (a, b))
            # lem:family-ci (exponent bookkeeping of the two binomials)
            j_ = a % 2
            w_ = (a - j_) // 2
            z_ = (2 * b - a - j_) // 2
            ck.check("lem:family-ci binomials homogeneous and in the ideal",
                     j_ + w_ + z_ == b and a == j_ + 2 * w_ and 2 * b - a == j_ + 2 * z_, (a, b))
            ck.equal("lem:family-ci c_1 = (a-1)(b-1), c_2 = (b-1)(2b-a-1)",
                     (conductor([a, b, 2 * b]), conductor([2 * b, 2 * b - a, b])),
                     ((a - 1) * (b - 1), (b - 1) * (2 * b - a - 1)), (a, b))
            ck.equal("lem:family-ci sigma(A) = ceil((c_1+c_2)/d) = b - 1",
                     -((-((a - 1) * (b - 1) + (b - 1) * (2 * b - a - 1))) // d), b - 1, (a, b))
            # lem:family-residue
            u, v = uv(A0, b)
            for s in range(b):
                i_s = ((-s) * pow(a, -1, b)) % b if b > 1 else 0
                ck.equal("lem:family-residue u_s = ceil(i_s a / b), v_s = 2 i_s - u_s",
                         (u[s], v[s]), (-((-i_s * a) // b), 2 * i_s + ((-i_s * a) // b)), (a, b, s))
            ck.check("lem:family-residue similarly ordered at m = b", similarly_ordered(u, v), (a, b))
            ck.equal("lem:family-residue U = a, V = 2b - a - 2", (max(u), max(v)), (a, max(2 * b - a - 2, 0)), (a, b))
            beta, gamma, lam = spectrum(u, v)
            expect = sorted([comb(b, 2)] + [ph[j] for j in range(2, b) for _ in (0, 1)], reverse=True)
            ck.equal("thm:family-spectrum spectrum at m = b", lam, [x for x in expect if x > 0], (a, b))
            # engine-checked part
            if b <= BENG:
                H = b + 2
                series = engine_a.sumset_series(A0, H)
                ck.check("lem:family-ci CM (engine B)", engine_b.is_cm(A0), (a, b))
                ck.equal("cor:family-mass mass by definition", M4v, two_generator_mass(a, b) + two_generator_mass(b, 2 * b - a), (a, b))
                tots = [0]
                for h in range(1, H + 1):
                    bits = series[h]
                    size = bin(bits).count("1")
                    ck.equal("lem:family-ci |h A_0|", size, (h + 1) ** 2 if h <= b - 1 else 2 * b * h + 1 - (b - 1) ** 2, (a, b, h))
                    D = engine_a.undershoot_scan(bits, 0, h * d)
                    tots.append(sum(D))
                    ck.equal("eq:family-tot / Tot = M4 for h >= b-1", tots[-1], Tot_formula(a, b, h), (a, b, h))
                    runs = runs_outside(bits, h * d)
                    expect = {}
                    for j in range(1, min(h, b - 1) + 1 if h < b - 1 else b):
                        times = 1 if j == 1 else 2
                        lm = level_multiset(b, a, j)
                        for k in lm:
                            expect[k] = expect.get(k, 0) + times * lm[k]
                    if h < b - 1:
                        lm = level_multiset(b, a, h + 1)
                        for k in lm:
                            expect[k] = expect.get(k, 0) + lm[k]
                    ck.equal("thm:family-runs run multiset", runs, expect, (a, b, h))
                    Dmirror = engine_a.undershoot_scan(engine_a.sumset_bits([0, b - a, b, 2 * b], h), 0, h * d) if b - a != a else D
                    ck.equal("thm:family-runs distribution of D_h unchanged under a -> b - a", sorted(D), sorted(Dmirror), (a, b, h))
                ck.check("lem:family-residue Tot nondecreasing and concave",
                         all(tots[i] <= tots[i + 1] for i in range(H)) and
                         all(tots[i + 1] - tots[i] <= tots[i] - tots[i - 1] for i in range(1, H)), (a, b))
                ck.equal("lem:family-residue Tot = M4 exactly for h >= b - 1",
                         [h for h in range(1, H + 1) if tots[h] == M4v], list(range(max(b - 1, 1), H + 1)) if b > 1 else [], (a, b))
    # thm:three-distance and cor:quadratic-pieces
    ngap = 0
    for b in range(2, BTD + 1):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            ph = phi_family(a, b)
            cs, q, eps = convergent_data(a, b)
            K = len(cs)
            for h in range(2, b + 1):
                sp = {}
                for g in spacings(b, a, h):
                    sp[g] = sp.get(g, 0) + 1
                ck.equal("thm:three-distance gap multiset", sp, three_distance_gaps(a, b, h), (a, b, h))
                ngap += 1
                k, m, r = cf_index(a, b, h)
                ck.check("eq:cf-index (k, m) != (K-1, c_K)", (k, m) != (K - 1, cs[K - 1]), (a, b, h))
                val = ((h - q[k]) * comb(eps[k], 2) + (r + 1) * comb(eps[k - 1] - m * eps[k], 2)
                       + (q[k] - r - 1) * comb(eps[k - 1] - (m - 1) * eps[k], 2))
                ck.equal("thm:three-distance phi formula", val, ph[h], (a, b, h))
            if b <= BFORM:
                tot = [Tot_formula(a, b, h) for h in range(0, b)]
                covered = []
                for k in range(0, K):
                    for m in range(1, cs[k] + 1):
                        if (k, m) == (K - 1, cs[K - 1]):
                            continue
                        Hs = m * q[k] + q[k - 1] + 1
                        dl = eps[k] * (eps[k - 1] - m * eps[k])
                        for j in range(0, q[k] + 1):
                            hh = Hs - 2 + j
                            if hh <= b - 1:
                                ck.equal("cor:quadratic-pieces", tot[hh], tot[Hs - 2] + 2 * ph[Hs - 1] * j - dl * j * j, (a, b, k, m, j))
                        covered.append((Hs - 2, Hs + q[k] - 2))
                covered.sort()
                ck.check("cor:quadratic-pieces intervals tile [0, b-1] with shared endpoints",
                         covered[0][0] == 0 and covered[-1][1] == b - 1 and
                         all(covered[i][1] == covered[i + 1][0] for i in range(len(covered) - 1)), (a, b))
    ck.record("three_distance_gap_multisets", ngap)
    # a = 2, b = 9
    ck.stated("a=2, b=9: partial quotients of 2/9", convergent_data(2, 9)[0], [4, 2])
    tots = [Tot_formula(2, 9, h) for h in range(0, 10)]
    ck.stated("a=2, b=9: Tot(4..9)", tots[4:10], [120, 127, 132, 135, 136, 136])
    t6 = tots[6:10]
    ck.stated("a=2, b=9: third difference of Tot(6..9)", t6[3] - 3 * t6[2] + 3 * t6[1] - t6[0], 1)
    ck.stated("a=2, b=9: last block of phi is h = 6..9", [h for h in range(2, 10) if cf_index(2, 9, h)[:2] == cf_index(2, 9, 9)[:2]], [6, 7, 8, 9])
    A0 = [0, 2, 9, 18]
    series = engine_a.sumset_series(A0, 9)
    ck.stated("a=2, b=9: Tot(4..9) by engine A", [engine_a.total_undershoot(series[h], 0, 18 * h) for h in range(4, 10)],
              [120, 127, 132, 135, 136, 136])
    ck.finish()


if __name__ == "__main__":
    main()

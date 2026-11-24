# Computations behind the paper

This table maps every computational statement of the paper to the script that reproduces it,
the range checked at the default bounds, and the output file (all under `outputs/`). Counts are
the numbers of passed instances in the default run (`sh run_all.sh`); every output also records
its bounds. "Stated" means the exact value printed in the paper is asserted. Results are named by
descriptive names with LaTeX labels; theorem numbers are not used.

Conventions: "reduced pairs p <= P" means all 1 <= r < p <= P with gcd(p, r) = 1. "Sets with
d <= D" means all A_0 = {0 < a_1 < ... < d} with the stated number of elements, gcd 1, d <= D.
All arithmetic is exact (integers and rationals).

## Default and quick bounds, and runtime

| Script | Default bounds | Quick bounds | Default runtime |
|---|---|---|---|
| `checks/check_sumsets.py` | three- and four-element sets d <= 20, five-element d <= 12, levels h <= 8 (extended to the threshold for the threshold checks); census d <= 30 (four), d <= 13 (five); families d <= 40 | 9, 7, 5; 14, 9; 14 | 15.4 s |
| `checks/check_family.py` | engine-checked b <= 25; formula-checked b <= 60; three distances b <= 80 | 12; 30; 30 | 6.5 s |
| `checks/check_noncm.py` | four-element d <= 23, five-element d <= 14 | 12, 8 | 1.8 s |
| `checks/check_windows.py` | overshoot form and corrections d <= 11 (five-element d <= 8), h <= 5; hinge (three, four, five elements) d <= 10, h <= 4; gcd reduction on five menus, h <= 5 | 7, 3; 7, 3 | 19.5 s |
| **Total** (`sh run_all.sh`) | | | **about 43 s** (quick: about 2 s) |

Runtimes: one core of a 2-core Linux machine, Python 3.13.

## Sections 3 and 4: sumsets and four-element sets

| Label | Statement | Range checked | Script | Output |
|---|---|---|---|---|
| thm:cm | CM (engine B: Apery set of size d) implies h A_0 = I_h (engines A and C); the converse (a non-CM set has a failing level) is checked in `check_noncm.py` | 939 CM sets (three, four, five elements), h <= 8 | `check_sumsets.py` | `sumsets.json` |
| engines | engine B rebuilds h A_0 and the undershoot of every target from the homogenized Apery set alone; both equal engine A | every set of the main ranges (CM or not), h <= 8, every target in [0, hd] | same | same |
| lem:residue | n in I_h iff u_s <= k <= N - v_s; also h A_0 for CM sets | every modulus of every CM set (540744 targets) | same | same |
| prop:counting | (a) every target; (b) Tot(h) = sum f(H_k); (c) three windows per level; (a) also for I_h at every level of non-CM sets | 487980 + 304344 targets; 8592 levels; 25776 windows | same | same |
| lem:spectrum-mass, lem:merge, thm:spectrum, prop:strict-levels | as stated, at every modulus | 1074 (set, modulus) pairs, 8592 levels | same | same |
| ex:spectrum (a) | {0,2,5,10}, m = 10: positive spectrum (14,6,2,2), bound 22 and Tot(3) = 23 at h = 3; similarly ordered at m = 5 but not at m = 10 (stated); u, v recorded | | same | same |
| ex:spectrum (b) | {0,2,8,11} is CM; M(S_1)+M(S_2) = 18; Tot(h) = 19, 19, 18, 18, ...; (u_6,v_6) = (2,0), (u_1,v_1) = (1,1), interval [1,1]; five inversions, union of their intervals [1,2] (stated) | CM tested by engine B and h <= 15; Tot for h <= 12 | same | same |
| opening of the introduction | {0,5,7,14}: CM, similarly ordered at m = 7 and at m = 14, positive spectrum at m = 14 (32, 16, 8, 5, 3, 1), M(S_1)+M(S_2) = 65, phi_{7,5}(1..7) = 21, 11, 5, 3, 2, 1, 0 from circular gaps, Tot(h) = 32, 48, 56, 61, 64, 65, 65, ... (h = 1..10); {0,4,9,13}: not CM, failing levels 4 to 7, M(S_1)+M(S_2) = 40, Tot(h) = 22, 36, 45, 48, 49, 46, 44, 40, ... (h = 1..12), conductors recorded (stated); {0,2,8,11} as in ex:spectrum (b), positive spectrum at m = 11 (10, 5, 3) (stated) | | same | same |
| opening of Section 3 | {0,2,3,5}: CM, only modulus 5, Tot(h) = 2 for h = 1..8, pairs (1,1) and (2,0) occur (stated) | | same | same |
| eq:reg-cm, prop:threshold | reg = max(u+v) = U+V; U+V-1 = ceil((c_1+c_2)/d) = sigma(A) (sigma computed from its definition); reg = sigma + 1; Tot(h) = mass iff h >= sigma; at every modulus m >= 2 with similar ordering: c_1 + c_2 = (U+V-2)m + s_2 - s_1 + 2 with 2 <= s_2 - s_1 + 2 <= m, and Tot(h) = mass iff h >= sigma; at every modulus m >= 2 with an index attaining U and V: U+V-1 = ceil((c_1+c_2)/m) and reg = sigma + 1 (reg from eq:reg-cm at m = d) | 487 similarly ordered CM sets at m = d and 110 similarly ordered (set, m < d) pairs; sigma = ceil((c_1+c_2)/d) for all 939 CM sets | same | same |
| remark after prop:threshold | an index attaining both U and V (at m = d) gives U+V-1 = sigma and reg = sigma + 1, with or without similar ordering | every CM set with such an index (360 of them not similarly ordered) | same | same |
| remarks after eq:reg-cm and after prop:threshold | {0,1,2,5}: s = 1, 2 give (1,1), Apery points (4,6), (3,7); {0,3,4,6}: CM, similarly ordered, (c_1,c_2) = (6,2), decomposition with empty central interval at h = 1 < ceil(8/6) = 2 (stated) | | same | same |
| lem:diagram, thm:similar, prop:shape | CM iff integral diagram; order ideal of d points; eq:uv-diagram (every residue, r = 0 included); (1) iff (2) iff (3); rectangle or L-shape, a_1 W + a_2 Z = 0 mod d, both notches >= 2 forces not similarly ordered | four-element sets d <= 20 (997 sets) and the 406 sets {0,a,a+1,d} with d <= 30 (1403 diagrams, 1011 CM) | same | same |
| end of Section 4 (after prop:families) | {0,2,5,10}: diagram [0,4] x [0,1], layer 4 values 8 and 11 (stated) | | same | same |
| fig:diagrams | Apery diagrams of {0,5,7,14} (rectangle [0,6] x [0,1], no crossing layer) and {0,2,8,11} (L-shape [0,2] x [0,2] plus (3,0), (3,1); layers 2 and 3 cross 11) (stated) | | same | same |
| prop:families | (a) {0,1,b,d} CM iff e = 0 or q + e >= b, then similarly ordered; (b) A_0 CM and similarly ordered iff d - A_0 is, on every CM four-element set, and every CM {0,a,d-1,d} similarly ordered; (c) {0,a,d-a,d} similarly ordered iff d = 3, pairs of {0,1,2,3} (stated); (d) {0,a,a+1,d} CM, criterion | d <= 40 ((d) diagram criterion for d <= 30) | same | same |
| eq:lambda-test | a four-element set is CM iff len(w + kd) - len(w) >= k for every Apery element w of S_1 modulo d (checked for k <= 2d) | four-element sets d <= 20 | same | same |

## Section 5: the family {0, a, b, 2b}

| Label | Statement | Range checked | Script | Output |
|---|---|---|---|---|
| lem:family-ci | binomial exponents; CM; |h A_0|; c_1, c_2; sigma = b - 1 | coprime a < b <= 60 (engine checks b <= 25, h <= b+2) | `check_family.py` | `family.json` |
| lem:family-ci (computer algebra, optional) | the ideal of the curve, by Groebner-basis elimination, equals the ideal of the two binomials, over the rationals, GF(2) and GF(3); two negative controls detected | coprime a < b <= 12 (45 pairs per field) | `cas/check_ideal_sympy.py` (sympy; not in `run_all.sh`) | `cas/output/ideal.json` |
| lem:family-residue | u_s, v_s modulo b; similar ordering; U = a, V = 2b-a-2; Tot nondecreasing, concave, = M_4 exactly for h >= b-1 | b <= 60 (engine b <= 25) | same | same |
| thm:family-runs | run multiset of [0,2hb] outside h A_0; invariance under a -> b - a | b <= 25, h <= b+2 (3812 levels) | same | same |
| thm:family-spectrum, eq:family-tot | spectrum at m = b; Tot formula; phi nonincreasing | b <= 60 (spectrum), b <= 25 (Tot against engine A) | same | same |
| cor:family-mass | M_4 = 2 Phi(b,a) - binom(b,2) = M(a,b) + M(b,2b-a); phi and Phi also computed from their definition (circular spacings) | b <= 60 (1101 pairs) | same | same |
| thm:three-distance | gap multiset and phi formula, 2 <= h <= b | b <= 80 (103423 gap multisets) | same | same |
| cor:quadratic-pieces | Tot(H-2+j) = Tot(H-2) + 2 phi(H-1) j - delta j^2; tiling of [0, b-1] | b <= 60 (57728 values) | same | same |
| after cor:quadratic-pieces | a = 2, b = 9: 2/9 = [0;4,2], last block of phi is h = 6..9, Tot(4..9) = 120,127,132,135,136,136, third difference of Tot(6..9) is 1 (stated, also by engine A) | | same | same |

## Section 6: non-Cohen-Macaulay sets

| Label | Statement | Range checked | Script | Output |
|---|---|---|---|---|
| prop:failing-interval | per class the failing levels are the corner interval; kappa_1 = u_s, mu_t = v_s; tau - 1 = max(kappa_{j+1} + mu_j - 2) | four-element d <= 23 (1567 sets, 632 non-CM), five-element d <= 14 | `check_noncm.py` | `noncm.json` |
| prop:overshoot | total of I_h = M(S_1) + M(S_2) when hd >= c_1 + c_2; Tot(h) >= mass, equality iff not failing | 13469 levels | same | same |
| before cor:four-overshoot (cited) | (tau - 1) d >= c_1 + c_2 + d for non-CM four-element sets: finite evidence only, the proof is cited | 632 non-CM sets, d <= 23 | same | same |
| cor:four-overshoot | Tot(tau - 1) > M(S_1) + M(S_2) | same 632 sets | same | same |
| ex:five-element | {0,5,9,11,20}: c_1 = 18, c_2 = 44, mass 74, only failing level 3, tau = 4, Tot(1..3) = 53, 69, 74, Tot = 74 for h >= 4, 60 < 62 (stated) | | same | same |

## Section 7: corrections and windows

| Label | Statement | Range checked | Script | Output |
|---|---|---|---|---|
| lem:elementary (ii) | physical undershoot = (T mod g) + g D_h(floor(T/g)) | five physical menus, h <= 5 | `check_windows.py` | `windows.json` |
| prop:overshoot-form | (a)-(d); (b) at every anchor x >= c_1 in I_h | all three-, four-element sets d <= 11 and five-element d <= 8, h <= 5 | same | same |
| prop:correction | sum_W D_h = sum_W D_I + Psi_h(W); Psi_h([0,hd]) > 0 iff failing | same, four windows per level (5280) | same | same |
| lem:hinge | H_C(h; l_0, q_0) equals the direct count of the union, for every union of corners over consecutive positions | three-, four-, five-element sets d <= 10, h <= 4, every 0 <= l_0 <= h + 1 and l_0 - 1 <= q_0 <= h + 1 | same | same |
| thm:hinge | window totals and Psi_h from corners; comparison criterion for adding c versus raising the cap | three-, four-, five-element sets d <= 10, h <= 4, five windows per level, every added c | same | same |

## Statistics for the computation section (recorded, not stated in the text)

| Statistic | Value | Script | Output |
|---|---|---|---|
| CM four-element sets [count, similarly ordered at m = d] | d <= 16: [311, 155]; d <= 22: [798, 311]; d <= 30: [1966, 577] | `check_sumsets.py` | `sumsets.json` |
| CM five-element sets | d <= 13: [293, 139] | same | same |
| non-CM four-element sets, d <= 23 | 632 of 1567; 4977 failing levels in all | `check_noncm.py` | `noncm.json` |

Not reproduced here (not stated in the text): the coverage of similarly ordered sets by the
families of prop:families, the head comparison and non-CM correction statistics.


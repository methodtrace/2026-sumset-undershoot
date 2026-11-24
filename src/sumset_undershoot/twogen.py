"""Two-generator theory: circular spacings phi_{n,a}(h), Phi(n,a), the insertion numbers
lambda = L_j R_j, Euclid's algorithm and the quadratic form, the block formula for the mass,
gap-run multisets, Moree's count and canonical continued fractions.

Notation follows the paper: <x>_n is the least nonnegative residue, phi_{n,a}(h) is the sum of
binomial(g, 2) over the circular spacings g of {0, a, ..., (h-1)a} in Z/nZ (zero for h > n),
and Phi(n, a) = sum_{h=1}^{n} phi_{n,a}(h).
"""
from math import comb, gcd


# ---------------------------------------------------------------- circular spacings

def spacings(n, a, h):
    """Circular spacings of the h points 0, a, ..., (h-1)a in Z/nZ (definition)."""
    pts = sorted({(i * a) % n for i in range(h)})
    assert len(pts) == h, "points must be distinct (h <= n, gcd(n,a) = 1)"
    out = [pts[i + 1] - pts[i] for i in range(h - 1)]
    out.append(n - pts[-1] + pts[0])
    return out


def phi_definition(n, a, h):
    if h > n:
        return 0
    return sum(comb(g, 2) for g in spacings(n, a, h))


def Phi_definition(n, a):
    return sum(phi_definition(n, a, h) for h in range(1, n + 1))


def insertion_numbers(n, a):
    """[(L_j, R_j) for j = 1..n-1] with L_j = min_{d<=j} <da>_n, R_j = min_{d<=j} <-da>_n."""
    out = []
    L = R = n
    for j in range(1, n):
        L = min(L, (j * a) % n)
        R = min(R, (-j * a) % n)
        out.append((L, R))
    return out


def phi_sequence(n, a):
    """[phi_{n,a}(1), ..., phi_{n,a}(n)] from the insertion identity lambda(j) = L_j R_j."""
    seq = [comb(n, 2)]
    for (L, R) in insertion_numbers(n, a):
        seq.append(seq[-1] - L * R)
    return seq


def Phi_insertion(n, a):
    """Phi(n, a) = sum_j j * lambda(j)."""
    return sum(j * L * R for j, (L, R) in enumerate(insertion_numbers(n, a), start=1))


def Phi_rectangles(n, a):
    """Lemma lem:rect: Phi(n,a) = sum_{W,H >= 1} (n - |a[0,W) + [0,H)|) in Z/nZ."""
    total = 0
    for W in range(1, n):
        X = {(i * a) % n for i in range(W)}
        for H in range(1, n):
            cover = {(x + y) % n for x in X for y in range(H)}
            total += n - len(cover)
    return total


# ---------------------------------------------------------------- masses

def delta(h, q, p):
    return (h * q) // p - ((h - 1) * q) // p


def mass_block(p, q):
    """Lemma lem:block, first form: M(p,q) = sum_{h=1}^{p-1} delta_h(q) phi_{p,q}(h)."""
    if p == 1:
        return 0
    ph = phi_sequence(p, q % p)
    return sum(delta(h, q, p) * ph[h - 1] for h in range(1, p))


def mass_block_lambda(p, q):
    """Lemma lem:block, second form: sum_j floor(jq/p) lambda_{p,q}(j)."""
    if p == 1:
        return 0
    return sum(((j * q) // p) * L * R
               for j, (L, R) in enumerate(insertion_numbers(p, q % p), start=1))


# ---------------------------------------------------------------- Euclid and the quadratic form

def P_sum(n, s, c):
    return sum(comb(j * n + s, 2) for j in range(1, c + 1))


def P_closed(n, s, c):
    num = 2 * n * n * c ** 3 + 3 * n * (n + 2 * s - 1) * c * c + (n * n + 6 * n * s - 3 * n + 6 * s * s - 6 * s) * c
    assert num % 12 == 0
    return num // 12


def euclid_chain(q, p):
    """Theorem thm:closed-form data for coprime q > p >= 1 (or q >= 2, p = 1).

    Gives lists n (n_0..n_{K+1}), c (c_1..c_K), u (u_0..u_{K+1}), v (v_0..v_{K+1}),
    P (P_1..P_K), all indexed from their natural start (c[0] is c_1, P[0] is P_1)."""
    n = [q, p]
    c = []
    while n[-1] != 0:
        cq, rem = divmod(n[-2], n[-1])
        c.append(cq)
        n.append(rem)
    K = len(c)
    u = [0, 1]
    v = [1, 0]
    for i in range(1, K + 1):
        u.append(c[i - 1] * u[i] + u[i - 1])
        v.append(c[i - 1] * v[i] + v[i - 1])
    P = [P_closed(n[i], n[i + 1], c[i - 1]) for i in range(1, K + 1)]
    return {"n": n, "c": c, "u": u, "v": v, "P": P, "K": K}


def quadratic_form(q, p):
    """(Phi(q,p), M(p,q), Phi(p,q)) for coprime q > p >= 1 by Theorem thm:closed-form."""
    ch = euclid_chain(q, p)
    K, P, u, v = ch["K"], ch["P"], ch["u"], ch["v"]
    A = sum(P[i - 1] * u[i] * u[i] for i in range(1, K + 1))
    B = sum(P[i - 1] * u[i] * v[i] for i in range(1, K + 1))
    C = sum(P[i - 1] * v[i] * v[i] for i in range(1, K + 1))
    return A, B, C


def Phi(n, a):
    """Phi(n, a) for coprime n >= 1 and any a (depends on a mod n), via the quadratic form."""
    if n == 1:
        return 0
    r = a % n
    return quadratic_form(n, r)[0]


def mass(p, q):
    """M(p, q) for coprime p, q >= 1 via the quadratic form."""
    if p == 1 or q == 1:
        return 0
    lo, hi = min(p, q), max(p, q)
    return quadratic_form(hi, lo)[1]


def cf_digits(p, r):
    """Canonical digits a_1..a_K of r/p = [0; a_1, ..., a_K] (a_K >= 2 unless p = 2... r = 1)."""
    out = []
    while r:
        a, s = divmod(p, r)
        out.append(a)
        p, r = r, s
    return out


def digit_sum(p, r):
    return sum(cf_digits(p, r))


def units(p):
    return [r for r in range(1, p) if gcd(r, p) == 1]


def inverse(r, p):
    return pow(r, -1, p)


# ---------------------------------------------------------------- run multisets

def level_runs(n, a, h):
    """Multiset {g-1 : g >= 2 a circular spacing at level h}, as {length: count}."""
    out = {}
    for g in spacings(n, a, h):
        if g >= 2:
            out[g - 1] = out.get(g - 1, 0) + 1
    return out


def add_multiset(target, src, times=1):
    for k in src:
        if times * src[k]:
            target[k] = target.get(k, 0) + times * src[k]
    return target


def G_multiset(n, a):
    out = {}
    if n >= 2:
        for h in range(1, n):
            add_multiset(out, level_runs(n, a, h))
    return out


def runs_by_levels(p, q):
    """Theorem thm:runs, first formula."""
    out = {}
    if p >= 2:
        for h in range(1, p):
            add_multiset(out, level_runs(p, q % p, h), delta(h, q, p))
    return out


def runs_by_chain(p, q):
    """Theorem thm:runs, third formula: R_{p,q} = sum_i c_i G_{n_i, n_{i+1}} (q > p)."""
    ch = euclid_chain(q, p)
    out = {}
    for i in range(1, ch["K"] + 1):
        ni, ni1 = ch["n"][i], ch["n"][i + 1]
        if ni >= 2:
            add_multiset(out, G_multiset(ni, ni1), ch["c"][i - 1])
    return out


def moree_count(p, q):
    """Moree's count nu = rho*sigma - 1 with rho p + sigma q = 1 + pq."""
    sigma = (pow(q, -1, p)) % p
    sigma = sigma if sigma else p
    rho = (1 + p * q - sigma * q) // p
    assert rho * p + sigma * q == 1 + p * q and 1 <= rho <= q - 1 and 1 <= sigma <= p - 1
    return rho * sigma - 1


def bachman_lengths(p, r):
    """Set of run lengths of <p, q>, q = r mod p, from Bachman's continued-fraction description."""
    rs = [p, r]
    Z = []
    while rs[-1] != 0:
        Zi, rem = divmod(rs[-2], rs[-1])
        Z.append(Zi)
        rs.append(rem)
    t = len(Z)  # r_t = 1, r_{t+1} = 0
    out = set()
    for i in range(1, t):
        for z in range(Z[i - 1]):
            out.add(rs[i - 1] - z * rs[i] - 1)
    for z in range(Z[t - 1] - 1):
        out.add(rs[t - 1] - z - 1)
    return out

"""Engine A (direct sumsets).

h A_0 is computed by repeated shifted OR of a Python integer used as a bitset (bit n set
means n is a sum of h elements of A_0); the undershoot is read off by a linear scan.
"""


def sumset_bits(A0, h):
    """Bitset of h A_0 = {sums of exactly h elements of A_0}; A_0 must contain 0."""
    cur = 1
    for _ in range(h):
        nxt = 0
        for a in A0:
            nxt |= cur << a
        cur = nxt
    return cur


def sumset_series(A0, H):
    """List [bits(0 A_0), bits(1 A_0), ..., bits(H A_0)]."""
    out = [1]
    cur = 1
    for _ in range(H):
        nxt = 0
        for a in A0:
            nxt |= cur << a
        cur = nxt
        out.append(cur)
    return out


def elements(bits):
    s = bin(bits)[:1:-1]
    return [i for i, c in enumerate(s) if c == "1"]


def undershoot_scan(bits, L, U):
    """[D(T) for T in L..U] with D(T) = T - max{s in set : s <= T}; needs an element <= L."""
    below = bits & ((1 << (L + 1)) - 1)
    if below == 0:
        raise ValueError("no element <= L")
    last = below.bit_length() - 1
    out = []
    for T in range(L, U + 1):
        if (bits >> T) & 1:
            last = T
        out.append(T - last)
    return out


def total_undershoot(bits, L, U):
    return sum(undershoot_scan(bits, L, U))

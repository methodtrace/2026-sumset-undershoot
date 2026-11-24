"""Three independent undershoot engines used as oracles.

engine A (direct):        h A_0 by repeated shifted OR of an integer bitset.
engine B (homogenized):   minimum factorization lengths, the Apery set of the homogenized
                          semigroup with respect to (1,0) and (1,d), and h A_0 rebuilt from it.
engine C (intersection):  the gap sets of the two numerical semigroups S_1 and S_2 and the
                          intersection form I_h = {n in [0,hd] : n in S_1, hd - n in S_2}.

The three modules share no code with each other or with the rest of the package.
"""

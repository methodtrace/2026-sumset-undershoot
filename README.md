# Undershoot of h-fold sumsets: code and data

Code and data for the paper

> Prashant Chaudhary and Kapil Khandelwal, *Undershoot of h-fold sumsets: Apery spectra, similar
> ordering and the three-distance theorem*.

The repository recomputes every computational statement made or relied on in the paper: the
worked examples and stated numbers, and finite checks of the identities and inequalities on
explicit ranges against three independent engines. Every check uses exact integer or rational
arithmetic. No result of the paper depends on computation.

## Requirements and running

Python 3.9 or later, standard library only. Nothing to install. One optional check, `cas/check_ideal_sympy.py`, needs the computer-algebra library sympy and is not part of `run_all.sh` (see below).

```sh
sh run_all.sh            # default bounds: under a minute on one core
sh run_all.sh --quick    # small bounds: about 2 seconds
```

`run_all.sh` deletes and rewrites `outputs/`, runs every check, writes
`outputs/MANIFEST.sha256` (SHA-256 of every output file, sorted by path) and compares it with
`expected/MANIFEST.sha256` (or `expected/MANIFEST-quick.sha256` with `--quick`). It exits with
status 0 only if every check passes and the manifest matches. Each script can also be run alone,
for example `python3 checks/check_family.py --quick`; each writes one JSON summary to `outputs/`
and exits with a nonzero status on any failure.

Observed runtime on one core: see `COMPUTATIONS.md`.

## Determinism

Every output is byte-for-byte reproducible: integer and exact rational arithmetic only; JSON
with sorted keys, fixed separators and a trailing newline; no timestamps, timings, paths, host
names or floating point numbers; no randomness (sampled rows and pairs are chosen by fixed
arithmetic rules). Running the suite twice, from different directories or with different
Python versions, produces identical bytes; the committed manifests in `expected/` record them.

## Layout

```
README.md            this file
COMPUTATIONS.md      table: each computational statement, range, script and output file
LICENSE              Apache-2.0 (code)
LICENSE-data.md      CC BY 4.0 (outputs and expected manifests)
CITATION.cff
run_all.sh           runs everything, writes and compares the manifest
src/sumset_undershoot/  the package
  engines/           three independent undershoot engines used as oracles
    direct.py          h A_0 by repeated shifted OR of an integer bitset
    homogenized.py     Apery set of the homogenized semigroup, h A_0 rebuilt from it
    intersection.py    the two numerical semigroups and the intersection form I_h
  semigroup.py       Apery sets, membership, conductor, gaps, mass (definition and Apery formula)
  twogen.py          circular gaps phi, Phi and the two-generator mass, used by the family
  sumsets.py         residue intervals, spectrum, diagrams, corners, corrections, hinge formulas
  family.py          the family {0, a, b, 2b}
  report.py          deterministic JSON reports
checks/              one script per topic (see below)
tools/write_manifest.py
cas/                 optional computer-algebra check of lem:family-ci (needs sympy), with its expected outputs
expected/            committed manifests of the default and quick runs
outputs/             written by the run (JSON summaries and the manifest)
```

## Checks, by topic

Results are named by short descriptive names with the LaTeX label in parentheses (theorem numbers
are not stable). `COMPUTATIONS.md` gives the ranges and counts for each statement.

| Script | Results checked |
|---|---|
| `checks/check_sumsets.py` | Cohen-Macaulay test in three engines (thm:cm); h A_0 and the undershoot rebuilt by the homogenized engine and compared with direct sumsets; residue intervals (lem:residue); exact counting (prop:counting); spectrum and mass (lem:spectrum-mass); merge (lem:merge); spectrum bounds and monotone totals (thm:spectrum); strict levels (prop:strict-levels); examples (ex:spectrum), including {0,2,8,11}; regularity (eq:reg-cm); threshold (prop:threshold) at every similarly ordered modulus, the joint-maximum remark and the {0,1,2,5}, {0,3,4,6} remarks; common diagram (lem:diagram); similar ordering (thm:similar); shape (prop:shape); families (prop:families) and the test (eq:lambda-test); census of CM four- and five-element sets |
| `checks/check_family.py` | structure (lem:family-ci); residue intervals modulo b (lem:family-residue); run lengths (thm:family-runs); spectrum and totals (thm:family-spectrum, eq:family-tot); mass of the family (cor:family-mass), with phi and Phi also computed from their definition; three distances (thm:three-distance); quadratic pieces (cor:quadratic-pieces) and the example a = 2, b = 9 |
| `checks/check_noncm.py` | failing levels form intervals, and the corner formula for tau (prop:failing-interval); overshoot of the mass (prop:overshoot); the cited inequality (tau - 1) d >= c_1 + c_2 + d (finite evidence only); four-element overshoot (cor:four-overshoot); five-element example (ex:five-element) |
| `checks/check_windows.py` | reduction by the gcd (lem:elementary (ii)); head and reflected tail (prop:overshoot-form); correction (prop:correction); hinge formulas (lem:hinge against the direct count of each union; thm:hinge on three-, four- and five-element sets) |

## Optional computer-algebra check

`cas/check_ideal_sympy.py` checks lem:family-ci with Groebner bases in sympy, independently of
`src/`: for every coprime 0 < a < b <= 12, the ideal of the monomial curve of {0, a, b, 2b},
computed by elimination, equals the ideal generated by the two binomials of the lemma, over the
rationals and the fields with 2 and 3 elements; two negative controls must be detected. It takes
about a minute (`--quick`: b <= 6, a few seconds). Its output holds counts and booleans only and
is compared byte for byte with `cas/expected-ideal.json`.

```sh
python3 -m venv .venv && .venv/bin/pip install sympy
.venv/bin/python cas/check_ideal_sympy.py            # ends with "output: identical to expected-ideal.json"
```

## Imported results

- The four-element overshoot (cor:four-overshoot) uses an inequality proved in a companion
  paper; `checks/check_noncm.py` reports finite evidence for it only.

## License

Code: Apache License 2.0 (`LICENSE`). Outputs and expected manifests: CC BY 4.0 (`LICENSE-data.md`).
Please cite the release using `CITATION.cff` (doi:10.5281/zenodo.22952563).

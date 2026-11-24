"""Deterministic reporting for the check scripts.

Every check script builds a Checker, records named assertions and stated values, and
writes one JSON summary. JSON is written with sorted keys, fixed separators and a
trailing newline; it never contains timings, paths, host names or floating point numbers.
"""
import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAX_LOGGED_FAILURES = 20


def to_jsonable(x):
    from fractions import Fraction
    if isinstance(x, bool) or x is None or isinstance(x, (int, str)):
        return x
    if isinstance(x, Fraction):
        return str(x)
    if isinstance(x, float):
        raise TypeError("floating point values are not allowed in outputs")
    if isinstance(x, dict):
        return {str(k): to_jsonable(x[k]) for k in x}
    if isinstance(x, (list, tuple)):
        return [to_jsonable(v) for v in x]
    if isinstance(x, (set, frozenset)):
        return sorted(to_jsonable(v) for v in x)
    if hasattr(x, "c") and isinstance(getattr(x, "c"), dict):
        return [[i, j, str(x.c[(i, j)])] for (i, j) in sorted(x.c)]
    raise TypeError("cannot serialise %r" % type(x))


def dumps(obj):
    return json.dumps(to_jsonable(obj), sort_keys=True, indent=1, separators=(",", ": "),
                      ensure_ascii=True) + "\n"


def write_json(relpath, obj):
    path = os.path.join(ROOT, "outputs", relpath)
    text = dumps(obj)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write(text)


def standard_parser(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--quick", action="store_true",
                        help="use the small bounds (a few seconds) instead of the default bounds")
    return parser


class Checker:
    """Collects named exact checks. Each name accumulates a count of passed instances."""

    def __init__(self, name, description, bounds):
        self.name = name
        self.description = description
        self.bounds = dict(bounds)
        self.counts = {}
        self.failures = []
        self.nfail = 0
        self.values = {}

    def check(self, label, condition, detail=None):
        if condition:
            self.counts[label] = self.counts.get(label, 0) + 1
        else:
            self.nfail += 1
            self.counts.setdefault(label, 0)
            if len(self.failures) < MAX_LOGGED_FAILURES:
                self.failures.append({"check": label, "detail": detail})
        return bool(condition)

    def equal(self, label, got, expected, detail=None):
        ok = got == expected
        if not ok:
            detail = {"got": got, "expected": expected, "where": detail}
        return self.check(label, ok, detail)

    def stated(self, label, got, expected):
        """Assert a value printed in the paper and record it in the summary."""
        self.equal(label, got, expected)
        self.values[label] = got

    def record(self, label, value):
        """Record a computed value that the paper does not state (informational)."""
        self.values[label] = value

    def finish(self):
        summary = {
            "check": self.name,
            "description": self.description,
            "bounds": self.bounds,
            "passed_instances": self.counts,
            "total_instances": sum(self.counts[k] for k in self.counts),
            "failures": self.failures,
            "failure_count": self.nfail,
            "values": self.values,
            "status": "pass" if self.nfail == 0 else "FAIL",
        }
        write_json(self.name + ".json", summary)
        line = "%-34s %s  (%d instances, %d failures)" % (
            self.name, summary["status"], summary["total_instances"], self.nfail)
        print(line)
        sys.exit(0 if self.nfail == 0 else 1)

#!/usr/bin/env python3
"""Write outputs/MANIFEST.sha256 (SHA-256 of every output file, sorted by path) and, if an
expected manifest is given and exists, compare the two. Exit status 1 on any difference."""
import hashlib
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")


def manifest_lines():
    paths = []
    for base, dirs, files in os.walk(OUT):
        dirs.sort()
        for name in files:
            rel = os.path.relpath(os.path.join(base, name), OUT).replace(os.sep, "/")
            if rel != "MANIFEST.sha256":
                paths.append(rel)
    lines = []
    for rel in sorted(paths):
        with open(os.path.join(OUT, rel), "rb") as fh:
            lines.append("%s  %s\n" % (hashlib.sha256(fh.read()).hexdigest(), rel))
    return lines


def main():
    lines = manifest_lines()
    with open(os.path.join(OUT, "MANIFEST.sha256"), "w", encoding="ascii", newline="\n") as fh:
        fh.writelines(lines)
    print("manifest: %d output files" % len(lines))
    if len(sys.argv) > 1:
        expected = os.path.join(ROOT, sys.argv[1])
        if not os.path.exists(expected):
            print("manifest: MISSING expected manifest %s (no comparison made)" % sys.argv[1])
            return 1
        with open(expected, encoding="ascii") as fh:
            want = fh.readlines()
        if want == lines:
            print("manifest: identical to %s" % sys.argv[1])
            return 0
        got_map = {l.split("  ", 1)[1].strip(): l.split("  ", 1)[0] for l in lines}
        want_map = {l.split("  ", 1)[1].strip(): l.split("  ", 1)[0] for l in want}
        for key in sorted(set(got_map) | set(want_map)):
            if got_map.get(key) != want_map.get(key):
                print("manifest: DIFFERS %s" % key)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/bin/sh
# Run every check, write outputs/MANIFEST.sha256 and compare it with
# the committed expected manifest. Usage: sh run_all.sh [--quick]
set -u
cd "$(dirname "$0")" || exit 2
PY="${PYTHON:-python3}"
PYTHONDONTWRITEBYTECODE=1
export PYTHONDONTWRITEBYTECODE
MODE=""
EXPECTED="expected/MANIFEST.sha256"
if [ "${1:-}" = "--quick" ]; then
    MODE="--quick"
    EXPECTED="expected/MANIFEST-quick.sha256"
fi
rm -rf outputs
mkdir outputs
status=0
start=$(date +%s)
for script in \
    checks/check_sumsets.py \
    checks/check_family.py \
    checks/check_noncm.py \
    checks/check_windows.py
do
    "$PY" "$script" $MODE || status=1
done
"$PY" tools/write_manifest.py "$EXPECTED" || status=1
end=$(date +%s)
echo "elapsed: $((end - start)) s"
if [ "$status" -eq 0 ]; then
    echo "ALL CHECKS PASSED"
else
    echo "SOME CHECKS FAILED"
fi
exit "$status"

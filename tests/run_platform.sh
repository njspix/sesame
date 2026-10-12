#!/bin/sh
# A named --platform and a bead count the registry does not know
# (zhou-lab/sesame#2). With --platform the sample is read anyway and one
# warning says how many of the platform's probes have no address in this IDAT
# (they come out NA); without it, a count the registry does not know -- or a
# second platform in the batch -- is still refused.
#
# Real IDATs: EPIC+ (1,055,583 beads) is a revision of EPIC (1,051,815) and
# loses only a handful of probes; HM450 read as EPIC loses most of them, which
# the warning has to say. Needs the store and the test IDATs; SKIPs without.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
[ -f "$yhome/EPIC/EPIC.ordering.tsv.gz" ] || {
    echo "SKIP platform: no EPIC ordering in $yhome"; exit 0; }
[ -d "$idats" ] || { echo "SKIP platform: $idats not found"; exit 0; }

epic="$idats/EPIC/GSM2995280_201868590258_R01C01"
plus="$idats/EPIC+/202761280074_R08C01"
hm450="$idats/HM450/3999492009_R01C01"

PASS=0; FAIL=0
## label, expected exit status, expected stderr fragment ("" = no warning),
## then the preprocess arguments
check() {
    lbl=$1; want_rc=$2; want=$3; shift 3
    rc=0
    YAME_DATA_HOME="$yhome" "$bin" preprocess --prep "" --output beta \
        --out "$work/$lbl" "$@" >/dev/null 2>"$work/$lbl.err" || rc=$?
    if [ "$rc" -ne "$want_rc" ]; then
        echo "FAIL $lbl: exit $rc, wanted $want_rc"
        sed 's/^/    /' "$work/$lbl.err" | head -3; FAIL=$((FAIL+1)); return
    fi
    if [ -z "$want" ]; then
        if grep -q "warning" "$work/$lbl.err"; then
            echo "FAIL $lbl: warned"; sed 's/^/    /' "$work/$lbl.err" | head -3
            FAIL=$((FAIL+1)); return
        fi
    elif ! grep -qF -- "$want" "$work/$lbl.err"; then
        echo "FAIL $lbl: stderr lacks '$want'"
        sed 's/^/    /' "$work/$lbl.err" | head -3; FAIL=$((FAIL+1)); return
    fi
    echo "ok   $lbl"; PASS=$((PASS+1))
}

## the registry's own count: no warning
check epic-as-epic 0 "" --platform EPIC "$epic"
## EPIC+ named as EPIC: read, with the count of probes it lacks
check epicplus-as-epic 0 \
    "has 1055583 beads, not EPIC's 1051815; 2 of 866553 EPIC probes" \
    --platform EPIC "$plus"
## HM450 named as EPIC: read, and the warning shows most probes are missing
check hm450-as-epic 0 "has 622399 beads, not EPIC's 1051815" \
    --platform EPIC "$hm450"
pct=$(sed -n 's/.*EPIC probes (\([0-9.e+-]*\)%).*/\1/p' "$work/hm450-as-epic.err")
if awk -v p="$pct" 'BEGIN { exit !(p > 40) }'; then
    echo "ok   hm450-as-epic reports $pct% NA"; PASS=$((PASS+1))
else
    echo "FAIL hm450-as-epic: reported '$pct'% NA, expected most probes"
    FAIL=$((FAIL+1))
fi
## no --platform: an unknown count is still refused ...
check epicplus-detect 1 "cannot identify platform from 1055583 beads" "$plus"
## ... and so is a second platform in a detected batch (that sample fails)
check mixed-detect 1 "not the batch platform EPIC" "$epic" "$plus"

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]

#!/bin/bash
# Run the FactorialMomentsTask chain (Run 3 MC, LHC25f3 / 545312) on every
# AO2D listed in data/545312/file.list.
#
# Usage:  ./run_mc.sh [AO2D file ...]      (default: all files from file.list)
# Output: results/<AO2D_NAME>/AnalysisResults.root + run.log
#
# Notes:
#  * --configuration must be given to every device (single JSON, all devices).
#  * --aod-file must be given to every device, otherwise the reader sees "".
#  * Device process switches (processMonteCarlo / processCentralityRun3) are
#    set in the JSON, NOT on the CLI: a CLI bool flag always ends up "true".
set -eo pipefail
cd "$(dirname "$0")"
# no `set -u` here: the O2Physics init.sh sourced by ./enable is not -u clean
source ./enable >/dev/null 2>&1
export ALICEO2_CCDB_NOTOKENCHECK=1

JSON="--configuration json://$PWD/mc_chain_config.json"
# default DPL shm segment is far too small for the MC tables (McParticles
# message ~55 MB) -> give every device the same, large segment.
SHM="--shm-segment-size 8589934592"   # 8 GiB, the option is in BYTES

if [ "$#" -gt 0 ]; then
  FILES=("$@")
else
  mapfile -t FILES < data/545312/file.list
fi

for f in "${FILES[@]}"; do
  f=$(readlink -f "$f")
  tag=$(basename "$f" .root)
  out="results/$tag"
  mkdir -p "$out"
  AOD="--aod-file $f"
  echo "=== $tag -> $out ==="
  rc=0
  (
    cd "$out"
    echo | \
    o2-analysis-event-selection-service  -b ${JSON} ${AOD} ${SHM} --evselOpts.isMC 1 | \
    o2-analysis-propagationservice       -b ${JSON} ${AOD} ${SHM} | \
    o2-analysis-trackselection           -b ${JSON} ${AOD} ${SHM} | \
    o2-analysis-multcenttable            -b ${JSON} ${AOD} ${SHM} | \
    o2-analysis-cf-factorial-moments     -b ${JSON} ${AOD} ${SHM} > run.log 2>&1
  ) || rc=$?
  if [ "$rc" -ne 0 ]; then
    echo "  FAILED rc=$rc (see $out/run.log)"
    grep -m5 -E "FATAL|First error|Couldn.t open file|no input " "$out/run.log" || true
    exit "$rc"
  fi
  ls -la "$out/AnalysisResults.root"
done

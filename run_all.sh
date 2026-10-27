#!/usr/bin/env bash
# =====================================================================
# run_all.sh  -  run the whole pipeline in order.
# The workshop is meant to be run STEP BY STEP (see README). This is
# only here for a quick full pass once you understand the steps.
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"
for s in steps/0*.sh steps/1*.sh; do
  echo; echo "##################  $s  ##################"
  bash "$s"
done
echo; echo "All steps complete. See results/ ."

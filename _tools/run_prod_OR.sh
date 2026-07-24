#!/bin/bash
# Full production (extract -> render) for the OR model (model0 = Mu10ORDoubleMu).
# All plots: tables + mass pages, no test cap, no skip_mass.
set -o pipefail

BASE=/vols/cms/jtafoyav/parking/bdt/icenet
TOOLS=$BASE/_tools
LOG=$TOOLS/prod_run_OR.log

cd "$BASE" || exit 1
source /home/hep/jtafoyav/vols/dependencies/anaconda3/etc/profile.d/conda.sh
conda activate icenet
source setenv.sh
source setproxy.sh

cd "$TOOLS" || exit 1

echo "==================================================================" | tee -a "$LOG"
echo "[$(date)] START full production, model0 (Mu10ORDoubleMu)"          | tee -a "$LOG"
echo "==================================================================" | tee -a "$LOG"

echo "[$(date)] EXTRACT begin"                                            | tee -a "$LOG"
root -l -b -q 'getSignalEff_extract.C("tables model0")' >> "$LOG" 2>&1
rc=$?
echo "[$(date)] EXTRACT done (rc=$rc)"                                    | tee -a "$LOG"
if [ $rc -ne 0 ]; then
    echo "[$(date)] EXTRACT FAILED -> aborting, not rendering"           | tee -a "$LOG"
    exit $rc
fi

echo "[$(date)] RENDER begin"                                             | tee -a "$LOG"
root -l -b -q 'getSignalEff_render.C("tables model0")' >> "$LOG" 2>&1
rc=$?
echo "[$(date)] RENDER done (rc=$rc)"                                     | tee -a "$LOG"

echo "[$(date)] ALL DONE (rc=$rc)"                                        | tee -a "$LOG"
exit $rc

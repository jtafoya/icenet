import pickle
import sys
sys.path.insert(0, "/home/hep/jtafoyav/vols/parking/bdt/icenet")
from icenet.tools import aux

run_id = "2026-03-03_20-54-49_lxdgpu00"
#type="XGB"
type="XGB-NOJETS"

MODELDIR  = f"../checkpoint/dqcd/config__tune0_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2018_firstTry/{run_id}/{type}"
OUTFILE   = type.replace('-','_')+f"_2018_scenarioA_Reproduction.model"  # or .json, .ubj

print(f"Searching model in directory: {MODELDIR}")

# Find best epoch automatically (minimum val loss)
pkl_path = aux.create_model_filename(path=MODELDIR, label=type, epoch=-1, filetype=".pkl")

with open(pkl_path, "rb") as f:
    data = pickle.load(f)

data["model"].save_model(OUTFILE)
print(f"Saved to {OUTFILE}")


######
# -> Passing model to dqcd analysis framework:
#
# Models are stored in dqcd/nanoaod_base_analysis/data/cmssw/CMSSW_13_0_13/src/DQCD/Modules/data/
#
# The model is imported by the dqcd-modules repository: https://github.com/ic-dqcd/dqcd-modules/blob/main/python/BDTULinferenceFinal.py
# Additional variables (e.g. MET_*) used in the training are also added to BDTULinferenceFinal.py
# N.B. the full path to such file looks like e.g. /home/hep/jtafoyav/vols/parking/nanoaod_base_analysis_13/dqcd/nanoaod_base_analysis/data/cmssw/CMSSW_13_0_13/src/DQCD/Modules/python
#
# Once you have this BDT inference module set up in the CMSSW in nanoaod_base_analysis, you can call it using a module in your dqcd directory like this: https://github.com/ic-dqcd/dqcd/blob/main/config/modules_nojet.yaml#L84
#
#
#####

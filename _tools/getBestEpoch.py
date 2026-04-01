import pickle
import sys
sys.path.insert(0, "/home/hep/jtafoyav/vols/parking/bdt/icenet")
from icenet.tools import aux



##### 2018, with MET variables

##type="XGB"
#type="XGB-NOJETS"

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2018_firstTry/2026-03-03_20-54-49_lxdgpu00/{type}"	#v0
#OUTFILE   = type.replace('-','_')+f"_2018_scenarioA_Reproduction.model"  # or .json, .ubj
##


#MODELDIR  = f"../checkpoint/dqcd//{type}"	#v_index


##### 2024, without MET variables

##type="XGB-NOMET"
#type="XGB-NOJETS-NOMET"

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_noMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_noMET/2026-03-16_09-15-23_lxbgpu00/{type}"	#v0
#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_noMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_noMET/2026-03-18_11-48-13_lxdgpu01/{type}"	#v1
#OUTFILE   = type.replace('-','_')+f"_2024_Mu10orDoubleMu_noMET.model"  # or .json, .ubj
## best model at epoch [144] with validation loss = 0.2281

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_Mu10_noMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_noMET/2026-03-16_09-15-45_lxdgpu00/{type}"	#v0
#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_Mu10_noMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_noMET/2026-03-18_11-26-53_lxdgpu01/{type}"	#v1
#OUTFILE   = type.replace('-','_')+f"_2024_Mu10_noMET.model"  # or .json, .ubj
## best model at epoch [83] with validation loss = 0.1724

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_DoubleMu_noMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_noMET/2026-03-16_09-15-23_lxbgpu01/{type}"	#v0
#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_DoubleMu_noMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_noMET/2026-03-18_11-03-16_lxdgpu01/{type}"	#v1
#OUTFILE   = type.replace('-','_')+f"_2024_DoubleMu_noMET.model"  # or .json, .ubj
## best model at epoch [141] with validation loss = 0.2313


##### 2024, with MET variables

##type="XGB"
type="XGB-NOJETS"

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_withMET/2026-03-15_17-26-06_lxbgpu00/{type}"	#v0
#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_withMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_withMET/2026-03-18_11-52-51_lxcgpu02/{type}"	#v1
#OUTFILE   = type.replace('-','_')+f"_2024_Mu10orDoubleMu_withMET.model"  # or .json, .ubj
## best model at epoch [147] with validation loss = 0.2298

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_Mu10_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10/2026-03-15_17-26-42_lxbgpu01/{type}"	#v0
#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_Mu10_withMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_withMET/2026-03-18_11-34-00_lxcgpu02/{type}"	#v1
#OUTFILE   = type.replace('-','_')+f"_2024_Mu10_withMET.model"  # or .json, .ubj
## best model at epoch [85] with validation loss = 0.1698

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_DoubleMu_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu/2026-03-15_17-17-30_lxdgpu00/{type}"	#v0
#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_DoubleMu_withMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_withMET/2026-03-18_11-02-10_lxcgpu03/{type}"	#v1
#OUTFILE   = type.replace('-','_')+f"_2024_DoubleMu_withMET.model"  # or .json, .ubj
## best model at epoch [124] with validation loss = 0.2313


##### 2024, with MET variables, v2 (with muonSV_phi - MET_phi)

##type="XGB"
type="XGB-NOJETS"

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_withMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_withMET_v2/2026-03-29_22-07-42_lxcgpu02/{type}"	#v2
#OUTFILE   = type.replace('-','_')+f"_2024_Mu10orDoubleMu_withMET_v2.model"  # or .json, .ubj
## best model at epoch [129] with validation loss = 0.2096

#MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_Mu10_withMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_withMET_v2/2026-03-29_22-08-14_lxdgpu01/{type}"	#v2
#OUTFILE   = type.replace('-','_')+f"_2024_Mu10_withMET_v2.model"  # or .json, .ubj
## best model at epoch [80] with validation loss = 0.1622

MODELDIR  = f"../checkpoint/dqcd/config__tune0_2024_DoubleMu_withMET_new.yml/modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_withMET_v2/2026-03-29_22-14-40_lxbgpu00/{type}"	#v2
OUTFILE   = type.replace('-','_')+f"_2024_DoubleMu_withMET_v2.model"  # or .json, .ubj
## best model at epoch [105] with validation loss = 0.2176

######

print(f"Searching model in directory: {MODELDIR}")

# Find best epoch automatically (minimum val loss)
pkl_path = aux.create_model_filename(path=MODELDIR, label=type, epoch=-1, filetype=".pkl")

with open(pkl_path, "rb") as f:
    data = pickle.load(f)

data["model"].save_model(f"./models/{OUTFILE}")
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

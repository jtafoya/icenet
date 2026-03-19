import sys
sys.path.insert(0, '/vols/cms/jtafoyav/parking/bdt/icenet')
import pickle, numpy as np

### 2018
#PATH_TO_EVAL="/vols/cms/jtafoyav/parking/bdt/icenet/figs/dqcd/config__tune0_new.yml/inputmap__mc_map__scenarioA_all.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2018_firstTry/2026-03-03_20-54-49_lxdgpu00/eval/eval_results.pkl"

### 2024, with MET
#PATH_TO_EVAL="/vols/cms/jtafoyav/parking/bdt/icenet/figs/dqcd/config__tune0_2024_withMET_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_withMET/2026-03-18_11-52-51_lxcgpu02/eval/eval_results.pkl"
#PATH_TO_EVAL="/vols/cms/jtafoyav/parking/bdt/icenet/figs/dqcd/config__tune0_2024_Mu10_withMET_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_withMET/2026-03-18_11-34-00_lxcgpu02/eval/eval_results.pkl"
#PATH_TO_EVAL="/vols/cms/jtafoyav/parking/bdt/icenet/figs/dqcd/config__tune0_2024_DoubleMu_withMET_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_withMET/2026-03-18_11-02-10_lxcgpu03/eval/eval_results.pkl"

### 2024, no MET
#PATH_TO_EVAL="/vols/cms/jtafoyav/parking/bdt/icenet/figs/dqcd/config__tune0_2024_noMET_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_noMET/2026-03-18_11-48-13_lxdgpu01/eval/eval_results.pkl"
#PATH_TO_EVAL="/vols/cms/jtafoyav/parking/bdt/icenet/figs/dqcd/config__tune0_2024_Mu10_noMET_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_noMET/2026-03-18_11-26-53_lxdgpu01/eval/eval_results.pkl"
PATH_TO_EVAL="/vols/cms/jtafoyav/parking/bdt/icenet/figs/dqcd/config__tune0_2024_DoubleMu_noMET_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_noMET/2026-03-18_11-03-16_lxdgpu01/eval/eval_results.pkl"



print(f"Analysing {PATH_TO_EVAL}")
#with open('/PATH/TO/2026-03-18_11-26-53_lxdgpu01/eval/eval_results.pkl', 'rb') as f:
with open(PATH_TO_EVAL, 'rb') as f:
  r = pickle.load(f)

# Two models: index 0 = XGB(-NOMET), index 1 = XGB-NOJETS(-NOMET)
# Use 'inclusive' for a global working point (all signal points combined)
# or e.g. 'F[0]: mpi = 4 & mA = 1.33 & ctau = 1' for a specific point

metric = r['roc_mstats']['inclusive'][1]   # model index 1 = XGB-NOJETS(-NOMET)

fpr        = metric.fpr         # background efficiency array
tpr        = metric.tpr         # signal efficiency array
thresholds = metric.thresholds  # BDT score thresholds

# Find threshold at fpr = 1e-4 (interpolating)
idx = np.searchsorted(fpr, 1e-4)
print(f"Threshold at fpr=1e-4: {thresholds[idx]:.4f}")
print(f"  -> signal eff (tpr): {tpr[idx]:.4f}")

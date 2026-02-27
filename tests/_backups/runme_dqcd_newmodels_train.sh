#!/bin/sh
#
# Execute training and evaluation for the DQCD analysis
#
# Remember to execute first: runme_dqcd_newmodels_init_yaml.sh (only once, and just once)
__conda_setup="$('/vols/cms/khl216/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/vols/cms/khl216/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/vols/cms/khl216/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/vols/cms/khl216/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate icenet

ICEPATH="/vols/cms/khl216/icenet_Mikael/icenet"
cd $ICEPATH
echo "$(pwd)"
source $ICEPATH/setenv.sh

CONFIG="tune0_new.yml"
DATAPATH="/vols/cms/khl216"

CONDITIONAL=0
MAX=5000000    # Tune according to maximum CPU RAM available

python analysis/dqcd.py --runmode genesis  --maxevents $MAX --inputmap mc_map__scenarioA_all.yml --config $CONFIG --datapath $DATAPATH
python analysis/dqcd.py --runmode train    --maxevents $MAX --inputmap mc_map__scenarioA_all.yml --modeltag scenarioA_all_no_DA_old_BDT_with_dRmSV_final --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL
python analysis/dqcd.py --runmode eval     --maxevents $MAX --inputmap mc_map__scenarioA_all.yml --modeltag scenarioA_all_no_DA_old_BDT_with_dRmSV_final --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL
python analysis/dqcd.py --runmode optimize --maxevents $MAX --inputmap mc_map__scenarioA_all.yml --modeltag scenarioA_all_no_DA_old_BDT_with_dRmSV_final --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL

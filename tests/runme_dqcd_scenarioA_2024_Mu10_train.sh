#!/bin/sh
#
# Execute training and evaluation for the DQCD analysis
#
# Remember to execute first: runme_dqcd_newmodels_init_yaml.sh (only once, and just once)

__conda_setup="$('/home/hep/jtafoyav/vols/dependencies/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/hep/jtafoyav/vols/dependencies/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/hep/jtafoyav/vols/dependencies/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/hep/jtafoyav/vols/dependencies/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate icenet

ICEPATH="/home/hep/jtafoyav/vols/parking/bdt/icenet"
cd $ICEPATH
echo "$(pwd)"
source $ICEPATH/setenv.sh
source $ICEPATH/setproxy.sh

CONFIG="tune0_2024_Mu10_new.yml"
DATAPATH="/vols/cms/khl216"

CONDITIONAL=0
MAX=5000000    # Tune according to maximum CPU RAM available


MODELTAG="scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10"

python analysis/dqcd_2024.py --runmode genesis  --maxevents $MAX --inputmap mc_map__scenarioA_all_2024.yml --config $CONFIG --datapath $DATAPATH
python analysis/dqcd_2024.py --runmode train    --maxevents $MAX --inputmap mc_map__scenarioA_all_2024.yml --modeltag $MODELTAG --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL
python analysis/dqcd_2024.py --runmode eval     --maxevents $MAX --inputmap mc_map__scenarioA_all_2024.yml --modeltag $MODELTAG --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL
python analysis/dqcd_2024.py --runmode optimize --maxevents $MAX --inputmap mc_map__scenarioA_all_2024.yml --modeltag $MODELTAG --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL

#!/bin/sh
#
# Execute distributed deployment for the DQCD analysis
#
# Remember to execute first: runme_dqcd_vector_init_yaml.sh (only once, and just once)

export HTC_PROCESS_ID=$1
export HTC_QUEUE_SIZE=$2

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

export GRID_ID=$HTC_PROCESS_ID
export GRID_NODES=$HTC_QUEUE_SIZE

CONFIG="tune0_2024_new.yml"
DATAPATH="/vols/cms/khl216"

CONDITIONAL=0


MODELTAG="scenarioA_all_no_DA_old_BDT_2024"

python analysis/dqcd_deploy_2024.py --runmode deploy --use_conditional $CONDITIONAL --inputmap 'include/scenarioA_all_model_points_2024_deploy.yml' --modeltag $MODELTAG --grid_id $GRID_ID --grid_nodes $GRID_NODES --config $CONFIG --datapath $DATAPATH
#python analysis/dqcd_deploy_2024.py --runmode deploy --use_conditional $CONDITIONAL --inputmap 'include/QCD_2024_deploy.yml' --modeltag $MODELTAG --grid_id $GRID_ID --grid_nodes $GRID_NODES --config $CONFIG --datapath $DATAPATH

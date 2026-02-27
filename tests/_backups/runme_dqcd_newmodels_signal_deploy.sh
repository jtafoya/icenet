#!/bin/sh
#
# Execute distributed deployment for the DQCD analysis
# Run with: source runme.sh

# Remember to execute first: runme_dqcd_newmodels_init_yaml.sh (only once, and just once)

export HTC_PROCESS_ID=$1
export HTC_QUEUE_SIZE=$2

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
#CONFIG="tune0_new_DA.yml"
DATAPATH="/vols/cms/khl216"
CONDITIONAL=0

python analysis/dqcd_deploy.py --runmode deploy --use_conditional $CONDITIONAL --inputmap 'include/scenarioA_all_model_points_deploy.yml' --modeltag scenarioA_all_no_DA_old_BDT_final_variable_set2 --grid_id $GRID_ID --grid_nodes $GRID_NODES --config $CONFIG --datapath $DATAPATH
#python analysis/dqcd_deploy.py --runmode deploy --use_conditional $CONDITIONAL --inputmap 'include/scenarioA_all_model_points_deploy.yml' --modeltag scenarioA_all_DA_new --grid_id $GRID_ID --grid_nodes $GRID_NODES --config $CONFIG --datapath $DATAPATH
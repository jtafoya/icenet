executable	= tests/runme_dqcd_scenarioA_2024_train.sh
output		= iceqsub/output/outputfile.$(CLUSTER)
error		= iceqsub/error/errorfile.$(CLUSTER)
log		= iceqsub/log/dqcd_scenarioA_train_2024.job.$(CLUSTER).log

#Resource request
#request_gpus	= 1   # combined moved to high-mem CPU (needs ~90GB, GPU nodes cap at 125GB and are scarce)
request_memory	= 150000
request_cpus	= 4
+lxfw = true
#+MaxRuntime	= 86400
+MaxRuntime	= 86400
periodic_release = (HoldReasonCode == 34) && (HoldReasonSubCode == 0)

#Notification
notify_user	= j.tafoya.vargas@cern.ch
notification	= Complete

queue

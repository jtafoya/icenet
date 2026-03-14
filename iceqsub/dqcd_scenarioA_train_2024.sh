executable	= tests/runme_dqcd_scenarioA_2024_train.sh
output		= iceqsub/output/outputfile.$(CLUSTER)
error		= iceqsub/error/errorfile.$(CLUSTER)
log		= iceqsub/log/dqcd_scenarioA_train_2024.job.$(CLUSTER).log

#Resource request
request_gpus	= 1
request_memory	= 100G
+MaxRuntime	= 86400
periodic_release = (HoldReasonCode == 34) && (HoldReasonSubCode == 0)

#Notification
notify_user	= j.tafoya.vargas@cern.ch
notification	= Complete

queue 

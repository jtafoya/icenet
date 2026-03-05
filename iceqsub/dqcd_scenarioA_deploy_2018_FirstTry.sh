executable	= tests/runme_dqcd_scenarioA_deploy.sh
arguments	= "$(PROCESS) 1000"
output		= iceqsub/output/outputfile.$(CLUSTER)
error		= iceqsub/error/errorfile.$(CLUSTER)
log		= iceqsub/log/dqcd_scenarioA_train_2018_FirstTry.job.$(CLUSTER).log

#Resource request
+MaxRuntime = 6500
periodic_release = (HoldReasonCode == 34) && (HoldReasonSubCode == 0)

queue 1000

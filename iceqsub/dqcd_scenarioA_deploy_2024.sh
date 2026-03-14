executable	= tests/runme_dqcd_scenarioA_2024_deploy.sh
arguments	= "$(PROCESS) 100"
#arguments	= "$(PROCESS) 1000"
output		= iceqsub/output/outputfile.$(CLUSTER)
error		= iceqsub/error/errorfile.$(CLUSTER)
log		= iceqsub/log/dqcd_scenarioA_deploy_2024.job.$(CLUSTER).log

#Resource request
+MaxRuntime = 6500
#+MaxRuntime = 40000
periodic_release = (HoldReasonCode == 34) && (HoldReasonSubCode == 0)

queue 100
#queue 1000

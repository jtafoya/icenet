#Signal
#NJobs	= 100
#MaxRuntimeVal	= 6500

#QCD
NJobs	= 1000
#NJobs	= 4000
MaxRuntimeVal	= 12000
#MaxRuntimeVal	= 40000

executable	= tests/runme_dqcd_scenarioA_2024_DoubleMu_deploy_signal.sh
arguments	= "$(PROCESS) $(NJobs)"
output		= iceqsub/output/outputfile.$(CLUSTER)
error		= iceqsub/error/errorfile.$(CLUSTER)
log		= iceqsub/log/dqcd_scenarioA_deploy_2024_DoubleMu_signal.job.$(CLUSTER).log

#Resource request
+MaxRuntime = $(MaxRuntimeVal)
periodic_release = ((HoldReasonCode == 34) && (HoldReasonSubCode == 0)) || (HoldReasonCode == 26)

queue $(NJobs)

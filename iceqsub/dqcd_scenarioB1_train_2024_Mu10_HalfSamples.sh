executable	= tests/runme_dqcd_scenarioB1_2024_Mu10_train_HalfSamples.sh
output		= iceqsub/output/outputfile.$(CLUSTER)
error		= iceqsub/error/errorfile.$(CLUSTER)
log		= iceqsub/log/dqcd_scenarioB1_train_2024_Mu10_HalfSamples.job.$(CLUSTER).log

#Resource request
#request_gpus	= 1   # HalfSamples: reads full 50% QCD (~430M evts) -> ~165GB, GPU caps at 125GB
request_memory	= 185000
request_cpus	= 4
+lxfw = true
#+MaxRuntime	= 86400
+MaxRuntime	= 86400
periodic_release = (HoldReasonCode == 34) && (HoldReasonSubCode == 0)

#Notification
notify_user	= j.tafoya.vargas@cern.ch
notification	= Complete

queue

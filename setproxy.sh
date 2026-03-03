#!/bin/sh
#

# Set grid certificate proxy
# (use after loading setenv.sh)

BASE_DIR=/home/hep/jtafoyav/vols/parking/bdt/icenet

source /cvmfs/cms.cern.ch/cmsset_default.sh
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
source /cvmfs/grid.cern.ch/alma9-ui-current/etc/profile.d/setup-alma9-test.sh

echo "Be sure to create a valid GRID certificate before submitting jobs!"
echo "DO: voms-proxy-init --rfc -voms cms --valid 192:00"
#voms-proxy-init --rfc -voms cms --valid 192:00




export CMT_ON_HTCONDOR="1"
#export X509_USER_PROXY="/vols/cms/jtafoyav/parking/nanoaod_base_analysis_13/dqcd/nanoaod_base_analysis/x509up"
export X509_USER_PROXY="$BASE_DIR/x509up"

#    source /home/hep/jtafoyav/.bashrc
#    source /home/hep/jtafoyav/.bash_profile

#    # Jaime uses this
#    source /cvmfs/grid.cern.ch/umd-c7ui-latest/etc/profile.d/setup-c7-ui-example.sh

# Display proxy info
echo "Using proxy: $X509_USER_PROXY"
ls -l $X509_USER_PROXY || echo "Warning: proxy not found"
#echo "Checking proxy info:"
#which voms-proxy-init
echo " "

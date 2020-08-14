#!/bin/bash
source #cmssetsh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/py2-xgboost/0.80-ikaegh/etc/profile.d/init.sh

hostname

cd #cmsswdir; eval `scramv1 runtime -sh`; cd #writedir;
echo "Begin script..."
root -b -q #macrofile
echo "End script!"

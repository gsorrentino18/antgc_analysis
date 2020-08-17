#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/py2-xgboost/0.80-ikaegh/etc/profile.d/init.sh

hostname

cd /local/cms/user/wadud/aNTGCmet/CMSSW_10_2_23/src/ggAnalysis/ggNtuplizer/test/; eval `scramv1 runtime -sh`; cd /data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis_skim/phoID/2020_08_13_Prompt//WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8/;
echo "Begin script..."
root -b -q /data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/batch/jobs/2020_08_13_Propmt//WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8//WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8_006.C
echo "End script!"

#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /local/cms/user/wadud/aNTGCmet/CMSSW_10_2_23/src/ggAnalysis/ggNtuplizer/test/; eval `scramv1 runtime -sh`; cd /data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis_skim/phoID/2020_08_06/promptDRPtMatching//GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8/;
echo "Begin script..."
root -b -q /data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/batch/jobs/2020_8_6/promptDRPtMatching//GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8//GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8_006.C
echo "End script!"

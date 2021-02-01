#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh

hostname

cd /local/cms/user/wadud/aNTGCmet/CMSSW_10_2_23/src/ggAnalysis/ggNtuplizer/test/; eval `scramv1 runtime -sh`; cd /data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataE_01FebV2//SinglePhotonRun2017E31Mar2018v1MINIAOD/;
echo "Begin script..."
root -b -q /data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataE_01FebV2//SinglePhotonRun2017E31Mar2018v1MINIAOD//SinglePhotonRun2017E31Mar2018v1MINIAOD_004.C
echo "End script!"

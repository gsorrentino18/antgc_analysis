#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh

hostname

cd /local/cms/user/wadud/aNTGCmet/CMSSW_10_2_23/src/ggAnalysis/ggNtuplizer/test/; eval `scramv1 runtime -sh`; cd /data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/wz_28Jan//WZTuneCP513TeVpythia8/;
echo "Begin script..."
root -b -q /data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/wz_28Jan//WZTuneCP513TeVpythia8//WZTuneCP513TeVpythia8_001.C
echo "End script!"

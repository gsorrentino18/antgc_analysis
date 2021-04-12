#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh

hostname

cd /local/cms/user/wadud/aNTGCmet/CMSSW_10_2_23/src/ggAnalysis/ggNtuplizer/test/; eval `scramv1 runtime -sh`; cd /data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_20GeV_noID//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph/;
echo "Begin script..."
root -b -q /data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_20GeV_noID//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_009.C
echo "End script!"

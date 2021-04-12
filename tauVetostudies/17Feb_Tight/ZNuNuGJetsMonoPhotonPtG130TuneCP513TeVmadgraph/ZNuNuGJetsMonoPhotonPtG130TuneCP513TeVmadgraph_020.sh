#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh

hostname

cd /local/cms/user/wadud/aNTGCmet/CMSSW_10_2_23/src/ggAnalysis/ggNtuplizer/test/; eval `scramv1 runtime -sh`; cd /data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/17Feb_4//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph/;
echo "Begin script..."
root -b -q /data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/17Feb_4//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_020.C
echo "End script!"

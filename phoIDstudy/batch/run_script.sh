#!/bin/bash
source #cmssetsh

hostname

cd #cmsswdir; eval `scramv1 runtime -sh`; cd #writedir;
echo "Begin script..."
root -b -q #macrofile
echo "End script!"

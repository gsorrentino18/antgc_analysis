#!/bin/bash

ahadd=/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/scripts/ahadd.py
cmsswDir=/local/cms/user/wadud/aNTGCmet/CMSSW_10_2_23/src/ggAnalysis/ggNtuplizer/test/

cd ${cmsswDir}; eval `scramv1 runtime -sh`; cd -;

function mergeFilesInDir(){
	_dir=$1
	_outDir=$2
	echo "Merging root files in " ${_dir}
	_outFile=${_outDir}/$(basename -- "$_dir")".root"
	if [ $( find ${_dir} -name "*.root" | wc -l )  == "0" ];
	then
		echo "Error! No root files in directory "${_dir}
		return
	fi

	if [ -f ${_outFile} ]; then
		echo "Outfile exists! ->"${_outFile}
		return
	fi
	python ${ahadd} -f ${_outFile} $(find ${_dir}/ -name "*.root")
}


searchPath=$1
searchPath=$(readlink -f ${searchPath}/)

outDir=$2
if [ -z "$outDir" ]; then
	outDir=${searchPath}/merged/
	# rm -rf ${outDir}
fi
outDir=$(readlink -f ${outDir}/)


mkdir -p ${outDir}
for directory in $(find "${searchPath}" -maxdepth 1 -mindepth 1 -not -path "*merged*" -type d);
do
	mergeFilesInDir ${directory} ${outDir}
done



# if [ ! -f ${outDir}/SinglePhotonRun2017.root ]; then
# 	python ${ahadd} -f ${outDir}/SinglePhotonRun2017.root $(find ${outDir} -name "*SinglePhotonRun2017*Nov2017v1MINIAOD.root" -not -name "*merged*" -printf "%p ")
# fi
# python ${ahadd} ${outDir}/GJets.root ${outDir}/GJetsHT*TuneCP513TeVmadgraphMLMpythia8.root
# python ${ahadd} ${outDir}/aNTGC_0p0003_0p0000004_0p_0p.root ${outDir}/aNTGC0p00030p00000040p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p0003_0p_0p_0p.root ${outDir}/aNTGC0p00030p0p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p0005_0p0000008_0p_0p.root ${outDir}/aNTGC0p00050p00000080p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p0005_0p_0p_0p.root ${outDir}/aNTGC0p00050p0p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p0008_0p000001_0p_0p.root ${outDir}/aNTGC0p00080p0000010p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p0008_0p_0p_0p.root ${outDir}/aNTGC0p00080p0p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p0010_0p000005_0p_0p.root ${outDir}/aNTGC0p00100p0000050p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p0010_0p_0p_0p.root ${outDir}/aNTGC0p00100p0p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p_0p0000004_0p_0p.root ${outDir}/aNTGC0p0p00000040p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p_0p0000008_0p_0p.root ${outDir}/aNTGC0p0p00000080p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p_0p000001_0p_0p.root ${outDir}/aNTGC0p0p0000010p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p_0p000005_0p_0p.root ${outDir}/aNTGC0p0p0000050p0p*.root
# python ${ahadd} ${outDir}/aNTGC_0p_0p_0p_0p.root ${outDir}/aNTGC0p0p0p0p*.root


# rm -f ${outDir}/SinglePhotonRun2017*Nov2017v1MINIAOD.root
# rm -f ${outDir}/aNTGC*200500.root
# rm -f ${outDir}/aNTGC*5001200.root
# rm -f ${outDir}/GJetsHT*TuneCP513TeVmadgraphMLMpythia8.root
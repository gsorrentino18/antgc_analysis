#!/bin/bash

### syntax:  bash csub.sh [path_to_ntuple_lists]  [write_directory]
set +H
cmsswfile="/cvmfs/cms.cern.ch/cmsset_default.sh"
cmsswDir=/local/cms/user/wadud/aNTGCmet/CMSSW_10_2_23/src/ggAnalysis/ggNtuplizer/test/
workDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"/

#sampleSheet=/local/cms/user/gsorrent/antgc_analysis/data/METv5Ntuples_MC_AND_DATASETS.csv
sampleSheet=/local/cms/user/gsorrent/antgc_analysis/tauVetostudies/data/METv5Ntuples_MC_AND_DATASETS.csv
jobListFile=/local/cms/user/gsorrent/antgc_analysis/tauVetostudies/batch/jobList.txt
# jobListFile=$(readlink -m ${jobListFile})
writeDir=$1
writeDir=$(readlink -m ${writeDir})/
jobsDir=$1
jobsDir=$(readlink -m ${jobsDir})/

# espresso     = 20 minutes
# microcentury = 1 hour
# longlunch    = 2 hours
# workday      = 8 hours
# tomorrow     = 1 day
# testmatch    = 3 days
# nextweek     = 1 week
jobflavor=workday
splitfiles=2


className="tauVeto"
ccfilepath="/local/cms/user/gsorrent/antgc_analysis/tauVetostudies/tauVetoV3.cc"

macroTemplate=${workDir}/macroTemplateV2.C
runScriptTemplate=${workDir}/run_script.sh
condorCFGtemplate=${workDir}/condor_job.sh

bdtPath="/hdfs/cms/user/wadud/anTGC/BDTdata/optimizedV0/aNTGC_photon_BDT.model"
dataPUfile="/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root"
CHEffArea="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt"
WCHEffArea="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfWorstChargedHadrons_70percentBased.txt"
PhoEffArea="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt"
NHEffArea="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt"
############################################################################

current_date_time=$(date +%Y-%m-%d_%H-%M-%S)

echo $current_date_time;
echo -e '\nWork directory = '${workDir}''
echo 'Write directory = '${writeDir}''
echo 'Job directory = '${jobsDir}''
echo 'Sample sheet = '${sampleSheet}''


function preSelectDataset(){

	#### grab inputs ####
	fileListPath=$(echo $1 | tr -d '\040\011\012\015')
	jobName=$(echo $2 | tr -d '\040\011\012\015')
	jobDir=$(echo $3 | tr -d '\040\011\012\015')
	writeOutDir=$(echo $4 | tr -d '\040\011\012\015')
	xSec=$(echo $5 | tr -d '\040\011\012\015')
	mcPUdist=$(echo $6 | tr -d '\040\011\012\015')
	ccfilecopy=$7


	#### check that file list exists ####
	if [ ! -f ${fileListPath} ]; then
		echo "\tError! Input file list not found! Offending file: " ${fileListPath}
		exit
	fi

	logDir=${jobDir}/log/
	
	if [ ! -d "${logDir}" ]; then
		echo -e "\tError! Job directory has not been created! "${jobDir}
		exit
	fi

	if [ ! -d "${writeOutDir}" ]; then
		echo -e "\tError! Write directory has not been created! "${jobDir}
		exit
	fi

	#### Show inputs ####
	echo	-e		"\tCreating job for " ${fileListPath}
	echo 	-e		"\t\tJob name = "${jobName}
	echo 	-e		"\t\tMC pileup file = "${mcPUdist}
	echo 	-e		"\t\tData pileup file = "${dataPUfile}
	echo 	-e		"\t\tCross section = "${xSec}
	echo 	-e		"\t\tJob directory = "${jobDir}
	echo 	-e		"\t\tOuput directory = "${writeOutDir}
	echo 	-e		"\t\tLog directory = "${logDir}


	nFiles=$(sed -n '=' ${fileListPath} | wc -l)
	echo	-e	"\t\t# of files = " ${nFiles}

	#### output root file ###
	outFile=${writeOutDir}/${jobName}.root
	echo	-e	"\t\tOutput file = "${outFile}

	### prepare root macro ###
	rootMacro=${jobDir}/${jobName}.C
	cp ${macroTemplate} ${rootMacro}
	
	sed -i 's|#className|'${className}'|g' ${rootMacro}
	sed -i 's|#macroname|'${jobName}'|g' ${rootMacro}
	sed -i 's|#fileList|'${fileListPath}'|g' ${rootMacro}
	sed -i 's|#outfilepath|'${outFile}'|g' ${rootMacro}
	sed -i 's|#ccfilepath|'${ccfilecopy}'|g' ${rootMacro}
	sed -i 's|#mcPU|'${mcPUdist}'|g' ${rootMacro}
	sed -i 's|#dataPU|'${dataPUfile}'|g' ${rootMacro}
	sed -i 's|#xSec|'${xSec}'|g' ${rootMacro}
	sed -i 's|#CHEffArea|'${CHEffArea}'|g' ${rootMacro}
	sed -i 's|#WCHEffArea|'${WCHEffArea}'|g' ${rootMacro}
	sed -i 's|#PhoEffArea|'${PhoEffArea}'|g' ${rootMacro}
	sed -i 's|#NHEffArea|'${NHEffArea}'|g' ${rootMacro}
	sed -i 's|#BDTpath|'${bdtPath}'|g' ${rootMacro}

	### prepare run script ###
	runScript=${jobDir}/${jobName}.sh
	cp ${runScriptTemplate} ${runScript}
	sed -i 's|#cmsswdir|'${cmsswDir}'|g' ${runScript}
	sed -i 's|#macrofile|'${rootMacro}'|g' ${runScript}
	sed -i 's|#cmssetsh|'${cmsswfile}'|g' ${runScript}
	sed -i 's|#writedir|'${writeOutDir}'|g' ${runScript}
	chmod +x ${runScript}

	### prepare condor script ###
	condorCFG=${jobDir}/condor_${jobName}.sh
	cp ${condorCFGtemplate} ${condorCFG}
	sed -i 's|#script|'${runScript}'|g' ${condorCFG}
	sed -i 's|#logDir|'${logDir}'|g' ${condorCFG}
	sed -i 's|#jobName|'${jobName}'|g' ${condorCFG}
	sed -i 's|#jobflavour|'${jobflavor}'|g' ${condorCFG}

	chmod +x ${condorCFG}
	
	cd ${jobDir}
	condor_submit ${condorCFG}
	cd ${workDir}
	echo -e "\n"
}


[ ! -f $sampleSheet ] && { echo "$sampleSheet : Sample sheet not found"; exit 99; }
[ ! -f $jobListFile ] && { echo "$jobListFile : Job list not found"; exit 99; }


readarray -t jobList < $jobListFile


{
	read
	while IFS=, read -r shortName dataset xSec xSecUnc singleJobFileList mcPUfile Nevents SumW SumW2 Neff lumi
	do
		if [[ ! " ${jobList[@]} " =~ " ${shortName} " ]]; then
                        #echo " ${jobList[@]} "
                        #echo " ${shortName} "
                        continue
		fi

		echo -e "\n***********************************************************************************************************************************************"
		echo -e "Preparing job:\n\tDataset = ${dataset}\n\tNtuple List = ${singleJobFileList}\n\txSec = ${xSec}\n\tmc pileup file = ${mcPUfile}"

		jobBaseName=$(basename "${singleJobFileList}")
		jobBaseName="${jobBaseName%.*}"
		# jobBaseName=$(echo ${jobBaseName} | tr -cd [:alnum:])

		jobDir=${jobsDir}/${jobBaseName}/
		jobOutDir=${writeDir}/${jobBaseName}/

		if [ -d "${jobDir}" ]; then
			echo -e "\tError! Job directory already exists "${jobDir} "Skipping!"
			continue
		fi
		
		## create job directory
		mkdir -p "${jobDir}/log/"
		mkdir -p ${jobOutDir}

		nFiles=$(sed -n '=' ${singleJobFileList} | wc -l)
		echo	-e	"\t# of files in base job = " ${nFiles} "\n\n"

		split -d -a 3 -l ${splitfiles} ${singleJobFileList} "${jobDir}/${jobBaseName}_"

		for subJobList in $(find "${jobDir}" -name "${jobBaseName}_*");
		do

			jobName=$(basename ${subJobList})
			jobName="${jobName%.*}"

			ccfilecopy=${writeDir}/$(basename $ccfilepath)

			if [ -f $ccfilecopy ]; then
				echo -e "\tUsing existing cc file: "$ccfilecopy
			else
				echo -e "\tCopied CC file "$ccfilepath "to" $ccfilecopy
				cp $ccfilepath $ccfilecopy
			fi

			preSelectDataset ${subJobList} ${jobName} ${jobDir} ${jobOutDir} ${xSec} ${mcPUfile} ${ccfilecopy}

		done
		
		echo -e "***********************************************************************************************************************************************\n"

	done
} < ${sampleSheet}


echo "Submission complete!"
current_date_time=$(date +%Y-%m-%d_%H-%M-%S)
echo $current_date_time;

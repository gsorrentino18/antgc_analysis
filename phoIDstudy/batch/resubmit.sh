#!/bin/bash

searchDir=$1

if [ -z "$searchDir" ]; then
	exit
fi

FILES=($(find ${searchDir} -name "*.err" 2>/dev/null)) 

declare finished_jobs
declare running_jobs

echo -e "Resubmitting Failed Tasks: "
FILES=($(find ${searchDir} -name "*.log" 2>/dev/null)) 

function resubmit(){
	file=$1
	fileName=$(basename $file)
	pathName=$(dirname $file)
	pathName=${pathName%log}
	condorFile=${pathName}/condor_${fileName}
	condorFile=${condorFile%.*}.sh
	

	if [ ! -f $condorFile ]; then
		return 0
	fi

	# echo $condorFile
	# condor_submit $condorFile
	jobFiles=${pathName}/${fileName%.*}
	echo $jobFiles
}

for file in "${FILES[@]}"; do
	file=$(readlink -e $file)
	fileName=$(basename $file)

	# if [[ $fileName == *"QCDHT100to200"* ]]; then
	# 	continue
	# fi
	
	normTermLine=$(tac $file | grep -m 1 "Normal termination")
	killLine=$(tac $file | grep -m 1 "Job was aborted by the user")
	updateLine=$(tac $file | grep -m 1 "Image size of job updated")
	submittedN=$(tac $file | grep -n -m 1 "Job submitted from host" | cut -d : -f 1)
	hostLine=$(tac $file | grep -m 1 "Job executing on host")	
	hostLineN=$(tac $file | grep -n -m 1 "Job executing on host" | cut -d : -f 1)

	jobID=${hostLine#*\(}
	jobID=${jobID%\)*}
	jobID=${jobID%00.000*}	

	hostName=${hostLine#*<}
	hostName=${hostName%:*}

	# if ( [ -z "$updateLine" ] ); then
	# 	echo -e "No update: " $file
	# 	continue
	# fi
	

	updateN=$(tac $file | grep -n -m 1 "Image size of job updated" | cut -d : -f 1)
	
	if grep -q "Normal termination (return value 0)" "$file"; then
		finished_jobs+="\t"
		finished_jobs+=$fileName
		finished_jobs+="\n"
		continue;	
	fi

	## skip if idle
	# if ( [ "$submittedN" -lt "$updateN" ] && [ -z "$killLine" ]); then
	# 	continue
	# fi

	## kill if running on scorpion 6
	# if ( [ "$hostLine" == *"128.101.221.176"* ] && [ "$hostLineN" -lt "$submittedN" ] ); then
	# if ( [ "$hostLine" == *"128.101.221.178"* ] && [ "$hostLineN" -lt "$submittedN" ] ); then #scorpion8
	#scorpion8
	# if [[ "$hostLine" == *"128.101.221.178"* ]]; then
	# 	echo "Killing "$jobID " running on scorpion8: " $(basename $file) $hostLineN
	# 	condor_rm $jobID
	# fi

	# if [[ "$hostLine" == *"128.101.221.176"* ]]; then
	# 	echo "Killing "$jobID " running on scorpion6: " $(basename $file) $hostLineN
	# 	condor_rm $jobID
	# fi

	# if [[ "$hostLine" == *"128.101.221.176"* ]]; then
	# 	echo "Killing "$jobID " running on scorpion6: " $(basename $file) $hostLineN
	# 	condor_rm $jobID
	# fi


	# if ( [ -z "$normTermLine" ] ); then
		
	# 	continue
	# fi


	if ( [ ! -z "$normTermLine" ] && [ ! -z "$killLine" ] ); then
		normTermLineN=$(tac $file | grep -n -m 1 "Normal termination" | cut -d : -f 1)
		killLineN=$(tac $file | grep -n -m 1 "Job was aborted by the user" | cut -d : -f 1)
		
		if ( ( [ "$updateN" -lt "$killLineN" ]  && [ "$updateN" -lt "$normTermLineN" ] ) || ( [ "$submittedN" -lt "$killLineN" ] && [ "$submittedN" -lt "$normTermLineN" ] ) );then
			continue
		fi

		# echo 1: $updateN $normTermLineN $killLineN

		resubmit $file
	elif ( [ ! -z "$normTermLine" ] ); then
		normTermLineN=$(tac $file | grep -n -m 1 "Normal termination" | cut -d : -f 1)

		if ( [ "$updateN" -lt "$normTermLineN" ] || [ "$submittedN" -lt "$normTermLineN" ] );then
			continue
		fi

		# echo 2: $submittedN $updateN $normTermLineN	
		resubmit $file
	elif ( [ ! -z "$killLine" ] ); then
		killLineN=$(tac $file | grep -n -m 1 "Job was aborted by the user" | cut -d : -f 1)
		updateN=$(tac $file | grep -n -m 1 "Image size of job updated" | cut -d : -f 1)

		if ( [ "$updateN" -lt "$killLineN" ] || [ "$submittedN" -lt "$killLineN" ] );then
			continue
		fi

		# echo 3: $submittedN $updateN $killLineN
		resubmit $file
	fi
done


# nRunning=$(echo -e "$running_jobs" | wc -l)
# echo -e "\n\n\n\nRunning jobs:" $nRunning
# echo -e "$running_jobs" | sort

nFinished=$(echo -e "$finished_jobs" | wc -l)
echo -e "\n\n\n\nFinished jobs:" $nFinished
echo -e "$finished_jobs" | sort
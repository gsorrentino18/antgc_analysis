#!/bin/bash

workDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"/

function getXsec(){
	_logFile=$1
	_SampleName=$2

	xsec=$(grep "After filter: final cross section =" ${_logFile})

	if [ -z "$xsec" ]
	then
		continue
	fi

	xsec=$(echo $xsec | sed "s/After filter: final cross section =//g")
	xsec=$(echo $xsec | sed "s/+-/,/g")
	xsec=$(echo $xsec | sed "s/ //g")
	xsec=$(echo $xsec | sed "s/nb/,nb/g")
	xsec=$(echo $xsec | sed "s/pb/,pb/g")
	xsec=$(echo $xsec | sed "s/fb/,fb/g")

	echo ${_SampleName},${miniaodFile},${xsec}
}



searchPath=$1
writeFile=$2

for jobDir in $(find "${searchPath}" -maxdepth 1 -mindepth 1 -type d); do	
	_sampleName=$(basename -- "$jobDir")
	_sampleName=$(echo ${_sampleName} | tr -cd [:alnum:])
	
	for logFile in $(find $jobDir -mindepth 1 -type f -name "*.err"); do
		getXsec ${logFile} ${_sampleName} | tee -a ${writeFile}
	done
done
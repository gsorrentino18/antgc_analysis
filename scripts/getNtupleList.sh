#!/bin/bash

searchPath=$(readlink -e $1)
writeDir=$(readlink -e $2)

function listNtuples(){
	SearchPath=$1
	outfile=$(basename ${SearchPath})
	outfile=${writeDir}/${outfile}.txt
	echo $outfile
	find ${SearchPath} -name "*.root" -type f | tee ${outfile}
	echo -e "\n\n\n"
}

mkdir -p $writeDir

for directory in $(find "${searchPath}" -maxdepth 1 -mindepth 1 -type d -not -path '*/\.*'); do
	listNtuples ${directory}
done

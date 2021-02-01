#!/bin/bash

searchPath=$1
searchPath=$(readlink -f ${searchPath}/)


for directory in $(find "${searchPath}" -maxdepth 1 -mindepth 1 -not -name "merged*" -type d | sort);
do
	echo "Deleting "$directory
	rm -rf $directory
done
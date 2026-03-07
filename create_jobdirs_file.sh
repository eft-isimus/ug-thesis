#!/bin/bash

path="$1"
cd "$path" 
find . -type f -name "job.pbs" -printf "%h\n" | sort > jobdirs.txt

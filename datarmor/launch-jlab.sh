#!/bin/bash
set -e

# Usage:
#   bash
#   source activate pangeo
#   ./launch-jlab.sh 

source activate pangeon

echo "Launching job"
s=`qsub jlab.pbs`
sjob=${s%.*}
echo ${s}

#qstat -f ${sjob}


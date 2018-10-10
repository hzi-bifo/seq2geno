#!/bin/bash
python -V
echo '==in set_env=='
ENV=$1
CONFIG=$2
SEQ2GENO_HOME=$3
echo $ENV
echo $CONFIG
echo $SEQ2GENO_HOME
source activate $ENV

python door.internal.py --env $ENV --config $CONFIG --home $SEQ2GENO_HOME

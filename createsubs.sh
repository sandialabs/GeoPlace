# Copyright 2016 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of
# this work by or on behalf of the U.S. Government. Export of this
# program may require a license from the United States Government.

#e.g. ./createsubs.sh 4 problem.txt output.txt
export NUM_SUBFOOTPRINTS=$1
export PROBLEM_FILE=$2
export OUTPUT_FILE=$3

echo "
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"

if [ -z "$NUM_SUBFOOTPRINTS" ]
then
  export NUM_SUBFOOTPRINTS=4
  echo "Num subfootprints  is unset, so setting it to '$NUM_SUBFOOTPRINTS'"
else
  echo "var is set to '$NUM_SUBFOOTPRINTS'"
fi

if [ -z "$PROBLEM_FILE" ]
then
  export PROBLEM_FILE="problem_instance.txt"
  echo "Problem file is unset, so setting it to '$PROBLEM_FILE'"
else
  echo "Problem file is set to '$PROBLEM_FILE'"
fi

if [ -z "$OUTPUT_FILE" ]
then
   export OUTPUT_FILE="output.txt"
   echo "Output file  is unset, so we set it to '$OUTPUT_FILE'"
else
   echo "Output file is set to '$OUTPUT_FILE'"
fi


#export SOLVER=cbc
#export SOLVER=cplex
export SOLVER="--solver=cplex"
#export SOLVER=gurobi

#export OPTIONS="--solver-options=threads=4"
#export OPTIONS="--report-timing --solver-options=\"mip_strategy_rinsheur=100\""
export OPTIONS="--report-timing"
#export OPTIONS="--report-timing --verbose"

#export MODEL="geo_sub_alt.py"
export MODEL="geo_sub_priority.py"

export DATA=

#export POST="--postprocess cover_priority.py"

export STREAM="--stream-output"

export OUT_FILE=$OUTPUT_FILE

#export SAVE_MODEL="--save-model=geo_sub_model.lp"
export SAVE_MODEL=

#pyomo --solver=${SOLVER} ${OPTIONS} ${MODEL} ${DATA} ${POST} ${STREAM} | tee ${OUT_FILE}
#pyomo solve ${SAVE_MODEL} --solver=${SOLVER} ${OPTIONS} ${MODEL} ${DATA} ${POST} ${STREAM} | tee ${OUT_FILE}
pyomo solve ${SAVE_MODEL} ${SOLVER} ${OPTIONS} ${POST} ${STREAM} ${MODEL} ${DATA} ${POST} ${STREAM} | tee ${OUT_FILE}

python post_process.py

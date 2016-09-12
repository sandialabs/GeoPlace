# Copyright 2016 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of
# this work by or on behalf of the U.S. Government. Export of this
# program may require a license from the United States Government.

#export SOLVER=--solver=cbc,gurobi,glpk
export SOLVER=--solver=cplex

#export OPTIONS="--solver-options=threads=4 --report-timing --solver-options=\"mip_strategy_rinsheur=100\" --verbose"
export OPTIONS="--report-timing"

export MODEL=geo_cover.py

export DATA=

#export POST=geo_cover.py

export STREAM=--stream-output

export OUT_FILE=output.txt

echo "
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"

pyomo solve ${SOLVER} ${OPTIONS} ${MODEL} ${DATA} ${POST} ${STREAM} | tee ${OUT_FILE}

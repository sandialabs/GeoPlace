#!/bin/bash

# Copyright 2016 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of
# this work by or on behalf of the U.S. Government. Export of this
# program may require a license from the United States Government.

max=9
for i in `seq 5 $max`
do
  output_name="multi_output_$i.txt"
  source createsubs.sh $i problem_instance.txt $output_name
done

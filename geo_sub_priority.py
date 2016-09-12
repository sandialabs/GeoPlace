"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import pyomo.environ as pyomo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv

import os # Imported for grabbing set environment variables
try:
	import cPickle as pickle
except ImportError:
	import pickle

###############################################################################
#####################Begin global variable initialization######################
###############################################################################
offset = [1,1]    			# Where in pixelspace the chips start. 
                                    # This is anywhere within [1,8] x [1,8]
#TODO: Generate this offset randomly, using a variable domain.
#TODO: Add weights for the pixels of interest (so that some can be more 
#      important than others.

#num_footprints = 4         # Number of footprints to be placed
num_subfoots_str = os.environ.get("NUM_SUBFOOTPRINTS",4)         # Number of footprints to be placed
if not num_subfoots_str:
   num_subfoots_str = '4'

num_subfoots = int(num_subfoots_str)

# Get the problem text file
problem_file = os.environ.get("PROBLEM_FILE","problem_simple.txt")
if not problem_file:
   problem_file = "problem_simple.txt"

B = 16                      #Bandwidth

sub_footprints = {}         # The subfootprints to be placed; assumed to be
                            # rectangles
#foot_corners = {}          # The upper-left and lower-right coordinates
                            # for each of the subfootprints. Assumed to be
                            # rectangles. In chip space. TODO: next version,
                            # as we optimize the offset, will be in pixel space
# sub_footprint = {'1', [[ll_x,ll_y],[ur_x,ur_y]], '2', [[ll_x,ll_y],[ur_x,ur_y]], ...}

chip_x_len = 12             # Number of chip columns
chip_y_len = 12             # Number of chip rows

chip_domain = {}            # The space of possible 2-dimensional 
                            # chip coordinates

chip_corners = {}           # The lower-left and upper-right coordinates
                            # for each of the chips. Assumed to be
                            # rectangles. In chip space. TODO: next version,
                            # as we optimize the offset, will be in pixel space

chips_to_cover = []

#chip_priorities = []

pixels_in_chip_x_len = 100   # Number of pixels in the x direciton of a chip

pixels_in_chip_y_len = 100   # Number of pixels in the y direciton of a chip

epsilon = 0.0001

###############################################################################

#
# Determine chip corner coordinates based on the offset
#

#def chip_corners_from_offset(off=(1,1), x_len=100, y_len=100, size=8):
#    """
#    This method populates the chip_corners dictionary based on the selected
#    offset.
#    TODO: Consider the effect of the offset on the chips
#    TODO: Properly implement (using modulus?) the wrapping of chips for
#          defining corners
#    """
#    for c in range(0,num_chips):
#        chip_corners[c][0] = off + 
#        chip_corners[c][1] = off + chip_size + 
###############################################################################

###############################################################################
def list_chips_of_interest():
    """
    This method builds the list of chips of interest.
    """
    with open(problem_file ,'rb') as file:
        reader = csv.reader(file, delimiter=',')
        for row in reader:
            if reader.line_num == 1:
                chip_x_len = int(row[0])
                chip_y_len = int(row[1])
                chip_priorities = [[0 for i in range(chip_x_len)] for j in range(chip_y_len)]
            elif reader.line_num == 2:
                max_priority = float(row[0])
            else:
                y = chip_y_len - reader.line_num + 2
                for x,value in enumerate(row):
                    if value is ' ':
                        chip_priorities[x][y] = 0.0
                    else:
                        chip_priorities[x][y] = float(value)
                        chips_to_cover.append([x,y])
    return chip_priorities

            
###############################################################################

#
# Create the Pyomo model
#

def pyomo_create_model(options=None, model_options=None):
    """
    This function gets called after the pyomo preprocessing step when 
    this file is passed to the 'pyomo solve' command. This is the core of
    the file, akin to main(). First the global variables are initialized.
    Next, the Pyomo model is setup and initialized: ConcreteModel, Sets,
    Params, Vars, Constraints and then the Objective.
    """

    import math

    chip_priorities = list_chips_of_interest()

    M = pyomo.ConcreteModel()
    M.num_subfoots = pyomo.Param(within=pyomo.PositiveIntegers,
                                 initialize=num_subfoots)

    foot_compares = {}
    for i in range(num_subfoots):
        for j in range(num_subfoots):
            if i < j:
                foot_compares[i,j] = 1
            else:
                foot_compares[i,j] = 0

    # A Pyomo Set representing the sub-footprints
    M.F = pyomo.RangeSet(0, num_subfoots-1, ordered=True)

    linear_chips_to_cover = []
    for chip in chips_to_cover:
        linear_chips_to_cover.append(chip[0]*chip_y_len + chip[1])

    linear_chip_priorities = [chip_priorities[x][y] for x in range(chip_x_len) for y in range(chip_y_len)]

    # A Pyomo Set representing the active chips
#    M.Ca = pyomo.RangeSet(0, len(chips_to_cover)-1, ordered=True)
    M.Ca = pyomo.Set(initialize=linear_chips_to_cover)
    chips = [(x,y) for x in range(chip_x_len) for y in range(chip_y_len)]
    M.C = pyomo.RangeSet(0, len(chips)-1, ordered=True)

    # Variable
    # Each subfootprint's korner coordinates
    # with the V's the way they are, it is 40% slower to enforce these as integer
#    M.Fkllx = pyomo.Var(M.F, within=pyomo.NonNegativeIntegers, bounds=(0, chip_x_len-1), initialize=0)
#    M.Fklly = pyomo.Var(M.F, within=pyomo.NonNegativeIntegers, bounds=(0, chip_y_len-1), initialize=0)
#    M.Fkurx = pyomo.Var(M.F, within=pyomo.NonNegativeIntegers, bounds=(0, chip_x_len-1), initialize=chip_x_len-1)
#    M.Fkury = pyomo.Var(M.F, within=pyomo.NonNegativeIntegers, bounds=(0, chip_y_len-1), initialize=chip_y_len-1)
    M.Fkllx = pyomo.Var(M.F, within=pyomo.NonNegativeReals, bounds=(0, chip_x_len-1), initialize=0)
    M.Fklly = pyomo.Var(M.F, within=pyomo.NonNegativeReals, bounds=(0, chip_y_len-1), initialize=0)
    M.Fkurx = pyomo.Var(M.F, within=pyomo.NonNegativeReals, bounds=(0, chip_x_len-1), initialize=chip_x_len-1)
    M.Fkury = pyomo.Var(M.F, within=pyomo.NonNegativeReals, bounds=(0, chip_y_len-1), initialize=chip_y_len-1)

    # Variable
    # X1 indicates whether that chip location is covered by that footprint
    ## about 40% faster with X1 Binary, even though it increases the number of binary variables
    M.X1 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip being covered by subfootprint
#    M.X1 = pyomo.Var(M.F, M.C, within=pyomo.NonNegativeReals) # Indicator variable for chip being covered by subfootprint
    # define M.X = sum(M.X1[f]) and bounds=(0,1) or as binary - just slows down the solver 

    # Variable
    # V1..4 indicates whether that chip is (V1) to the right of the lower left corner, (V2) above ll y, 
    # (V3) left of ur x, (V4) below ur y, for that footprint
    ## Can we make these non-binary? Not easily... These appear essential, and force everyone else to be binary.
    M.V1 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip coverage (llx)
    M.V2 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip coverage (lly)
    M.V3 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip coverage (urx)
    M.V4 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip coverage (ury)

    # Variable
    # Area of each footprint. 
    # A useful value to have.
    M.A = pyomo.Var(M.F, within=pyomo.NonNegativeReals)
    # Variable
    # Maximum area over all footprints
    M.Amax = pyomo.Var(within=pyomo.NonNegativeReals)


    # Constraint
    # ll <= ur  korner coordinates
    ## I don't think this constraint is strictly necessary, because they won't cover anything if it is violated
    ## However, including this seems to speed up the solve marginally
    # For each sub-footprint, defined as a rectangle with corners
    # 0 and 1, 0.x <= 1.x AND 0.y <= 1.y
    #
    def c0x_rule(M, f):
        return M.Fkllx[f] <= M.Fkurx[f]
    M.c0x = pyomo.Constraint(M.F, rule=c0x_rule)

    def c0y_rule(M, f):
        return M.Fklly[f] <= M.Fkury[f]
    M.c0y = pyomo.Constraint(M.F, rule=c0y_rule)

   # Constraint
   ## Order the footprints by increasing lexicographic coordinates, x then y
   ## Enforces that 1.llx <= 2.llx 
   ## Tie-breaker speeds up the program about 60%
    def c1_rule(M, f):
        if f < M.num_subfoots - 1:
#            return M.Fkllx[f] <= M.Fkllx[f+1] # just x
            return M.Fkllx[f] + chip_x_len*M.Fklly[f] <= M.Fkllx[f+1] + chip_x_len*M.Fklly[f+1] + 1
        else:
            return pyomo.Constraint.Skip
    M.c1 = pyomo.Constraint(M.F, rule=c1_rule)

   
## Constraints  
## determine covered chips

    # M.X1[f,c] indicates whether or not footprint f covers active chip c
    ## M.X1 can be 1  only if  chip is within the footprints  rectangle, otherwise it must be zero (less than 1)
    ## that is, if the footprint doesn't cover the chip, then M.X1 = 0
    
    ## Each of the V's tests the sidedness of one of the corner coordinates wrt a chip coordinates

    ## chips outside the subfootprint rectangle have V = 0
    ## M.V? can be 1  only if  chip is inside the footprint interval, otherwise it must be zero (less than one)
    #Define the constraint of V1, lx
    def cvv1_rule(M, f, c):
        return (M.Fkllx[f] - (1 - M.V1[f,c])*chip_x_len <= chips[c][0])
    M.cvv1 = pyomo.Constraint(M.F, M.C, rule=cvv1_rule)

    #Define the constraint of V2, ly
    def cvv2_rule(M, f, c):
        return (M.Fklly[f] - (1 - M.V2[f,c])*chip_y_len <= chips[c][1]) 
    M.cvv2 = pyomo.Constraint(M.F, M.C, rule=cvv2_rule)

    #Define the constraint of V3, ux
    def cvv3_rule(M, f, c):
        return (M.Fkurx[f] + (1 - M.V3[f,c])*chip_x_len >= chips[c][0])
    M.cvv3 = pyomo.Constraint(M.F, M.C, rule=cvv3_rule)

    #Define the constraint of V4, uy
    def cvv4_rule(M, f, c):
        return (M.Fkury[f] + (1 - M.V4[f,c])*chip_y_len >= chips[c][1]) 
    M.cvv4 = pyomo.Constraint(M.F, M.C, rule=cvv4_rule)

    ## chips inside the subfootprint footprint rectangle have V = 1
    ## M.V? can be 0   only if  chip is outside the footprints interval, otherwise it must be one (more than zero)
    #Define the constraint of V1, lx
    def cv1_rule(M, f, c):
        return (chip_x_len+1)*M.V1[f,c] + M.Fkllx[f] >= chips[c][0] + 1 
    M.cv1 = pyomo.Constraint(M.F, M.C, rule=cv1_rule)

    #Define the constraint of V2, ly
    def cv2_rule(M, f, c):
        return (chip_y_len+1)*M.V2[f,c] + M.Fklly[f] >= chips[c][1] + 1 
    M.cv2 = pyomo.Constraint(M.F, M.C, rule=cv2_rule)

    #Define the constraint of V3, ux
    def cv3_rule(M, f, c):
        return M.Fkurx[f] - (chip_x_len+1)*M.V3[f,c] <= chips[c][0] - 1 
    M.cv3 = pyomo.Constraint(M.F, M.C, rule=cv3_rule)

    #Define the constraint of V4, uy
    def cv4_rule(M, f, c):
        return M.Fkury[f] - (chip_y_len+1)*M.V4[f,c] <= chips[c][1] - 1 
    M.cv4 = pyomo.Constraint(M.F, M.C, rule=cv4_rule)

    # This next constraint ensures active chips are within exactly one footprint
    # But we want more. We want any covered chip to be one

    # If a chip is inside the rectangle, all four v are 1, and x = 1
    ## M.X1 can be 0   only if  chip is outside the footprints rectangle, otherwise it must be one    
    ## If the sum of v is 4, then the chip is inside all intervals, i.e. inside the rectangle
    #Ensure that X1[f,c] is equal to one if V1 and V2 and V3 and V4 are all 1
    def cv_rule(M, f, c):
        return M.X1[f,c] >= M.V1[f,c] + M.V2[f,c] + M.V3[f,c] + M.V4[f,c] - 3
    M.cv = pyomo.Constraint(M.F, M.C, rule=cv_rule)

    # If a chip is outside the rectangle, at least one of the v's is 0, and x = 0
    ## X1 must be less than each of the V's, 
    ## i.e. if chip is uncovered, then X1 must be 0
    def cxv1_rule(M, f, c):
        return M.X1[f,c] <= M.V1[f,c];
    M.cxv1 = pyomo.Constraint(M.F, M.C, rule=cxv1_rule)
    def cvx2_rule(M, f, c):
        return M.X1[f,c] <= M.V2[f,c];
    M.cvx2 = pyomo.Constraint(M.F, M.C, rule=cvx2_rule)
    def cvx3_rule(M, f, c):
        return M.X1[f,c] <= M.V3[f,c];
    M.cvx3 = pyomo.Constraint(M.F, M.C, rule=cvx3_rule)
    def cvx4_rule(M, f, c):
        return M.X1[f,c] <= M.V4[f,c];
    M.cvx4 = pyomo.Constraint(M.F, M.C, rule=cvx4_rule)


    ## Given the above, M.X1 absolutely indicates whether a chip is covered (1) by a footprint, or not (0).
    
    # Constrain non-overlap by requiring each chip is in only one footprint:
    def cx_rule(M, c):
        return sum ( M.X1[f,c] for f in M.F ) <= 1
    M.cx = pyomo.Constraint(M.C, rule=cx_rule)

    # Objectives:
    #
    # Define a sub-footprint's area
    #
    def carea_rule(M, f):
        return M.A[f] == sum(M.X1[f,c] for c in M.C)
    M.carea = pyomo.Constraint(M.F, rule=carea_rule)

    #
    # Define the maximum sub-footprint area.
    #
    def camax_rule(M, f):
        return M.Amax >= M.A[f]
    M.camax = pyomo.Constraint(M.F, rule=camax_rule)

    # Constrain
    # Bandwidth is not exceeded
    def cbwidth_rule(M):
        return num_subfoots*M.Amax <= B
    M.cbwidth = pyomo.Constraint(rule=cbwidth_rule)

    # Objective is to minimize bandwidth if we are constrained to cover all of the active chips
    # Objective function: Maximize the Bandwidth over the area of the sub
    # -regions of interest (chips to cover)
    #
    # NOTE: without the X1, when using the area of the placed sub-footprints in the objective
    #       function, we get a (quadratic?) indefinite matrix.
#    M.obj = pyomo.Objective(expr=sum(M.X1[f,c] for f in M.F for c in M.C),
#                                     sense=pyomo.minimize)
#    M.obj = pyomo.Objective(expr=M.Amax, sense=pyomo.minimize)
#    M.obj = pyomo.Objective(expr=M.Amax + float(1)/(len(chips))*sum(M.A[f] for f in M.F), sense=pyomo.minimize)

    # Objective 
    # Maximize the weighted value of the covered active chips
    #   the penalty term prefers smaller footprints if it doesn't lose any active chips
    M.obj_priority = pyomo.Objective(expr=sum(M.X1[f,c]*(1-epsilon*f)*linear_chip_priorities[c] for f in M.F for c in M.C) 
- sum(M.X1[f,c] for f in M.F for c in M.C)*0.01, sense=pyomo.maximize)

    # return the Pyomo model (this is required).	
    return M

################################################################################

#
# Print the generated footprint placement
#

def pyomo_postprocess(options, instance, results):
    foots = []
    total_area = 0
    print '***********************************************************' \
          '*****************'
    print '***********************************************************' \
          '*****************'
    print ""
    print "Number of footprints placed: %d" % (num_subfoots)
    print ""
    print '***********************************************************' \
          '*****************'
    print '***********************************************************' \
          '*****************'
    print ""
    print "Footprint placements: (ll_x,ll_y), (ur_x,ur_y):"
    out_file = open("geo_sub_solution.txt","w")
    for f in instance.F:
        ll_x = instance.Fkllx[f].value
        ll_y = instance.Fklly[f].value
        ur_x = instance.Fkurx[f].value
        ur_y = instance.Fkury[f].value
        print "Footprint %d placement: (%2d, %2d)  ,   (%2d, %2d)" % (f, ll_x, ll_y, ur_x, ur_y)
        print "Footprint %d area: %d" % (f, (instance.A[f].value))
        out_str = str(ll_x) + " " + str(ll_y) + " " + str(ur_x) + " " + str(ur_y)
        out_file.write(out_str + "\n")
    total_area = sum(instance.X1[f,c] for f in instance.F for c in instance.C)

    # weighted area
    chip_priorities = list_chips_of_interest()
    linear_chip_priorities = [chip_priorities[x][y] for x in range(chip_x_len) for y in range(chip_y_len)]
    weighted_area = sum(linear_chip_priorities[c]*instance.X1[f,c] for f in instance.F for c in instance.C)

    # graphically plot the subfootprints, depending on these flags 
    dox = 1;
    dov = 0;

    chips = [(x,y) for x in range(chip_x_len) for y in range(chip_y_len)]
    print " number of subfootprint %d. " % (len( instance.F ))
    print " max area %d. " % (instance.Amax)

    for f in instance.F:

 	print " Coordinate directions: "	
        print "  y -> "
        print "x"
        print "|"
        print "v"
            
        print "Footprint %d area %f" % (f,instance.A[f])

	if dox>0: 
 	# X1
	    print "Foot %d X1: " % (f)
 	    oldx = 0;
	    for c in instance.C:
		    x = chips[c][0]
		    y = chips[c][1]
		    if x>oldx:
		  	print " "
			oldx = x
		    xx = instance.X1[f,c].value
		    #print "Foot %d chip %d %d x = %d" % (f, x, y, xx)
		    print "%d" % (xx),
	    print " "

        if dov>0:
        # V1
	    print "Foot %d V1: " % (f)
 	    oldx = 0;
	    for c in instance.C:
		x = chips[c][0]
		y = chips[c][1]
		if x>oldx:
		  	print " "
			oldx = x
		xx = instance.V1[f,c].value
		#print "Foot %d chip %d %d v1 = %d" % (f, x, y, xx)
		print "%d" % (xx),
	    print " "

        # V2
	    print "Foot %d V2: " % (f)
 	    oldx = 0;
	    for c in instance.C:
		x = chips[c][0]
		y = chips[c][1]
		if x>oldx:
		  	print " "
			oldx = x
		xx = instance.V2[f,c].value
		#print "Foot %d chip %d %d v1 = %d" % (f, x, y, xx)
		print "%d" % (xx),
	    print " "

        # V3
	    print "Foot %d V3: " % (f)
 	    oldx = 0;
	    for c in instance.C:
		x = chips[c][0]
		y = chips[c][1]
		if x>oldx:
		  	print " "
			oldx = x
		xx = instance.V3[f,c].value
		#print "Foot %d chip %d %d v1 = %d" % (f, x, y, xx)
		print "%d" % (xx),
	    print " "

        # V4
	    print "Foot %d V4: " % (f)
 	    oldx = 0;
	    for c in instance.C:
		x = chips[c][0]
		y = chips[c][1]
		if x>oldx:
		  	print " "
			oldx = x
		xx = instance.V4[f,c].value
		#print "Foot %d chip %d %d v1 = %d" % (f, x, y, xx)
		print "%d" % (xx),
	    print " "

    out_file.close()
    print "Weighted area of placed sub-footprints (obj):  %d" % (weighted_area)
    print "The objective function value is:               %d" % (instance.obj_priority)
    print "Number (area) of chips in footprints:          %d" % (total_area)
    print "Number (area) of chips wanting to be covered:  %d" % (len(chips_to_cover))
    print ""
    print '***********************************************************' \
          '*****************'
    print '***********************************************************' \
          '*****************'
################################################################################

def main():
	print("This script is not meant to be run at the command line." \
		  "Please run createsubs.sh instead.")
	return 0

################################################################################

if __name__=='__main__':
  exit(main())


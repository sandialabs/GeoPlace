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
num_subfoots = 4         # Number of footprints to be placed

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

pixel_x_len = 100           # Number of pixels in the x direciton

pixel_y_len = 100           # Number of pixels in the y direciton

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
    #"""
    chips_to_cover.append([0,0])
    chips_to_cover.append([1,0])
    chips_to_cover.append([0,1])
    chips_to_cover.append([1,1])
    chips_to_cover.append([0,2])
    chips_to_cover.append([1,2])
    chips_to_cover.append([0,3])
    chips_to_cover.append([1,3])
    chips_to_cover.append([0,4])
    chips_to_cover.append([1,4])
    chips_to_cover.append([0,5])
    chips_to_cover.append([1,5])
    chips_to_cover.append([0,6])
    chips_to_cover.append([1,6])
    chips_to_cover.append([0,7])
    chips_to_cover.append([1,7])
    chips_to_cover.append([0,8])
    chips_to_cover.append([1,8])


    chips_to_cover.append([3,8])
    chips_to_cover.append([4,8])
    chips_to_cover.append([4,7])
    chips_to_cover.append([5,7])

    chips_to_cover.append([8,3])
    chips_to_cover.append([8,4])
    chips_to_cover.append([8,5])
    chips_to_cover.append([9,3])
    chips_to_cover.append([9,4])
    chips_to_cover.append([9,5])
    #"""
    #chips_to_cover.append([0,0])    
    #chips_to_cover.append([4,4])    
    #chips_to_cover.append([8,8])    
    #chips_to_cover.append([3,7])    
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

    list_chips_of_interest()
    print chips_to_cover
    print len(chips_to_cover)

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

    # A Pyomo Set representing the active chips
#    M.Ca = pyomo.RangeSet(0, len(chips_to_cover)-1, ordered=True)
    M.Ca = pyomo.Set(initialize=linear_chips_to_cover)
    chips = [(x,y) for x in range(chip_x_len) for y in range(chip_y_len)]
    M.C = pyomo.RangeSet(0, len(chips)-1, ordered=True)

    # A Pyomo Variable representing each placed footprint's korner.
    M.Fkllx = pyomo.Var(M.F, within=pyomo.NonNegativeIntegers, bounds=(0, chip_x_len-1), initialize=0)
    M.Fklly = pyomo.Var(M.F, within=pyomo.NonNegativeIntegers, bounds=(0, chip_y_len-1), initialize=0)
    M.Fkurx = pyomo.Var(M.F, within=pyomo.NonNegativeIntegers, bounds=(0, chip_x_len-1), initialize=chip_x_len-1)
    M.Fkury = pyomo.Var(M.F, within=pyomo.NonNegativeIntegers, bounds=(0, chip_y_len-1), initialize=chip_y_len-1)

    M.covered = pyomo.Var(pyomo.RangeSet(0,chip_x_len-1), pyomo.RangeSet(0,chip_y_len-1))

    M.Yff = pyomo.Set(initialize=foot_compares.keys())
    M.Y1 = pyomo.Var(M.Yff, within=pyomo.Binary) # Indicator variable for disjunctive constraints
    M.Y2 = pyomo.Var(M.Yff, within=pyomo.Binary) # Indicator variable for disjunctive constraints
    M.Y3 = pyomo.Var(M.Yff, within=pyomo.Binary) # Indicator variable for disjunctive constraints

    M.X1 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for disjunctive constraints

    M.V1 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip coverage (llx)
    M.V2 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip coverage (lly)
    M.V3 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip coverage (urx)
    M.V4 = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for chip coverage (ury)
    M.V = pyomo.Var(M.F, M.C, within=pyomo.Binary) # Indicator variable for overall chip coverage

    # A Pyomo convenience variable representing placed footprint
    # areas.

    M.Fax = pyomo.Var(M.F, within=pyomo.PositiveIntegers, initialize=chip_x_len)
    M.Fay = pyomo.Var(M.F, within=pyomo.PositiveIntegers, initialize=chip_y_len)

    # Area of each footprint. A useful value to have.
    M.A = pyomo.Var(M.F, within=pyomo.PositiveIntegers)
    # Variable representing the maximum footprint's area.
    M.Amax = pyomo.Var(within=pyomo.PositiveIntegers)

    #
    # For each sub-footprint, defined as a rectangle with corners
    # 0 and 1, 0.x <= 1.x AND 0.y <= 1.y
    #
    def c0x_rule(M, f):
        return M.Fkllx[f] <= M.Fkurx[f] 
    M.c0x = pyomo.Constraint(M.F, rule=c0x_rule)

    def c0y_rule(M, f):
        return M.Fklly[f] <= M.Fkury[f] 
    M.c0y = pyomo.Constraint(M.F, rule=c0y_rule)

    #
    # This constraint and the subsequent constraint are to ensure that
    # sub-footprints (rectangles) cannot overlap.
    #
    def c1_rule(M, f):
        if f < M.num_subfoots - 1:
            return M.Fkllx[f] <= M.Fkllx[f+1]
        else:
            return pyomo.Constraint.Skip
    M.c1 = pyomo.Constraint(M.F, rule=c1_rule)

    #
    # disjunctive constraints over every pair
    # ((M.Fkurx[f1] < M.Fkllx[f2]) or (M.Fklly[f1] > M.Fkury[f2]) or (M.Fkury[f1] < M.Fklly[f2]))
    #

    def c2a_rule(M, f1, f2):
        if f1 < f2:
            return (M.Fkllx[f2] - M.Fkurx[f1] + chip_x_len*(1 - M.Y1[f1,f2]) >= 1)
        else:
            return pyomo.Constraint.Skip
    M.c2a = pyomo.Constraint(M.F, M.F, rule=c2a_rule)

    def c2b_rule(M, f1, f2):
        if f1 < f2:
            return (M.Fklly[f1] - M.Fkury[f2] + chip_y_len*(1 - M.Y2[f1,f2]) >= 1)
        else:
            return pyomo.Constraint.Skip
    M.c2b = pyomo.Constraint(M.F, M.F, rule=c2b_rule)

    def c2c_rule(M, f1, f2):
        if f1 < f2:
            return (M.Fklly[f2] - M.Fkury[f1] + chip_y_len*(1 - M.Y3[f1,f2]) >= 1)
        else:
            return pyomo.Constraint.Skip
    M.c2c = pyomo.Constraint(M.F, M.F, rule=c2c_rule)

    def c2d_rule(M, f1, f2):
        if f1 < f2:
            return M.Y1[f1,f2] + M.Y2[f1,f2] + M.Y3[f1,f2] >= 1
        else:
            return pyomo.Constraint.Skip
    M.c2d = pyomo.Constraint(M.F, M.F, rule=c2d_rule)

    #
    # Each of the chips containing pixels to be covered need to be covered by
    # a subfootprint
    # ((M.Fkllx[f] <= M.Ck2cx[c]) and (M.Fkurx[f] >= M.Ck2cx[c]) and (M.Fklly[f] <= M.Ck2cy[c]) and (M.Fkury[f] >= M.Ck2cy[c])):
    #

    # M.X1[f,c] indicates whether or not footprint f covers active chip c
    def c3a_rule(M, f, c):
        return (M.Fkllx[f] - (1 - M.X1[f,c])*chip_x_len <= chips[c][0])
#        return (M.Fkllx[f] + (1 - M.X1[f,c])*chip_x_len >= chips[c][0])
    M.c3a = pyomo.Constraint(M.F, M.C, rule=c3a_rule)

    def c3b_rule(M, f, c):
        return (M.Fkurx[f] + (1 - M.X1[f,c])*chip_x_len >= chips[c][0])
    M.c3b = pyomo.Constraint(M.F, M.C, rule=c3b_rule)

    def c3c_rule(M, f, c):
        return (M.Fklly[f] - (1 - M.X1[f,c])*chip_y_len <= chips[c][1]) 
    M.c3c = pyomo.Constraint(M.F, M.C, rule=c3c_rule)

    def c3d_rule(M, f, c):
        return (M.Fkury[f] + (1 - M.X1[f,c])*chip_y_len >= chips[c][1]) 
    M.c3d = pyomo.Constraint(M.F, M.C, rule=c3d_rule)

    # This next constraint ensures active chips are within exactly one footprint
    # But we want more. We want any covered chip to be one
    def c3j_rule(M, c):
        return sum(M.X1[f,c] for f in M.F) == 1
    M.c3j = pyomo.Constraint(M.Ca, rule=c3j_rule)


    #Define the constraint of V1
    def cv1_rule(M, f, c):
        return (chip_x_len+1)*M.V1[f,c] + M.Fkllx[f] >= chips[c][0] + 1 
    M.cv1 = pyomo.Constraint(M.F, M.C, rule=cv1_rule)

    #Define the constraint of V2
    def cv2_rule(M, f, c):
        return (chip_y_len+1)*M.V2[f,c] + M.Fklly[f] >= chips[c][1] + 1 
    M.cv2 = pyomo.Constraint(M.F, M.C, rule=cv2_rule)

    #Define the constraint of V3
    def cv3_rule(M, f, c):
        return M.Fkurx[f] - (chip_x_len+1)*M.V3[f,c] <= chips[c][0] - 1 
    M.cv3 = pyomo.Constraint(M.F, M.C, rule=cv3_rule)

    #Define the constraint of V4
    def cv4_rule(M, f, c):
        return M.Fkury[f] - (chip_y_len+1)*M.V4[f,c] <= chips[c][1] - 1 
    M.cv4 = pyomo.Constraint(M.F, M.C, rule=cv4_rule)

    #Ensure that X1[f,c] is equal to one if V1 and V2 and V3 and V4 are all 1
    def cv_rule(M, f, c):
        return M.X1[f,c] >= M.V1[f,c] + M.V2[f,c] + M.V3[f,c] + M.V4[f,c] - 3
    M.cv = pyomo.Constraint(M.F, M.C, rule=cv_rule)

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

    #
    # Define the sub-footprint area through two Pyomo.Constraints 
    #

#    def c4_rule(M, f):
#        return M.Fax[f] == (M.Fkurx[f] - M.Fkllx[f]) + 1
#    M.c4 = pyomo.Constraint(M.F, rule=c4_rule)
#
#    def c5_rule(M, f):
#        return M.Fay[f] == (M.Fkury[f] - M.Fklly[f]) + 1
#    M.c5 = pyomo.Constraint(M.F, rule=c5_rule)

    #
    # Objective function: Maximize the Bandwidth over the area of the sub
    # -regions of interest (chips to cover)
    #
    # NOTE: When using the area of the placed sub-footprints in the objective
    #       function, we get a (quadratic?) indefinite matrix.
#    M.obj = pyomo.Objective(expr=sum((M.Fax[f]+M.Fay[f])/len(chips_to_cover) for f in M.F),
#                                     sense=pyomo.minimize)
#    M.obj = pyomo.Objective(expr=sum(M.X1[f,c] for f in M.F for c in M.C),
#                                     sense=pyomo.minimize)
#    M.obj = pyomo.Objective(expr=M.Amax, sense=pyomo.minimize)

    M.obj = pyomo.Objective(expr=M.Amax + float(1)/(len(chips))*sum(M.A[f] for f in M.F), sense=pyomo.minimize)

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
#        print "Footprint %d area: %d" % (f, (instance.Fax[f].value)*(instance.Fay[f].value))
        print "Footprint %d area: %d" % (f, (instance.Fax[f].value)*(instance.Fay[f].value))
        #total_area += (instance.Fax[f].value)*(instance.Fay[f].value)
        out_str = str(ll_x) + " " + str(ll_y) + " " + str(ur_x) + " " + str(ur_y)
        out_file.write(out_str + "\n")
    total_area = sum(instance.X1[f,c] for f in instance.F for c in instance.C)
    out_file.close()
    print "Total area of placed sub-footprints :    %d" % (total_area)
    print "Number (area) of chips to be covered:    %d" % (len(chips_to_cover))
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


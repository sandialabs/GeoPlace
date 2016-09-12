"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import post_process_multi as pp
import pyomo.environ as pyomo
import numpy as np
import csv

try:
    import cPickle as pickle
except ImportError:
    import pickle

chip_x_len = 12  # Number of chip columns
chip_y_len = 12  # Number of chip rows

epsilon  = 0.001
epsilon2 = 0.00001

def list_chips_of_interest(problem_file):
    """
    This method builds the list of chips of interest.
    """
    chips_to_cover = []
    with open(problem_file, 'rb') as file:
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
                for x, value in enumerate(row):
                    if value is ' ':
                        chip_priorities[x][y] = 0.0
                    else:
                        chip_priorities[x][y] = float(value)
                        chips_to_cover.append([x, y])
    return chip_priorities, chips_to_cover


###############################################################################

#
# Create the Pyomo model
#

def pyomo_create_model(options=None, model_options=None, chip_priorities=[], chips_to_cover=[], num_subfoots=0, B=0,
                       chip_x_len=0, chip_y_len=0):
    import math
    if chip_priorities == []:
        # No data passed in, so call this
        problem_file = 'problem_instance.txt'
        chip_priorities, chips_to_cover = list_chips_of_interest(problem_file)
        B = 16
        num_subfoots = 4
        chip_x_len = 12
        chip_y_len = 12

    M = pyomo.ConcreteModel()
    M.num_subfoots = pyomo.Param(within=pyomo.PositiveIntegers,
                                 initialize=num_subfoots)

    foot_compares = {}
    for i in range(num_subfoots):
        for j in range(num_subfoots):
            if i < j:
                foot_compares[i, j] = 1
            else:
                foot_compares[i, j] = 0

    # A Pyomo Set representing the sub-footprints
    M.F = pyomo.RangeSet(0, num_subfoots - 1, ordered=True)

    linear_chips_to_cover = []
    for chip in chips_to_cover:
        linear_chips_to_cover.append(chip[0] * chip_y_len + chip[1])

    linear_chip_priorities = [chip_priorities[x][y] for x in range(chip_x_len) for y in range(chip_y_len)]

    # A Pyomo Set representing the active chips
    M.Ca = pyomo.Set(initialize=linear_chips_to_cover)
    chips = [(x, y) for x in range(chip_x_len) for y in range(chip_y_len)]
    M.C = pyomo.RangeSet(0, len(chips) - 1, ordered=True)

    # Variable
    # Each subfootprint's korner coordinates
    # with the V's the way they are, it is 40% slower to enforce these as integer

    M.Fkllx = pyomo.Var(M.F, within=pyomo.NonNegativeReals, bounds=(0, chip_x_len - 1), initialize=0)
    M.Fklly = pyomo.Var(M.F, within=pyomo.NonNegativeReals, bounds=(0, chip_y_len - 1), initialize=0)
    M.Fkurx = pyomo.Var(M.F, within=pyomo.NonNegativeReals, bounds=(0, chip_x_len - 1), initialize=chip_x_len - 1)
    M.Fkury = pyomo.Var(M.F, within=pyomo.NonNegativeReals, bounds=(0, chip_y_len - 1), initialize=chip_y_len - 1)

    # Variable
    # X1 indicates whether that chip location is covered by that footprint
    M.X1 = pyomo.Var(M.F, M.C, within=pyomo.Binary)  # Indicator variable for chip being covered by subfootprint

    # Variable
    # V1..4 indicates whether that chip is (V1) to the right of the lower left corner, (V2) above ll y,
    # (V3) left of ur x, (V4) below ur y, for that footprint
    ## Can we make these non-binary? Not easily... These appear essential, and force everyone else to be binary.
    M.V1 = pyomo.Var(M.F, M.C, within=pyomo.Binary)  # Indicator variable for chip coverage (llx)
    M.V2 = pyomo.Var(M.F, M.C, within=pyomo.Binary)  # Indicator variable for chip coverage (lly)
    M.V3 = pyomo.Var(M.F, M.C, within=pyomo.Binary)  # Indicator variable for chip coverage (urx)
    M.V4 = pyomo.Var(M.F, M.C, within=pyomo.Binary)  # Indicator variable for chip coverage (ury)

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
            return M.Fkllx[f] + chip_x_len * M.Fklly[f] <= M.Fkllx[f + 1] + chip_x_len * M.Fklly[f + 1] + 1
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
    # Define the constraint of V1, lx
    def cvv1_rule(M, f, c):
        return (M.Fkllx[f] - (1 - M.V1[f, c]) * chip_x_len <= chips[c][0])

    M.cvv1 = pyomo.Constraint(M.F, M.C, rule=cvv1_rule)

    # Define the constraint of V2, ly
    def cvv2_rule(M, f, c):
        return (M.Fklly[f] - (1 - M.V2[f, c]) * chip_y_len <= chips[c][1])

    M.cvv2 = pyomo.Constraint(M.F, M.C, rule=cvv2_rule)

    # Define the constraint of V3, ux
    def cvv3_rule(M, f, c):
        return (M.Fkurx[f] + (1 - M.V3[f, c]) * chip_x_len >= chips[c][0])

    M.cvv3 = pyomo.Constraint(M.F, M.C, rule=cvv3_rule)

    # Define the constraint of V4, uy
    def cvv4_rule(M, f, c):
        return (M.Fkury[f] + (1 - M.V4[f, c]) * chip_y_len >= chips[c][1])

    M.cvv4 = pyomo.Constraint(M.F, M.C, rule=cvv4_rule)

    ## chips inside the subfootprint footprint rectangle have V = 1
    ## M.V? can be 0   only if  chip is outside the footprints interval, otherwise it must be one (more than zero)
    # Define the constraint of V1, lx
    def cv1_rule(M, f, c):
        return (chip_x_len + 1) * M.V1[f, c] + M.Fkllx[f] >= chips[c][0] + 1

    M.cv1 = pyomo.Constraint(M.F, M.C, rule=cv1_rule)

    # Define the constraint of V2, ly
    def cv2_rule(M, f, c):
        return (chip_y_len + 1) * M.V2[f, c] + M.Fklly[f] >= chips[c][1] + 1

    M.cv2 = pyomo.Constraint(M.F, M.C, rule=cv2_rule)

    # Define the constraint of V3, ux
    def cv3_rule(M, f, c):
        return M.Fkurx[f] - (chip_x_len + 1) * M.V3[f, c] <= chips[c][0] - 1

    M.cv3 = pyomo.Constraint(M.F, M.C, rule=cv3_rule)

    # Define the constraint of V4, uy
    def cv4_rule(M, f, c):
        return M.Fkury[f] - (chip_y_len + 1) * M.V4[f, c] <= chips[c][1] - 1

    M.cv4 = pyomo.Constraint(M.F, M.C, rule=cv4_rule)

    # This next constraint ensures active chips are within exactly one footprint
    # But we want more. We want any covered chip to be one

    # If a chip is inside the rectangle, all four v are 1, and x = 1
    ## M.X1 can be 0   only if  chip is outside the footprints rectangle, otherwise it must be one
    ## If the sum of v is 4, then the chip is inside all intervals, i.e. inside the rectangle
    # Ensure that X1[f,c] is equal to one if V1 and V2 and V3 and V4 are all 1
    def cv_rule(M, f, c):
        return M.X1[f, c] >= M.V1[f, c] + M.V2[f, c] + M.V3[f, c] + M.V4[f, c] - 3

    M.cv = pyomo.Constraint(M.F, M.C, rule=cv_rule)

    # If a chip is outside the rectangle, at least one of the v's is 0, and x = 0
    ## X1 must be less than each of the V's,
    ## i.e. if chip is uncovered, then X1 must be 0
    def cxv1_rule(M, f, c):
        return M.X1[f, c] <= M.V1[f, c]

    M.cxv1 = pyomo.Constraint(M.F, M.C, rule=cxv1_rule)

    def cvx2_rule(M, f, c):
        return M.X1[f, c] <= M.V2[f, c]

    M.cvx2 = pyomo.Constraint(M.F, M.C, rule=cvx2_rule)

    def cvx3_rule(M, f, c):
        return M.X1[f, c] <= M.V3[f, c]

    M.cvx3 = pyomo.Constraint(M.F, M.C, rule=cvx3_rule)

    def cvx4_rule(M, f, c):
        return M.X1[f, c] <= M.V4[f, c]

    M.cvx4 = pyomo.Constraint(M.F, M.C, rule=cvx4_rule)

    ## Given the above, M.X1 absolutely indicates whether a chip is covered (1) by a footprint, or not (0).

    # Constrain non-overlap by requiring each chip is in only one footprint:
    def cx_rule(M, c):
        return sum(M.X1[f, c] for f in M.F) <= 1

    M.cx = pyomo.Constraint(M.C, rule=cx_rule)

    # Objectives:
    #
    # Define a sub-footprint's area
    #
    def carea_rule(M, f):
#        return M.A[f] >= sum(M.X1[f, c] for c in M.C)
    # samitch-change above doesn't seem to work
        return M.A[f] == sum(M.X1[f, c] for c in M.C)

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
        return num_subfoots * M.Amax <= B

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
    # Note on epsilon regularization term:
    #   We get more value for covering a chip with a smaller-indexed subfootprint.
    #   We don't think this results in better quality solutions, but it improves solution time 
    #   because it reduces symmetry in the objective

##    M.obj_priority = pyomo.Objective(
##        expr=sum(M.X1[f, c] * (1 - epsilon * f) * linear_chip_priorities[c] for f in M.F for c in M.C)
##             - sum(M.X1[f, c] for f in M.F for c in M.C) * 0.01, sense=pyomo.maximize)

    # TODO: Consider alternative regularization techniques, ones that actually capture something we care about.

    # Prefer high value chips in small footprints:
    # For two solutions that cover the same chips, prefer the high priority chips to lie 
    #   in the small subfootprints, because then they get more bandwidth.
    #   Thus the larger subfootprints cover the lower priority chips, when possible.
    # Here is a start:
    # Recall M.A[f] = the size of footprint f, the number of chips it contains, its area
    #   M.A[f] = sum(M.X1[f,c] for f in M.F) 
    # f gets credit for covering chip c, but reduce this by some epsilon times the size of f
    # shorthand: p = linear_chip_priorities
    # sinc p[c] is constant, the following is still linear
    #   obj: + ( M.X1[f,c] - epsilon * M.A[f] ) * p[c]
    #
    # However, this isn't good enough because the size-penalty applies regardless of M.X1, so is meaningless.
    # Rather, we'd like
    #   obj: + z[f,c] * p[c]
    #   z[f,c] - s[f,c] = ( M.X1[f,c] - epsilon * M.A[f] ), z and s in [0,1], z and s real 
    # where s is a slack variable. The trick to getting z and s to take on the values we want, 
    # is to ensure z is 0 whenever M.X1[f,c] is 0, and s is zero whenever M.X1[f,c] is 1
    #   z[f,c] <=     M.X1[f,c]
    #   s[f,c] <= 1 - M.X1[f,c]
    # and we must replace the equalities by inequality. Since the objective encourages z to be large, we can use
    #   z[f,c] - s[f,c] <= ( M.X1[f,c] - epsilon * M.A[f] ), z and s in [0,1], z and s real 
    # Since a larger M.A results in a lower objective function value, we can use 
    #   M.A[f] >= sum(M.X1[f,c] for f in M.F) 
    # 
    # The downside to all this is it adds a lot of variables, two new variables for each f and c combination. 
    # At least they are reals and not binary...
 
    # Here is its implementation:
    
    # z in [0,1], real, z < x1
    M.z = pyomo.Var(M.F, M.C, within=pyomo.NonNegativeReals, bounds=(0, 1), initialize=0)
    def zx_rule(M, f, c):
        return M.z[f,c] <= M.X1[f,c]
    M.zx = pyomo.Constraint(M.F, M.C, rule=zx_rule)

    # s in [0,1], real, z < x1
    M.s = pyomo.Var(M.F, M.C, within=pyomo.NonNegativeReals, bounds=(0, 1), initialize=1)
    def sx_rule(M, f, c):
        return M.s[f,c] <= 1 - M.X1[f,c]
    M.sx = pyomo.Constraint(M.F, M.C, rule=sx_rule)

    # z - s = x - epsilon A
    def zsxa_rule(M, f, c):
        return M.z[f,c] - M.s[f,c] <= M.X1[f,c] - epsilon * M.A[f]  
    M.zsxa = pyomo.Constraint(M.F, M.C, rule=zsxa_rule)

    """
    Scott's elucidation on the slack variable machinery:
    To summarize
     
                 - sum(M.X1[f, c] for f in M.F for c in M.C) * epsilon2, sense=pyomo.maximize)
    This is meant to penalize covering zero-priority chips, so the subfootprints are only as big as they need to be to cover the valuable chips.
     
    The idea for the rest is to prefer that high value chips end up in small subfootprints rather than big subfootprints. 
    This is accomplished by setting up a multiplication of the subfootprint size times the chip priority in the objective.
    The goal here is for z to be almost X, but penalized by the size of the subfootprint. 
    This then gets multiplied by the chip priority in the objective, so smaller subfootprints have a small penalty for the high-value chips.
    Slack s is just machinery to turn off the penalty, and allow z to be equal to x, in the case that x = 0.
    (Not being able to turn off the penalty is what is wrong with the simpler objective commented out in lines 301-303.)
     
    The description in the code is hard for me to follow, and I wrote it!  
    So let's try an example.
    Suppose there are three active chips abc, with priorities 5,3,3, and we have two subfootprints that can cover any two of these chips.
    We'd prefer solution0 subfootA = (a) covering (5) and subfootB = (b,c) covering (3,3) rather than solution1 subfootA = (a,b) covering (5,3) and subfootB = (c) covering (3).
     
    For both solutions, the first two constraints are about the same, just assigning the non-zero z to the right subfootprint
    M.z(f, c) <= 1, for the chips c in footprint f, and z is 0 for other chip-footprint combinations
    M.s(f,c)  = 0  for the chips c in footprint f, and >= 0 for other chips-footprint combinations
     
    The third constraint is 
    z - s = x - epsilon Area
    In solution0, z[A,a] = 1 - 1 epsilon, z[B,b] = 1 - 2 epsilon, z[B,c] = 1 - 2 epsilon
    In solution1, z[A,a] = 1 - 2 epsilon, z[A,b] = 1 - 2 epsilon, z[B,c] = 1 - 1 epsilon
    And in the objective each of these gets multiplied by the chip priority, so 
    Solution0 obj = const  (1 * 5  + 2 *3 + 2 *3) eps = const  17*eps
    Solution1 obj = const  (2 * 5 +  2* 3 + 1 *3) eps = const19 *eps 
    Hence solution0 is preferred, as we wanted!
    """




    # obj z * p 
    # note this makes the penalty for large footprints redundant, except if a subfootprint has only zero-priority chips
#        expr=sum(M.z[f, c] * linear_chip_priorities[c] for f in M.F for c in M.C), sense=pyomo.maximize)
    M.obj_priority = pyomo.Objective(
        expr=sum(M.z[f, c] * linear_chip_priorities[c] for f in M.F for c in M.C) 
             - sum(M.X1[f, c] for f in M.F for c in M.C) * epsilon2, sense=pyomo.maximize)

    # return the Pyomo model (this is required).
    return M


################################################################################

#
# Print the generated footprint placement
#

def pyomo_postprocess(options, instance, chip_priorities=[], chips_to_cover=[], num_subfoots=0, chip_layout=[], pixel_coordinates=[],pixel_priorities=[],pixel_row_len=[],pixel_col_len=[],title=""):
    if chips_to_cover == []:
        chip_priorities, chips_to_cover = list_chips_of_interest('problem_instance.txt')
        num_subfoots = 4

    chip_x_len = np.shape(chip_priorities)[0]
    chip_y_len = np.shape(chip_priorities)[1]
    chip_row_len = chip_x_len
    chip_col_len = chip_y_len
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
    out_file = open("geo_sub_solution.txt", "w")
    for f in instance.F:
        ll_x = instance.Fkllx[f].value
        ll_y = instance.Fklly[f].value
        ur_x = instance.Fkurx[f].value
        ur_y = instance.Fkury[f].value
        print "Footprint %d placement: (%2d, %2d)  ,   (%2d, %2d)" % (f, ll_x, ll_y, ur_x, ur_y)
        out_str = str(ll_x) + " " + str(ll_y) + " " + str(ur_x) + " " + str(ur_y)
        out_file.write(out_str + "\n")
    total_area = sum(instance.X1[f, c] for f in instance.F for c in instance.C)

    # weighted area
    linear_chip_priorities = [chip_priorities[x][y] for x in range(chip_x_len) for y in range(chip_y_len)]
    weighted_area = sum(linear_chip_priorities[c] * instance.X1[f, c] for f in instance.F for c in instance.C)

    # graphically plot the subfootprints, depending on these flags
    dox = 1
    dov = 0

    chips = [(x, y) for x in range(chip_x_len) for y in range(chip_y_len)]
    for f in instance.F:

        print " Coordinate directions: "
        print "  y -> "
        print "x"
        print "|"
        print "v"

        if dox > 0:
            # X1
            print "Foot %d X1: " % (f)
            oldx = 0
            for c in instance.C:
                x = chips[c][0]
                y = chips[c][1]
                if x > oldx:
                    print " "
                    oldx = x
                xx = instance.X1[f, c].value
                # print "Foot %d chip %d %d x = %d" % (f, x, y, xx)
                print "%d" % (xx),
            print " "

        if dov > 0:
            # V1
            print "Foot %d V1: " % (f)
            oldx = 0
            for c in instance.C:
                x = chips[c][0]
                y = chips[c][1]
                if x > oldx:
                    print " "
                    oldx = x
                xx = instance.V1[f, c].value
                # print "Foot %d chip %d %d v1 = %d" % (f, x, y, xx)
                print "%d" % (xx),
            print " "

            # V2
            print "Foot %d V2: " % (f)
            oldx = 0
            for c in instance.C:
                x = chips[c][0]
                y = chips[c][1]
                if x > oldx:
                    print " "
                    oldx = x
                xx = instance.V2[f, c].value
                # print "Foot %d chip %d %d v1 = %d" % (f, x, y, xx)
                print "%d" % (xx),
            print " "

            # V3
            print "Foot %d V3: " % (f)
            oldx = 0
            for c in instance.C:
                x = chips[c][0]
                y = chips[c][1]
                if x > oldx:
                    print " "
                    oldx = x
                xx = instance.V3[f, c].value
                # print "Foot %d chip %d %d v1 = %d" % (f, x, y, xx)
                print "%d" % (xx),
            print " "

            # V4
            print "Foot %d V4: " % (f)
            oldx = 0
            for c in instance.C:
                x = chips[c][0]
                y = chips[c][1]
                if x > oldx:
                    print " "
                    oldx = x
                xx = instance.V4[f, c].value
                # print "Foot %d chip %d %d v1 = %d" % (f, x, y, xx)
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
    pp.PlotChipResults(chips_to_cover,chip_priorities,num_subfoots,chip_row_len,chip_col_len, chip_layout,pixel_coordinates,pixel_priorities,pixel_row_len,pixel_col_len,title)


################################################################################

def main():
    print("This script is not meant to be run at the command line." \
          "Please run createsubs.sh instead.")
    return 0


################################################################################

if __name__ == '__main__':
    exit(main())

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
T = 0   							# number of placement points
t_points = []               		# enumerated list of placement points
C = 0     							# num. of coverage points/constraints
c_points = []               		# enumerated list of timesteps

c_t_dict = {}                       # hash table of t_points that cover a
                                    # given coverage point
t_n_dict = {}                       # hash table of number of coverage 
                                    # points covered by a footprint
overlaps = {}                       # hash table of pairs of t_points and
                                    # corresponding normalized overlap
sigma = 0                           # normalizing constant equal to one
                                    # over the sum of pairwise distances
###############################################################################

#
# Initialize the dictionary for defining which footprints cover a given
# cover point
#
def initialize_dictionary():
    """
    This function initializes the hash table mapping coverage points 
    covered by a given footprint placement. In otherwords, when a 
    footprint is placed at placement points t1, t2, ..., tn coverage
    point ci is covered.
    """

    global T
    global t_points
    global C
    global c_points

    f = open('SimpleMPS/run/covered_by_footprints.txt', 'r')

    flag = True;
    for line in f:
        del_str = line.rstrip('\n').split()
        if flag:
            T = int(del_str[1])
            C = int(del_str[0])
            print "T = %d" % (T)
            print "C = %d" % (C)
            flag = False
        elif len(del_str):
            c_t_dict[int(del_str[0])] = map(int,del_str[1:])
    f.close()

    t_points = list(range(1,T+1))
    c_points = list(range(1,C+1))

    f = open('SimpleMPS/run/footprint_covers.txt', 'r')

    flag = True;
    for line in f:
        del_str = line.rstrip('\n').split()
        if flag:
            flag = False
        elif len(del_str):
            t_n_dict[int(del_str[0])] = len(del_str[1:])
    f.close()

###############################################################################

#
# Initialize the dictionary for defining which footprint pairs to
# consider.
#
def initialize_pairs():
    """
    This function initializes the hash table mapping placement pairs
    with corresponding distance between the footprint pair.
    """

    global sigma

    f = open('SimpleMPS/run/footprint_overlaps.txt', 'r')

    flag = True
    for line in f:
        del_str = line.rstrip('\n').split()
        count = 1
        if flag:
            flag = False
            continue
        elif len(del_str):
            while count < len(del_str):
                #if a key with swapped values exists, skip
                if (not overlaps.has_key(
                    (
                        int(del_str[count]),
                        int(del_str[0]))
                    )
                ):
                    overlaps[int(del_str[0]),int(del_str[count])] = float(del_str[count+1])
                count += 2
    f.close()

# This code will remain for future consideration. The initial goal was to
# define the sigma parameter(s) based on what the maximum potential overlap
# could be for a given problem instance.

#    run_sum = 0
#    for i in overlaps:
#        run_sum += overlaps[i]
#    sigma = 1/run_sum
    sigma = 0.001

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

    initialize_dictionary()
    initialize_pairs()

    M = pyomo.ConcreteModel()

    # A Pyomo Set representing the set of potential footprint placement
    # points.
    M.T = pyomo.Set(initialize=t_points)

    # A Pyomo Set representing the set of coverage points to be covered
    # by footprints
    M.C = pyomo.Set(initialize=c_points)

    # A Pyomo Set representing pairwise footprints. The goal is to only
    # consider pairs that correspond to footprints close enough to overlap.
    M.d = pyomo.RangeSet(1,len(overlaps))
    
    # A Pyomo Variable representing placed footprints.
    M.dT = pyomo.Var(M.T, within=pyomo.Binary)

    # A Pyomo Variable representing pairs of placed footprints. Again, the
    # goal is to only have pairs for footprints close enough To overlap.
    M.p = pyomo.Var(M.d, within=pyomo.Binary)

    #
    # For each coverage point, the sum of footprint placements covering
    # this point must be greater than or equal to one. Over all cover
    # points this ensures the entire domain is covered by footprints.
    #
    def c0_rule(M, c):
        dT_sum = 0
        for i in c_t_dict.get(c):
            dT_sum += M.dT[i]
        return dT_sum >= 1
    M.c0 = pyomo.Constraint(M.C, rule=c0_rule)

    #
    # Only one footprint can be placed at a placement point
    #
    def c1_rule(M, t):
        return M.dT[t] <= 1
    M.c1 = pyomo.Constraint(M.T, rule=c1_rule)

    #
    # A footprint pair is one only if both placement points of the pair
    # have footprints. The next three constraints accomplish this. This
    # is equivalent to a logical AND (p[i,j] = 1 iff dT[i] AND dT[j] = 1).
    #
    def c2_rule(M, d):
        return M.p[d] >= M.dT[overlaps.keys()[d-1][0]] + M.dT[overlaps.keys()[d-1][1]] - 1
    M.c2 = pyomo.Constraint(M.d, rule=c2_rule)

    def c3_rule(M, d):
        return M.p[d] <= M.dT[overlaps.keys()[d-1][0]]
    M.c3 = pyomo.Constraint(M.d, rule=c3_rule)

    def c4_rule(M, d):
        return M.p[d] <= M.dT[overlaps.keys()[d-1][1]]
    M.c4 = pyomo.Constraint(M.d, rule=c4_rule)

    #
    # Objective function: Minimize the number of footprints used to
    # satisfy the constraint of covering all coverage points while also
    # minimizing the overlap between placed footprints.
    #
    if 1:
        M.obj = pyomo.Objective(expr=sum(M.dT[t] for t in M.T) + sigma*
                            sum(overlaps.values()[i-1]*M.p[i] for i in M.d),
                            sense=pyomo.minimize)
    #
    # Objective function: Minimize the number of footprints used to
    # satisfy the constraint of covering all coverage points while also
    # minimizing the overlap between placed footprints and maximizing the
    # number of coverage points contained in each placed footprint.
    #
    else:
        M.obj = pyomo.Objective(expr=sum(M.dT[t] for t in M.T) + sigma*
                            sum(overlaps.values()[i-1]*M.p[i] for i in M.d) - 
                            sigma*sum(t_n_dict.values()[i-1]*M.dT[i] for i in
                            M.T), sense=pyomo.minimize)

    # return the Pyomo model (this is required).	
    return M

################################################################################

#
# Print the generated footprint placement
#

def pyomo_postprocess(options, instance, results):
    foots = []
    print '***********************************************************' \
          '*****************'
    print '***********************************************************' \
          '*****************'
    print ""
    print "Number of footprints placed: %d" % (sum(instance.dT[t].value for t
           in instance.T))
    print ""
    print '***********************************************************' \
          '*****************'
    print '***********************************************************' \
          '*****************'
    print ""
    print "Footprint placements:"
    for t in instance.T:
        if round(float(instance.dT[t].value)) == 1:
            print t
            foots.append(t)
    print ""
    print '***********************************************************' \
          '*****************'
    print '***********************************************************' \
          '*****************'

    out_fn = "SimpleMPS/run/solution.txt"
    out_file = open(out_fn, 'w')
    
    print >>out_file, ("%d" % (len(foots)))
    for i in foots:
        print >>out_file, ("%d footA.txt" % (i))

################################################################################

def main():
	print("This script is not meant to be run at the command line." \
		  "Please run createplacements.sh instead.")
	return 0

################################################################################

if __name__=='__main__':
  exit(main())


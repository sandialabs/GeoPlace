"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import test_generator
import pixel_shift
import numpy as np
import pyomo.environ as pyomo
import geo_sub_priority_multi as geosub

"""
Our objectives are the following:
1: Read a bitmap image and get a corresponding pixel_list
2: Find all shifts of the pixel_list
3: Generate chip_layouts per each
4: Loop over the chip layouts, applying either the MIP or a Heuristic or some other solver
"""

print "\n"
print "Copyright 2016 Sandia Corporation. Under the terms of Contract"
print "DE-AC04-94AL85000, there is a non-exclusive license for use of"
print "this work by or on behalf of the U.S. Government. Export of this"
print "program may require a license from the United States Government."
print "\n"

# Step 1: Read in the pixel list from a bitmap

# filename = "test_pixel_list.bmp"
#filename = "examples/mississippi/cities20.bmp"
#filename = "examples/mississippi/junctions40.bmp"
#filename = "examples/mississippi/JunctionsCities40.bmp"
#filename = "examples/mississippi/riversJunctions40.bmp"
#filename = "examples/mississippi/rivers40.bmp"
filename = "examples/mississippi/rivers80.bmp"

# machine generated, not read from file
# pixel_list, pixel_priorities, pixel_shape, chip_shape = test_generator.GenerateSimplePixelLayout()

# define number of chips in the footprint, in each direction
# TODO: This is hard-coded at the moment.
# these are assumed to divide evenly into the image shape
# e.g. 15 x 20 pixels = 3 x 4 chips of 5 x 5 pixels each
#
# for 20 = 15 x 20
#num_chips_per_row = 5 
#num_chips_per_col = 5
# for 40 = 30 x 40
#num_chips_per_row = 6
#num_chips_per_col = 8
# for 80 = 60 x 80
#num_chips_per_row = 10
#num_chips_per_col = 10 
# for 80 = 60 x 80, or more
num_chips_per_row = 15
num_chips_per_col = 20 

chip_shape = (num_chips_per_row, num_chips_per_col)

# read bitmap, its dimensions are pixel_shape, its list of non-zeros is pixel_list and those values are pixel_priorities 
pixel_list, pixel_priorities, pixel_shape = pixel_shift.BitmapToPixelList(filename)

# Step 2 and 3: Get all possible chip layouts from shifted pixels positions.

pixels_per_chip_row = pixel_shape[0] / chip_shape[0]
pixels_per_chip_col = pixel_shape[1] / chip_shape[1]

print "domain is", pixel_shape[0], "x", pixel_shape[1], "pixels"
print "footprint is", num_chips_per_row*pixels_per_chip_row, "x", num_chips_per_col*pixels_per_chip_col, "pixels"
print "footprint is", num_chips_per_row, "x", num_chips_per_col, "chips"
print "chip is", pixels_per_chip_row, "x", pixels_per_chip_col, "pixels"

# todo: modify pixel list and chip shape before passing in to Find, by adding padding of empty pixels.
# on the postprocess, remove that extra padding.
all_chip_layouts, pixel_configuration_list = pixel_shift.FindAllChipLayouts(pixel_list, pixels_per_chip_row,
                                                                            pixels_per_chip_col, chip_shape[0],
                                                                            chip_shape[1])

opt = pyomo.SolverFactory('gurobi')
#opt = pyomo.SolverFactory('cplex')

# Step 4: Loop over chip layouts

for index, chip_layout in enumerate(all_chip_layouts):

    print "chip layout %d of %d" % (index, len(all_chip_layouts)-1 )

    chip_priorities = pixel_shift.GenerateChipPriorities(chip_layout, pixel_configuration_list[index], pixel_priorities,
                                                         pixels_per_chip_row, pixels_per_chip_col)

    chips_to_cover = pixel_shift.ChipsToCover(chip_layout)

    M = geosub.pyomo_create_model(None, None, chip_priorities, chips_to_cover, 4, 16, chip_shape[0], chip_shape[1])

    instance = M.create()

    results = opt.solve(instance)
    instance.load(results)
    geosub.pyomo_postprocess(None, instance, chip_priorities, chips_to_cover, 4, chip_layout, pixel_configuration_list[index],pixel_priorities,pixel_shape[0],pixel_shape[1])

# for list_of_chips in all_chip_layouts:
# DO SOMETHING
#    pass

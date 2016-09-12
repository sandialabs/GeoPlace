"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import math
import numpy as np


def greedy(ff, bandwidth, chip_shape, chips_to_cover, chip_priorities):
    sz = int(math.floor(bandwidth / ff))
    subfoot_zoo = generate_subfoot_zoo(sz)

    covered_chip_matrix = np.zeros(chip_shape)
    potential_greedy_placements = []
    placement_value = []
    # for all subfootprint shapes 
    for subfoot_shape in subfoot_zoo:
        # for all upper-left corners of the subfootprint
        for col in range(chip_shape[1] - subfoot_shape[1] + 1):
            for row in range(chip_shape[0] - subfoot_shape[0] + 1):
                # calculate value of covered chips
                subfoot_priority = 0
                for chip_col in range(subfoot_shape[1]):
                    for chip_row in range(subfoot_shape[0]):
                        subfoot_priority += chip_priorities[row + chip_row, col + chip_col]
                # save value and shape of the placement
                potential_greedy_placements.append((row, col, subfoot_priority, subfoot_shape))
                placement_value.append(subfoot_priority)

    placement_index = [i[0] for i in sorted(enumerate(placement_value), key=lambda x: x[1])]

    covered_chip_matrix = np.zeros(chip_shape)

    footprint_choices = []
    f = 0
    p = len(placement_index) - 1
    while f < ff and p >= 0:
        footprint_choice = potential_greedy_placements[placement_index[p]]
        uncovered = True
        for chip_col in range(footprint_choice[3][1]):
            for chip_row in range(footprint_choice[3][0]):
                if covered_chip_matrix[footprint_choice[0] + chip_row, footprint_choice[1] + chip_col] == 1:
                    uncovered = False
        if uncovered:
            f = f + 1
            footprint_choices.append(footprint_choice)
            for chip_col in range(footprint_choice[3][1]):
                for chip_row in range(footprint_choice[3][0]):
                    covered_chip_matrix[footprint_choice[0] + chip_row, footprint_choice[1] + chip_col] = 1
        p = p - 1
    return footprint_choices


# sz = maximum number of chips in a single subfootprint
def generate_subfoot_zoo(sz):
    subfoot_zoo = []
    for i in range(1, sz + 1):
        for j in range(1, int(math.floor(sz / i)) + 1):
            subfoot_zoo.append((i, j))
    return subfoot_zoo

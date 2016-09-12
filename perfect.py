"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import csv
import numpy as np


def list_chips_of_interest(problem_file):
    """
    This method builds the list of chips of interest.
    """
    with open(problem_file, 'rb') as file:
        reader = csv.reader(file, delimiter=',')
        chips_to_cover_priorities = []
        chips_to_cover = []
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
                        chips_to_cover_priorities.append(float(value))
    return chip_priorities, chips_to_cover_priorities, chips_to_cover


def perfect(chip_priorities, bw):
    chip_pris_array = np.array(chip_priorities)
    value = sum(sorted(chip_pris_array.ravel())[-bw:])
    sorted_index_list = [i[0] for i in sorted(enumerate(chip_pris_array.ravel()), key=lambda x: x[1])][-bw:]
    indices = unravel(sorted_index_list, chip_pris_array.shape)
    return value, indices


def unravel(index_list, shape):
    num_cols = shape[1]

    indices = []
    for index in index_list:
        indices.append((index / num_cols, np.mod(index, num_cols)))

    return indices

"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import csv
import numpy as np

def generate_pixel_layout():
    chip_priorities, chips_to_cover = list_chips_of_interest('problem_instance.txt')
    chip_shape_in_pixels = (2,3)
    num_rows = chip_shape_in_pixels[0]*np.shape(chip_priorities)[0]
    num_cols = chip_shape_in_pixels[1]*np.shape(chip_priorities)[1]
    pixel_array = np.zeros((num_rows,num_cols))

    for chip in chips_to_cover:
        priority = 25*chip_priorities[chip[0]][chip[1]]
        pixel_row = chip_shape_in_pixels[0]*chip[0]
        pixel_col = chip_shape_in_pixels[1]*chip[1]
        pixel_array[pixel_row][pixel_col] = priority
        
    return pixel_array

def GenerateSimplePixelLayout():
    pixel_list = []
    pixel_list.append((0, 0))
    pixel_list.append((2, 2))
    pixel_list.append((1, 4))
    pixel_list.append((1, 7))
    pixel_list.append((1,1))

    pixel_priorities = [1,3,5,7,9]

    image_shape = (4,10) #4 rows of pixels, 10 columns of pixels
    chip_shape = (2,5) #2 rows of chips, 5 columns of chips
    return pixel_list, pixel_priorities, image_shape, chip_shape

def list_chips_of_interest(problem_file):
    """
    This method builds the list of chips of interest.
    """
    chips_to_cover = []
    with open(problem_file, 'rb') as file:
        reader = csv.reader(file, delimiter=',')
        for line in reader:
            if reader.line_num == 1:
                chip_row_len = int(line[0])
                chip_col_len = int(line[1])
                chip_priorities = [[0 for i in range(chip_row_len)] for j in range(chip_col_len)]
            elif reader.line_num == 2:
                max_priority = float(line[0])
            else:
#                row = chip_row_len - reader.line_num + 2
                row = reader.line_num - 3
                for col, value in enumerate(line):
                    if value is ' ':
                        chip_priorities[row][col] = 0.0
                    else:
                        chip_priorities[row][col] = float(value)
                        chips_to_cover.append([row, col])
    return chip_priorities, chips_to_cover


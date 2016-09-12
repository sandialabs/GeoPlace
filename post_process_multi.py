"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from geo_sub import list_chips_of_interest
import pickle
import numpy as np
import re
import csv

import os # Used to grab environment variables
def ReadSolution():

  # Read in the geo_sub_solution.txt file.
  sub_feet = []
  with open("geo_sub_solution.txt", "r") as f:
     for line in f:
        line = line.replace("\n","") #Remove newline character
        split_line = line.split(" ") #Split digits
        chip = [] 
        for item in split_line:
           int_value = int(round(float(item)))
           chip.append(int_value) # Cast as a number, then to an int.
        sub_feet.append(chip)
  return sub_feet

def GetChipsToCover():
    chips_to_cover = [] #Initialize empty list
    problem_file = os.environ.get("PROBLEM_FILE", "problem_instance.txt")
    if not problem_file:
       problem_file = "problem_simple.txt"
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
    return chips_to_cover

def GetNumSubfootprints():
   num_subfootprints = os.environ.get("NUM_SUBFOOTPRINTS", "4")
   if not num_subfootprints:
      num_subfootprints = 4
   return int(num_subfootprints) 


# def PlotChipResults(chips_to_cover,chip_priorities,num_subfoots,chip_row_len,chip_col_len):
def PlotPerfectChipResults(index_list,chip_row_len,chip_col_len,title=""):
    import re
    fig = plt.figure()
    if title:
        plt.suptitle(title)


    # chips and their priorities
    plt.subplot(132)
    plt.title("Chips to Cover")
    cm = plt.get_cmap('Greys')
    ax = fig.gca()
    ax.set_xticks(np.arange(0,chip_col_len+1,1))
    ax.set_yticks(np.arange(0,chip_row_len+1,1))
#   top left is 0,0
    plt.gca().invert_yaxis()

    # subfootprint layout
    plt.subplot(133)
    plt.title("Subfootprints Placed")
    ax = fig.gca()
    ax.set_xticks(np.arange(0,chip_col_len+1,1))
    ax.set_yticks(np.arange(0,chip_row_len+1,1))
#   top left is 0,0
    plt.gca().invert_yaxis()
    color_dictionary = {0: "red", 1: "teal", 2: "blue", 3: "green", 4: "black", 5: "red", 6: "purple", 7: "pink", 8: "lavender"}
    for count, foot in enumerate(index_list):
       index = count % len(color_dictionary)  # Repeat colors
       color = color_dictionary[count]
       #Draw rectangle....
       delta_x = 1
       delta_y = 1
       ax.add_patch(Rectangle((foot[1],foot[0]),delta_x,delta_y,facecolor=color))
    print("I'm inside the plotter!")
    plt.grid()
#    plt.show(block=False)   
    if title:
        plt.savefig(re.sub(' ','_',title) + ".png")   

def PlotPixels(pixel_row_len, pixel_col_len, chip_row_len, chip_col_len, pixel_coordinates, pixel_priorities):
    
    # pixels in chips
    fig = plt.figure()    
    plt.title("Pixels in Chips")
    cm = plt.get_cmap('Greys')
    ax = fig.gca()
    pixels_per_chip_row = pixel_row_len / chip_row_len
    pixels_per_chip_col = pixel_col_len / chip_col_len
    ax.set_xticks(np.arange(0,pixel_col_len+1,pixels_per_chip_row))
    ax.set_yticks(np.arange(0,pixel_row_len+1,pixels_per_chip_col))
    #ax.set_xticklabels(np.array([0,4,8]))
#   top left is 0,0

    plt.gca().invert_yaxis()
    max_pri = np.max(pixel_priorities)
    
    for index, pixel in enumerate(pixel_coordinates):
       pr = pixel[0] # row
       pc = pixel[1] # col
       px = pc
       py = pr
       priority = pixel_priorities[index]/ float(max_pri)
#      avoid almost-white squares by a shift
       priority_color = 0.1 + 0.9 * (priority)
       ax.add_patch(Rectangle((px,py),1,1,facecolor=cm(priority_color)))
    
    # add in chip lines
    plt.grid()
    plt.show()
    pickle.dump(ax, file('pixels_in_chips.pkl', 'w'))
    

def PlotChipResults(chips_to_cover,chip_priorities,num_subfoots,chip_row_len,chip_col_len, chip_layout,pixel_coordinates,pixel_priorities,pixel_row_len,pixel_col_len, title="", solution=[]):
    import re
    fig = plt.figure()
    if title:
        plt.suptitle(title)

    # chips and their priorities
    plt.subplot(121)
    plt.title("Chips to Cover")
    cm = plt.get_cmap('Greys')
    ax = fig.gca()
    ax.set_xticks(np.arange(0,chip_col_len+1,1))
    ax.set_yticks(np.arange(0,chip_row_len+1,1))
#   top left is 0,0
    plt.gca().invert_yaxis()
    #Color over these chips....
    max_pri = np.max(chip_priorities)
    for index, chip in enumerate(chips_to_cover):
       cr = chip[0] # row
       cc = chip[1] # col
       chipx = cc
       chipy = cr
       priority = chip_priorities[cr][cc]/max_pri
#      avoid almost-white squares by a shift
       priority_color = 0.1 + 0.9 * (priority)
       ax.add_patch(Rectangle((chipx,chipy),1,1,facecolor=cm(priority_color)))
    plt.grid()




    plt.subplot(122)
    plt.title("Subfootprints Placed")
    ax = fig.gca()
    ax.set_xticks(np.arange(0,chip_col_len+1,1))
    ax.set_yticks(np.arange(0,chip_row_len+1,1))
#   top left is 0,0
    plt.gca().invert_yaxis()
    color_dictionary = {0: "red", 1: "teal", 2: "blue", 3: "green", 4: "black", 5: "purple", 6: "lavendar", 7: "pink", 8: "yellow"}
    for count, foot in enumerate(solution):
       index = count % len(color_dictionary)  # Repeat colors
       color = color_dictionary[count]
       #Draw rectangle....
       delta_x = foot[3]-foot[1]+1
       delta_y = foot[2]-foot[0]+1
       ax.add_patch(Rectangle((foot[1],foot[0]),delta_x,delta_y,facecolor=color))
   
    plt.grid()
#    plt.show(block=False)   
    if title:
        plt.savefig(re.sub(' ','_',title) + ".png")   

    pickle.dump(ax, file(title + '.pkl', 'w'))

def GreedySolutionConvert(greedy_solution):
    # Greedy solution is of the form 
    # [(5, 10, 4405.0, (4, 1)),
    # (2, 2, 4202.0, (1, 4)),
    # (8, 4, 4142.0, (1, 4)),
    # (10, 10, 3732.0, (4, 1))]
    footprints = []
    for footprint in greedy_solution:
        print("The footprint in greedy solution is", footprint)
        corners = GenerateFootprintCorners(footprint)
        print("The footprint is", corners)
        footprints.append(corners)
    return footprints

def GenerateFootprintCorners(footprint):
    row_upper_left = footprint[0]
    col_upper_left = footprint[1]
    footprint_shape = footprint[3]
    num_rows_to_add = footprint_shape[0]-1 
    num_cols_to_add = footprint_shape[1]-1
    row_bottom_right = row_upper_left + num_rows_to_add
    col_bottom_right = col_upper_left + num_cols_to_add
    return [row_upper_left, col_upper_left, row_bottom_right, col_bottom_right]

def ReloadImage(filename):
    return pickle.load(file(filename))

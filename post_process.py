"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from geo_sub import list_chips_of_interest
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
#Fake initial work
num_subfootprints = GetNumSubfootprints()
chip_x_len = 12
chip_y_len = 12

chips_to_cover = GetChipsToCover()


fig = plt.figure(1)
plt.subplot(121)
ax = fig.gca()
ax.set_xticks(np.arange(0,chip_x_len+1,1))
ax.set_yticks(np.arange(0,chip_y_len+1,1))
plt.title("Chips to Cover")
#Color over these chips....
for chip in chips_to_cover:
   chipx = chip[0]
   chipy = chip[1]
   ax.add_patch(Rectangle((chipx,chipy),1,1,facecolor="grey"))
plt.grid()


sub_feet = ReadSolution()
plt.subplot(122)
plt.title("Subfootprints Placed")
ax = fig.gca()
ax.set_xticks(np.arange(0,chip_x_len+1,1))
ax.set_yticks(np.arange(0,chip_y_len+1,1))
color_dictionary = {0: "red", 1: "teal", 2: "blue", 3: "green", 4: "black", 5: "red", 6: "purple", 7: "pink", 8: "lavender"}
for count, foot in enumerate(sub_feet):
   index = count % len(color_dictionary)  # Repeat colors
   color = color_dictionary[count]
   #Draw rectangle....
   delta_x = foot[2]-foot[0]+1
   delta_y = foot[3]-foot[1]+1
   ax.add_patch(Rectangle((foot[0],foot[1]),delta_x,delta_y,facecolor=color))
   
plt.grid()
plt.show()   

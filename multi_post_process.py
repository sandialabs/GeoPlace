"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import glob

def ReadOutputFile(filename):
   f = open(filename, 'r')
   lines = f.readlines()
   for line in lines:
      if "objective function value" in line:
         print(line)

def GetListOfFilenames():
    # What if the file is old...? If we ran 2...5 and then 2..3, 4 and 5 would still exist
    filename_list = []
    for name in glob.glob('multi_output*.txt'):
       filename_list.append(name)
    return filename_list

filename_list = GetListOfFilenames()
print("The filename list is", filename_list)
for filename in filename_list:
   ReadOutputFile(filename)

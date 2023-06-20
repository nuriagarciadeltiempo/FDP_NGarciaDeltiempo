
import pandas as pd
import sys
import re
import csv

# read the output1 file
# is it necessary to remove childs from the output1 file
# create a dictionary with the Individual ID as a key and the gender as the value.

with open(r"pedigree_GP.csv") as output1_file:
    ID2gender = {}
    for line in output1_file:
        line = line.split("\t")
        if line[2] == '1\n':
            ID2gender[line[1]] = line[2].replace("1", "MALE")
        elif line[2] == '2\n':
            ID2gender[line[1]] = line[2].replace("2", "FEMALE")
        
        #print(ID2gender)

with open("ID2gender.txt", "w") as output1_file:
    for key, value in ID2gender.items():
        output1_file.write('%s:%s\n' % (key, value))






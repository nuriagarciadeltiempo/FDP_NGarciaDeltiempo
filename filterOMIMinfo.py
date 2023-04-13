
import pandas as pd


#import csv omim file as data frame
OMIM_file = pd.read_csv(r'OMIM_X-linked_non_adult_onset.csv', sep=';')

#remove columns 1,2 and 3 from OMIM_file
OMIM_file = OMIM_file.drop(OMIM_file.columns[[1,2,3]], axis=1)

#store OMIM_file into a new csv
OMIM_file.to_csv('OMIM_X-linked_non_adult_onset_removedcol.csv', sep=';')

#!/usr/bin/python

import sys
import os
import pandas


import sys
import pandas
counter = 0

for i in sys.argv[1:]:
    sample = i.split("/")[-1].split("_HV.count")[0]
    counter += 1
    if counter == 1:
        pandasdf = pandas.read_csv(i, sep = "\t", index_col=0,skiprows=[0], names = ["Geneid",sample], usecols=[0,6])
    else:
        pandasdf2 = pandas.read_csv(i, sep = "\t", index_col=0,skiprows=[0], names = ["Geneid",sample], usecols=[0,6])
        #pandasdf = pandas.merge(pandasdf, pandasdf2, left_index=True, right_index=True, how = 'outer').fillna(0) # If you want integers, set astype(int) in the end
        pandasdf = pandas.merge(pandasdf, pandasdf2, left_index=True, right_index=True, how = 'outer').fillna(0).astype(int)


pandasdf.to_csv("MergedTable_RawCounts.txt", sep = "\t",na_rep='NA')

#!/usr/bin/python

# script to harvest empirical likelihoods from random trees and print
# out % confidence boundaries

import sys
import pdb
import glob
import math

# input args
emp_tree_d = sys.argv[1]
write_d = sys.argv[2]
thresh = float(sys.argv[3])

# find the files
trees = glob.glob(emp_tree_d + "*")

# where to write the results
out_f = open(write_d + "/thresh.file",'w')
lik_dict = {}

for tree in trees:

    # open the likelihood file
    lik_file = open(tree + "/lik.file",'r')
    lik_lines = lik_file.readlines()

    for line in lik_lines:
        parts = line.strip().split()
        leaves = parts[0] + " " + parts[1]
        if leaves not in lik_dict:
            lik_dict[leaves] = [float(parts[2])]
        else:
            lik_dict[leaves].append(float(parts[2]))

    lik_file.close()

# sort the likelihoods
for leaves in lik_dict:
    lik_dict[leaves].sort()
    perc = lik_dict[leaves][int(len(lik_dict[leaves])*thresh)]
    # write out the results
    out_f.write(leaves + " " + str(perc) + "\n")

out_f.close()

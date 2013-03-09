#!/usr/bin/python
# AdaptML

import os
import sys
import pdb
import time
import random
import glob

from scipy.io import write_array
from scipy.io import read_array
from numpy.linalg import *
from numpy.core import *
from numpy.lib import *
from numpy import *

import rand_multitree as multitree
import rand_ML as ML

start_time = time.time()
sys.setrecursionlimit(10000)

# potential variables:
tree_fn = None
outgroup = None
migration_fn = None
mu_fn = None
thresh_fn = None
ultratree_fn = None
color_map_fn = None
write_dir = None
cdist = False
to_truncate = False
obs_states = None
iter_limit = 10000

# load inputs #
inputs = sys.argv
for ind in range(1,len(inputs)):
    arg_parts = inputs[ind].split('=')
    code = arg_parts[0]
    arg = arg_parts[1]
    if code == 'tree':
        tree_fn = arg
    elif code == 'truncate':
        to_truncate = True
    elif code == 'outgroup':
        outgroup = [arg]
    elif code == 'habitats':
        migration_fn = arg
    elif code == 'mu':
        mu_fn = arg
    elif code == 'color':
        color_fn = arg
    elif code == 'write':
        write_dir = arg
    elif code =='thresh':
        thresh_fn = arg
    elif code == 'cdist':
        cdist = True
    elif code == 'ultra':
        ultratree_fn = arg
    elif code == 'iters':
        iter_limit = int(arg)

# load the parameters
migration_file = open(migration_fn,'r')
migration_matrix = eval(migration_file.read())
mu_file = open(mu_fn,'r')
mu = float(mu_file.read().strip())

# build the tree #
tree_file = open(tree_fn,"r")
tree_string = tree_file.read().strip()
tree = multitree.multitree()
tree.build(tree_string)
tree_file.close()

# root it and remove zero branch lengths
min_len = min([b.length for b in tree.branch_list if b.length > 0.0])
for b in tree.branch_list:
    if b.length <= 0.0:
        b.length = min_len
    names_1 = b.ends[0].name_dict[b.ends[1]]
    names_2 = b.ends[1].name_dict[b.ends[0]]
    if names_1 == outgroup or names_2 == outgroup:
        tree.rootify(b)

# wipe the likelihoods off of the tree (just to be sure)
ML.TreeWipe(tree.root)

# learn the likelihoods
root = tree.root
kids = root.GetKids()
if kids[0].name in outgroup:
    true_root = kids[1]
else:
    true_root = kids[0]
ML.LearnLiks(tree,mu,migration_matrix,true_root)

# estimate the states
ML.EstimateStates(true_root)

# begin to generate the random data
count = len(glob.glob(write_dir + "/*"))
if count > 0:
    print "write directory not empty; exiting"
    sys.exit(1)

while count < iter_limit:

    # politely make a new directory
    dir_name = write_dir + "/" + str(count) + "/"

    os.mkdir(dir_name)
    count = len(glob.glob(write_dir + "/*"))

    ML.TreeWipe(tree.root)

    # shuffle the tree's leaves
    tree.LeafShuffle()
    
    # relearn the likelihoods using the new leaves and the old states
    ML.LearnShuffleLiks(tree,mu,migration_matrix,true_root)

    # write out the likelihoods at each internal node
    lik_file = open(dir_name + "lik.file",'w')
    files = {}
    files['lik'] = lik_file
    true_root.PieCharts(files)
    lik_file.close()
    
    print float(count) / iter_limit

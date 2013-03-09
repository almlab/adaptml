#!/usr/bin/python
# AdaptML

import os
import sys
import pdb
import time
import random

from numpy.linalg import *
from numpy.core import *
from numpy.lib import *
from numpy import *

import multitree
import ML

start_time = time.time()
sys.setrecursionlimit(25000)

# load inputs #
tree_filename = None
hab_num = 16
outgroup = None
write_dir = './'
rateopt = 'avg'
mu = 1.00000000001
habitat_thresh = 0.10
converge_thresh = 0.001

inputs = sys.argv
for ind in range(1,len(inputs)):
    arg_parts = inputs[ind].split('=')
    code = arg_parts[0]
    arg = arg_parts[1]
    if code == 'tree':
        tree_filename = arg
    elif code == 'init_hab_num':
        hab_num = int(arg)
    elif code == 'outgroup':
        outgroup = [arg]
    elif code == 'converge_thresh':
        converge_thresh = float(arg)
    elif code == 'write_dir':
        write_dir = arg
    elif code == 'mu':
        mu = float(arg)
    elif code == 'rateopt':
        rateopt = arg
    elif code == 'collapse_thresh':
        habitat_thresh = float(arg)

# track stats
stats_file = open(write_dir + '/stats.file','w')

# build the tree #
print "\nBuilding Tree"
tree_file = open(tree_filename,"r")
tree_string = tree_file.read().strip()
tree = multitree.multitree()
tree.build(tree_string)
tree_file.close()

# remove zero branches
min_len = min([b.length for b in tree.branch_list if b.length > 0.0])
for b in tree.branch_list:
    if b.length <= 0.0:
        b.length = min_len

################################
# build an initial rate matrix #
################################

# how many species are there?
species_dict = tree.species_count
total_leaves = sum(species_dict.values())
habitat_list = []
filter_list = species_dict.keys()
for i in range(hab_num):
    habitat_list.append("habitat " + str(i))

# create O(n^3) habitat matrix
print "Instantiating Habitat Matrix"
habitat_matrix = {}
for habitat in habitat_list:
    habitat_matrix[habitat] = {}
    for filt in filter_list:
        habitat_matrix[habitat][filt] = random.random()

    # normalize
    scale = sum(habitat_matrix[habitat].values())
    for filt in filter_list:
        habitat_matrix[habitat][filt] /= scale

score = -9999.99999999
diff = 1.0
old_diff = 1.0

print "Learning Habitats:"
while 1:
    counter = 0

    print "\t" + str(len(habitat_matrix)) + " habitats"
    print "\tRefinement Steps [d(Habitat Score)]: "
    stats_str = ""
    stats_str += "counter\t"
    stats_str += "habs\t"        
    stats_str += "ML score\t"
    stats_str += "mu\t"        
    stats_str += "\thabitat dist diff\t"        
    stats_file.write(stats_str + "\n")

    while 1:
        stats_str = ""
        stats_str += str(counter) + "\t"
        stats_str += str(len(habitat_matrix)) + "\t"        
        stats_str += str(score) + "\t"
        stats_str += str(mu) + "\t"        
        stats_str += str(diff) + "\t"        
        stats_file.write(stats_str + "\n")
        stats_file.flush()
        print "\t\t" + str(counter) + "\t" + str(diff)

        # wipe the likelihoods off of the tree
        ML.TreeWipe(tree)

        # learn the likelihoods
        ML.LearnLiks(tree,mu,habitat_matrix)

        # estimate the states (by making each node trifurcating...)
        ML.EstimateStates(tree.a_node,habitat_matrix)
        
        # upgrade guesses for mu and habitat matrix
        this_migrate = habitat_matrix
        mu, habitat_matrix = ML.LearnRates(tree,mu,habitat_matrix,rateopt)
        new_migrate = habitat_matrix

        # stop?
        old_diff = diff
        score, diff = ML.CheckConverge(tree,new_migrate,this_migrate)

        if diff < converge_thresh:
            break

        # this should break the loop if you end up bouncing back and
        # forth between the same values
        sig_figs = 8
        diff1 = math.floor(diff*math.pow(10,sig_figs))
        diff2 = math.floor(old_diff*math.pow(10,sig_figs))        
        if diff1 > 0:
            if diff1 == diff2:
                break
        if counter > 500:
            break
        counter += 1

    #########################
    # remove similar groups #
    #########################
    print "Removing Redundant Habitats"
    new_habitats = {}
    for habitat_1 in habitat_matrix:
        old_habitat = habitat_matrix[habitat_1]
        add_habitat = True
        for habitat_2 in new_habitats:
            new_habitat = new_habitats[habitat_2]
            score = 0
            for this_filter in old_habitat:
                diff = old_habitat[this_filter] - new_habitat[this_filter]
                score += math.pow(diff,2)
            if score < habitat_thresh:
                add_habitat = False
        if add_habitat:
            new_habitats[habitat_1] = habitat_matrix[habitat_1]
    if len(new_habitats) == len(habitat_matrix):
        break
    habitat_matrix = new_habitats
    if len(habitat_matrix) < 2:
        break
print "Learned " + str(len(habitat_matrix)) + " habitats",
print "in " + str(time.clock()) + " seconds"
stats_file.write("\nEnd Of Run\n")

############################
# take the best parameters #
############################

# find the branch you're interested in by looking at the partition on
# leaves the branch defines
for b in tree.branch_list:
    names_1 = b.ends[0].name_dict[b.ends[1]]
    names_2 = b.ends[1].name_dict[b.ends[0]]
    if names_1 == outgroup or names_2 == outgroup:
        tree.rootify(b)

# write out the results
mu_file = open(write_dir + '/mu.val','w')
mu_file.write(str(mu))
mu_file.close()

habitat_file = open(write_dir + '/habitat.matrix','w')
habitat_file.write(str(habitat_matrix))
habitat_file.close()
stats_file.close()


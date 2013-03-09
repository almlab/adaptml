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

# potential variables:
tree_fn = None
outgroup = None
migration_fn = None
mu_fn = None
thresh_fn = None
color_map_fn = None
write_dir = None
cdist = False
to_truncate = False
obs_states = None

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
        
# read in data files
params = {}
params['cdist'] = cdist

# migration
migration_f = open(migration_fn,'r')
migration_matrix = eval(migration_f.read())

# mu 
mu_f = open(mu_fn,'r')
mu = float(mu_f.read().strip())
mu_f.close()

# colors
color_f = open(color_fn,'r')
color_hash = {}
for line in color_f:
    parts = line.strip().split(' ')
    ring_number = parts[0]
    ring_state = parts[1]
    hex_vals = [hex(int(part)) for part in parts[2:]]
    hex_code = "#"
    if ring_number not in color_hash:
        color_hash[ring_number] = {}
    for hex_val in hex_vals:
        hex_color = hex_val[2:]
        if len(hex_color) < 2:
            hex_color = "0" + hex_color
        hex_code += hex_color
    color_hash[ring_number][ring_state] = hex_code
color_f.close()

# if specified, grab threshold file
thresh_dict = None
if thresh_fn is not None:
    # grab all of the thresholds and toss into a hash
    thresh_dict = {}
    thresh_f = open(thresh_fn,'r')
    thresh_lines = thresh_f.readlines()
    for line in thresh_lines:
        line_parts = line.strip().split()
        thresh = float(line_parts[2])
        if line_parts[0] not in thresh_dict:
            thresh_dict[line_parts[0]] = {}
        thresh_dict[line_parts[0]][line_parts[1]] = thresh
    thresh_f.close()

    # open up files for writing
    bars_fn = write_dir + "/bars.file"
    bars_f =  open(bars_fn,"w")
    bar_header = "\t"
    obs_states = migration_matrix[migration_matrix.keys()[0]].keys()
    for state in obs_states:
        bar_header += state + "\t"
    bars_f.write(bar_header + "\n")

    cluster_fn = write_dir + "/cluster.file"
    cluster_f =  open(cluster_fn,"w")
    prune_fn = write_dir + "/prune.file"
    prune_f =  open(prune_fn,"w")
    cdist_fn = write_dir + "/cdist.file"
    cdist_f =  open(cdist_fn,"w")

lik_fn = write_dir + "/lik.file"
lik_f =  open(lik_fn,"w")

# build the tree #
print "Building tree"
tree_f = open(tree_fn,"r")
tree_string = tree_f.read().strip()
tree = multitree.multitree()
tree.build(tree_string)
tree_f.close()

# root it and remove zero branch lengths
min_len = min([b.length for b in tree.branch_list if b.length > 0.0])
for b in tree.branch_list:
    if b.length <= 0.0:
        b.length = min_len
    names_1 = b.ends[0].name_dict[b.ends[1]]
    names_2 = b.ends[1].name_dict[b.ends[0]]
    if names_1 == outgroup or names_2 == outgroup:
        tree.rootify(b)

# write the data files
full_fn = write_dir + "/full.file"
full_f =  open(full_fn,"w")
migration_fn = write_dir + "/habitat.file"
migration_f =  open(migration_fn,"w")

# write the labels
habitats = migration_matrix.keys()
label_line = "LABELS"
rings = color_hash.keys()
rings.sort()
for ring in rings:
    symbols = color_hash[ring].keys()
    symbols.sort()
    for symbol in symbols:
        label_line += "," + symbol
color_line = "COLORS"
for ring in rings:
    symbols = color_hash[ring].keys()
    symbols.sort()
    for symbol in symbols:
        color_line += "," + color_hash[ring][symbol]
full_f.write(label_line + "\n")
full_f.write(color_line + "\n")

if thresh_fn is not None:
    prune_f.write(label_line + "\n")
    prune_f.write(color_line + "\n")
    cluster_f.write(label_line + "\n")
    cluster_f.write(color_line + "\n")

# embed important variables into tree
tree.migration_matrix = migration_matrix
tree.mu = mu
tree.states = obs_states
tree.color_hash = color_hash
tree.thresh_dict = thresh_dict
tree.color_hash = color_hash

# wipe the likelihoods off of the tree (just to be sure)
ML.TreeWipe(tree.a_node)

# learn the likelihoods
print "Learn habitat assignments"
ML.LearnLiks(tree,mu,migration_matrix,outgroup)

# estimate the states (by making each node trifurcating...)
root = tree.root
kids = root.GetKids()
if kids[0].name in outgroup:
    true_root = kids[1]
else:
    true_root = kids[0]
lik_score = ML.EstimateStates(true_root)
k = 1 + len(migration_matrix)*(len(migration_matrix.values()[0])-1)
aic_score = 2*k - 2*lik_score
lik_f.write("likelihood: " + str(lik_score) + "\n")
lik_f.write("aic: " + str(aic_score) + "\n")

# write out the habitat assignment of each leaf
for leaf in true_root.leaf_nodes:
    this_str = str(leaf) + "\t" + leaf.habitat
    migration_f.write(this_str + "\n")

# draw out the tree
files = {}
files['full'] = full_f

if thresh_fn is not None:
    files['cluster'] = cluster_f
    files['bar'] = bars_f
    files['prune'] = prune_f
    files['cdist'] = cdist_f

# figure out how far each leaf is from the outgroup:
if to_truncate:
    min_dist = 0.001
    limit_dist = 0.42
    leaf_dist_fn = write_dir + "/leaves.dist"
    leaf_dist_f =  open(leaf_dist_fn,"w")
    for leaf in tree.leaf_node_list:
        this_dist = leaf.DistTo(true_root,None,0)[1]
        leaf_dist_f.write(str(this_dist) + "\n")
        # truncate branches that are too long:
        if this_dist > limit_dist:
            leaf.TruncateDist(true_root,this_dist,min_dist,limit_dist)
    leaf_dist_f.close()

    # write out the tree (cheating a little to make the root branches
    # shorter)
    for branch in tree.root.child_branches:
        branch.length = min_dist

itol_fn = write_dir + "/itol.tree"
itol_f =  open(itol_fn,"w")
itol_f.write(true_root.treePrint("") + ";")
itol_f.close()

# print out the full file
print "Write out results"
true_root.FulliTol(files)

if thresh_fn is not None:

    strain_fn = write_dir + "/strain.names"
    strain_f = open(strain_fn,"w")

    # find the divergence points
    divergers = true_root.GetDivergencePoints()
    true_root.divergers = divergers

    true_root.ClusterTest(files,params)
    true_root.DrawSubclusters(None,cluster_f)
    tree.DrawLeaves(files)

    # print lists of constituents in each cluster
    cluster_roots = filter(lambda a: a.cluster_root,tree.node_dict.values())
    for cluster_root in cluster_roots:
        strains = cluster_root.leaf_nodes
        strain_f.write(str(len(strains)) + "\t" + str(cluster_root) +"\n")
        for strain in strains:
            strain_f.write("\t" + str(strain) + "\n")
        strain_f.write("\n")

    # retain significant clusters
    all_leaves = true_root.leaf_nodes
    for leaf in all_leaves:
        if not leaf.in_cluster:
            leaf.RemoveLeaf()
    true_root.UpPrune(files)

    # print pruned topology
    prune_topo_fn = write_dir + "/prune.tree"
    prune_topo_f =  open(prune_topo_fn,"w")
    prune_topo_f.write(true_root.treePrint("") + ";")
    prune_topo_f.close()

full_f.close()
migration_f.close()
if thresh_fn is not None:
    cluster_f.close()
    lik_f.close()
    prune_f.close()
    cdist_f.close()
    bars_f.close()
    strain_f.close()


# stick all the reconciliation-related code in here

import sys
import copy
import math
import pdb
from numpy import *
from numpy.matlib import *

def CheckConvergence(tree,old_score,params):

    thresh = float(params[1])

    total_prob = 0
    for i in tree.internal_node_list:
        for j in i.ML_probs:
            total_prob += sum(i.ML_probs[j].values())

    total_prob /= len(tree.internal_node_list)*len(tree.species_count)

    if abs(math.log(abs(old_score)) - math.log(abs(total_prob))) < thresh:
        return total_prob, 0
    else:
        return total_prob, 1
        
def EstimateStates(this_node):

    kids = this_node.GetKids()
    scores = {}
    habitats = this_node.ML_probs.values()[0].keys()

    for habitat in habitats:
        scores[habitat] = 0

    for kid in kids:
        for habitat in habitats:
            scores[habitat] += this_node.ML_probs[kid][habitat]

    max_score = -sys.maxint
    max_habitat = None

    for habitat in scores:
        if scores[habitat] > max_score:
            max_score = scores[habitat]
            max_habitat = habitat

    AssignState(this_node,habitats,max_habitat)

def AssignState(this_node,habitats,max_habitat):

    for habitat in habitats:
        if habitat is max_habitat:
            this_node.ML_state[habitat] = True
        else:
            this_node.ML_state[habitat] = False

    kids = this_node.GetKids()
    for kid in kids:
        kid_habitat = this_node.ML_path[kid][max_habitat]
        AssignState(kid,habitats,kid_habitat)

def LearnLiks(tree,mu,migration_matrix,true_root):

    true_kids = true_root.GetKids()
    for kid in true_kids:
        MLNode(true_root,kid,mu,migration_matrix)

def LearnShuffleLiks(tree,mu,migration_matrix,true_root):

    true_kids = true_root.GetKids()
    for kid in true_kids:
        MLShuffleNode(true_root,kid,mu,migration_matrix)

def TreeWipe(this_node):

    if len(this_node.true_state) < 1:
        this_node.true_state = this_node.ML_state.copy()

    this_node.ML_probs = {}
    this_node.ML_state = {}

    for kid in this_node.GetKids():
        TreeWipe(kid)

def MLNode(this_node,kid_node,mu,migration_matrix):

    # no need to recurse if you've already solved this problem
    if kid_node in this_node.ML_probs:
        return this_node.ML_probs[kid_node]

    # if not ...
    this_node.ML_probs[kid_node] = {}
    this_node.ML_path[kid_node] = {}

    # stop condition
    if kid_node is None:
        this_filter = this_node.species
        for this_habitat in migration_matrix:
            this_val = migration_matrix[this_habitat][this_filter]
            this_log = math.log(this_val)
            this_node.ML_probs[kid_node][this_habitat] = this_log
        return

    # if not a leaf, need to obtain probabilities of child values:
    child_nodes = filter(lambda i: i is not this_node, kid_node.other_nodes)
    if len(child_nodes) < 1:
        child_nodes = [None,None]
    MLNode(kid_node,child_nodes[0],mu,migration_matrix)
    MLNode(kid_node,child_nodes[1],mu,migration_matrix)    

    t = None
    for branch in this_node.branch_list:
        if kid_node in branch.ends:
            t = branch.length

    hab_num = float(len(migration_matrix))
    p_xx = 1.0/hab_num + (hab_num-1)/hab_num*math.exp(-hab_num*mu*t)
    p_xy = 1.0/hab_num*(1-math.exp(-hab_num*mu*t)) 

    habitats = migration_matrix.keys()

    # for each of the current nodes possible states
    for this_habitat in habitats:
        sum_prob = math.exp(-500)
        probs = {}
        
        for kid_habitat in habitats:
            
            probs[kid_habitat] = 0

            if this_habitat is kid_habitat:
                trans_prob = p_xx
            else:
                trans_prob = p_xy

            if trans_prob <= 0:
                pdb.set_trace()

            probs[kid_habitat] += math.log(trans_prob)

            # handle case at leaves
            if child_nodes == [None,None]:
                prob = kid_node.ML_probs[None][kid_habitat]
                probs[kid_habitat] += prob
            else:
                prob1 = kid_node.ML_probs[child_nodes[0]][kid_habitat]
                prob2 = kid_node.ML_probs[child_nodes[1]][kid_habitat]
                probs[kid_habitat] += prob1 + prob2

        high_prob = float(-sys.maxint)
        max_hab = None
        for a_habitat in probs:
            if probs[a_habitat] > high_prob:
                high_prob = probs[a_habitat]
                max_hab = a_habitat

        this_node.ML_probs[kid_node][this_habitat] = high_prob
        this_node.ML_path[kid_node][this_habitat] = max_hab

    return

def MLShuffleNode(this_node,kid_node,mu,migration_matrix):

    # no need to recurse if you've already solved this problem
    if kid_node in this_node.ML_probs:
        return this_node.ML_probs[kid_node]

    # if not ...
    this_node.ML_probs[kid_node] = {}
    this_node.ML_path[kid_node] = {}

    true_habitat = None
    for habitat in this_node.true_state:
        if this_node.true_state[habitat]:
            true_habitat = habitat

    # stop condition
    if kid_node is None:

        perturb_filter = this_node.perturb_species
        this_val = migration_matrix[true_habitat][perturb_filter]
        this_log = math.log(this_val)
        this_node.ML_probs[kid_node][true_habitat] = this_log

        return

    # if not a leaf, need to obtain probabilities of child values:
    child_nodes = filter(lambda i: i is not this_node, kid_node.other_nodes)
    if len(child_nodes) < 1:
        child_nodes = [None,None]
    MLShuffleNode(kid_node,child_nodes[0],mu,migration_matrix)
    MLShuffleNode(kid_node,child_nodes[1],mu,migration_matrix)    

    t = None
    for branch in this_node.branch_list:
        if kid_node in branch.ends:
            t = branch.length

    hab_num = float(len(migration_matrix))
    p_xx = 1.0/hab_num + (hab_num-1)/hab_num*math.exp(-hab_num*mu*t)
    p_xy = 1.0/hab_num*(1-math.exp(-hab_num*mu*t)) 

    # what's the kids habitat
    kid_habitat = None
    for habitat in kid_node.true_state:
        if kid_node.true_state[habitat]:
            kid_habitat = habitat

    # for each of the current nodes possible states
    prob = 0

    if kid_habitat is true_habitat:
        prob += math.log(p_xx)
    else:
        prob += math.log(p_xy)

    # if child
    if child_nodes == [None,None]:
        prob += kid_node.ML_probs[None][kid_habitat]
    else:
        prob += kid_node.ML_probs[child_nodes[0]][kid_habitat]
        prob += kid_node.ML_probs[child_nodes[1]][kid_habitat]

    this_node.ML_probs[kid_node][true_habitat] = prob

    return


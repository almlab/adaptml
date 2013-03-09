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
    return max_score


def AssignState(this_node,habitats,max_habitat):

    for habitat in habitats:
        if habitat is max_habitat:
            this_node.ML_state[habitat] = True
        else:
            this_node.ML_state[habitat] = False

    this_node.habitat = max_habitat
    this_node.orig_lik = this_node.GetLikelihood()
    this_node.orig_thresh = this_node.GetThresh()
    this_node.likelihood = this_node.GetLikelihood()
    this_node.threshold = this_node.GetThresh()    

    kids = this_node.GetKids()
    for kid in kids:
        kid_habitat = this_node.ML_path[kid][max_habitat]
        AssignState(kid,habitats,kid_habitat)

def LearnLiks(tree,mu,migration_matrix,outgroup):

    root = tree.root
    kids = root.GetKids()

    # figure out which kid isn't the outgroup
    if kids[0].name in outgroup:
        true_root = kids[1]
    else:
        true_root = kids[0]
        
    true_kids = true_root.GetKids()
    for kid in true_kids:
        MLNode(true_root,kid,mu,migration_matrix)

def TreeWipe(this_node):

    this_node.old_state = this_node.ML_state.copy()

    this_node.ML_probs = {}
    this_node.ML_state = {}

    for i in this_node.other_nodes:
        if len(i.ML_state) > 0:
            TreeWipe(i)

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

        probs = {}
        for kid_habitat in habitats:
            
            probs[kid_habitat] = 0

            if this_habitat is kid_habitat:
                trans_prob = p_xx
            else:
                trans_prob = p_xy

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

def MLNullNode(this_node,kid_node,eigs,MR_matrix):

    # no need to recurse if you've already solved this problem
    if kid_node in this_node.ML_null_probs:
        return this_node.ML_null_probs[kid_node]

    # if not ...
    this_node.ML_null_probs[kid_node] = {}
    species_key = eigs["key"]

    # stop condition
    if kid_node is None:

        spec_count = this_node.tree.species_count

        prob_dict = {}
        for species in spec_count:
            prob_dict[species] = 0
        
        for from_species in spec_count:
            scaling = float(spec_count[from_species])
            scaling /= sum(spec_count.values())
            for to_species in spec_count:
                s1 = species_key[from_species]
                s2 = species_key[to_species]
                scale_prob = scaling
                scale_prob *= MR_matrix[s1,s2]
                prob_dict[to_species] += scale_prob

        for species in prob_dict:
            this_prob = math.log(prob_dict[species])
            this_node.ML_null_probs[kid_node][species] = this_prob

        return this_node.ML_null_probs[kid_node]

    # if not a leaf, need to obtain probabilities of child values:
    child_nodes = filter(lambda i: i is not this_node, kid_node.other_nodes)
    if len(child_nodes) < 1:
        child_nodes = [None,None]
    MLNullNode(kid_node,child_nodes[0],eigs,MR_matrix)
    MLNullNode(kid_node,child_nodes[1],eigs,MR_matrix)    

    # now, work out transition probability matrix:
    if kid_node in this_node.branch_dict:
        branch_len = this_node.branch_dict[kid_node].length
    elif this_node in kid_node.branch_dict:
        branch_len = kid_node.branch_dict[this_node].length

    fst_matrix = eigs["vectors"]
    sec_matrix = diag(exp(eigs["values"]*branch_len))
    trd_matrix = eigs["inverse"]
    trans_matrix = real(fst_matrix*sec_matrix*trd_matrix)

    # for each of the current nodes possible states
    for this_species in this_node.tree.species_count:
        sum_prob = math.exp(-500)
        for kid_species in this_node.tree.species_count:

            prob = 0
            trans_prob = trans_matrix[species_key[kid_species],species_key[this_species]]
            if trans_prob < 0:
                trans_prob = math.exp(-500)
            
            prob += math.log(trans_prob)

            # handle case at leaves
            if child_nodes == [None,None]:
                prob +=  kid_node.ML_null_probs[None][kid_species]
            else:
                prob += kid_node.ML_null_probs[child_nodes[0]][kid_species]
                prob += kid_node.ML_null_probs[child_nodes[1]][kid_species]
            sum_prob += math.exp(prob)

        #########################################
        # use pseudocounts to prevent underflow #
        #########################################
        under_thresh = -150
        pseudocount = this_node.tree.species_count[this_species]
        sum_prob += pseudocount*exp(under_thresh)

        this_node.ML_null_probs[kid_node][this_species] = math.log(sum_prob)

    ############################
    ### SITE OF BIG PROBLEM? ###
    ############################

    """
    # rescale
    scale_num = sum(exp(array(this_node.ML_null_probs[kid_node].values())))
    for k in this_node.ML_null_probs[kid_node]:
        scale_prob = math.exp(this_node.ML_null_probs[kid_node][k])/scale_num
        this_node.ML_null_probs[kid_node][k] = math.log(scale_prob)
    """
    return

def GetTransitions(this_branch,parent,eigs,species_key):

    node1 = this_branch.ends[0]
    node2 = this_branch.ends[1]    
    
    # propogate
    if parent is node1:
        child = node2
    else:
        child = node1
    for i in child.child_branches:
        GetTransitions(i,child,eigs,species_key)

    # stop condition
    if node1.isLeaf() or node2.isLeaf():
        return
    if len(node1.ML_probs) < 1:
        return
    if len(node2.ML_probs) < 1:
        return

    # compute transition probability matrix
    branch_factor = 10
    branch_len = this_branch.length
    fst_matrix = eigs["vectors"]
    sec_matrix = diag(exp(eigs["values"]*branch_len*branch_factor))
    trd_matrix = eigs["inverse"]
    trans_matrix = real(fst_matrix*sec_matrix*trd_matrix)

    # now, get the priors
    keys = node1.ML_probs[node2].keys()
    p_child = {}
    p_paren = {}
    for i in keys:
        p_child[i] = 0
        p_paren[i] = 0        
    for nodes in child.ML_probs:
        if nodes is not parent:
            for state in keys:
                p_child[state] += child.ML_probs[nodes][state]
    for nodes in parent.ML_probs:
        if nodes is not child:
            for state in keys:
                p_paren[state] += parent.ML_probs[nodes][state]

    prob_matrix = matrix(ones((len(species_key),len(species_key))))
    for i in species_key:
        for j in species_key:
            child_ind = species_key[i]
            paren_ind = species_key[j]
            this_prob = 1
            
            this_prob *= trans_matrix[paren_ind,child_ind]
            this_prob *= math.exp(p_child[i])
            this_prob *= math.exp(p_paren[j])

            prob_matrix[child_ind,paren_ind] = this_prob

    """
    for i in range(len(species_key)):
        this_sum = sum(prob_matrix[i,:])
        for j in range(len(species_key)):
            prob_matrix[i,j] /= this_sum
    """

    keep_state = sum(diag(prob_matrix))
    swap_state = sum(prob_matrix) - keep_state

    child.confidence = math.log(keep_state) - math.log(swap_state)

    """
    foo = child.tree.node_dict["F_FAL492"]
    bar = child.tree.node_dict["F_FAL45"]
    if foo in parent.leaf_nodes and foo not in child.leaf_nodes:
        if bar in child.leaf_nodes:
            pdb.set_trace()
    """
    
    

# stick all the reconciliation-related code in here

import sys
import copy
import math
import pdb

from numpy import *
from numpy.matlib import *
from scipy.optimize import fminbound

def CheckConverge(tree,new_migrate,this_migrate):

    total_prob = 0
    for i in tree.internal_node_list:
        for j in i.ML_probs:
            total_prob += sum(i.ML_probs[j].values())

    total_prob /= len(tree.internal_node_list)*len(i.ML_probs)

    diff = 0
    for habitat in this_migrate:
        for filter in this_migrate[habitat]:
            p1 = this_migrate[habitat][filter]
            p2 = new_migrate[habitat][filter]            
            diff += abs(p1-p2)

    diff /= len(this_migrate)*len(this_migrate[habitat])

    return total_prob, diff
        
def EstimateStates(this_node,migration_matrix):

    habitats = migration_matrix.keys()

    if len(this_node.ML_probs) == 1:

        this_filter = this_node.species
        for habitat in habitats:
            this_prob = migration_matrix[habitat][this_filter]
            this_node.ML_state[habitat] = this_prob

        overall_lik = sum(this_node.ML_state.values())
        for habitat in this_node.ML_state:
            this_node.ML_state[habitat] /= overall_lik

    else:

        # establish the state dict
        states = {}
        for habitat in habitats:
            states[habitat] = 0

        # now, sum log likelihoods (consider trifurcating)
        for kid_node in this_node.ML_probs:
            for habitat in this_node.ML_probs[kid_node]:
                states[habitat] += this_node.ML_probs[kid_node][habitat]

        # scale log liks and convert to probs
        sum_prob = sum(exp(array(states.values())))

        min_float =  math.exp(-745)
        if sum_prob < min_float:
            sum_prob = min_float

        this_node.sum_lik = sum_prob

        for species in states:
            this_node.ML_state[species] = math.exp(states[species])/sum_prob

        # print
        # print this_node
        # print this_node.ML_state

        #if "1_1" in this_node.name:
        #    if "F_1" in this_node.name:
        #        pdb.set_trace()

    for kid_node in this_node.other_nodes:
        if len(kid_node.ML_state) < 1:
            EstimateStates(kid_node,migration_matrix)

def LearnMR(tree,migration_matrix,mu):

    new_migration_matrix = {}
    for i in migration_matrix:
        new_migration_matrix[i] = {}
        for j in migration_matrix[i]:
            new_migration_matrix[i][j] = math.exp(-24)

    hab_num = float(len(migration_matrix))
    
    for leaf in tree.leaf_node_list:

        # identify the parent
        parent_node = None
        for j in leaf.branch_list[0].ends:
            if j is not leaf:
                parent_node = j
        t = leaf.branch_list[0].length

        this_filter = leaf.species

        """
        for habitat in parent_node.ML_state:
            val = parent_node.ML_state[habitat]
            new_migration_matrix[habitat][this_filter] += val
        """

        p_xx = 1.0/hab_num + (hab_num-1)/hab_num*math.exp(-hab_num*mu*t)
        p_xy = 1.0/hab_num*(1-math.exp(-hab_num*mu*t)) 

        for parent_habitat in parent_node.ML_state:
            for leaf_habitat in leaf.ML_state:
                trans_prob = None
                if parent_habitat == leaf_habitat:
                    trans_prob = p_xx
                else:
                    trans_prob = p_xy
                prob = trans_prob*parent_node.ML_state[parent_habitat]
                new_migration_matrix[leaf_habitat][this_filter] += prob

    for this_habitat in new_migration_matrix:
        scale = sum(new_migration_matrix[this_habitat].values())
        for this_filter in new_migration_matrix[this_habitat]:
            new_migration_matrix[this_habitat][this_filter] /= scale

    return new_migration_matrix

def LearnRates(tree,mu,migration_matrix,rateopt):

    # update the migration matrix
    migration_matrix_opt = LearnMR(tree,migration_matrix,mu)

    # update the transition rate
    mu_0 = NodeRate(tree)

    if rateopt == 'avg':
        return mu_0, migration_matrix_opt
    elif rateopt == 'num':
        # optimize rate numerically
        mu_opt = fminbound(TestMu,0,mu_0*3,args=(tree,migration_matrix),xtol=1e-2)
        if mu_opt > mu_0*2.75:
            print "look into fminbound bounds"
            sys.exit(1)

        # put things back the way they used to be (need to do, since
        # there's a treewipe inside of testmu)
        TreeWipe(tree)
        LearnLiks(tree,mu,migration_matrix)
        EstimateStates(tree.a_node,migration_matrix)

        return mu_opt, migration_matrix_opt

def TestMu(mu,tree,migration_matrix):

    TreeWipe(tree)
    LearnLiks(tree,mu,migration_matrix)
    score, diff = CheckConverge(tree,migration_matrix,migration_matrix)

    # since you're trying to minimize, return inverse ...

    return -score

def LearnRootedRates(tree,eigs,Q_matrix):

    species_key = eigs["key"]

    Q_zero = matrix(zeros((len(species_key),len(species_key))))
    Q, total_array = NodeRootedRate(tree.root,eigs)

    for i in range(len(total_array)):
        Q[i,:] /= total_array[i]

    return Q

def NodeRate(tree):

    # think of a really small number
    eps = math.exp(-10)
    same_count = 0
    total_dist = 0
    counter = 0

    for branch in tree.branch_list:

        # maybe we shouldn't include leaf nodes (they probably make it
        # look like there are more changes than there ought to be)

        node1 = branch.ends[0]
        node2 = branch.ends[1]

        if len(node1.ML_probs) == 1:
            continue
        if len(node2.ML_probs) == 1:
            continue

        m1 = node1.ML_state
        m2 = node2.ML_state

        max1 = 0
        max1_state = None
        max2 = 0
        max2_state = None

        for i in m1:
            if m1[i] > max1:
                max1_state = i
                max1 = m1[i]
        
        for i in m2:
            if m2[i] > max2:
                max2_state = i
                max2 = m2[i]

        if max1_state == max2_state:
            same_count += 1

        counter += 1
        total_dist += branch.length

    hab_num = float(len(m1))
    t = total_dist/counter
    dist = 1 - float(same_count)/counter

    inner = 1 - dist*hab_num/(hab_num-1)

    this_mu = -math.log(inner) / (hab_num * t)

    return this_mu

def LearnLiks(tree,mu,migration_matrix):

    # compute likelihoods for all internal nodes, keyed by child nodes 
    for node in tree.internal_node_list:
        for kid in node.other_nodes:
            MLNode(node,kid,mu,migration_matrix)

def TreeWipe(tree):

    for i in tree.internal_node_list:
        i.ML_probs = {}
        i.ML_state = {}
    for i in tree.leaf_node_list:
        i.ML_probs = {}
        i.ML_state = {}

    """
    this_node.old_state = this_node.ML_state.copy()
    this_node.ML_probs = {}
    this_node.ML_state = {}
    for i in this_node.other_nodes:
        if len(i.ML_state) > 0:
            TreeWipe(i)
    """
        
def MLNode(this_node,kid_node,mu,migration_matrix):

    # no need to recurse if you've already solved this problem
    if kid_node in this_node.ML_probs:
        return this_node.ML_probs[kid_node]

    # if not ...
    this_node.ML_probs[kid_node] = {}

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

    # now, work out transition probabilities
    if kid_node in this_node.branch_dict:
        t = this_node.branch_dict[kid_node].length
    elif this_node in kid_node.branch_dict:
        t = kid_node.branch_dict[this_node].length

    hab_num = float(len(migration_matrix))
    p_xx = 1.0/hab_num + (hab_num-1)/hab_num*math.exp(-hab_num*mu*t)
    p_xy = 1.0/hab_num*(1-math.exp(-hab_num*mu*t)) 

    habitats = migration_matrix.keys()

    # for each of the current nodes possible states
    for this_habitat in habitats:
        sum_prob = math.exp(-500)
        for kid_habitat in habitats:

            if this_habitat is kid_habitat:
                trans_prob = p_xx
            else:
                trans_prob = p_xy

            if trans_prob <= 0:
                pdb.set_trace()

            prob = 0
            prob += math.log(trans_prob)

            # handle case at leaves
            if child_nodes == [None,None]:
                prob +=  kid_node.ML_probs[None][kid_habitat]
            else:
                prob += kid_node.ML_probs[child_nodes[0]][kid_habitat]
                prob += kid_node.ML_probs[child_nodes[1]][kid_habitat]
            sum_prob += math.exp(prob)

        this_node.ML_probs[kid_node][this_habitat] = math.log(sum_prob)

    return

    
    


import random
import re
import copy
import sys
import pdb

from numpy import *
from sets import Set

import multitree
import branch

class node:

    def __init__(self,name,arbre):

        self.name_backup = name # argh, i need to refactor
        self.name = name

        # assign species names by stripping '.X' or '_X'
        match_str = re.compile(r'[\._]\w*-')
        name = match_str.sub('-',name,0)
        # catch trailers
        match_str = re.compile(r'[\._]\w*')
        species = match_str.sub('',name,0)

        self.species = species
        self.perturb_species = species
	self.range = []
	self.range.append(self.species)
        self.tree = arbre
	self.lca_scores = dict()        # keep track of lca scores
	self.event_str = dict()         # keep track of events for each lca
	    
	# clear all of the other variables
        self.myBranch = None            # used to construct unrooted trees
	self.branch_list = []           # used to construct unrooted trees
	self.child_branches = []        # list of branches to kids
	self.parent_branch = []         # should have length 1
        
	self.subnodes = dict()          # all nodes below this
	self.leaves = None              # leaves below this
        self.leaf_nodes = []
        
        # attach to each node a serial number, to tell nodes apart
	#self.serial = random.randint(-sys.maxint,sys.maxint)
        self.serial = hash(self.name)
        self.newick_string = ""
        self.counts_dict = {}
        self.ancestors = []

	self.visited = False            # useful for rooting trees
        self.dup_count = 0              # number of duplications at this node

        # things for GetNodeLinkDict
        self.link_dict_visited = False

        # things necessary for UnrootedLeaving
        self.leaf_dict = {}
        self.name_dict = {}        
        self.branch_dict = {}
        self.unrooted_leaving_visited = False
        self.other_nodes = []

        # stuff for global reconcile:
        self.lca_lookups = {}  # store here 3 lookup tables, 1 for
                               # each of the possible parents of this node

        # ML things
        self.ML_probs = {}
        self.ML_path = {}
        self.ML_null_probs = {}        
        self.ML_state = {}
        self.old_state = array([0,0,0])
        self.sum_lik = 0

        # cluster things
        self.confidence = math.exp(100)
        self.cluster_root = False
        self.in_cluster = False
        self.is_rep = False
        self.prune_leaves = []
        self.max_filter = None
        self.habitat = None
        self.lik_diff = None
        self.thresh_diff = None
        self.orig_lik = None
        self.orig_thresh = None
        self.likelihood = None
        self.threshold = None
        self.divergers = None
        self.div_tested = False

    def __repr__(self):
	#if self.myBranch == None:
	#    return self.name
	#else:	
	#    return self.name + ":" + str(self.myBranch.length)
        return self.name 
    
    def __eq__(self,other):
	if self is None or other is None:
	    return False
	elif self.serial == other.serial:
	    return True
	else:
	    return False

    def __ne__(self,other):

	if self is None or other is None:
	    return True
	elif self.serial == other.serial:
	    return False
	else:
	    return True

    def __gt__(self,other):
	if self.name > other.name:
	    return True
	else:
	    return False

    def __lt__(self,other):
	if self.name < other.name:
	    return True
	else:
	    return False
    def __hash__(self):
        return self.serial

    def isLeaf(self):
	if len(self.child_branches) < 1:
	    return True
	else:
	    return False

    def GetKids(self):
        kids = []
        for i in self.child_branches:
            for j in i.ends:
                if j is not self:
                    kids.append(j)
        return kids
	
    def PrintStates(this_node):

        if len(this_node.ML_state) > 0:
            if max(this_node.ML_state.values()) > 0:

                print "node:",
                node_name = this_node.treePrint("",0)
                if node_name[0] is not "(":
                    print "(" + node_name + ")"
                else:
                    print node_name
                print "prob:",
                for state in this_node.ML_state:
                    print state + ":", round(this_node.ML_state[state],3), 
                print

                print "lik:", this_node.sum_lik
                
                #for state in this_node.ML_probs:
                #    print state + ":", round(this_node.ML_probs[state],3), 
                print

        for branch in this_node.child_branches:
            for kid_node in branch.ends:
                if kid_node is not this_node:
                    kid_node.PrintStates()

    def GetParent(this_node):

        parent = None

        if this_node.tree.root is this_node:
            return parent
        
        if len(this_node.parent_branch) > 0:
            for i in this_node.parent_branch[0].ends:
                if i is not this_node:
                    parent = i
        return parent
	
    # recursively print out nodes
    def BootPrint(self,newickString,verbose=1):

	count = 0

	for kid_branches in self.child_branches:
	    for kid_nodes in kid_branches.ends:
		if kid_nodes is not self:
		    if count == 0:
			newickString += "("
			newickString = kid_nodes.BootPrint(newickString,verbose)
		    else:
			newickString += ","
			newickString = kid_nodes.BootPrint(newickString,verbose)
			newickString += ")"
		    count += 1

        # add ability to print state probabilities
        boot_str = ""
        if len(self.child_branches) > 1:
            if len(self.ML_state) > 0:
                for i in self.ML_state:
                    boot_str += "|" + i + "|"
                    boot_str += "-" + str(round(self.ML_state[i],3))

        newickString += boot_str

        addString = ""
        if verbose:
            if len(self.parent_branch) > 0 and len(self.child_branches) == 0:
                addString = self.name + ":" + str(self.parent_branch[0].length) 
            elif len(self.parent_branch) > 0: 
                addString = ":" + str(self.parent_branch[0].length)
        else:
            if len(self.parent_branch) > 0 and len(self.child_branches) == 0:
                addString = self.species + ":" + str(self.parent_branch[0].length) 
            elif len(self.parent_branch) > 0: 
                addString = ":" + str(self.parent_branch[0].length)

	return newickString + addString
                

    def addBranch(self,length):

	self.myBranch = branch.branch(length)
	self.myBranch.addNode(self)

    # recursively print out nodes
    def treePrint(self,newickString,verbose=1):

	count = 0

	for kid_branches in self.child_branches:
	    for kid_nodes in kid_branches.ends:
		if kid_nodes is not self:
		    if count == 0:
			newickString += "("
			newickString = kid_nodes.treePrint(newickString,verbose)
		    else:
			newickString += ","
			newickString = kid_nodes.treePrint(newickString,verbose)
			newickString += ")"
		    count += 1

        addString = ""

        if verbose:
            if len(self.parent_branch) > 0 and len(self.child_branches) == 0:
                addString = self.name + ":" + str(self.parent_branch[0].length) 
            elif len(self.parent_branch) > 0: 
                addString = ":" + str(self.parent_branch[0].length)
        else:
            if len(self.parent_branch) > 0 and len(self.child_branches) == 0:
                addString = self.name 
	return newickString + addString

    # join two nodes in a central node
    def unite(self,node2):

	# create a new node
	new_node_name = self.name + "-" + node2.name
	center_node = node(new_node_name,self.tree)
	self.myBranch.addNode(center_node)
	node2.myBranch.addNode(center_node)
	return center_node

    # recursively impose a hierarchy
    def imposeHierarchy(self):

	# find all the child nodes that have yet to be visited and
	# assign their children
	for kid_branch in self.child_branches:
	    for node in kid_branch.ends:
		if not node.visited and node is not self:
		    node.visited = True
		    node.parent_branch.append(kid_branch)
		    for child_branch in node.branch_list:
			if child_branch is not kid_branch:
			    node.child_branches.append(child_branch)
			    node.imposeHierarchy()

    # recursively label the roots of subtrees w/ the leaves contained
    # below
    def subtreeLabel(self):

	leaf_vec = []
	
	if self.leaves is not None:
	    return self.leaves

	# recurse
        kids = self.GetKids()
        for kid in kids:
            child_leaves = kid.subtreeLabel()
            leaf_vec.extend(child_leaves)

	# once you hit a leaf
	if len(self.child_branches) == 0:
	    leaf_vec.append(self.species)
            self.leaf_nodes.append(self)
            self.prune_leaves.append(self)
        else:
            for kid in kids:
                self.leaf_nodes.extend(kid.leaf_nodes)
                self.prune_leaves.extend(kid.prune_leaves)
                
	# non-duplicates:
	leaf_vec = dict(map(lambda i: (i,1),leaf_vec))
	self.leaves = leaf_vec
	return leaf_vec.keys()

    # figure out which species tree node a gene tree node maps to
    def subtreeMap(self,species_node):
	
	# do the children of the current species node possess the
	# relevant genes?  if so, follow that child.  if not, return
	# the current node

	foundSet = 0
	for kid_branches in species_node.child_branches:
	    for kid_nodes in kid_branches.ends:
		if kid_nodes is not species_node:

		    s_leaves = kid_nodes.leaves
		    g_leaves = self.leaves

		    # count how many elements of the gene set of
		    # leaves are not in the species set of leaves

		    n = 0
		    for i in g_leaves.keys():
			if not i in s_leaves:
			    n += 1
			    break
		    
		    # if you can see somewhere to descend, follow it
		    if n == 0:
			return self.subtreeMap(kid_nodes)

	# once you've run out of places to descend, just return
	# wherever you've ended up:
	return species_node

    # recursively store at each node a hash of all nodes below
    def Find_Subnodes(self):

	for kid_branches in self.child_branches:
	    for kid_nodes in kid_branches.ends:
		if kid_nodes is not self:
		    new_dict = kid_nodes.Find_Subnodes()
		    for k in new_dict.keys():
			if not k in self.subnodes:
			    self.subnodes[k] = 1
		    # self.subnodes.extend(kid_nodes.Find_Subnodes())

	# self.subnodes.append(self)
	self.subnodes[self.species] = 1
	return self.subnodes

    # find the last common ancestor of two nodes.  will descend from
    # the subroot (self) looking to see which children possess both
    # nodes.  if neither children possess both nodes, return current
    # node

    def Find_LCA(self,node1,node2):

        for kid_branches in self.child_branches:
            for kid_nodes in kid_branches.ends:
                if kid_nodes is not self:
                    if node1.species in kid_nodes.subnodes:
                        if node2.species in kid_nodes.subnodes:
                            return kid_nodes.Find_LCA(node1,node2)
        return self

    # get list of all vertices in a tree
    def Fill_Node_Dict(vertex):
	
	# vertex.tree.node_dict[vertex.name] = vertex
        vertex.tree.node_dict[vertex.name] = vertex

	for kid_branches in vertex.child_branches:
	    for kid_nodes in kid_branches.ends:
		if kid_nodes is not vertex:
		    kid_nodes.Fill_Node_Dict()

    # initialize hash tables for lca scores for each node
    def Init_LCA_Scores(vertex,node_dict):

	for tnode in node_dict:
	    vertex.lca_scores[tnode.name] = node.maxInt

	for kid_branches in vertex.child_branches:
	    for kid_nodes in kid_branches.ends:
		if kid_nodes is not vertex:
		    kid_nodes.Init_LCA_Scores(node_dict)

    # determine if one node is descended from the other ...
    def Are_Related(node1,node2):

	verdict = False
	for e in node1.leaves.keys():
	    if e in node2.leaves:
		verdict = True

	return verdict

    # method to list all the leaves of this node, keyed by the parent
    def UnrootedLeaving(this_node):

        for branches in this_node.branch_list:
            if branches not in this_node.tree.branch_list:
                this_node.tree.branch_list.append(branches)

        this_node.tree.internal_node_list.append(this_node)
        other_nodes = this_node.GetOtherNodes()

        # get leaves in each direction
        for parent_node in other_nodes:

            child_nodes = list(Set(other_nodes).difference(Set([parent_node])))
            this_node.GetChildLeaves(parent_node,child_nodes)

        # when done, move on to neighboring nodes:
        for node in other_nodes:

            if len(node.leaf_dict) is not len(node.branch_list):
                node.UnrootedLeaving()
                
    def GetLeaves(this_node,parent_node):

        # if the leaf dict has already been defined:
        if parent_node in this_node.leaf_dict:
            return this_node.leaf_dict[parent_node], this_node.name_dict[parent_node]

        # if is a leaf
        if len(this_node.branch_list) < 2:
            
            this_node.tree.leaf_node_list.append(this_node)

            if this_node.species in this_node.tree.species_count:
                this_node.tree.species_count[this_node.species] += 1
            else:
                this_node.tree.species_count[this_node.species] = 1

            this_node.leaf_dict[parent_node] = [this_node.species]
            this_node.name_dict[parent_node] = [this_node.name]          
            return this_node.leaf_dict[parent_node], this_node.name_dict[parent_node]

        other_nodes = this_node.GetOtherNodes()
        child_nodes = list(Set(other_nodes).difference(Set([parent_node])))
        this_node.GetChildLeaves(parent_node,child_nodes)
            
        return this_node.leaf_dict[parent_node], this_node.name_dict[parent_node]

    def GetChildLeaves(this_node,parent_node,other_nodes):

        merge_list = []
        name_list = []
        for kid_node in other_nodes:

            if parent_node in this_node.leaf_dict:
                if len(this_node.leaf_dict[parent_node]) > 1:
                    continue

            kid_leaves, kid_names = kid_node.GetLeaves(this_node)
            
            if type(kid_leaves[0]) is not type(''):
                merge_list.extend(kid_leaves[0])
                name_list.extend(kid_names[0])
            else:
                merge_list.extend(kid_leaves)
                name_list.extend(kid_names)

            # this is a total hack
            this_node.leaf_dict[parent_node] = []
            this_node.name_dict[parent_node] = []

        merge_list.sort()
        name_list.sort()
        this_node.leaf_dict[parent_node].extend(merge_list)
        this_node.name_dict[parent_node].extend(name_list)

    def GetNodeLinkDict(this_node,node_link_dict):

        this_node.link_dict_visited = True
        leaf_dict = this_node.leaf_dict

        for key in leaf_dict:
            subleaves = leaf_dict[key]

            if type(subleaves[0]) == type([]):
                subleaves = subleaves[0]
            sub_str = repr(subleaves)
            
            if sub_str in node_link_dict:
                node_link_dict[sub_str].append((this_node,key))
            else:
                node_link_dict[sub_str] = [(this_node,key)]
            # node_link_dict[sub_str] = (this_node,key)

        # after that's been done, decide where to recurse
        for relative in leaf_dict:

            relative_keys = relative.leaf_dict.keys()
            subleaves = relative.leaf_dict[relative_keys[0]]
            subleaves.sort()

            if not relative.link_dict_visited:
                relative.GetNodeLinkDict(node_link_dict)

    # return a list of all bordering nodes
    def GetOtherNodes(this_node):

        # fill out the leaf dictionary
        other_nodes = []
        for branch in this_node.branch_list:
            for node in branch.ends:
                if node is not this_node:
                    other_nodes.append(node)
                    node.branch_dict[this_node] = branch
        this_node.other_nodes = other_nodes

        return other_nodes

    # get the distance between two nodes:
    # call as follows:
    # found, dist = node1.DistTo(node2)
    def DistTo(this_node,that_node,this_branch=None,dist=0):

        # what to do at the end
        if this_node is that_node:
            return True, dist

        # if not ...
        was_found = False
        found_dist = dist
        for i in this_node.branch_list:
            if i is not this_branch:
                for j in i.ends:
                    if j is not this_node:
                        found, new_dist = j.DistTo(that_node,i,dist+i.immutable_length)
                        if found is True:
                            was_found = found
                            found_dist = new_dist

        return was_found, found_dist
        
    def NodeWipe(this_node):

        for kid_branches in this_node.child_branches:
            for kid_nodes in kid_branches.ends:
                if kid_nodes is not this_node:
                    kid_nodes.NodeWipe()
                        
        this_node.child_branches = []
        this_node.parent_branch = []
        this_node.subnodes = {}
        this_node.leaves = None
        this_node.ML_state = {}
        this_node.ML_probs = {}

    def LearnMR(this_node,MR_matrix,species_key):

        if len(this_node.branches) < 1:
            pdb.set_trace()
            return
        for i in this_node.child_branches:
            for j in i.ends:
                if j is not this_node:
                    j.LearnMR(MR_matrix,species_key)
        
    def GetALeaf(this_node,branch_binary_choice):

        if len(this_node.child_branches) < 1:
            return this_node

        else:
            the_branch = this_node.child_branches[branch_binary_choice]
            for i in the_branch.ends:
                if i is not this_node:
                    return i.GetALeaf(branch_binary_choice)

    def TruncateDist(this_node,root,this_dist,min_dist,limit_dist):

        parent_dist = this_node.parent_branch[0].length

        # how much distance is necessary?
        nec_dist = this_dist - limit_dist

        # if all you need to do is shorten this branch
        if nec_dist < parent_dist:
            new_dist = parent_dist - nec_dist
            this_node.parent_branch[0].length = new_dist
        else:
            new_dist = min_dist
            this_dist -= parent_dist + min_dist
            this_node.parent_branch[0].length = new_dist
            parent = this_node.GetParent()
            parent.TruncateDist(root,this_dist,min_dist,limit_dist)
        return
        
    def RemoveLeaves(this_node,cluster_node):

        for leaf in cluster_node.leaf_nodes:
            if leaf not in this_node.leaf_nodes:
                pdb.set_trace()
            
            this_node.leaf_nodes.remove(leaf)
        parent = this_node.GetParent()
        if parent is not None:
            parent.RemoveLeaves(cluster_node)

    def GetSibling(this_node):

        sibling = None
        parent = this_node.GetParent()
        for kid in parent.GetKids():
            if kid is not this_node:
                sibling = kid
        return sibling

    def GetLikelihood(this_node):

        lik = 0
        hab = this_node.habitat
        if hab is None:
            return None
        for kid in this_node.GetKids():
            lik += this_node.ML_probs[kid][hab]
        return lik

    def GetThresh(this_node):

        thresh_dict = this_node.tree.thresh_dict
        leaf1 = this_node.GetALeaf(0).name
        leaf2 = this_node.GetALeaf(1).name

        # you're at a leaf
        if leaf1 == leaf2:
            return 0

        if thresh_dict is None:
            return None
        else:
            return thresh_dict[leaf1][leaf2]

    def SubtractLikelihood(this_node,lik_diff,thresh_diff):

        if this_node is this_node.tree.root:
            return

        # adjust likelihoods and thresholds
        this_node.likelihood -= lik_diff
        this_node.threshold -= thresh_diff
        parent_node = this_node.GetParent()
        if parent_node is not this_node.tree.root:
            parent_node.SubtractLikelihood(lik_diff,thresh_diff)

    # remove a leaf from a tree
    def RemoveLeaf(this_node):

        tree = this_node.tree

        if not this_node.isLeaf():
            print "this node isn't a leaf."
            sys.exit(1)

        parent_node = this_node.GetParent()
        kids = parent_node.GetKids()
        sib = filter(lambda i: i is not this_node,kids)[0]

        # what to do if parent is the root
        if tree.root is parent_node:
            tree.root = sib
            return

        # if not tied to root, slightly more complex ...
        grandparent = parent_node.GetParent()

        # build a new connecting branch
        branch_len = 0
        branch_len += sib.parent_branch[0].length
        branch_len += parent_node.parent_branch[0].length
        new_branch = branch.branch(branch_len)
        new_branch.addNode(sib)
        new_branch.addNode(grandparent)
        
        sib.parent_branch[0] = new_branch
        gp_branch = parent_node.parent_branch[0]
        
        for i in range(len(grandparent.child_branches)):
            if grandparent.child_branches[i] is gp_branch:
                grandparent.child_branches[i] = new_branch

        # remove the nodes from the tree's node-dictionary
        #del tree.node_dict[this_node.name]
        #del tree.node_dict[parent_node.name]

    def KeepOneLeaf(this_node):

        keep_leaves = []
        leaf_dict = {}

        tree = this_node.tree
        leaves = this_node.leaf_nodes
        f_leaves = filter(lambda i: not i.is_rep,leaves)
        k_leaves = filter(lambda i: i.is_rep,leaves)

        # get filter counts
        for filt in this_node.tree.species_count.keys():
            leaf_dict[filt] = 0
        for leaf in f_leaves:
            leaf_dict[leaf.species] += 1

        if len(k_leaves) > 0:
            keep_leaves.append(k_leaves[0])
        
        # keep the eldest leaf (that isn't a representative of a
        # sub-cluster)
        keep_leaf = f_leaves[0]

        for leaf_ind in range(1,len(f_leaves)):
            this_leaf = f_leaves[leaf_ind]
            this_leaf.RemoveLeaf()

        keep_leaf.is_rep = True
        keep_leaves.append(keep_leaf)
        parent = this_node.GetParent()
        parent.AddLeaf(keep_leaf)

        return keep_leaves, leaf_dict

    def DrawSubclusters(this_node,cluster_habitat,cluster_file):

        if this_node.isLeaf():
            return

        this_habitat = this_node.habitat
        this_filter = this_node.max_filter
        kids = this_node.GetKids()
        
        if this_node.cluster_root:
            for kid in kids:
                kid.DrawSubclusters(this_habitat,cluster_file)
            return

        if this_habitat is cluster_habitat:

            internal_line = this_node.PrintiTolInternal(10)
            cluster_file.write(internal_line + "\n")
            
        for kid in kids:
            kid.DrawSubclusters(cluster_habitat,cluster_file)

    def GetLikDiffs(this_node):

        lik_diff = 0
        thresh_diff = 0

        sis = this_node.GetSibling()
        dad = this_node.GetParent()
        uncle = dad.GetSibling()
        gdad = dad.GetParent()

        if gdad is this_node.tree.root:
            return 0, 0

        old_gdad_lik = gdad.likelihood
        new_gdad_lik = gdad.ML_probs[uncle][gdad.habitat]
        new_gdad_lik += gdad.ML_probs[dad][gdad.habitat]
        new_gdad_lik -= dad.ML_probs[this_node][dad.habitat]

        old_gdad_thresh = gdad.threshold
        new_gdad_thresh = uncle.threshold
        new_gdad_thresh += sis.threshold
        new_gdad_thresh += math.log(dad.GetTransProb(gdad))
        new_gdad_thresh += math.log(sis.GetTransProb(dad))
        new_gdad_thresh += math.log(uncle.GetTransProb(gdad))

        lik_diff = old_gdad_lik - new_gdad_lik
        thresh_diff = old_gdad_thresh - new_gdad_thresh

        # while you're at it, amend the likelihood of the parent
        dad.likelihood = dad.ML_probs[sis][dad.habitat]
        dad.threshold = sis.threshold + math.log(sis.GetTransProb(dad))

        return lik_diff, thresh_diff
        
    def GetTransProb(this_node,parent_node):

        if this_node.GetParent() is not parent_node:
            print "not parent passed"
            sys.exit(1)

        tree = this_node.tree
        mu = tree.mu
        h_num = float(len(tree.migration_matrix))
        t = this_node.parent_branch[0].length

        if this_node.habitat == parent_node.habitat:
            tp = 1.0/h_num + (h_num-1)/h_num*math.exp(-h_num*mu*t)
        else:
            tp = 1.0/h_num*(1-math.exp(-h_num*mu*t)) 

        return tp

    def PrintiTolLeaf(this_node):

        this_str = ""
        this_str += str(this_node)

        # states = this_node.tree.hidden_states
        color_hash = this_node.tree.color_hash
        this_state = this_node.species

        rings = color_hash.keys()
        rings.sort()
        for ring in rings:
            symbols = color_hash[ring].keys()
            symbols.sort()
            for symbol in symbols:
                this_str += ","
                if ring == "H":
                    this_str += str(0.0)
                else:
                    if symbol == this_state[int(ring)-1]:
                        this_str += str(10.0)
                    else:
                        this_str += str(0.0)
        return this_str

    def ParsePredictions(this_node):

        if this_node is this_node.tree.root:
            return

        habitat = this_node.habitat
        states = this_node.tree.states
        s_dict = this_node.tree.migration_matrix[habitat]
        rs_dict = dict(map(lambda i: (i[1],i[0]), s_dict.items()))
        keys = rs_dict.keys()
        keys.sort()
        keys.reverse()

        state1 = rs_dict[keys[0]]
        state2 = rs_dict[keys[1]]

        states = [state1,state2]
        vals = [keys[0], keys[1]]

        return states, keys

    def PrintiTolInternal(this_node,R_size,nodes=None):
        
        internal_line = ""
        if this_node is this_node.tree.root:
            return internal_line

        if nodes is None:
            # find nodes to define this subtree
            node1 = this_node.GetALeaf(0)
            node2 = this_node.GetALeaf(1)                
        else:
            node1 = nodes[0]
            if len(nodes) < 2:
                node2 = None
            else:
                node2 = nodes[1]

        color_hash = this_node.tree.color_hash
        habitat_choice = this_node.habitat.split(' ')[-1]

        # write out the results
        if node2 is not None:
            internal_line += str(node1) + "|"
            internal_line += str(node2)
        else:
            internal_line += str(node1)
        internal_line += ",R" + str(R_size)

        # add zeros for all the states
        rings = color_hash.keys()
        rings.sort()
        for ring in rings:
            symbols = color_hash[ring].keys()
            symbols.sort()
            for symbol in symbols:
                internal_line += ","
                if ring == 'H':
                    if symbol == habitat_choice:
                        internal_line += str(10.0)
                    else:
                        internal_line += str(0.0)
                else:
                    internal_line += str(0.0)
        return internal_line

    def GetDivergencePoints(this_node):

        divergers = []
        if this_node.isLeaf():
            return divergers

        this_hab = this_node.habitat
        for kid in this_node.GetKids():
            if kid.habitat != this_hab:
                divergers.append(kid)
            divergers.extend(kid.GetDivergencePoints())

        if this_hab != this_node.GetParent().habitat:
            this_node.divergers = divergers

        return divergers

    def ClusterTest(this_node,files,params):
        
        if this_node.isLeaf():
            this_node.div_tested = True
            return

        # make sure this node is a divergence point
        if this_node.divergers is None:
            print "no divergers in clustertest!"
            sys.exit(1)

        # test that you've checked out all divergence nodes below
        divergers = this_node.divergers

        for div_node in divergers:
            if not div_node.div_tested:
                div_node.ClusterTest(files,params)
                div_node.div_tested = True

        # see where the clusters are
        this_node.ClusterFind(this_node.habitat,files,params)

        return

    def ClusterFind(this_node,habitat,files,params):

        # return if you hit a leaf or another divergence node
        if len(this_node.leaf_nodes) < 1:
            return
        if this_node.habitat != habitat:
            return

        is_cluster = False

        # test first for homogeneity
        same_count = 0.0
        for leaf in this_node.leaf_nodes:
            if leaf.ML_state[this_node.habitat]:
                same_count += 1

        frac = same_count / len(this_node.leaf_nodes)
        if frac > 0.9:
            # now, test likelihood:
            lik = this_node.likelihood
            thresh = this_node.threshold

            # precision is about 8 points due to earlier
            # conversion to string
            lik_big = int(lik*math.pow(10,8))
            thresh_big = int(thresh*math.pow(10,8))
            if lik_big > thresh_big:
                is_cluster = True

        # if is a cluster
        if is_cluster:
            cluster_file = files['cluster']
            cdist_file = files['cdist']

            parent = this_node.GetParent()
            this_node.cluster_root = True

            if params['cdist']:
                # get the average distance between leaves in the cluster
                avg_dist = this_node.GetAverageClusterDist()
                cdist_file.write(str(avg_dist) + "\n")
            
            # remove leaf nodes from all upstream nodes
            parent.RemoveLeaves(this_node)

            # subtract likelihood
            gdad = parent.GetParent()
            if gdad is not None:
                lik_diff, thresh_diff = this_node.GetLikDiffs()
                gdad.SubtractLikelihood(lik_diff,thresh_diff)
            internal_line = this_node.PrintiTolInternal(50)
            cluster_file.write(internal_line + "\n")
            
        else:
            for kid in this_node.GetKids():
                kid.ClusterFind(habitat,files,params)

        return
            
    def UpPrune(this_node,files):

        tree = this_node.tree
        colors = tree.states
        prune_file = files['prune']
        bars_file = files['bar']
        
        # post-fix
        kids = this_node.GetKids()
        for kid in kids:
            kid.UpPrune(files)

        # find clusters
        if not this_node.cluster_root:
            return

        # remove all but one of the kids
        keep_leaves, leaf_dict = this_node.KeepOneLeaf()
        internal_line = this_node.PrintiTolInternal(10,keep_leaves)
        prune_file.write(internal_line + "\n")

        bar_str = ""
        for leaf in keep_leaves:
            bar_str += str(leaf) + "|"

        for state in this_node.tree.states:
            if state in leaf_dict:
                bar_str += "\t" + str(leaf_dict[state])
        bars_file.write(bar_str + "\n")

    # add a leaf back into a node's leaf nodes
    def AddLeaf(this_node,leaf):
        this_node.leaf_nodes.append(leaf)
        parent = this_node.GetParent()
        if parent is not None:
            parent.AddLeaf(leaf)
    
    # find the average distance between leaves in this cluster
    def GetAverageClusterDist(this_node):

        leaves = this_node.leaf_nodes
        total_dist = 0
        pair_count = 0
        for ind1 in range(len(leaves)):
            for ind2 in range(ind1+1,len(leaves)):
                dist = leaves[ind1].DistTo(leaves[ind2])[1]
                total_dist += dist
                pair_count += 1

        avg_dist = total_dist / pair_count
        return avg_dist

    def FulliTol(this_node,files):

        full_file = files['full']
        kids = this_node.GetKids()
        if len(kids) < 1:
            this_str = this_node.PrintiTolLeaf()
        else:
            this_str = this_node.PrintiTolInternal(10)
        full_file.write(this_str + "\n")

        for kid in kids:
            kid.FulliTol(files)

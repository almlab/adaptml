# object class for construction of rooted and unrooted binary trees.
#
# lawrence david - ldavid@mit.edu - 2006.

import re
import copy
import random
import math
import sys
import pdb

from numpy import *

import node
import branch

class multitree:
    
    # recursive module for building unrooted trees from newick-notated
    # strings.  
    def __init__(self):
	self.build_queue = []
        self.a_node = None
	self.a_branch = None
	self.root = None
	self.root_branch = None
	self.num_dups = 0
	self.num_losses = 0
	self.num_xfers = 0
	self.event_queue = []           # keep running list of events
	self.node_dict = dict()
	self.mapping_hash = dict()      # store node mappings
        self.loss_dict = {}             # cache res from countlosses
        self.lca_dict = {}
        self.rec_dict = {}
        self.leaf_count = 0             # how many leaves in tree
        self.species_count = {}         # dict w/ counts of species
        self.internal_node_list = []
        self.leaf_node_list = []
        self.branch_list = []
        self.liks = None
        self.migration_matrix = None
        self.mu = None        
        self.states = None
        self.thresh_dict = {}
        self.base_colors = None
        
    # print the tree (intended for rooted trees)
    def __repr__(self):
	if self.root == None:
	    return "not rooted"
	else:
	    newickString = ""
	    return self.root.treePrint(newickString)

    # get a list of species
    def GetSpeciesDict(self):
        species_count = {}

        for i in self.node_dict.values():
            if i.isLeaf():
                if i.species not in species_count:
                    species_count[i.species] = 1
                else:
                    species_count[i.species] += 1
        self.species_count = species_count
        return self.species_count

    # hijack old build, so as to give opportunity to remove bootstraps
    # ... 
    def build(self,newick):

        # remove any bootstraps from newick
        p = re.compile('\)[\d\.]+:')
        newick = p.sub('):',newick)
        self.RecursiveBuild(newick)  

    def RecursiveBuild(self,newick):

	# handle the case when you reach a trifurcating node (so
	# unrooted), by looking for when you've got a stretch of 2
	# colons unseparated by a parenthesis
	colon_match1 = re.findall(':',newick)
	colon_match2 = re.findall(':[^\)]*:[^\)]*:',newick)

	# if you reach a root
	if len(colon_match1) == 1:

	    root_node = self.build_queue[0]
	    
	    for kid_branches in root_node.branch_list:
		root_node.child_branches.append(kid_branches)

            # do rooted stuff
	    self.root = root_node
	    self.root.imposeHierarchy()
            self.labelSubtrees()
            self.Get_Subnodes()
            self.root.Fill_Node_Dict()

        # handle case of a tree with only 2 nodes:
        elif len(colon_match1) == 2:

	    # find an innermost set of parentheses:
            regex = re.search('\([^\(\)]*\)',newick).span()
	    clade = newick[regex[0]+1:regex[1]-1]

	    # split up the clade
	    children = re.split(',',clade)
	    son = children[0]
	    daughter = children[1]

            # build node1
            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)

	    # unite the two nodes
	    center_node = node1.unite(node2)
            
            # create a fake node:
            dummy_node = self.BuildNode("dummy_node:0.1")
            dummy_node.branch_list = []
            self.leaf_count -= 1
            dummy_node.addBranch(0.01)
            dummy_node.myBranch.addNode(center_node)
            self.a_branch = dummy_node.myBranch
            self.a_node = dummy_node
            center_node.UnrootedLeaving()

	# handle the trifurcating node of an unrooted tree
	elif len(colon_match1) == 3 and len(colon_match2) > 0:

            # note that this updated block can handle trees with leaf
            # count >= 3

            # pull out the first node
            regex = re.search('\([^,]*,',newick).span()
            first_leaf = newick[regex[0]+1:regex[1]-1]

            ####newick = re.sub(first_leaf + ',','',newick)
            newick = newick.replace(first_leaf + ',','')

            first_node = self.BuildNode(first_leaf)

            # first_node.branch_list = []
            # don't need to worry about clearing, since we'll be
            # adding branches anyway later down in this block.
            self.build_queue.append(first_node)

	    # join the last 2 nodes:
            regex = re.search('\(.*,[^\)]*\);',newick).span()
            clade = newick[regex[0]+1:regex[1]-2]

            # split up the clade
            children = re.split(',',clade)
            son = children[0]
            daughter = children[1]

            # build two nodes and join
            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)
            center_node = node1.unite(node2)

            # finally, join this new node with the remaining node:
            last_node = self.build_queue[0]
            regex = re.search(':[^,]*,',newick).span()
            last_node.myBranch.addNode(center_node)

            self.a_branch = last_node.myBranch

            # do some unrooted stuff, making sure you don't begin on a
            # leaf. 
            if len(self.a_branch.ends[0].branch_list) == 3:
                self.a_node = self.a_branch.ends[0]
                self.a_branch.ends[0].UnrootedLeaving()
            else:
                self.a_node = self.a_branch.ends[1]
                self.a_branch.ends[1].UnrootedLeaving()

	else:

	    # find an innermost set of parentheses:
	    regex = re.search('\([^\(\)]*\):',newick).span()
	    clade = newick[regex[0]+1:regex[1]-2]

	    # split up the clade
	    children = re.split(',',clade)
	    son = children[0]
	    daughter = children[1]

            # build nodes and unite
            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)
	    center_node = node1.unite(node2)

	    # replace the matched string with the new node:
	    new_newick = newick[0:regex[0]]
	    new_newick += center_node.name
	    new_newick += newick[regex[1]-1:]

	    # add the new node to the queue of extant nodes
	    self.build_queue.append(center_node)

	    # recurse
	    self.RecursiveBuild(new_newick)

    # create a new node:
    def BuildNode(this_tree,node_string):
        splitSon = re.split(':',node_string)

        # see if anyone in the queue has the same name
        node1match = False
        for i in range(len(this_tree.build_queue)):
            if this_tree.build_queue[i].name == splitSon[0]:
                node1 = this_tree.build_queue[i]
                this_tree.build_queue.remove(node1)
                node1match = True
                break
        if node1match is not True:
            node1 = node.node(splitSon[0],this_tree)
            this_tree.leaf_count += 1   # another leaf added to tree 

        node1.addBranch(splitSon[1])

        return node1

    # label the subtrees
    def labelSubtrees(self):

	if self.root == None:
	    print "you're trying to label an unrooted tree"
	else:
	    self.root.subtreeLabel()


    # get all subnodes (useful for finding LCAs)
    def Get_Subnodes(self):
	self.root.Find_Subnodes()


    def rootify(self,center_branch):

	# remember the old branch:
	self.root_branch = center_branch

	# create a new node
	center_name = center_branch.ends[0].name
	center_name += "-" + center_branch.ends[1].name
	center_node = node.node(center_name,self)
	self.root = center_node
	
	# give it children branches
	child1 = branch.branch(center_branch.length/2)
	child1.addNode(center_node)
	child1.addNode(center_branch.ends[0])
	child2 = branch.branch(center_branch.length/2)
	child2.addNode(center_node)
	child2.addNode(center_branch.ends[1])
	center_node.child_branches.append(child1)
	center_node.child_branches.append(child2)

	# erase the original branch from the child nodes branch_list
	for kids in center_branch.ends:
	    kids.branch_list.remove(center_branch)

	# impose a hierarchy from the root
	center_node.imposeHierarchy()
        self.labelSubtrees()
        self.Get_Subnodes()
        self.root.Fill_Node_Dict()

    def MatchNodes(this_tree,that_tree):

        match_dict = {}

        this_node_list = this_tree.node_dict.values()
        that_node_list = that_tree.node_dict.values()

        for this_node in this_node_list:
            for that_node in that_node_list:

                # do these two nodes have the same number of leaves?
                these_leaves = map(lambda i:i.name,this_node.leaf_nodes)
                those_leaves = map(lambda i:i.name,that_node.leaf_nodes)

                if len(these_leaves) != len(those_leaves):
                    continue

                leaf_match = True
                for this_leaf in these_leaves:
                    if this_leaf not in those_leaves:
                        leaf_match = False
                        break

                if leaf_match:
                    match_dict[this_node] = that_node
                    that_node_list.remove(that_node)
                
        return match_dict
        
    def DrawLeaves(this_tree,files):

        cluster_file = files['cluster']

        for leaf in this_tree.leaf_node_list:
            this_str = leaf.PrintiTolLeaf()
            cluster_file.write(this_str + "\n")


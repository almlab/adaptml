import copy;

class branch:

    def __init__(self,length):
        self.ends = []                  # nodes connected to
        self.length = float(length)	# length (can change during rooting)
        self.immutable_length = self.length # don't ever change this
        self.visited = False            # used for traversing the tree

    def __repr__(self):

	if len(self.ends) == 2:
	    print_string = "(" + self.ends[0].name + "," 
	    print_string += self.ends[1].name + "):" + str(self.immutable_length)
	else:
	    print_string = ":" + str(self.immutable_length)
	return print_string
	
    def addNode(self,node):
	self.ends.append(node)
	node.branch_list.append(self)

    # recursion for finding all of the branches in an unrooted tree
    def findBranches(self,all_branches):
	
	all_branches.append(self)
	self.visited = True

	for node in self.ends:
	    for brch in node.branch_list:
		if not brch.visited:
		    all_branches = brch.findBranches(all_branches)

	return all_branches

    # recusion to dump all the branches of an unrooted tree into the
    # provided dictionary
    def FillBranchDict(this_branch,branch_dict):

        for parent_node in this_branch.ends:
            for branch in parent_node.branch_list:
                if branch is not this_branch:
                    for child_node in branch.ends:
                        if child_node is not parent_node:
                            print child_node
                            print child_node.leaves
            

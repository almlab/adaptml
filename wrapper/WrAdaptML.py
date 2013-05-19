#!/usr/bin/python
#
# adaptML wrapper

import os
import pdb
import sys
import shutil
import subprocess as sub

print "\nwelcome to WrAdaptML, the AdaptML wrapper."

if len(sys.argv) <= 1:
    print 'Arguments are:'
    print 'tree'
    print 'init_hab_num'
    print 'outgroup'
    print 'converge_thresh'
    print 'write_dir'
    print 'rateopt'
    print 'collapse_thresh'
    print 'rand'

# read in the variables
inputs = sys.argv
tree_fn = None
rand_bool = False
rand_iter_num = str(100)
init_hab_num = "16"
outgroup = None
write_d = './'
rateopt = 'avg'
mu = 1.00000000001
collapse_thresh = 0.10
converge_thresh = 0.001
color_fn = None

for ind in range(1,len(inputs)):
    arg_parts = inputs[ind].split('=')
    code = arg_parts[0]
    arg = arg_parts[1]
    if code == 'tree':
        tree_fn = arg
    elif code == 'init_hab_num':
        init_hab_num = arg
    elif code == 'outgroup':
        outgroup = arg
    elif code == 'converge_thresh':
        converge_thresh = arg
    elif code == 'write_dir':
        write_d = arg
    elif code == 'rateopt':
        rateopt = arg
    elif code == 'collapse_thresh':
        collapse_thresh = arg
    elif code == 'rand':
        rand_bool = True
        rand_iter_num = arg
        
###############
# run AdaptML #
###############

print "Running AdaptML ..."
wr_path_s = os.path.abspath(sys.argv[0])
home_path = '/'.join(wr_path_s.split('/')[:5])
adaptml_x = home_path + "/habitats/trunk/AdaptML.py"
adaptml_l = []
adaptml_l.append("/usr/local/bin/python2.5")
adaptml_l.append(adaptml_x)
adaptml_l.append("tree=" + tree_fn)
adaptml_l.append("init_hab_num=" + init_hab_num)
adaptml_l.append("outgroup=" + outgroup)
adaptml_l.append("write_dir=" + write_d)
adaptml_l.append("collapse_thresh=" + collapse_thresh)
adaptml_l.append("converge_thresh=" + converge_thresh)
adaptml_l.append("rateopt=" + rateopt)
proc = sub.Popen(adaptml_l,stdout=sub.PIPE)

# wait until finished
stdout,stderr = proc.communicate()

##############
# run RandML #
##############

print "Running RandML ..."
emp_d = write_d + "/emp_trees/"
if os.path.exists(emp_d):
    shutil.rmtree(emp_d)
os.mkdir(emp_d)

randml_x = home_path + "/clusters/getstats/rand_JointML.py"
randml_l = []
randml_l.append("/usr/local/bin/python2.5")
randml_l.append(randml_x)
randml_l.append("tree=" + tree_fn)
randml_l.append("outgroup=" + outgroup)
randml_l.append("write=" + emp_d)
randml_l.append("habitats=" + write_d + "/habitat.matrix")
randml_l.append("mu=" + write_d + "/mu.val")
randml_l.append("iters=" + rand_iter_num)
proc = sub.Popen(randml_l,stdout=sub.PIPE)
# wait until finished
stdout,stderr = proc.communicate()


#!/usr/bin/python
#
# adaptML wrapper

import os
import pdb
import sys
import shutil
import subprocess as sub

if len(sys.argv) != 2:
    print "Please supply a single input file to read AdaptML settings from"

# Read input file
file_location = sys.argv[1]
file_handle = open(file_location, 'r')
run_configs = file_handle.read()
file_handle.close()

# Parse input file into dictionary
file_data = run_configs.split("\n")
run_config = {}
for line in file_data:
    line = line.split('=')
    if len(line) == 1:
        continue
    if len(line) != 2:
        raise ValueError("Cannot parse line "+str(line))
    key = line[0].strip()
    value = line[1].strip()
    run_config[key] = value

# Create variable defaults
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

# Set Variables
for code, arg in run_config.items():
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
base_path = os.path.dirname(os.path.realpath(__file__))+'/../'
adaptml_x = base_path + "/habitats/trunk/AdaptML.py"
adaptml_l = []
adaptml_l.append("python")
adaptml_l.append(adaptml_x)
adaptml_l.append("tree=" + tree_fn)
adaptml_l.append("init_hab_num=" + init_hab_num)
adaptml_l.append("outgroup=" + outgroup)
adaptml_l.append("write_dir=" + write_d)
adaptml_l.append("collapse_thresh=" + collapse_thresh)
adaptml_l.append("converge_thresh=" + converge_thresh)
adaptml_l.append("rateopt=" + rateopt)
proc = sub.Popen(adaptml_l,stdout=sub.PIPE, stderr=sub.PIPE)

# wait until finished
stdout,stderr = proc.communicate()

if stderr != '':
    print stdout
    print 'Errors:'
    print stderr
    sys.exit()

##############
# run RandML #
##############

print "Running RandML ..."
emp_d = write_d + "/emp_trees/"
if os.path.exists(emp_d):
    shutil.rmtree(emp_d)
os.mkdir(emp_d)

randml_x = base_path + "/clusters/getstats/rand_JointML.py"
randml_l = []
randml_l.append("python")
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


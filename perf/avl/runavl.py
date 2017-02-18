# -*- encoding: utf-8 -*-
#
# Nemo equilibrium and performance evaluation
#
# Python script to run avl with some pre-defined
# parameters on several speeds
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Aug. 2016

import os
import subprocess
from jinja2 import Environment, FileSystemLoader

# User defined variables

XFOIL_DIR = "./" # CHANGE TO YOUR XFOIL PATH
CONFIG_NAME = "nemo"
TEMPLATE = "avl_input.j2"

# Auxiliary functions

def populate_vars(n_cases,**kwargs):
    vars = dict()
    vars.update(kwargs)
    vars['runcases'] = list()
    for i in range(n_cases):
        vars['runcases'].append({'id':i+1})
    return vars
    
# Command line interface

if __name__=="__main__":
    # Create environment
    env = Environment(loader=FileSystemLoader('./'),trim_blocks=True)
    template = env.get_template(TEMPLATE)
    # Create input file
    with open("input","w") as fin:
        fin.write(template.render(
            vars=populate_vars(13, iw=2.5)))
    #cat = subprocess.Popen(('cat', 'input'),stdout=subprocess.PIPE)
    #output = subprocess.check_output(('../../../tools/AVL/avl nemo'), stdin=cat.stdout)
    #cat.wait()
    #print cat, output

# Change log

# 2016-08-27 - Gustavo Violato
# ----------------------------
# First release

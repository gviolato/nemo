#!/usr/bin/python
# -*- encoding: utf-8 -*-
#
# Nemo airfoil selection tools
#
# Python script to run xfoil with some pre-defined
# parameters on several airfoils in batch mode
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Nov. 2015

import os
import subprocess
from jinja2 import Environment, FileSystemLoader

# User defined variables

XFOIL_DIR = "../../tools/Xfoil" # CHANGE TO YOUR XFOIL PATH
COORD_DIR = "../../nemo/dbfiles/subsistemas/airfoils/uiuc/coordinates"
WING_DIR = "../../nemo/wings"
AIRFOIL_FILE = "test.txt"
TEMPLATE = "input.j2"

N_CRIT = 7.5
RE_SQRTCL = 356535
N_ITER = 50
CL_MIN = 0.1
CL_MAX = 1.1
D_CL = 0.05

# Auxiliary functions

def populate_vars(foil_name):
    vars = dict()
    vars['airfoil_path'] = COORD_DIR + '/' + foil_name + '.dat'
    vars['ncrit'] = "{:.1f}".format(N_CRIT)
    vars['re_sqrtcl'] = "{:d}".format(RE_SQRTCL)
    vars['niter'] = "{:d}".format(N_ITER)
    vars['cl_min'] = "{:.2f}".format(CL_MIN)
    vars['cl_max'] = "{:.2f}".format(CL_MAX)
    vars['d_cl'] = "{:.2f}".format(D_CL)
    vars['polar_file'] = WING_DIR + "/" + foil_name + '.pol'
    return vars
    
# Command line interface

if __name__=="__main__":
    # Create environment
    env = Environment(loader=FileSystemLoader('./'))
    template = env.get_template(TEMPLATE)
    # For each airfoil, run input on xfoil
    with open(AIRFOIL_FILE,'r') as fid:
        print "bla"
        for line in fid.readlines():
            print line
            with open(XFOIL_DIR+"/input","w") as fin:
                fin.write(template.render(
                    vars=populate_vars(line.strip())))
            cat = subprocess.Popen(('cat', XFOIL_DIR+'/input'),stdout=subprocess.PIPE)
            output = subprocess.check_output((XFOIL_DIR+'/xfoil'), stdin=cat.stdout)
            cat.wait()
            print cat, output

            

# Change log

# 2015-11-16 - Gustavo Violato
# ----------------------------
# First release

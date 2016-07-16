#!/usr/bin/python
# -*- encoding: utf-8 -*-
#
# Nemo parametric hull shapes
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Jun. 2016

# Dependencies

import numpy as np
from jinja2 import Environment, FileSystemLoader
import itertools as it
import subprocess
import os
import time

import hull_geom as geo
import hull_dynamics as dyn

MLT_HOME = '/home/gustavo/prjs/hph/tools/mlt933'
TYPE = 'krecov';

L   = 2.3;         # hull length in m
alt  = 0.06;
Prof = 0.23;       # Draft
Aber = 0.20/2;     # Beam
w = 1.50           # Separation

mass = 87.+14
cgZ = np.array([0.15,0.,0.2])

class BuoyShape():
    def __init__(self,length,depth,width,height,pwr,root):
        self.length = length
        self.keel   = lambda x: geo.normal_rectangle(x,depth/length)
        self.cap    = lambda x: -1*geo.normal_pol4(x,height/length)
        self.waterline   = lambda x: -1*geo.normal_pol4(x,width/length)
        self.bottom   = lambda x: geo.oval(x,pwr,root)
        self.top   = lambda x: -1*geo.normal_pol4(x,1)

shape = BuoyShape(L,Prof,Aber,alt,2,2)

env = Environment(loader=FileSystemLoader('./templates'))
sh_template = env.get_template('autorun.j2')
sh_context = dict()
sh_context['hull_offsets'] = 'nemohull_' + TYPE + '.txt'
mlt_template = env.get_template('nemo_'+TYPE+'.j2')
mlt_context = dict()

(phi, theta, dz) = dyn.equilibrium(mass,cgZ,shape,w,verbose=1)
print "Found equilibirum position. Equilibrium variables are:"
print "Phi={:.1f}".format(phi)+"deg"
print "Theta={:.1f}".format(theta)+"deg"
print "Dz={:.0f}".format(dz*1000)+"mm"

drafts = np.linspace(0,-0.99,20)
thetas = np.linspace(-1.1*theta,1.1*theta,10)

for (d,t) in it.product(drafts,thetas):
    res_fname = 'drag_' + TYPE + '_d{:.0f}t{:.2f}.dat'.format(abs(100*d),t)
    if os.path.exists('./results/'+res_fname):
        print 'Skipping file ' + res_fname
        continue
    print 'Running mlt to generate ' + res_fname
    sh_context['result_file'] = res_fname
    with open(MLT_HOME + '/autorun.sh','w') as fid:
        fid.write(sh_template.render(sh_context))
    mlt_context['trim']='{:.4f}'.format(t)
    mlt_context['sinkage']='{:.4f}'.format(d)
    with open(MLT_HOME + '/in.mlt','w') as fid:
        fid.write(mlt_template.render(mlt_context))
    time.sleep(1)
    subprocess.call(['./autorun.sh','&>/dev/null'],cwd=MLT_HOME)

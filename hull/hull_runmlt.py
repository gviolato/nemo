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
TYPE = 'assim_trapz';

L      = 2.30;       # hull length (LOA) in m
alt    = 0.05;       # Height of top
pontal = 0.18;       # Draft
mboca  = 0.29/2      # Half-beam
w      = 1.35        # Separation

mass = 87.+14
cgZ = np.array([0.074,0.,0.2])

class BuoyShape():
    def __init__(self,length,depth,width,height):
        self.length = length
        self.keel   = lambda x: geo.trapz(x,depth/length,perc=0.92)
        self.waterline   = lambda x: -geo.normal_pol(x,width/length,4)+0.01*np.sin(np.pi*x)
        self.bottom   = lambda x: geo.wigley_cross(x)
        self.cap    = lambda x: -1*geo.normal_pol(x,height/length,8)
        self.top   = lambda x: geo.oval(x,6,2)

shape = BuoyShape(L,pontal,mboca,alt)

mlt_pnts = geo.geom_to_mlt(shape, 0, 71, 61)
np.savetxt(MLT_HOME + '/examples/nemohull_' + TYPE + '.txt',
           -1*mlt_pnts,fmt='%.6f',delimiter=',')

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
thetas = np.arange(-6,6.1,0.5)

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

#!/usr/bin/python
# -*- encoding: utf-8 -*-
#
# Script to prepare hull Buoancy Actions look up table
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Oct. 2016

# Dependencies

import numpy as np

import hull_geom as geo
import hull_dynamics as dyn

PRINT_B = True
PRINT_M = True

GRAV  = 9.80665
RHO_W = 1000
LT2M3 = 1000

# user defined variables
L      = 2.30;       # hull length (LOA) in m
alt    = 0.05;       # Height of top
pontal = 0.18;       # Draft
mboca  = 0.29/2      # Half-beam
w      = 1.35        # Separation

cgZ = np.array([0., 0., 0.2])
DZS = np.arange(0.25,-0.051,-0.01)
THETAS = np.arange(-6,6.1,0.5)

class BuoyShape():
    def __init__(self,length,depth,width,height):
        self.length = length
        self.keel   = lambda x: geo.trapz(x,depth/length,perc=0.92)
        self.waterline   = lambda x: -geo.normal_pol(x,width/length,4)+0.01*np.sin(np.pi*x)
        self.bottom   = lambda x: geo.wigley_cross(x)
        self.cap    = lambda x: -1*geo.normal_pol(x,height/length,8)
        self.top   = lambda x: geo.oval(x,6,2)

def pprint(mat,fmt='{:.3f},'):
    dim = mat.ndim
    retstr = '['
    if dim==1:
        for e in mat:
            retstr += fmt.format(e)
        retstr = retstr[:-1] + ']'
    else:
        for submat in mat:
            retstr += pprint(submat) + ','
    retstr = retstr[:-1] + ']'
    return retstr

if __name__=="__main__":
    shape = BuoyShape(L,pontal,mboca,alt)
    B = np.zeros((len(DZS),len(THETAS)))
    M = np.zeros((len(DZS),len(THETAS)))
    for i,dz in enumerate(DZS):
        for j,theta in enumerate(THETAS):
            la = dyn.long_actions(0., cgZ, np.array([0.,theta,0.]),
                                  np.array([0.,0.,dz]),shape,w)
            B[i,j] = la[0]
            M[i,j] = la[1]

    print pprint(-DZS)
    if PRINT_B:
        print pprint(B.T/(GRAV*RHO_W)*LT2M3)

    if PRINT_M:
        print pprint(M.T/(GRAV*RHO_W)*LT2M3)


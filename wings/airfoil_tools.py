#!/usr/bin/python
# -*- encoding: utf-8 -*-

# Nemo airfoil selection tools
#
# Python script to read, sort, plot and classify airfoils
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Nov. 2015

import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

# User defined variables

AIRFOIL_DB = os.environ['NEMO_ROOT']+"/dbfiles/airfoils"
TCK_MIN = 0.10 # Minimum acceptable thickness
TCK_MAX = 0.15 # Maximum acceptable thickness

# Auxiliary functions

def ReadAirfoil(filepath):
    foil = dict()
    with open(filepath,"r") as fid:
        foil['name'] = fid.readline()
    foil['coords'] = np.loadtxt(filepath,skiprows=1)
    return foil

def SplitUpperLower(airfoil):
    le_idx = np.argmin(airfoil['coords'][:,0])
    if airfoil['coords'][le_idx,1]>0:
        up = np.vstack((airfoil['coords'][0:le_idx+1,:],[0,0]))
        lo = np.vstack(([0,0],airfoil['coords'][le_idx+1:,:]))
    else:
        up = np.vstack((airfoil['coords'][0:le_idx,:],[0,0]))
        lo = np.vstack(([0,0],airfoil['coords'][le_idx:,:]))
    return (up,lo)

def GetFoilThickness(airfoil):
    xc = np.linspace(0,1)
    xc = xc.reshape((len(xc),1))
    (up,lo) = SplitUpperLower(airfoil)
    idx_up = np.argsort(up[:,0])
    up = up[idx_up,:]
    up_i = np.interp(xc,up[:,0],up[:,1])
    lo_i = np.interp(xc,lo[:,0],lo[:,1])
    return np.hstack((xc, up_i-lo_i))

def GetFoilCamberline(airfoil):
    xc = np.linspace(0,1)
    xc = xc.reshape((len(xc),1))
    tck = GetFoilThickness(airfoil)
    (up,lo) = SplitUpperLower(airfoil)
    lo = np.interp(xc,lo[:,0],lo[:,1])
    return np.hstack((xc, lo+tck[:,1:]/2.))
       
def GetFoilProps(airfoil):
    foilprops=dict()
    tck = GetFoilThickness(airfoil)
    camberl = GetFoilCamberline(airfoil)
    foilprops['t/c']=np.max(tck[:,1])
    foilprops['h/c']=np.max(camberl[:,1])
    return foilprops

def PlotFoil(airfoil,plt_chord=True,plt_camber=True):
    plt.plot(airfoil['coords'][:,0],airfoil['coords'][:,1])
    if plt_chord:
        plt.plot([0,1],[0,0],'k')
    if plt_camber:
        camberl = GetFoilCamberline(airfoil)
        plt.plot(camberl[:,0],camberl[:,1],'k--')
    plt.axis('equal')
    plt.axis('off')
    plt.show()
    return 0

# Command line interface

if __name__=="__main__":
    foilfiles = glob(AIRFOIL_DB + "/uiuc/coordinates/*.dat")
    foilfiles.sort()
    for fpath in foilfiles:
        try:
            airfoil = ReadAirfoil(fpath)
            props = GetFoilProps(airfoil)
            if props['t/c']>0.1 and props['t/c']<0.15:
                print os.path.basename(fpath)[:-4]
        except:
            print "Erro! " + fpath


# Change log

# 2015-11-08 - Gustavo Violato
# ----------------------------
# First release

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

AIRFOIL_DB = os.environ['NEMO_ROOT']+"/dbfiles/subsistemas/airfoils"
TCK_MIN = 0.06 # Minimum acceptable thickness
TCK_MAX = 0.08 # Maximum acceptable thickness
CBR_MIN = 0.035 # Minimum acceptable camber
CBR_MAX = 0.045 # Maximum acceptable camber

POLAR_HDR_NLINES = 12

# Auxiliary functions

def WriteAirfoil(filepath,airfoil):
    np.savetxt(filepath, airfoil['coords'], fmt='%10.6f',
               delimiter=' ', header=airfoil['name'].strip(),
               comments='')
    
def ReadAirfoil(filepath):
    foil = dict()
    with open(filepath,"r") as fid:
        foil['name'] = fid.readline().strip()
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
    a = np.linspace(0,np.pi,70)
    xc = (1-np.cos(a))/2
    xc = xc.reshape((len(xc),1))
    (up,lo) = SplitUpperLower(airfoil)
    idx_up = np.argsort(up[:,0])
    up = up[idx_up,:]
    up_i = np.interp(xc,up[:,0],up[:,1])
    lo_i = np.interp(xc,lo[:,0],lo[:,1])
    return np.hstack((xc, up_i-lo_i))

def GetFoilCamberline(airfoil):
    a = np.linspace(0,np.pi,70)
    xc = (1-np.cos(a))/2
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
            if props['t/c']>TCK_MIN and props['t/c']<TCK_MAX:
                if props['h/c']>CBR_MIN and props['h/c']<CBR_MAX:
                    print os.path.basename(fpath)[:-4]
        except:
            print "Error! " + fpath


# Change log

# 2015-11-08 - Gustavo Violato
# ----------------------------
# First release

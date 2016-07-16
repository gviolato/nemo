#!/usr/bin/python
# -*- encoding: utf-8 -*-
#
# Nemo parametric hull shapes
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Mar. 2016

# Dependencies

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

import hull_geom as geo
import hull_dynamics as dyn
import hull_plot as plot

# User define variables

PLOTS = 1
NOT_CONTOURS=1
SAVE  = 0
EQUIL = 1

HULL_NAME = 'nemohull_krecov.txt' # nemohull_wig.txt

L   = 2.3;         # hull length in m
alt  = 0.06;
Prof = 0.24;       # Draft
Aber = 0.20/2;     # Beam
w = 1.50           # Separation
    
V_TO = [0.1,0.2,0.5,0.7,1.0,1.2,1.5,1.7,2.0,2.2,2.5,2.7,3.0]

# Define the total mass and CG position
mass = 87.+14
cgZ = np.array([0.1,0.25,0.2])

### Shape-defining functions


# class BuoyShape():
#     def __init__(self,length,depth,width,height,pwr,root):
#         self.length = length
#         self.keel   = lambda x: geo.wigley_keel(x,depth/length)
#         self.cap    = lambda x: -1*geo.normal_pol4(x,height/length)
#         self.waterline   = lambda x: geo.wigley_wl(x,width/length)
#         self.bottom   = lambda x: geo.wigley_cross(x)
#         self.top   = lambda x: -1*geo.normal_pol4(x,1)

class BuoyShape():
    def __init__(self,length,depth,width,height,pwr,root):
        self.length = length
        self.keel   = lambda x: geo.normal_pol4(x,depth/length)
        self.cap    = lambda x: -1*geo.normal_pol4(x,height/length)
        self.waterline   = lambda x: -1*geo.normal_pol4(x,width/length)
        self.bottom   = lambda x: geo.oval(x,pwr,root)
        self.top   = lambda x: -1*geo.normal_pol4(x,1)
        
if __name__=="__main__":
    # Define the hull shape and dimensions
    shape = BuoyShape(L,Prof,Aber,alt,4,2)
    # Equilibrium calculation
    if EQUIL:
        (phi, theta, dz) = dyn.equilibrium(mass,cgZ,shape,w,verbose=1)
        print "Found equilibirum position. Equilibrium variables are:"
        print "Phi={:.1f}".format(phi)+"deg"
        print "Theta={:.1f}".format(theta)+"deg"
        print "Dz={:.0f}".format(dz*1000)+"mm"
    else:
        for v in V_TO:
            mass = 102-0.5*1000*v**2*0.222*1.001/9.80665
            (_,_,dz) = dyn.equilibrium(mass,cgZ,shape,w)
            print "{:.1f}\t{:.3f}".format(v,-1*dz/Prof)

    if PLOTS:
        # Plot equilibrium scene
        angles = np.array([phi,theta,0])
        position = np.array([0.,0.,dz])
        cgPos_e = cgZ+position
        plot.plot_scene(shape,w,angles,position,cgPos_e)
        # ax = plt.gca()
        #ax.elev = 90
        #ax.azim = 0
        if not NOT_CONTOURS:
            xs = np.linspace(-1,1,51)
            ys = np.linspace(-1,1,51)
            ycs = np.zeros((51,51))
            zcs = np.zeros((51,51))
            for i,x in enumerate(xs):
                cs = geo.cross_section(shape.length/2*(x+1),
                                       shape.length,shape.keel,
                                       shape.waterline,shape.bottom, 26)
                ycs[i,:] = np.append(-1*cs[0,:0:-1],cs[0,:])
                zcs[i,:] = np.append(cs[1,:0:-1],cs[1,:])
                plt.figure()
                wls = plt.contour(xs,ys,zcs, [-dz])
            plt.figure()
            for wlc in wls.collections:
                wl = wlc.get_paths()[0]
                v = wl.vertices
                xnorm = v[:,0]
                ynorm = v[:,1]
                fy_spl = interpolate.RectBivariateSpline(xs,ys,ycs)
                x_wl = shape.length/2*(xnorm+1)
                y_wl = fy_spl.ev(xnorm,ynorm)
                plt.plot(x_wl,y_wl,'b')
            L_s = max(x_wl)-min(x_wl)
            print 'L* =', L_s
            plt.axis('equal')

    if SAVE:
        mlt_pnts = geo.geom_to_mlt(shape, dz, 51, 41)
        np.savetxt(HULL_NAME,-1*mlt_pnts,fmt='%.6f',delimiter=',')
        
    if PLOTS:
        plt.show()

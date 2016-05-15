#!/usr/bin/python
# -*- encoding: utf-8 -*-
#
# Nemo parametric hull shapes
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Mar. 2016

# Dependencies

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

import hull_geom as geo
import hull_dynamics as dyn

def getLimits(left_bot,right_bot,left_top,right_top):
    limits = list()
    for ax in [0,1]:
        limits.append(max([left_bot[ax,:].max(),right_bot[ax,:].max(),
                           left_top[ax,:].max(),right_top[ax,:].max()]))
        limits.append(min([left_bot[ax,:].min(),right_bot[ax,:].min(),
                           left_top[ax,:].min(),right_top[ax,:].min()]))
    return tuple(limits)

def plot_scene(shape,w,angles,position,cgPos_e):
    # Create figure and axes to hold plot
    fig= plt.figure()
    ax = fig.gca(projection='3d')
    # Calculate buoy shape coordinates in its intrinsic coordinate
    # system
    csB_coords_b = geo.getShapePoints(shape,'bottom',25)
    csT_coords_b = geo.getShapePoints(shape,'top',25)

    # Calculate rotated and translated coordinates for left and right
    # buoys
    bl_pos = np.array([0.,
                       -w/2*np.cos(np.radians(angles[0])),
                       -w/2*np.sin(np.radians(angles[0]))])
    csLB_e = dyn.buoy_to_earth(csB_coords_b,angles,position+bl_pos)
    csRB_e = dyn.buoy_to_earth(csB_coords_b,angles,position-bl_pos)
    csLT_e = dyn.buoy_to_earth(csT_coords_b,angles,position+bl_pos)
    csRT_e = dyn.buoy_to_earth(csT_coords_b,angles,position-bl_pos)
    # Plot water surface
    (xM,xm,yM,ym) = getLimits(csLB_e,csRB_e,csLT_e,csRT_e)
    ax.plot_trisurf([xM+1,xM+1,xm-1,xm-1],
                    [yM+1,ym-1,ym-1,yM+1],
                    [0.,0.,0.,0.,0.],color=(0.,0.,1.,0.5))
    # Plot buoys
    ax.plot_trisurf(csLB_e[0,:],csLB_e[1,:],csLB_e[2,:],color='c',linewidth=0.)
    ax.plot_trisurf(csLT_e[0,:],csLT_e[1,:],csLT_e[2,:],color='c',linewidth=0.)
    ax.plot_trisurf(csRB_e[0,:],csRB_e[1,:],csRB_e[2,:],color='c',linewidth=0.)
    ax.plot_trisurf(csRT_e[0,:],csRT_e[1,:],csRT_e[2,:],color='c',linewidth=0.)
    # Plot CG
    ax.scatter(cgPos_e[0],cgPos_e[1],cgPos_e[2],s=25,c='r',marker='o')
    # Scale and present plot
    L = shape.length
    ax.auto_scale_xyz([-L,L],[-L,L],[-L,L])
    plt.axis('off')
    

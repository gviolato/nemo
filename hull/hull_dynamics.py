#!/usr/bin/python
# -*- encoding: utf-8 -*-
#
# Nemo parametric hull shapes
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Mar. 2016

# Dependencies

import hull_geom as geo
import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optim
import Polygon, Polygon.IO

# User define variables

X_AX = 0
Y_AX = 1
Z_AX = 2

D2R = np.pi/180.
RHO_W = 1000      # Water density [kg/m^3]
GRAV = 9.80665    # Gravity [m/s^2]


def to_polygon(array):
    pnts = list()
    for (x,y) in zip(array[0,:].tolist(), array[1,:].tolist()):
        pnts.append([x,y])
    return Polygon.Polygon(pnts)

def to_array(poly):
    return np.array([[p[0],p[1]] for p in poly[0]])

def underwater_part(cs):
    offset = 0.1;
    lowest = min(cs[1,:])
    if lowest < 0:
        p_cs = to_polygon(cs)
        zm = lowest-offset
        yM = max(cs[0,:])+offset
        ym = min(cs[0,:])-offset
        wb_pnts = ((ym,0.),(ym,zm),(yM,zm),(yM,0.))
        p_wb = Polygon.Polygon(wb_pnts)
        p_uw = p_wb & p_cs
        if len(p_uw):
            return to_array(p_uw)
        else:
            return None
    else:
        return None

def RotMatrix(angles):
    # Initialization
    phi   = D2R*angles[X_AX]
    theta = D2R*angles[Y_AX]
    psi   = D2R*angles[Z_AX]
    R = np.zeros((3,3))
    # Preparing terms
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    # Matrix calculation
    R[0,0] = ctheta*cpsi
    R[0,1] = ctheta*spsi
    R[0,2] = -1*stheta
    R[1,0] = sphi*stheta*cpsi-ctheta*spsi
    R[1,1] = sphi*stheta*spsi+cphi*cpsi
    R[1,2] = sphi*ctheta
    R[2,0] = cphi*stheta*cpsi+sphi*spsi
    R[2,1] = cphi*stheta*spsi-sphi*cpsi
    R[2,2] = cphi*ctheta
    return R

def earth_to_buoy(e_coords, angles, translate):
    return np.dot(RotMatrix(angles),e_coords-translate[:,np.newaxis])

def buoy_to_earth(b_coords, angles, translate):        
    return np.dot(RotMatrix(angles).T,b_coords) + translate[:,np.newaxis]
    
def buoyancy_actions(shape,angles,translate,cg_pos):
    """ Returns buoyancy actions: force and moment with respec to cg
    as a function of buoy orientation and vertical displacement """
    # TODO: A diagram explaining reference axis and origin conventions

    l           = shape.length
    waterline_f = shape.waterline
    
    # Auxiliar variables
    buoy_orig = np.array([l/2,0,0])

    # Cross section types
    parts = ['top','bottom']
    # Stations
    stt = np.linspace(0,l,100)

    # Force and Moment Initialization
    F_z = 0;
    M_x = 0;
    M_y = 0;
    for p in parts:
        if p=='top':
            center_f = shape.cap
            cross_f  = shape.top
        elif p=='bottom':
            center_f = shape.keel
            cross_f  = shape.bottom
        areas = list()
        ct_pos = np.empty((3,0))
        for s in stt:
            # Get the cross section
            cs = geo.cross_section(s,l,center_f,waterline_f,cross_f)
            cs = np.hstack((np.vstack((cs[0,:0:-1],cs[1,:0:-1])),
                            np.vstack((-1*cs[0,:],cs[1,:]))))
            csB = np.vstack((s*np.ones((1,cs.shape[1])),cs)) - \
                  buoy_orig[:,np.newaxis]*np.ones((1,cs.shape[1]))
            # Get the coordinates on the earth reference system
            csE = buoy_to_earth(csB,angles,translate)
            # Only get the points with Ze<0 (underwater points)
            csE_cs_uw = underwater_part(csE[1:,:])
            if csE_cs_uw is None:
                areas.append(0)
                ct_pos = np.hstack((ct_pos,np.array([[0.],[0.],[0.]])))
                continue
            x_uw = np.interp(csE_cs_uw[:,0],csE[1,:],csE[0,:])
            csE_uw = np.vstack((x_uw,csE_cs_uw.T))
            # Transform underwater points back to buoy reference system
            # and calculate underwater area and centroid position
            csB_uw = earth_to_buoy(csE_uw,angles,translate)
            (uw_area_i, Ct_B) = geo.centroid_and_area(csB_uw.T)
            Ct_E = buoy_to_earth(Ct_B[:,np.newaxis],angles,translate)
            # Append to the list of areas and centroid positions on the
            # earth reference system
            areas.append(uw_area_i)
            ct_pos = np.hstack((ct_pos,Ct_E))

        areas = np.array(areas)
        ct_pos = ct_pos.T
        # For the force, just integrate the submerged areas over the station
        # positions and multiply by water density and gravity
        F_z += integrate.simps(areas,stt)*RHO_W*GRAV

        # For the moment, split into two components (Mx and My) by integrating
        # the product of forces by their relative distances to the CG
        cgE = buoy_to_earth(cg_pos[:,np.newaxis],angles,translate)
        dY_vec  = ct_pos[:,Y_AX] - cgE[Y_AX]
        dX_vec  = ct_pos[:,X_AX] - cgE[X_AX]
        M_x += integrate.simps(areas*dY_vec.T,stt)*RHO_W*GRAV
        M_y += -1*integrate.simps(areas*dX_vec.T,stt)*RHO_W*GRAV

    return (F_z, M_x, M_y)


def actions(Mass,cgZ,angles,position,shape,w):
    # Initialize return variables
    F_z = -Mass*GRAV
    M_x = 0
    M_y = 0
    # Calculate buoyancy for both buoys
    bl_pos = np.array([0.,
                       -w/2*np.cos(np.radians(angles[X_AX])),
                       -w/2*np.sin(np.radians(angles[X_AX]))])
    F_z_l, M_x_l, M_y_l = buoyancy_actions(shape,angles,
                                           position+bl_pos,
                                           cgZ+np.array([0.,w/2,0.]))
    br_pos = np.array([0.,
                       w/2*np.cos(np.radians(angles[X_AX])),
                       w/2*np.sin(np.radians(angles[X_AX]))])
    F_z_r, M_x_r, M_y_r = buoyancy_actions(shape,angles,
                                           position+br_pos,
                                           cgZ+np.array([0.,-w/2,0.]))
    # Resultant actions
    F_z += F_z_l + F_z_r
    M_x += M_x_l + M_x_r
    M_y += M_y_l + M_y_r
    
    return np.array([F_z, M_x, M_y])

def callbackfun(x,f):
    print 'Solution:', x
    print 'Residual:', f


def long_actions(Mass,cgZ,angles,position,shape,w):
    (F_z, _, M_y) = actions(Mass,cgZ,angles,position,shape,w)
    return np.array([F_z, M_y])
    
def equilibrium(mass,cgZ,shape,w,verbose=None):

    # !! Comment modafoca.
    obj_fun = lambda x: actions(mass,cgZ,np.array([x[0],x[1],0.]),
                                np.array([0.,0.,x[2]]),shape,w)

    if verbose:
        sol = optim.root(obj_fun,np.array([0.,0.,0.02]),method='krylov',
                         callback=callbackfun)
    else:
        sol = optim.root(obj_fun,np.array([0.,0.,0.02]),method='krylov')
        
    return tuple(sol.x)

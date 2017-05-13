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
from scipy import interpolate, integrate
from scipy.optimize import fsolve

# User defined variables

NC  = 30;   # Number of points used to define cross sections

# Functions

def _wigRalpha_(alpha, b, c, r):
    return 2*np.tan(alpha)*np.sqrt(1-r*np.sin(alpha)/b)+c/b

def normal_pol(ksi, frac, deg=4):
    return frac*(np.power(ksi,deg)-1)

def wigley_wl(ksi, frac, f0=1.):
    x = (ksi+1.)/2.
    X = 4.*x*(1.-x)
    return frac*X**f0

def wigley_keel(ksi, frac, f2=0.):
    x = (ksi+1.)/2.
    X = 4.*x*(1.-x)
    return -1*frac*X**f2

def wigley_cross(x, b, c, f1=1.):
    xf = x/b
    return c*np.sqrt(1-xf**(1./f1))

def wigley1_cross_R(x, b, c, rb=0.03):
    r = b*rb
    alpha = fsolve(_wigRalpha_, np.radians(45), args=(b,c,r))
    return np.where(x<=r*np.sin(alpha),
                    c*np.sqrt(1-r*np.sin(alpha)/b) \
                    + r*np.cos(alpha) - np.sqrt(r**2 - x**2),
                    c*np.sqrt(1-x/b))

def oval(x,b,c,pwr,root):
    xf = x/b
    return c*np.power(1-np.power(xf,pwr),1./root)

def normal_rectangle(ksi, frac):
    try:
        return -np.ones((len(ksi),))*frac
    except:
        return -1*frac

def trapz(ksi, frac, perc=0.9):
    ret_val = list()
    if type(ksi)!=type(np.array([])):
        ksi = [ksi]
    for k in ksi:
        if abs(k)>perc:
            ret_val.append(frac*((abs(k)-perc)/(1-perc)-1))
        else:
            ret_val.append(-1*frac)
    if len(ret_val)==1:
        return ret_val[0]
    else:
        return np.array(ret_val)

def line_coord(s, l, nf):
    return l*nf( 2*s/l-1 )

def cross_section(s, l, nf_c, nf_b, cross_fun, npnts=NC+1):
    b  = line_coord(s,l,nf_b)
    c  = line_coord(s,l,nf_c)
    bN = b/2*(1+np.cos(np.linspace(-1,0,npnts)*np.pi))
    y  = cross_fun(bN,b,c)
    return np.vstack((bN,y))

def getShapePoints(shape,part,nstt,npnts=NC+1):
    if part=='top':
        center_f = shape.cap
        cross_f  = shape.top
    elif part=='bottom':
        center_f = shape.keel
        cross_f  = shape.bottom
    sts = np.linspace(0,shape.length,nstt)
    xcs = np.empty((0))
    ycs = np.empty((0))
    zcs = np.empty((0))
    for s in sts:
        cs = cross_section(s,shape.length,center_f,shape.waterline,
                           cross_f, npnts=npnts)
        xcs = np.append(xcs,s*np.ones((2*cs.shape[1]-1,))-shape.length/2)
        ycs = np.append(ycs,np.append(-1*cs[0,:0:-1],cs[0,:]))
        zcs = np.append(zcs,np.append(cs[1,:0:-1],cs[1,:]))
    return np.vstack((xcs, ycs, zcs))

def geom_to_mlt(shape, dz, n_stt, n_wl, l_star=None):
    mlt_pnts = np.empty((n_stt, n_wl))
    if l_star is None:
        sts = np.linspace(0,shape.length,n_stt)
    else:
        l_bow   = (shape.length-l_star)/2.
        l_stern = shape.length-l_bow
        sts = np.linspace(l_bow,l_stern,n_stt)
    depth = -1*shape.keel(0)*shape.length
    zs  = np.linspace(dz,depth,n_wl)
    for i,s in enumerate(sts):
        cs = cross_section(s,shape.length,shape.keel,shape.waterline,
                           shape.bottom)
        z_u, idx_u = np.unique(cs[1,:],return_index=True)
        if len(idx_u)==1:
            mlt_pnts[i,:] = cs[0,0]*np.ones((1,len(zs)))
        else:
            f_wl = interpolate.interp1d(-1*cs[1,:],-1*cs[0,:],
                                        bounds_error=False, fill_value=0.)
            mlt_pnts[i,:] = f_wl(zs)[::-1]
    return mlt_pnts


def area_fun(wl_dist,keel_dist,cross_fun):
    area_vec = np.zeros(1)
    latN = np.linspace(0,1,NC+1)
    for (b,c) in zip(wl_dist[1:-1],keel_dist[1:-1]):
        prof = c*cross_fun(latN)
        lat = b*latN
        area_vec = np.append(area_vec,2*integrate.simps(prof,lat))
    return np.append(area_vec,0)

def volume(l, nf_c, nf_b, nf_cross):
    s = np.linspace(0,l,101)
    b = line_coord(s, l, nf_b)
    c = line_coord(s, l, nf_c)
    a = area_fun(b, c, nf_cross)
    return integrate.simps(a,s)

def centroid_and_area(coords):
    Ct_x = coords[0,0]
    cy = np.append(coords[:,1],coords[0,1])
    cz = np.append(coords[:,2],coords[0,2])
    aux_vec = cy[:-1]*cz[1:]-cy[1:]*cz[:-1]
    area = 0.5*np.sum(aux_vec)
    Ct_y = 1/(6*area)*np.dot(cy[:-1]+cy[1:],aux_vec)
    Ct_z = 1/(6*area)*np.dot(cz[:-1]+cz[1:],aux_vec)
    return (abs(area), np.array([Ct_x,Ct_y,Ct_z]))

# -*- encoding: utf-8 -*-
#
# Nemo basic longitudinal equilibrium/stability analysis
#
# Python script to calculate equilibrium and stability conditions,
# help on canard sizing, facilitate working with AVL and understanding
# pre-takeoff and flying conditions
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Aug. 2016

import numpy as np
from scipy.optimize import fsolve
from jinja2 import Environment, FileSystemLoader
import yaml

# Physical constants

GRAV  = 9.80665 # gravity [m/s^2]
RHO_W = 1020    # (sea) water density [kg/m^3]

CG_HEIGHT = True

def read_dims(dims_file):
    with open(dims_file,'r') as fid:
        dims = yaml.load(fid)
    return dims

def cma_s(surf):
    cr = np.sqrt(surf['S']/surf['AR'])*2/(1+surf['lambda'])
    ctip = cr*surf['lambda']
    return 2./3*(cr+ctip-cr*ctip/(cr+ctip))

def r_cg_s(surf, d_cg):
    return np.sqrt((surf['pos']+d_cg)**2 + surf['h_cg']**2)
    
def equil_coefs(dims, fp):
    wing = dims['wing']
    cnrd = dims['canard']
    keel = dims['keel']
    prop = dims['prop']
    
    q = 0.5*RHO_W*fp['V']**2

    cma  = cma_s(wing)
    r_w  = r_cg_s(wing,fp['delta_cg'])
    delta_w  = np.arctan((wing['pos']+fp['delta_cg'])/wing['h_cg'])
    
    C_factor = cnrd['S']/wing['S']
    c_c  = cma_s(cnrd)
    xc0  = r_cg_s(cnrd,-1*fp['delta_cg'])
    
    CLeq = fp['mass']*GRAV/(q*wing['S'])
    
    F = (C_factor*cnrd['Cla'] + wing['Cla'],
         C_factor*cnrd['Cla'],
         CLeq-C_factor*cnrd['Cl0']-wing['Cl0']-wing['Cla']*fp['iw'])

    M = (r_w/cma*wing['Cla'],
         C_factor*xc0/cma*cnrd['Cla']+
         r_w/cma*(wing['Cl0']+wing['Cla']*(fp['iw']-delta_w)),
         C_factor*xc0/cma*cnrd['Cla'],
         keel['S']/wing['S']*keel['h_cg']/cma*keel['Cd']-
         fp['T']*prop['h_cg']/(q*cma*wing['S'])-
         C_factor*c_c/cma*cnrd['Cm']-wing['Cm']-C_factor*xc0/cma*cnrd['Cl0']+
         r_w/cma*delta_w*(wing['Cl0']+wing['Cla']*fp['iw'])) 
    
    return F, M

def eq(p,F,M):
    alpha, dp = p
    eq1 = F[0]*alpha + F[1]*dp - F[2]
    eq2 = M[0]*alpha**2 + M[1]*alpha + M[2]*dp - M[3]
    return (eq1, eq2)

def solve_equilibrium(Vs, dims, fp, Thrust=None):
    guess = np.pi/180
    wing = dims['wing']
    cnrd = dims['canard']
    delta_w = np.arctan((wing['pos']+fp['delta_cg'])/wing['h_cg'])
    eq_matrix = np.zeros((len(Vs),10))
    fp['T'] = 0
    for i,v in enumerate(Vs):
        fp['V'] = v
        if Thrust is not None:
            fp['T'] = Thrust[i]
        q = 0.5*RHO_W*v**2
        CLeq = fp['mass']*GRAV/(q*wing['S'])
        a, dp = fsolve(eq,(guess,guess),args=equil_coefs(dims,fp))
        Cl_w  = wing['Cl0'] + wing['Cla']*(a+fp['iw'])
        L_w   = q*wing['S']*Cl_w/GRAV
        Cl_c  = cnrd['Cl0'] + cnrd['Cla']*(a+dp)
        L_c   = q*cnrd['S']*Cl_c/GRAV
        Xn    = (wing['pos']+
                 cnrd['pos'])/(1+
                               (wing['Cla']*wing['S']/(cnrd['S']*cnrd['Cla'])))
        Xcg = r_cg_s(wing,fp['delta_cg'])*np.sin(delta_w-a)
        eq_matrix[i,:] = np.array([v,a*180/np.pi,dp*180/np.pi,
                                   CLeq,Cl_w,L_w,Cl_c,L_c,Xn,Xcg])

    return eq_matrix

def write_runcase(runfilename, Vs, dims, fp,
                  template='run_template.j2', to=False,
                  alpha=0., elev=0., Thrust=None):
    # Template loading
    env = Environment(loader=FileSystemLoader('./'))
    template = env.get_template(template)
    # Used constants
    guess = np.pi/180
    # Wing properties
    wing = dims['wing']
    prop = dims['prop']
    cr = np.sqrt(wing['S']/wing['AR'])*2/(1+wing['lambda'])
    cma = cma_s(wing)
    delta_w = np.arctan((wing['pos']+fp['delta_cg'])/wing['h_cg'])
    # Runcases
    runcases = []
    fp['T'] = 0. # Default value for thrust
    for i,v in enumerate(Vs):
        runcases.append(dict())
        runcases[-1]['rho']   = RHO_W
        runcases[-1]['grav']  = GRAV
        runcases[-1]['mass'] = fp['mass']
        runcases[-1]['id'] = i+1
        runcases[-1]['name'] = 'V={:.2f}'.format(v)
        q = 0.5*RHO_W*v**2
        if Thrust is not None:
            fp['T'] = Thrust[i]
        if not to:
            runcases[-1]['alpha'] = 0.0
            CLeq = fp['mass']*GRAV/(q*wing['S'])
            runcases[-1]['cl'] = CLeq
            fp['V'] = v
            a, dp = fsolve(eq,(guess,guess),args=equil_coefs(dims,fp))
        else:
            if isinstance(alpha,np.ndarray):
                runcases[-1]['alpha'] = alpha[i]
            else:
                runcases[-1]['alpha'] = alpha
            runcases[-1]['dp'] = elev
            runcases[-1]['cl'] = 0.0
            a, dp = (alpha,elev)
        runcases[-1]['cm']  =  -1*fp['T']*prop['h_cg']/(q*wing['S']*cma)
        runcases[-1]['vel'] = v
        if not CG_HEIGHT:
            runcases[-1]['x_cg'] = -1*(r_cg_s(wing,fp['delta_cg'])*np.sin(delta_w-a)-cr*0.25)
            runcases[-1]['z_cg'] = 0.0
        else:
            runcases[-1]['x_cg'] = -1*(wing['pos']+fp['delta_cg']-cr*0.25)
            runcases[-1]['z_cg'] = wing['h_cg']
    with open(runfilename,'w') as fid:
        fid.write(template.render(runcases=runcases))
        
if __name__=="__main__":
    fp = {'V':4.,'iw':2.5*np.pi/180,'delta_cg':0.06,'mass':101}
    dims = read_dims('nemo_dims.yml')
    Vs = np.arange(3.3,7.0,0.3)
    Ts = np.array([66.81,62.09,60.08,60.00,61.33,63.79,
                   67.29,71.84,77.63,84.65,95.37,112.26,140.75])
    # print solve_equilibrium(Vs,dims,fp)
    write_runcase('nemo.run',Vs,dims,fp,template='to_template.j2',
                  to=True,elev=9.,alpha=np.arange(-3.5,3.6,0.55))#Thrust=Ts)
    #Vs = np.arange(1.,3.1,0.25)
    #write_runcase('nemo_to.run', Vs, dims, fp,
    #              template='to_template.j2', to=True,
    #              alpha=3.0, elev=8.0)

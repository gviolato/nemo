# -*- encoding: utf-8 -*-
#
# Nemo Hydrofoil Simulator (NHS)
#
# Script to simulate Nemo's Dynamics
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Sep. 2016

# Dependencies
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
#from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import RegularGridInterpolator
import yaml
import sys

# Program variables
TIME_TOL = 1E-10

# User defined variables
SIM_DEFS_FILE = './sim_defs.yml'

GRAV    = 9.80665
RHO_W   = 1020
RHO_AIR = 1.225
LT2M3   = 1000  # Liters to Cubic Meters

INC_HYDRO = True
INC_PROP  = True
INC_HULL  = True

DEBUG_FLAG = True

# Helper Classes
class StateVec():
    def __init__(self, channels, y=None):
        if 'Total' not in channels.keys():
            raise Exception('State channels must have a "Total" key')
        self.channels = channels
        for key,val in channels.iteritems():
            if key!='Total':
                setattr(self, key, None)
            else:
                setattr(self, key, val)
        if y is not None:
            self.populate(y)

    def populate(self,y):
        for k,i in self.channels.iteritems():
            if k!='Total':
                setattr(self,k,y[i])

    def vec(self):
        ret = np.empty((self.Total,1))
        for k,i in self.channels.iteritems():
            if k!='Total':
                ret[i] = getattr(self,k)
        return ret

    def labels_vec(self):
        ret = [None]*self.Total
        for k,i in self.channels.iteritems():
            if k!='Total':
                ret[i] = k
        return ret

# Auxiliar Functions
def _inlimits_(val,low,hi):
    if val>=low and val<=hi:
        return True
    else:
        return False

def _checkVar_(name,sim_defs):
    if name in sim_defs['state']:
        return 0
    elif name in sim_defs['control']:
        return 1
    else:
        raise Exception("Variable name {} is neither a state "
                        "or control variable. Check your defin"
                        "ition file".format(v['xname']))

def _getNVarFunc_(basefun,vartypes,names):
    return lambda x,u: basefun([getattr((x,u)[t],n)
                                for t,n in zip(vartypes,names)])

def _buildCoefDef_(sim_defs,coefgroup,hph='nemo'):
    CoefDef = dict()
    CoefDict = sim_defs[hph][coefgroup]
    for k,v in CoefDict.iteritems():
        if not isinstance(v,dict):
            CoefDef[k] = v
        else:
            n_vars = (len(v.keys())-1)/2
            points = tuple([v['x{}'.format(i+1)] for i in xrange(n_vars)])
            f = RegularGridInterpolator(points,v['vals'],method='linear',
                                        bounds_error=False)
            names = [v['namex{}'.format(i+1)] for i in xrange(n_vars)]
            vartypes = [_checkVar_(n,sim_defs) for n in names]
            CoefDef[k] = _getNVarFunc_(f,vartypes,names)

    return CoefDef

def _loadCoefs_(X,U,names,CoefDef):
    coefs = tuple()
    for n in names:
        try:
            c = CoefDef[n]
        except KeyError:
            raise Exception("Missing coeficient {} from model".format(n))
        if isinstance(c,float):
            coefs += (c,)
        else:
            coefs += (c(X,U)[0],)
    return coefs

# Simulator Functions
def dragBuildUp(Cd_basic,Cl,X,dims):
    """ Adds all sources of drag to the basic Cd """
    # Increasing induced drag from wing about 10% from Hoerner (XI-C, ~11-26)
    # "A single strut in the center increases drag due to lift by some 10%"
    Cdi_inc = 0.1*(Cd_basic - dims['aerocoefs']['cd_prof'])

    # Hydrofoil wave drag, based on Breslin, J. "Hydrofoil Wave Drag Theory"
    # Wing depth is supposed to be around a quarter of span ~1.8/4=0.45m
    F_c    = X.V/np.sqrt(GRAV*dims['cref'])
    Cd_h_w = Cl**2*(0.5/(np.pi*dims['w_AR']) + 0.13/F_c**2)

    # Interference drag, based on Hoerner VIII-4 (8-10)
    # Cdt = Drag/(qt^2) = 17*(t/c)^2-0.05
    ws_t = (dims['w_root_t'] + dims['strut_t'])/2
    rc_t = (dims['rudder_t'] + dims['canard_t'])/2

    Cd_int_ws = dims['ws_drag_factor']* \
                (17*(ws_t/dims['w_root_c'])**2-0.05)*ws_t**2/dims['Sref']
    Cd_int_rc = (17*(rc_t/dims['c_root_c'])**2-0.05)*rc_t**2/dims['Sref']

    # Spray drag, based on Hoerner X-C (10-13)
    # Cdx = Drag/(qx^2) = 0.24*(t/x)^2
    Cd_spray = 0.24*dims['strut_t']**2/dims['Sref']

    # Air resistance (Rider + Hulls), based on Hoerner (3-13) for Rider, and
    # Hoerner XIII-1 (13-1) for Hulls
    Cd_air = (0.54 + 0.1*2*dims['s_hull'])/dims['Sref']*RHO_AIR/RHO_W

    # Build-up
    Cd_total = Cd_basic + Cdi_inc + Cd_h_w + Cd_int_ws + Cd_int_rc + \
               Cd_spray + Cd_air

    # Arbitrary safety factor for other drag sources such as:
    # - Air resistance of frame and interferences
    # - Drag from surface follower
    # - Spray drag from rudder
    # - Propeller/Strut interaction
    SF = dims['drag_SF']

    return SF*Cd_total, (Cdi_inc, Cd_h_w, Cd_int_ws + Cd_int_rc,
                         Cd_spray, Cd_air)

def surface_follower(X):
    """ Returns canard angular deflection based on vehicle height """
    return min(np.pi/180*(35*X.H+21.5), np.pi/180*6)

def controls(t, X, pilot):
    """ Returns input vector for current state of control variables """
    #if t<10:
    #    dp = 3*np.pi/180
    #else:
    #    dp = min(5*np.pi/180,(t-7)*np.pi/180)
    dp  = np.pi/180*3. #surface_follower(X)
    pwr =pilot['pwr'] #+ 30.*np.sin(np.pi*t/10.) # min(1,t/pilot['t_pwr_ramp'])
    return (dp, pwr)

def nemo_hydro_actions(X,U,dims,AeroCoefs):
    """ Calculates hydrodynamic forces and moments """
    q = 0.5*RHO_W*X.V**2
    S = dims['Sref']
    c = dims['cref']
    
    (cl0, cla, cldp, clq) = _loadCoefs_(X,U,['cl0', 'cla', 'cldp', 'clq'],
                                        AeroCoefs)
    Cl_total = cl0 + cla*X.alpha + cldp*U.dp + clq*X.q
    L = q*S*Cl_total

    (cd0, k1, k) = _loadCoefs_(X,U,['cd0', 'k1', 'k'], AeroCoefs)
    Cd_basic = cd0 + k1*X.alpha + k*X.alpha**2
    Cd_total, _ = dragBuildUp(Cd_basic,Cl_total,X,dims)
    D = q*S*Cd_total
    
    (cm0, cma2, cmdp, cmq) = _loadCoefs_(X,U,['cm0', 'cma2', 'cmdp', 'cmq'],
                                         AeroCoefs)
    M = q*S*c*(cm0 + cma2*X.alpha**2 + cmdp*U.dp + cmq*X.q)
    return L, D, M

def nemo_prop_actions(X,U,dims,PropCoefs):
    """ Propeller forces and moments """
    TY = 0
    eff, = _loadCoefs_(X,U,['eff'],PropCoefs)
    T  = min(eff*U.pwr/X.V,100)
    M_T = T*dims['prop_h']
    return TY, T, M_T

def nemo_hull_actions(X,U,dims,HullCoefs):
    """ Forces and moments caused by the hulls """
    (dispV, dispM, dragF) = _loadCoefs_(X,U,['dispV','dispM','dragF'],
                                        HullCoefs)
    q = 0.5*RHO_W*X.V**2
    
    B_H = dispV*RHO_W*GRAV/LT2M3
    D_H = dragF*q/1e3
    M_H = dispM*RHO_W*GRAV/LT2M3 + B_H*dims['hull_cg_xoff']
    return B_H, D_H, M_H

def nemo_dynamics(t, y, sim_defs,
                  AeroCoefGen, PropCoefGen, HullCoefGen):
    """ Nemo rigid body dynamics """
    nemo = sim_defs['nemo'] #Nemo Parameters
    Mass = nemo['Mass']
    
    X = StateVec(sim_defs['state'],y)
    Xdot = StateVec(sim_defs['state'])
    U = StateVec(sim_defs['control'])
    
    U.dp, U.pwr = controls(t, X, sim_defs['pilot'])

    if INC_HYDRO and _inlimits_(X.alpha,-0.0611,0.0253) \
       and _inlimits_(X.V,1.,10.):
        L, D, M = nemo_hydro_actions(X, U, nemo, AeroCoefGen)
    else:
        L, D, M = (0, 0, 0)

    if INC_PROP:
        _, T, M_T = nemo_prop_actions(X, U, nemo, PropCoefGen)
    else:
        T, M_T = (0, 0)
        
    if INC_HULL and _inlimits_(X.H,-0.18,0.05):
        B_H, D_H, M_H = nemo_hull_actions(X, U, nemo, HullCoefGen)
        if X.V < 0.11:
            D_H = 0.
        #D_H = np.sign(X.V)*20.if abs(X.V)>0.05 else 0.
    else:
        B_H, D_H, M_H = (0, 0, 0)

    if DEBUG_FLAG:
        print "Time: ", t
        print L, D, M
        print T, M_T
        print B_H, D_H, M_H
        print X.vec().T
            
    
    gamma = X.theta - X.alpha
    
    Xdot.q = (M + M_T + M_H)/nemo['Iyy']
    Xdot.V = 1/Mass*(-D + T*np.cos(X.alpha) - D_H +
                     (B_H - Mass*GRAV)*np.sin(gamma))
    if X.V < 0.1:
        Xdot.alpha = X.q
    else:
        Xdot.alpha = X.q - 1/(Mass*X.V)*(L + T*np.sin(X.alpha) +
                                         (B_H - Mass*GRAV)*np.cos(gamma))
    Xdot.theta = X.q
    Xdot.X = X.V*np.cos(gamma)
    Xdot.H = -1*X.V*np.sin(gamma)
    
    return Xdot.vec()

# Program Execution
if __name__=="__main__":
    # Read simulation parameters
    with open(SIM_DEFS_FILE,'r') as fid:
        defs = yaml.load(fid.read())
    # Auxiliar variables definitions
    y0 = defs['init_state']
    ip = defs['integration'] # integration parameters
    t0, tf, dt = ip['t_start'], ip['t_end'], ip['t_step']
    # Setting up the integrator
    r = ode(nemo_dynamics).set_integrator('dopri5')
    r.set_initial_value(y0,t0)
    r.set_f_params(defs,
                   _buildCoefDef_(defs,'aerocoefs'),
                   _buildCoefDef_(defs,'propcoefs'),
                   _buildCoefDef_(defs,'hullcoefs'))
    # Setting up the output
    nsteps = int(np.ceil((tf-t0)/dt)+1)
    nchans = len(y0)+1
    yout = np.zeros((nsteps,nchans))
    yout[0,:] = np.append([t0],y0)
    # Integrating loop
    cnt = 1
    while r.successful() and r.t-tf<-TIME_TOL:
        if not cnt%10:
            sys.stdout.write('\rTime {:4.1f}'.format(r.t))
            sys.stdout.flush()
        r.integrate(r.t+dt)
        yout[cnt,:] = np.append([r.t],r.y)
        cnt += 1
    # Plotting results
    sys.stdout.write('\n')
    sys.stdout.flush()
    X = StateVec(defs['state'])
    lbls = X.labels_vec()
    f, axarr = plt.subplots(3,2, sharex='col')
    MultFac = [1.,1.,180/np.pi,180/np.pi,1.,180/np.pi]
    #print 'Max factor:{:.3f}'.format(
    #    abs(min(np.diff(yout[:,5]*np.sin(yout[:,4]))/dt)/GRAV)+1.)
    for i,ax in enumerate(axarr.flat):
        ax.plot(yout[:,0],yout[:,i+1]*MultFac[i])
        ax.set_title(lbls[i])
    plt.show()
    

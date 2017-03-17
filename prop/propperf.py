# -*- encoding: utf-8 -*-
#
# propperf.py: Calculator of propeller performance
#              according to adkin's analysis method
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# v0.1 - Mar. 2017

# Dependencies
import os, sys
import re
import itertools
import yaml
from glob import glob
import argparse

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.environ['NEMO_ROOT'],
                             '../tools/AirfoilPreppy/src'))

import airfoilprep as ap

# Program variables
TOL = 1e-4
ALPHA = 0.5
RPM2RPS = np.pi/30.

RUN_FILE  = 'runs.yml'

# Comand line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runfile",
                    help="YAML file with info on how "
                    "to run calculations", default = RUN_FILE)

# Auxiliary classes

class propfoil(ap.Airfoil):
    def __init__(self, pfs):
        self.Re = []
        polars = []
        for pf in pfs:
            polars.append(self._readpolarfile_(pf))
        ap.Airfoil.__init__(self,polars)
        self.pol_ext = self.extrapolate(1.3).getPolar(self.Re[0])
            

    def _readpolarfile_(self, pf):
        pol = dict()
        with open(pf,'r') as fid:
            lines = fid.readlines()
        pol['Re'] = float(lines[1].split(' ')[0])
        pdata = np.loadtxt(pf,delimiter='\t',skiprows=3)
        pol['alpha'] = pdata[:,0]
        pol['cd'] = pdata[:,1]
        pol['cl'] = pdata[:,2]
        pol['cm'] = pdata[:,3]
        self.Re.append(pol['Re'])
        return ap.Polar(**pol)

    def get_coefs(self,alpha):
        pol_ext = self.pol_ext
        cl = np.interp(alpha,pol_ext.alpha,pol_ext.cl)
        cd = np.interp(alpha,pol_ext.alpha,pol_ext.cd)
        cm = np.interp(alpha,pol_ext.alpha,pol_ext.cm)
        return cl, cd, cm
        

# Auxiliary functions
def get_propinfo(propfile):
    with open(propfile,'r') as fid:
        propinfo = yaml.load(fid.read())

    prop = dict()
    prop.update(propinfo['general'])
    prop.update(propinfo['geometry'])
    for i,foil in enumerate(prop['airfoils']):
        pol_dict = propinfo['polars']
        foil_name = prop['airfoils'][i]
        prop['airfoils'][i] = propfoil([pol_dict[foil_name]])

    return prop, propinfo['general']['name']

def prepare_propdir(propname, erase=False):
    if os.path.isdir(propname):
        if erase:
            for f in glob(os.path.join(propname,'*.out')):
                os.remove(f)
    else:
        os.mkdir(propname)
    
def lookup(filename,xvar):
    matrix = np.loadtxt(filename,delimiter=',')
    return np.interp(xvar,matrix[:,0],matrix[:,1])

def evalformula(formula, op):
    pat = re.compile('.*\$\((.*?)\)')
    m = pat.match(formula)
    g = m.groups()
    pat_sub = re.compile('\$\((.*?)\)')
    formula = re.sub(pat_sub,str(op[g[0]]),formula)
    return eval(formula)

def _total_local_speed_(op,perf_dist):
    local_V = op['V']*(1+perf_dist[:,3])
    local_T = op['RPM']*RPM2RPS*perf_dist[:,0]*(1-perf_dist[:,4])
    return np.sqrt(local_V**2 + local_T**2)

def _local_cy_(cl, cd, phi):
    return cl*np.cos(phi) - cd*np.sin(phi)

def _local_cx_(cl, cd, phi):
    return cl*np.sin(phi) + cd*np.cos(phi)


def propcalcs(xi, sigma, B,  airfoil, beta, tsr, phi):
    alpha = beta - np.degrees(phi)
        
    cl, cd, _ = airfoil.get_coefs(alpha)
    
    Cy = _local_cy_(cl, cd, phi)
    Cx = _local_cx_(cl, cd, phi)
        
    K  = Cy/(4*np.sin(phi)**2)
    Kl = Cx/(4*np.cos(phi)*np.sin(phi))

    f = (B/2)*(1-xi)/np.sin(np.arctan(tsr))
    F = (2/np.pi)*np.arccos(np.exp(-f))
    
    a  = np.min([sigma*K/(F - sigma*K), 0.7])
    al = np.min([sigma*Kl/(F + sigma*Kl), 0.7])

    return cl, cd, a, al

def adkins_iter(xi, sigma, B, airfoil, beta, tsr):
    """ 
    Implements Adkins iteration scheme for propeller performance 
    analysis
    """
    # Initialize
    err = 1.
    phiguess = np.arctan(tsr/xi)
    # Iterate until convergence
    while err>TOL:
        cl, cd, a, al = propcalcs(xi, sigma, B,
                                  airfoil, beta, tsr, phiguess)
        
        phi = np.arctan((1+a)/(1-al)/(tsr*xi))
        err = abs(phi-phiguess)
        phiguess = ALPHA*phi + (1-ALPHA)*phiguess

    # Re-calculate one more time
    cl, cd, a, al = propcalcs(xi, sigma, B, airfoil, beta, tsr, phiguess)
    
    return np.array([cl,cd,a,al, phiguess])

def calc_perf(prop, op):
    # Calculate prop tipspeed ratio
    tsr = op['RPM']*RPM2RPS*prop['D']/2./op['V']
    
    # Apply Adkins iterative procedure to each radius station
    # Perf dist will have 6 columns:
    # radius, cl, cd, a, a', phi (radians)
    perf_dist = np.zeros((len(prop['rs']),6),dtype=np.float)
    for i,r in enumerate(prop['rs']):
        beta = prop['betas'][i]
        airfoil = prop['airfoils'][i]
        c = prop['chords'][i]/1000.
        sigma = prop['B']*c/(2*np.pi*r)
        perf_dist[i,0]  = r
        perf_dist[i,1:] = adkins_iter(2.*r/prop['D'], sigma, prop['B'],
                                      airfoil, beta, tsr)

    # Calculate local performance parameters from
    # converged results
    Cy = _local_cy_(perf_dist[:,1], perf_dist[:,2], perf_dist[:,5])
    Cx = _local_cx_(perf_dist[:,1], perf_dist[:,2], perf_dist[:,5])
    
    W  = _total_local_speed_(op,perf_dist)

    T_dist = 0.5*op['RHO']*W**2*prop['B']*prop['chords']*Cy/1000.
    Q_dist = 0.5*op['RHO']*W**2*prop['B']*prop['chords']*Cx*prop['rs']/1000.

    # Integrate performance coefficients to get thrust
    # and torque
    rs = np.append(prop['rs'],prop['D']/2.)
    T_dist = np.append(T_dist,0.)
    Q_dist = np.append(Q_dist,0.)
    
    T = integrate.simps(T_dist,rs)
    Q = integrate.simps(Q_dist,rs)
    
    # Calculate efficiency
    eff = T*op['V']/(Q*op['RPM']*RPM2RPS)

    # Return performance
    return T, Q, eff, perf_dist


# Program execution
if __name__=="__main__":
    args = parser.parse_args()
    
    with open(args.runfile,'r') as fid:
        runinfo = yaml.load(fid.read())

    PROPSRAN = []
    for runname, op in runinfo.iteritems():
        
        prop, propname = get_propinfo(op['PROPFILE'])
        
        if propname not in PROPSRAN:
            prepare_propdir(propname, erase=True)
            PROPSRAN.append(propname)
        else:
            prepare_propdir(propname)
        
        resfile = os.path.join(propname,'{}.out'.format(runname))
        with open(resfile,'w') as fid:
            fid.write('# Results for prop {}:\n'.format(propname))
            fid.write('# Water density was {:.1f} kg*m^-3\n'.format(op['RHO']))
            fid.write('V,RPM,T,Q,pwr,eff,J,dist_file\n')

        iterargs = tuple()
        argsnames = []
        for name, args in op['_COMBINE'].iteritems():
            argsnames.append(name)
            iterargs += (args,)

        for idx, comb in enumerate(itertools.product(*iterargs)):
            for i,name in enumerate(argsnames):
                op[name] = comb[i]

            if '_OVERRIDE' in op.keys():
                for name, formula in op['_OVERRIDE'].iteritems():
                    op[name] = evalformula(formula, op)
                
            T, Q, eff, perf_dist = calc_perf(prop, op)

            omg = op['RPM']*RPM2RPS
            pwr = Q*omg
            J   = op['V']/(omg*prop['D']/2.)
        
            alpha_dist = prop['betas'] - np.degrees(perf_dist[:,-1])
            perf_dist = np.hstack((perf_dist,alpha_dist[:,np.newaxis]))

            dist_file_name = '{}_dist_{}.out'.format(runname,
                                                     idx)
            dist_file_path = os.path.join(propname,dist_file_name)
            np.savetxt(dist_file_path,perf_dist,delimiter=',',
                       header='r,cl,cd,a,al,phi,alpha')

            with open(resfile,'a') as fid:
                fid.write('{:.3f},{:.1f},{:.3f},{:.3f},'\
                          '{:.3f},{:.5f},{:.3f},{}\n'.format(op['V'],
                                                             op['RPM'],
                                                             T,Q,pwr,
                                                             eff,J,
                                                             dist_file_name))

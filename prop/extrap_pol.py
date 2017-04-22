# -*- encoding: utf-8 -*-
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# v0.1 - Mar. 2017

# Dependencies
import os, sys

import numpy as np
sys.path.append(os.path.join(os.environ['NEMO_ROOT'],
                             '../tools/AirfoilPreppy/src'))

import airfoilprep as ap

# Program variables
TOL = 1e-4
ALPHA = 0.5
RPM2RPS = np.pi/30.

# Auxiliary classes

class propfoil(ap.Airfoil):
    def __init__(self, pfs):
        self.Re = []
        polars = []
        for pf in pfs:
            polars.append(self._readpolarfile_(pf))
        ap.Airfoil.__init__(self,polars)
        self.pol_ext = self.extrapolate(1.8)
            

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
        

# Program execution
if __name__=="__main__":
    RE = 1e5
    foil = propfoil(['hs1712_Re5e4.pol'])
    alpha = np.concatenate((np.arange(-180,-14,5),
                            np.arange(-14,20,1),
                            np.arange(20,181,5)))
    foil2 = foil.pol_ext.interpToCommonAlpha(alpha)
    pol2 = foil2.getPolar(RE)
    with open('hs1712_ext_re5e4.pol','w') as polarfile:
        polarfile.write('# Reynolds Number\n{}\n'.format(RE))
        polarfile.write('# Polar data\n'
                        '# alpha, cl, cd\n')
        
        for a,cl,cd in zip(pol2.alpha,pol2.cl,pol2.cd):
            polarfile.write('{:.6f}\t{:.6f}\t{:.6f}\n'.format(a,cl,cd))


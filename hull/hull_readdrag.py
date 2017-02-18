#!/usr/bin/python
# -*- encoding: utf-8 -*-
#
# Nemo parametric hull shapes
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Jun. 2016

# Dependencies

import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import re
from scipy.interpolate import RegularGridInterpolator

PRINT_LUT = True

TYPE = 'assim_trapz'

DRAFT = 0.19
PAT = re.compile('.*_d(.*)t(.*).dat')

NV = 20
NT = 25
ND = 20

# Helper class to automatically add item to list when
# querying unexistent item's index
class IndexList(list):
    def index(self,val):
        try:
            return super(IndexList,self).index(val)
        except ValueError:
            self.append(val)
            return len(self)-1
        except Exception as e:
            raise e

def read_results(result_fmt,verbose=False):
    # Initializing data
    drag_data = np.empty((NV,NT,ND))
    draft_list = IndexList()
    theta_list = IndexList()
    speed_list = IndexList()
    # Looping through result files
    DRAG_FILES = glob(result_fmt)
    for f in DRAG_FILES:
        if verbose:
            print 'Processing: ',f
        m = PAT.search(f)
        draft_m  = (1-float(m.group(1))/100.)*DRAFT
        theta_deg = float(m.group(2))
        D = np.loadtxt(f,delimiter=',')
        total_drag = np.sum(D[:,1:],1)
        if not len(speed_list):
            speed_list = IndexList(D[:,0])
        drag_data[:,theta_list.index(theta_deg),
                  draft_list.index(draft_m)] = total_drag
    # Re-organizing vectors in order (glob is random)
    ts = list(np.argsort(theta_list))
    theta_list.sort()
    ds = list(np.argsort(draft_list))
    draft_list.sort()
    # Re-organizing final result
    drag_data = drag_data[:,:,ds]
    drag_data = drag_data[:,ts,:]
    return (drag_data, (speed_list, theta_list, draft_list))

def drag2CdS(drag_data, (ss,ts,ds), rho_mlt=1020):
    CdS = np.zeros(drag_data.shape)
    for i,s in enumerate(ss):
        q = 0.5*rho_mlt*s**2
        CdS[i,:,:] = drag_data[i,:,:]*1e6/q
    return CdS

def drag_fun(idxs,drag_data):
    return RegularGridInterpolator(idxs, drag_data)

def pprint(mat,fmt='{:.3f},'):
    dim = mat.ndim
    retstr = '['
    if dim==1:
        for e in mat:
            retstr += fmt.format(e)
        retstr = retstr[:-1] + ']'
    else:
        for submat in mat:
            retstr += pprint(submat) + ','
    retstr = retstr[:-1] + ']'
    return retstr

if __name__=='__main__':
    (drag_data, idxs) = read_results('./results/drag_'+TYPE+'_d*t*.dat',
                                     verbose=False)
    CdS = drag2CdS(drag_data, idxs)
    hull_drag = RegularGridInterpolator(idxs, drag_data)
    if PRINT_LUT:
        idx_names = ['V','Theta','Disp']
        factors = [1.,np.pi/180,1.]
        bias = [0.,0.,-0.190]
        for i,x in enumerate(idxs):
            print idx_names[i], ':'
            print pprint(np.array(x)*factors[i]+bias[i])
        print 'Values:'
        print pprint(CdS)


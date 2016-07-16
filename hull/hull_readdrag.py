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

DRAFT = 0.24
PAT = re.compile('.*_d(.*)t(.*).dat')

NV = 25
NT = 10
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

def drag_fun(idxs,drag_data):
    return RegularGridInterpolator(idxs, drag_data)

if __name__=='__main__':
    (drag_data, idxs) = read_results('./results/drag_krecov_d*t*.dat',
                                     verbose=True)
    hull_drag = RegularGridInterpolator(idxs, drag_data)

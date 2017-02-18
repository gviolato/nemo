# -*- encoding: utf-8 -*-
#
# Nemo equilibrium and performance evaluation
#
# Python script to read avl ouput and build 
# equilibrium parameter matrix
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Aug. 2016

# Dependencies
import re
import numpy as np
import os.path as pth
from matplotlib.pyplot import subplots, show
from scipy.optimize import curve_fit

# User defined variables

# Search Patterns
A     = re.compile(r'\s+Alpha.*?=\s*([\d\.\-\+E]+)')
DP    = re.compile(r'\s+elevator.*?=\s*([\d\.\-\+E]+)')
CL_EQ = re.compile(r'\s+CLtot.*?=\s*([\d\.\-\+E]+)')
CD    = re.compile(r'\s+CDtot.*?=\s*([\d\.\-\+E]+)')
CM    = re.compile(r'.*Cmtot.*?=\s*([\d\.\-\+E]+)')
XREF  = re.compile(r'\s+Xref.*?=\s*([\d\.\-\+E]+)')
CD_IND = re.compile(r'.*CDind.*?=\s*([\d\.\-\+E]+)')

VARS_FT = {'alpha':A,
           'dp':DP,
           'cl_eq':CL_EQ,
           'cd':CD,
           'cm':CM,
           'cdind': CD_IND,
           'xref':XREF}

FT_VAR_NAMES = ['alpha', 'dp', 'cl_eq', 'cd', 'cm', 'cdind', 'xref']

XN    = re.compile(r'\s*Neutral point  Xnp.*?=\s*([\d\.\-\+E]+)')
CLA   = re.compile(r'.*?CLa\s*?=\s*([\d\.\-\+E]+)')
CMA   = re.compile(r'.*?Cma\s*?=\s*([\d\.\-\+E]+)')
CLQ   = re.compile(r'.*?CLq\s*?=\s*([\d\.\-\+E]+)')
CMQ   = re.compile(r'.*?Cmq\s*?=\s*([\d\.\-\+E]+)')
CLDP  = re.compile(r'.*?CLd1\s*?=\s*([\d\.\-\+E]+)')
CMDP  = re.compile(r'.*?Cmd1\s*?=\s*([\d\.\-\+E]+)')

VARS_ST = {'cla': CLA,
           'cma':CMA,
           'clq':CLQ,
           'cmq':CMQ,
           'cldp':CLDP,
           'cmdp':CMDP}

ST_VAR_NAMES = ['cla','cma','clq','cmq','cldp','cmdp']

def lin(x,b,a):
    return b + a*x

def parab(x,c,b,a):
    return c + b*x + a*x**2

def jumplines(fid,n):
    for _ in xrange(n):
        fid.readline()

def read_FN(fn_file):
    with open(fn_file,'r') as fid:
        jumplines(fid,18)
        line = fid.readline()
        clw = float(line.strip().split()[3])
        jumplines(fid,2)
        line = fid.readline()
        clc = float(line.strip().split()[3])
    return (clw, clc)

def read_patterns(out_file, var_dict, var_names):
    res = dict()
    for varn in var_dict.keys():
       res[varn] = 0 
    with open(out_file,'r') as fid:
        for line in fid:
            for varn, pat in var_dict.iteritems():
                m = pat.match(line) 
                if m is not None:
                    res[varn] = float(m.group(1))

    ret = tuple()
    for v in var_names:
        ret += (res[v],)
                        
    return ret

def read_results(res_dir,VS):
    res_fmt= pth.join(res_dir,'case_{}_{}.out')
    eq_matrix = np.zeros((len(VS),10))
    st_matrix = np.zeros((len(VS),6))
    for i,v in enumerate(VS):
        (a, dp, cl_eq, cd,
         cm, cdind, xref) = read_patterns(res_fmt.format(i+1,'FT'),
                                          VARS_FT,
                                          FT_VAR_NAMES)
        st = read_patterns(res_fmt.format(i+1,'ST'),
                           VARS_ST,
                           ST_VAR_NAMES)
        (clw, clc) = read_FN(res_fmt.format(i+1,'FN'))
        eq_matrix[i,:] = np.array([v,a,dp,cl_eq,clw,clc,cd,cm,cdind,xref])
        st_matrix[i,:] = np.array(list(st))
    return eq_matrix, st_matrix

if __name__=='__main__':
    VS = np.arange(3.3,6.1,0.3)
    dps = [0,3,6,9]
    fmts = ['o-','r+-','gx-','ks-']
    eqd = dict()
    std = dict()
    cdfit = np.zeros((3,4))
    for i,dp in enumerate(dps):
        k = 'dp{}'.format(dp)
        eq, st = read_results('./results_'+k,VS)
        eqd[k] = eq
        std[k] = st
        popt, pcov = curve_fit(parab,eq[:,1]*np.pi/180,eq[:,6])
        cdfit[:,i] = popt.T
    alpha = eqd['dp6'][:,1]*np.pi/180
    cma = std['dp6'][:,1]
    popt, pcov = curve_fit(lin,alpha,cma)
    print 'Cma2 : {}'.format(popt[1]/2.)
    f, axarr = subplots(3,2, sharex='col')
    for i,ax in enumerate(axarr.flat):
        for j,dp in enumerate(dps):
            k = 'dp{}'.format(dp)
            ax.plot(eqd[k][:,1]*np.pi/180,std[k][:,i],fmts[j])
            if i==1 and j==1:
                a_s = np.arange(-4,3,0.1)
                ax.plot(a_s*np.pi/180,popt[0]+popt[1]*a_s*np.pi/180,'--')
        ax.set_title(ST_VAR_NAMES[i])
    c_tot_names = ['CL','CD','CM']
    c_tot_idxs  = [3,6,7]
    f2, axarr2 = subplots(3,1, sharex='col')
    for i,ax in enumerate(axarr2.flat):
        for j,dp in enumerate(dps):
            k = 'dp{}'.format(dp)
            ax.plot(eqd[k][:,1]*np.pi/180,eqd[k][:,c_tot_idxs[i]],fmts[j])
        ax.set_title(c_tot_names[i])
    f3, axarr3 = subplots(3,1, sharex='col')
    for i,ax in enumerate(axarr3.flat):
        ax.plot(dps,cdfit[i,:])
    show()
    
# Change log

# 2016-08-28 - Gustavo Violato
# ----------------------------
# First release

#!/usr/bin/python
# -*- encoding: utf-8 -*-
#
# Preliminary propeller calculations
#
#
# Diego Montero, Fernando Valentini, Gustavo Violato
# First Release: Apr. 2016

# External Modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# User defined vars
PWR_INPUT = 300*0.95 # Available power after drive train [W]
RHO_H20   = 999      # Water density [kg/m^3]
SPEED     = 5.       # Cruise Speed [m/s]
DES_EFF   = 0.99     # Desired efficiency

def D_from_Eff(eff):
    return np.sqrt(2*PWR_INPUT/(np.pi*RHO_H20*
                                np.multiply(1-eff,np.power(SPEED/eff,3))))

def Eff_Error(eff,spd,D):
    return np.power(2*PWR_INPUT/(np.pi*RHO_H20*D**2*(1-eff)),1./3.)*eff-spd

def Eff_from_Speed(spd,D):
    effs = np.empty((0,),np.float)
    for s in spd:
        effs = np.append(effs,fsolve(Eff_Error,0.7,(s,D),factor=0.1))
    return effs

if __name__=="__main__":    
    effs = np.linspace(0.8, 0.999, 50)
    Ds = D_from_Eff(effs)
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.plot(Ds,effs)
    ax1.plot(D_from_Eff(DES_EFF),DES_EFF,'ro')
    print "Diameter for {:.1f}%Eff@5m/s is {:.3f}m".format(100*DES_EFF,D_from_Eff(DES_EFF))
    ax2 = plt.subplot(212)
    spd_range = np.linspace(3.,8.,20)
    effs_spd = Eff_from_Speed(spd_range,D_from_Eff(DES_EFF))
    ax2.plot(spd_range,effs_spd)
    plt.show()

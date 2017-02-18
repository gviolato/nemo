# coding: utf-8

import numpy as np

x = (np.cos(np.linspace(0,np.pi))+1.)/2.
y = np.sqrt(0.5**2 - (x-0.5)**2)

top = np.hstack((x[:,None],y[:,None]))
bottom = np.hstack((top[::-1,0:1],-1*top[::-1,1:]))

circ = np.vstack((top,bottom[1:-1,:]))
np.savetxt('circ.dat',circ,delimiter=',',fmt='%.5f')

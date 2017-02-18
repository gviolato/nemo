# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

FILES = ['drag_K4Cov.dat','drag_wig.dat','drag_krecov.dat','drag_a_rec.dat','drag_a_trap1.dat','drag_a_trap2.dat']

results = dict()

f = plt.figure()

for fname in FILES:
    data = np.loadtxt(fname,delimiter=',')
    dragT = data[:,1] + data[:,2]
    pwr = np.multiply(dragT*1000,data[:,0])
    results[fname] = np.hstack((data[:,0][:,None],dragT[:,None]))
    plt.subplot(211)
    plt.plot(data[:,0],dragT)
    plt.legend(FILES,loc='best')
    plt.subplot(212)
    plt.plot(data[:,0],pwr)

plt.show()

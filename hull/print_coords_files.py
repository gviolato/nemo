# coding: utf-8

import numpy as np

import hull_geom as geo

L   = 2.3;         # hull length in m
alt  = 0.06;
Prof = 0.24;       # Draft
Aber = 0.20/2;     # Beam

class BuoyShape():
    def __init__(self,length,depth,width,height,pwr,root):
        self.length = length
        self.keel   = lambda x: geo.normal_pol4(x,depth/length)
        self.cap    = lambda x: -1*geo.normal_pol4(x,height/length)
        self.waterline   = lambda x: -1*geo.normal_pol4(x,width/length)
        self.bottom   = lambda x: geo.oval(x,pwr,root)
        self.top   = lambda x: -1*geo.normal_pol4(x,1)

if __name__=="__main__":
    shape = BuoyShape(L,Prof,Aber,alt,4,2)
    xs = np.linspace(-1,1,21)
    STTS = list()
    for i,x in enumerate(xs):
        offset = 0
        if i==0:
            offset = 0.001
        elif i==20:
            offset = -0.001
            station = shape.length/2*(x+1) + offset
        STTS.append(station)
        cs_b = geo.cross_section(station,
                                 shape.length,shape.keel,
                                 shape.waterline,shape.bottom, 50)
        cs_t = geo.cross_section(station,
                                 shape.length,shape.cap,
                                 shape.waterline,shape.top, 50)
        cs = np.hstack((cs_b,cs_t[:,-2::-1]))*1000
        np.savetxt('./results/cs_{:02d}.dat'.format(i),cs.T,
                   delimiter=',',fmt='%.5f',header='HULL SECTION')

    np.savetxt('./results/offsets.dat',np.array(STTS)*1000,fmt='%.3f')

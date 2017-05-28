# coding: utf-8

import numpy as np

import hull_geom as geo

NSTTS = 101
TOL = 1E-3

L      = 2.30;       # hull length (LOA) in m
alt    = 0.005;       # Height of top
pontal = 0.23;       # Draft
mboca  = 0.29/2      # Half-beam

class BuoyShape():
    def __init__(self,length,depth,width,height):
        self.length = length
        self.keel   = lambda x: geo.trapz(x,depth/length,perc=0.96)
        self.waterline   = lambda x: -geo.normal_pol(x,width/length,4)+0.01*np.sin(np.pi*x)
        self.bottom   = lambda x,b,c: geo.wigley1_cross_R(x,b,c, rb=0.04)
        self.cap   = lambda x: -geo.trapz(x,height/length,perc=0.995)
        self.top   = lambda x,b,c: -geo.trapz(x,height,perc=0.995)

if __name__=="__main__":
    shape = BuoyShape(L,pontal,mboca,alt)
    xs = np.cos(np.linspace(-1,0,NSTTS)*np.pi)
    STTS = list()
    for i,x in enumerate(xs):
        offset = 0
        if i==0:
            offset = 0.00025
        elif i==NSTTS-1:
            offset = -0.00025
        if abs(abs(x)-0.96)<TOL:
            x = np.sign(x)*0.96
        station = shape.length/2*(x+1) + offset
        STTS.append(station)
        cs_b = geo.cross_section(station,
                                 shape.length,shape.keel,
                                 shape.waterline,shape.bottom, 70)
        #cs_t = geo.cross_section(station,
        #                         shape.length,shape.cap,
        #                         shape.waterline,shape.top, 30)
        #cs = np.hstack((cs_b,cs_t[:,-2::-1]))*1000
        cs = cs_b*1000
        cs = np.vstack((cs,station*1000*np.ones((1,cs.shape[1]))))
        cs = cs[:,::-1]
        cs_mirror = np.vstack((-1*cs[0,-2::-1],cs[1,-2::-1],cs[2,-2::-1]))
        cs = np.hstack((cs,cs_mirror))
        np.savetxt('./cross_sections/cs_{:02d}.sldcrv'.format(i),cs.T,
                   delimiter=' ',fmt='%.5f')#,header='HULL SECTION {}'.format(i))

    np.savetxt('./cross_sections/offsets.dat',np.array(STTS)*1000,fmt='%.3f')

## Simples calculo de viga "I"para resistir ao momento

import numpy as np
from scipy.optimize import fsolve

## Dados

# Tensao de escoamento do aluminio 6061-0
# (o que o 6061-T6 vira nos pontos de solda)
SIGMA_Y_AL = 55E6   # Pa

# Carga de reacao equivalente em um dos botes
LOAD = 100*9.81     # N

# Braco da viga (metade da separacao dos botes)
D_ARM = 0.65        # m

# Relacao de largura e altura em diversas
# barras chatas de aluminio

AC_DATA_MM = [[12.7,3.17],[15.87,3.17],[25.4,3.17],[25.4,4.76],[31.75,3.17]] #mm


## Funcoes auxiliares

def p_equation(p, E):
    return p**2 - p/E + 1./12.

if __name__=="__main__":
    # Razao escoamento/momento
    YL = SIGMA_Y_AL/(LOAD*D_ARM)
    
    for i,ac in enumerate(AC_DATA_MM):
        A = ac[0]/1000.
        C = ac[1]/1000.
        E = 4*A**2*C*YL
        p = fsolve(p_equation, 4., args=E)
        print "Option {}:".format(i)
        print "A={:.1f}mm\nC={:.1f}mm\nH={:.1f}mm\n\n".format(A*1000,
                                                              C*1000,
                                                              p[0]*A*1000)


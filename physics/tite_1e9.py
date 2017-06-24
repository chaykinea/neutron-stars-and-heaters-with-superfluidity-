__author__ = 'maryhallow'

import numpy as np
from data import loaddata


def coef_init(rho, comp):

    global coeff, f_1, f_2, f_3, f_4, f_5

    if comp==1:

        f_1 = lambda Y: coeff[0] * np.power(Y, -coeff[1]) * (coeff[2] * np.power(Y, 2) + coeff[3] * np.power(Y, 4) - 1)
        f_2 = lambda Y: coeff[8] * (Y ** (coeff[10] - (coeff[9] * np.log10(Y) * np.log10(Y))))
        f_3 = lambda Y: coeff[11] * np.sqrt(1. / (Y * Y + coeff[12] * coeff[12])) * (1. - coeff[13] * Y * Y)
        f_4 = lambda Y: coeff[4] * np.power(Y, coeff[5]) * (1 + coeff[6] * np.power(Y, 2) - coeff[7] * np.power(Y, 4))
        f_5 = lambda Y: -0.4

        if rho==8:
            coeff = np.array([0.2420, 0.4844, 38.35, 0.8680, 5.184, 1.651, -0.04390, 0.001929, 3.462e4, 2.728, 4.120, 2.161, 2.065, 0.008442])
        elif rho==9:
            coeff = np.array([0.1929, 0.4239, 48.72, 1.423, 5.218, 1.652, 0.001037, 0.004236, 3.605e4, 2.119, 4.014, 1.943, 1.788, 0.01758])
        elif rho==10:
            coeff = np.array([0.1686, 0.3967, 55.94, 1.992, 5.208, 1.651, 0.03235, 0.005417, 3.652e4, 1.691, 3.930, 2.021, 1.848, 0.02567])
            
    elif comp==2:

        f_1 = lambda Y: coeff[0] * np.power(Y, coeff[1]*np.log10(Y) + coeff[2])
        f_2 = lambda Y: coeff[6] * np.power(Y, coeff[7]*np.log10(Y)*np.log10(Y) + coeff[8])
        f_3 = lambda Y: coeff[9] * np.sqrt(Y / (Y * Y + coeff[10] * coeff[10]))
        f_4 = lambda Y: coeff[3] * np.power(Y, coeff[4]*np.log10(Y) + coeff[5])
        f_5 = lambda Y: -0.2

        if rho == 8:
            coeff = np.array([5.161, 0.03319, 1.654, 3.614, 0.02933,  1.652, 1.061e5, 1.646,  3.707, 4.011, 1.153])
        elif rho == 9:
            coeff = np.array([5.296, 0.07402, 1.691, 3.774, 0.08210, 1.712, 1.057e5, 1.915, 3.679, 3.878, 1.110])
        elif rho == 10:
            coeff = np.array([5.386, 0.1027, 1.719, 3.872, 0.1344, 1.759, 1.056e5, 1.881, 3.680, 3.857, 1.102])


def T_b(Y,rho_star):

    temp = 1.e7 * f_4(Y) + 1.e7 * (f_1(Y) - f_4(Y)) * ((1 + (rho_star/f_2(Y)) ** f_3(Y)) ** f_5(Y))
    return temp

def tite_2(rho,g_s,rho_st, comp):

    tite = open('data/tite2.dat', 'w')

    g_s0 = 2.4271e14
    g_s *= 1e14

    coef_init(rho, comp)

    T_loc = np.logspace(5.2,6.71,50)
    if comp==1:
        p = [coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5],coeff[6],coeff[7],coeff[8],coeff[9],coeff[10],coeff[11],coeff[12],coeff[13]]
    else:
        p = [coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5],coeff[6],coeff[7],coeff[8],coeff[9],coeff[10]]

    tite.write(str(p) + '\n')
    tite.write('log(rho_b) = %2.1f g/cm3, g_s = %1.4e cm/sec2\n' % (rho,g_s))
    tite.write('------------------------\n')

    for i in range(0,len(T_loc)):

            Y = T_loc[i] * 1e-6 * np.power( g_s0/g_s , 1./4 )
            temp_T_b = T_b(Y,rho_st)
            tite.write('%6.3f %7.5f \n' % ((np.log10(temp_T_b)),np.log10(T_loc[i])))

    tite.close()

def init(rho_bound,rho_star, comp=1):
    tite_2(rho_bound, loaddata.g_surface(), rho_star, comp)




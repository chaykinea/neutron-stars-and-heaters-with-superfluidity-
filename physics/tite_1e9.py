__author__ = 'maryhallow'

import numpy
from scipy.interpolate import *
from data import loaddata
import matplotlib.pylab as plt

def coef_init():

    global p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10, p_11, p_12, p_13, p_14

    log_density = numpy.array([8.0,9.0,10.0])

    p_rho_8 = numpy.array([0.2420, 0.4844, 38.35, 0.8680, 5.184, 1.651,
                   -0.04390, 0.001929, 3.462e4, 2.728, 4.120, 2.161, 2.065, 0.008442])

    p_rho_9 = numpy.array([0.1929, 0.4239, 48.72, 1.423, 5.218, 1.652,
                   0.001037, 0.004236, 3.605e4, 2.119, 4.014, 1.943, 1.788, 0.01758])

    p_rho_10 = numpy.array([0.1686, 0.3967, 55.94, 1.992, 5.208, 1.651,
                    0.03235, 0.005417, 3.652e4, 1.691, 3.930, 2.021, 1.848, 0.02567])

    p_1 = interpolate.interp1d(log_density, numpy.array([p_rho_8[0],p_rho_9[0],p_rho_10[0]]), kind='linear')
    p_2 = interpolate.interp1d(log_density, numpy.array([p_rho_8[1],p_rho_9[1],p_rho_10[1]]), kind='linear')
    p_3 = interpolate.interp1d(log_density, numpy.array([p_rho_8[2],p_rho_9[2],p_rho_10[2]]), kind='linear')
    p_4 = interpolate.interp1d(log_density, numpy.array([p_rho_8[3],p_rho_9[3],p_rho_10[3]]), kind='linear')
    p_5 = interpolate.interp1d(log_density, numpy.array([p_rho_8[4],p_rho_9[4],p_rho_10[4]]), kind='linear')
    p_6 = interpolate.interp1d(log_density, numpy.array([p_rho_8[5],p_rho_9[5],p_rho_10[5]]), kind='linear')
    p_7 = interpolate.interp1d(log_density, numpy.array([p_rho_8[6],p_rho_9[6],p_rho_10[6]]), kind='linear')
    p_8 = interpolate.interp1d(log_density, numpy.array([p_rho_8[7],p_rho_9[7],p_rho_10[7]]), kind='linear')
    p_9 = interpolate.interp1d(log_density, numpy.array([p_rho_8[8],p_rho_9[8],p_rho_10[8]]), kind='linear')
    p_10 = interpolate.interp1d(log_density, numpy.array([p_rho_8[9],p_rho_9[9],p_rho_10[9]]), kind='linear')
    p_11 = interpolate.interp1d(log_density, numpy.array([p_rho_8[10],p_rho_9[10],p_rho_10[10]]), kind='linear')
    p_12 = interpolate.interp1d(log_density, numpy.array([p_rho_8[11],p_rho_9[11],p_rho_10[11]]), kind='linear')
    p_13 = interpolate.interp1d(log_density, numpy.array([p_rho_8[12],p_rho_9[12],p_rho_10[12]]), kind='linear')
    p_14 = interpolate.interp1d(log_density, numpy.array([p_rho_8[13],p_rho_9[13],p_rho_10[13]]), kind='linear')

def f_1(rho,Y):
    return p_1(rho)*numpy.power(Y,-p_2(rho)) * (p_3(rho)*numpy.power(Y,2) + p_4(rho)*numpy.power(Y,4) - 1)

def f_2(rho,Y):
    return p_9(rho) * (Y ** (p_11(rho) - (p_10(rho) * numpy.log10(Y) * numpy.log10(Y))))

def f_3(rho,Y):
    return p_12(rho) * numpy.sqrt(1./(Y*Y + p_13(rho)*p_13(rho))) * (1. - p_14(rho)*Y*Y)

def f_4(rho,Y):
    return p_5(rho)*numpy.power(Y,p_6(rho)) * (1 + p_7(rho)*numpy.power(Y,2) - p_8(rho)*numpy.power(Y,4))

def f_5(rho,Y):
    return -0.4

def T_b(rho,Y,rho_star):

    temp = 1.e7 * f_4(rho,Y) + 1.e7 * (f_1(rho,Y) - f_4(rho,Y)) * ((1 + (rho_star/f_2(rho,Y)) ** f_3(rho,Y)) ** f_5(rho,Y))
    return temp

def tite_2(rho,g_s,rho_st):

    tite = open('data/tite2.dat', 'w')

    g_s0 = 2.4271e14
    g_s *= 1e14

    coef_init()

    T_loc = numpy.logspace(5.2,6.71,50)

    p = [p_1(rho),p_2(rho),p_3(rho),p_4(rho),p_5(rho),p_6(rho),p_7(rho),p_8(rho),p_9(rho),p_10(rho),p_11(rho),p_12(rho),p_13(rho),p_14(rho)]
    tite.write(str(p) + '\n')
    tite.write('log(rho_b) = %2.1f g/cm3, g_s = %1.4e cm/sec2\n' % (rho,g_s))
    tite.write('------------------------\n')

    for i in range(0,len(T_loc)):

            Y = T_loc[i] * 1e-6 * numpy.power( g_s0/g_s , 1./4 )
            temp_T_b = T_b(rho,Y,rho_st)
            tite.write('%6.3f %7.5f \n' % ((numpy.log10(temp_T_b)),numpy.log10(T_loc[i])))

    tite.close()

def init(rho_bound,rho_star):
    tite_2(rho_bound,loaddata.g_surface(),rho_star)

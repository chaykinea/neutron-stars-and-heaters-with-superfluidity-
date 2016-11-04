#     Neutron star thermal evolution
# --------------------------------------
#            routines.py
# --------------------------------------
# This module contains different
# routines that is better to put
# alone for code consistency

from control.constants import *
from data import loaddata
import numpy


def pF(nd):                                                      # Fermi momentum for degenerate Fermi-gas
                                                                 # nd - number density of the gas in fm

    if nd > 0.:
        p=h*numpy.power( 3.*pi*pi*nd, 1./3. )/fm
    else:
        p=0.
    return p

def sqr(x):

    return x*x

def _Meff(n):                                                    # n in fm^-3

    kf = numpy.power(3.*pi*pi*n, 1./3.)                          # Fermi momentum in fm^-1 */
    f = loaddata.kf(kf)

    if(f<0.7):
        f=0.7

    if(f>1.):
        f=1

    return f


def Tri_diag_matrix_solver_init(N):

    global gamma, beta

    gamma =  numpy.zeros(N+1)
    beta  = numpy.zeros(N+1)


def Tri_diag_matrix_solver(a, b, c, d, N):

    gamma[0] = beta[0] = 0
    solution = numpy.zeros(N+1)

    for i in range (0,N):
        gamma[i+1] = -c[i]/(a[i]*gamma[i]+b[i])
        beta[i+1]  = (d[i]-a[i]*beta[i])/(a[i]*gamma[i]+b[i])

    solution[N] = 0

    for i in range(N,0,-1):
        solution[i-1] = gamma[i] * solution[i] + beta[i]

    return solution[:-1]





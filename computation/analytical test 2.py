__author__ = 'maryhallow'

import numpy
from control.constants import *
import matplotlib.pylab as plt
from other import routines


def Temperature_const(r):
    return numpy.sinh(k_0 * r)/r

def Temperature_transient(r,t):
    return numpy.sin(k_0 * r)/r * numpy.exp(-2*Q_0/C_0*t)

def Q(T):
    return Q_0*T

def init_const(number_of_zones=20):

    global Nzones,R,N

    R = 10

    Nzones = number_of_zones
    N = Nzones - 1

    global dr, dr_b, r, r_b, T_new

    r_b = numpy.linspace(0, R, Nzones + 1)                       # boundaries of cells in the mesh
    r = (r_b[1:] + r_b[:-1]) / 2                                 # cells in the mesh
    dr = r[1:] - r[:-1]                                          # distance between cells (between their middle points)
    dr_b = r_b[1:] - r_b[:-1]                                    # distance between boudaries
    T_new = numpy.ones(Nzones+1)

    global kappa_0,Q_0,k_0

    kappa_0 = 20
    Q_0 = 10
    k_0 = numpy.sqrt(Q_0/kappa_0)

    global A, B, kappa_coef, L_coef, C_coef

    A = numpy.zeros((Nzones,3))
    B = numpy.zeros(Nzones)
    C_coef = 1 / (4*pi*routines.sqr(r) * (dr_b))
    kappa_coef = 4*pi*routines.sqr(r_b[1:-1])/dr
    L_coef = (4*pi*kappa_0*(numpy.sinh(k_0 * R) - R*k_0*numpy.cosh(k_0*R)))

    routines.Tri_diag_matrix_solver_init(Nzones)                 # Initializing tri diagonal matrix solver

def init_transient(number_of_zones=20,delta_t=0.001,time_max=20):

    global Nzones,R,N

    R = 10

    Nzones = number_of_zones
    N = Nzones - 1

    global dr, dr_b, r, r_b,dt,t,t_max

    dt = delta_t
    t = 0.                                                       # initial time in seconds
    t_max =time_max

    r_b = numpy.linspace(0, R, Nzones + 1)                       # boundaries of cells in the mesh
    r = (r_b[1:] + r_b[:-1]) / 2                                 # cells in the mesh
    dr = r[1:] - r[:-1]                                          # distance between cells (between their middle points)
    dr_b = r_b[1:] - r_b[:-1]                                    # distance between boudaries

    global kappa_0,Q_0,k_0,C_0

    kappa_0 = 20
    Q_0 = 10
    C_0 = 1e3
    k_0 = numpy.sqrt(Q_0/kappa_0)

    global A, B, kappa_coef, L_coef, C_coef

    A = numpy.zeros((Nzones,3))
    B = numpy.zeros(Nzones)
    C_coef = dt / (4*pi*routines.sqr(r) * dr_b)
    kappa_coef = 4*pi*routines.sqr(r_b[1:-1])/dr
    L_coef = (4*pi*kappa_0*(numpy.sin(k_0 * R) - R*k_0*numpy.cos(k_0*R)))

    routines.Tri_diag_matrix_solver_init(Nzones)                 # Initializing tri diagonal matrix solver

def main_const():

    global T,T_analytic

    T = T_update_const()
    T_analytic = Temperature_const(r)                            # analytical calculation

def main_transient():

    global T,T_analytic,t

    T = Temperature_transient(r,0)

    while True:                                                  # numerical calculation loop

        t += dt

        if (t>t_max):
            break

        x = T_update_transient(T,T)
        T = T_update_transient(T,x)

    T_analytic = Temperature_transient(r,t-dt)                    # analytical calculation


def T_update_const():                                             # Implicit Euler scheme

    B[N] = -L_coef * C_coef[N]

    A[0,2] = C_coef[0] * kappa_coef[0] * kappa_0
    A[0,0] = 0.

    A[N,2] = 0.
    A[N,0] = C_coef[N] * kappa_coef[N-1] * kappa_0

    A[1:N,2] = C_coef[1:N] * kappa_coef[1:N] * kappa_0
    A[1:N,0] = C_coef[1:N] * kappa_coef[0:N-1] * kappa_0

    A[:,1] = (Q_0 + A[:,2] + A[:,0])

    return routines.Tri_diag_matrix_solver(-A[:,0], A[:,1], -A[:,2],B,Nzones)


def T_update_transient(T,T_func):                                  # Implicit Euler scheme

    B = T - Q(T_func)*dt/C_0
    B[N] -= L_coef * numpy.exp(-2*Q_0/C_0*t) * C_coef[N]/C_0

    A[0,2] = C_coef[0] * kappa_coef[0] * kappa_0 / C_0
    A[0,0] = 0.

    A[N,2] = 0.
    A[N,0] = C_coef[N] * kappa_coef[N-1] * kappa_0 / C_0

    A[1:N,2] = C_coef[1:N] * kappa_coef[1:N] * kappa_0 / C_0
    A[1:N,0] = C_coef[1:N] * kappa_coef[0:N-1] * kappa_0 / C_0

    A[:,1] = (1 + A[:,2] + A[:,0])

    return routines.Tri_diag_matrix_solver(-A[:,0], A[:,1], -A[:,2],B,Nzones)

def test_partition_error_constant():

    partition = numpy.array([16,32,64,128,256,512,1024,2048])
    error = numpy.zeros(len(partition))
    index = 0
    for i in partition:

        init_const(number_of_zones=i)
        main_const()
        error[index] = numpy.max(numpy.abs((T_analytic - T)))
        index += 1

        if(index==2):

            plt.figure(1)
            plt.plot(r,T_analytic,'b-',label='Analytical solution')
            plt.plot(r,T,'r^',label='Numerical solution')
            plt.xlabel('r')
            plt.ylabel('T(r)')
            plt.legend(loc='upper left')

    plt.figure(2)
    plt.grid()
    plt.plot(numpy.log10(partition),numpy.log10(error),'b-')
    plt.xlabel('log (dr)')
    plt.ylabel('log (Error)')

    print(error[:-1]/error[1:])
    print(error)

    plt.show()


def test_partition_error_transient():

    partition = numpy.array([16,32,64,128,256,512,1024])
    error = numpy.zeros(len(partition))
    index = 0

    for i in partition:

        init_transient(number_of_zones=i)
        main_transient()
        error[index] = numpy.max(numpy.abs((T_analytic - T)))
        index += 1
        print(i)

        if(index==2):
            plt.figure(1)
            plt.plot(r,T_analytic,'b-',label='Analytical solution')
            plt.plot(r,T,'r^',label='Numerical solution')
            plt.xlabel('r')
            plt.ylabel('T(r)')
            plt.legend(loc='upper left')

    plt.figure(2)
    plt.grid()
    plt.plot(numpy.log10(partition),numpy.log10(error),'b-')
    plt.xlabel('log (dr)')
    plt.ylabel('log (Error)')

    print(error[:-1]/error[1:])
    print(error)

    plt.show()


def test_time_error_transient():

    delta_t = numpy.array([0.0125,0.025,0.05,0.1,0.2,0.4,0.8,1.6])[::-1]
    error = numpy.zeros(len(delta_t))
    index = 0
    for i in delta_t:

        init_transient(number_of_zones=500,delta_t=i)
        main_transient()
        error[index] = numpy.max(numpy.abs((T_analytic - T)))
        index += 1

        if(index==2):
            plt.figure(1)
            plt.plot(r,T_analytic,'b-',label='Analytical solution')
            plt.plot(r,T,'r^',label='Numerical solution')
            plt.xlabel('r')
            plt.ylabel('T(r)')
            plt.legend(loc='upper left')

    plt.figure(2)
    plt.grid()
    plt.plot(numpy.log10(delta_t),numpy.log10(error),'b-')
    plt.xlabel('log (dt)')
    plt.ylabel('log (Error)')

    print(error[:-1]/error[1:])
    print(error)

    plt.show()

test_partition_error_constant()
#test_time_error_transient()
#test_partition_error_transient()



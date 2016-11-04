import numpy
from control.constants import *
import matplotlib.pylab as plt
from other import plot
from other import routines
import timeit

def C(a):                                                        # heat capacity function for analytical test
    return C_0 * a

def kappa(a):                                                    # thermal conductivity function for analytical test
    return kappa_0 * a

def Q(a):                                                        # neutrino emissivity function for analytical test
    return Q_0 * a * a

def Temperature(r, t):                                           # Analytical solution
    return T_0 * numpy.exp(-gamma * t) * numpy.sqrt((numpy.sin(k * r)/(k * r)))

def init(delta_t=1e6,number_of_zones=20,max_time=2e10):

    global Nzones, T_0, R, dt, t, t_max

    T_0 = 1e10                                                   # initial temperature in K
    R = 1.e6                                                     # star radius in cm
    dt = delta_t
    t = 0.                                                       # initial time in seconds
    t_max = max_time

    global Nzones, N, timestep_number

    Nzones = number_of_zones                                     # number of zones star will be logarithmically divided into
    N = Nzones - 1                                               # We introduce N to prevent us from writing 'Nzones-1' a lot of time in code below
    timestep_number = 10                                         # data of the simulation will be revealed timestep_number times

    global dr, dr_b, r, r_b

    r_b = numpy.logspace(-1, numpy.log10(R), Nzones + 1)         # boundaries of cells in the mesh
    r = (r_b[1:] + r_b[:-1]) / 2                                 # cells in the mesh
    dr = r[1:] - r[:-1]                                          # distance between cells (between their middle points)
    dr_b = r_b[1:] - r_b[:-1]                                    # distance between boudaries

    global k, C_0, kappa_0, Q_0, gamma

    k = pi/R                                                     # minimal wavenumber of a thermal wave
    C_0 = 1e12                                                   # heat capacity coefficient
    kappa_0 = 1e12                                               # thermal conductivity coefficient
    Q_0 = 10                                                     # neutrino emissivity coefficient
    gamma = (0.5 * kappa_0 * k * k + Q_0)/C_0                    # coefficient in exp. decay factor [exp(-gamma t)]

    global A, B, C_coef, kappa_coef, L_coef

    A = numpy.zeros((Nzones,3))
    B = numpy.zeros(Nzones)
    C_coef = dt / (4*pi*routines.sqr(r) * (dr_b))
    kappa_coef = 4*pi*routines.sqr(r_b[1:-1])/dr
    L_coef = 2*pi*kappa_0*T_0*T_0/k * (numpy.sin(k*r_b[-1])-k*r_b[-1]*numpy.cos(k*r_b[-1]))

    routines.Tri_diag_matrix_solver_init(Nzones)                 # Initializing tri diagonal matrix solver

def main():

    global t, T, T_analytic

    time_snapshot = numpy.logspace(numpy.log10(t+1), numpy.log10(t_max), timestep_number)
    T = Temperature(r, t)
    snapshot_counter = 0

    while True:                                                  # numerical calculation loop

        t += dt

        if (t>t_max):
            print(' ')
            break

        x = T_update(T,T)
        T = T_update(T,x)

        if t >= time_snapshot[snapshot_counter]:
            print('T difference after 1 iteration:', numpy.max(numpy.abs(T-x)))
            snapshot_counter += 1
            print('current time = %-6.2f years\n' % (t/yrtosec))

    T_analytic = Temperature(r,t-dt)                             # analytical calculation

def T_update(T,T_func):                                          # Implicit Euler scheme

    T_b = (T_func[1:] + T_func[:-1]) / 2
    B = T - dt*Q(T_func)/C(T_func)
    B[N] -= C_coef[N] * L_coef * numpy.exp(-2*gamma*t)/C(T_func[N])

    A[0,2] = C_coef[0] / C(T_func[0]) * kappa_coef[0] * kappa(T_b[0])
    A[0,0] = 0.

    A[N,2] = 0.
    A[N,0] = C_coef[N] / C(T_func[N]) * kappa_coef[N-1] * kappa(T_b[N-1])

    A[1:N,2] = C_coef[1:N] / C(T_func[1:N]) * kappa_coef[1:N] * kappa(T_b[1:N])
    A[1:N,0] = C_coef[1:N] / C(T_func[1:N]) * kappa_coef[0:N-1] * kappa(T_b[0:N-1])
    B[1:N] = T[1:N] - dt*Q(T_func[1:N])/C(T_func[1:N])

    A[:,1] = (1 + A[:,2] + A[:,0])

    return routines.Tri_diag_matrix_solver(-A[:,0], A[:,1], -A[:,2],B,Nzones)


def test_starter():                                              # use this function in order to launch test simulation

    print ('Analytical test starts... \n')
    init(number_of_zones=80,delta_t=1e6)
    print ('Initialization completed \n')
    print ('Initial parameters:')
    print ('----------------------')
    print ('Nzones: %d' % (Nzones))
    print ('t_step: %1.1e sec' % (dt))
    print ('t_max:  %1.1e years' % (t_max))
    print ('Radius:  %1.1e cm' % (R))
    print ('T_init: %1.1e K \n' % (T_0))
    print ('Implicit Euler scheme is used for numerical computation \n')

    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()

    print (" ")
    print ('calculation time: %3.5f sec' % (stop - start))
    print ('max error: %3.5f' % numpy.max(numpy.abs((T_analytic  - T)/T_analytic )))

    plot.double_plot(r, T_analytic , r, T, 'T(r)', 1,'analytically','numerically','r, cm', 'T, K', '-b','oy')

def test_partition_error():

    partition = numpy.array([4,8,16,32,64,128,256,512,1024])
    error = numpy.zeros(len(partition))
    index = 0
    for i in partition:

        init(number_of_zones=i)
        main()
        error[index] = numpy.max(numpy.abs((T_analytic - T)))
        index += 1
        print('Number of spherical layers = ', i,'\n')

        if(index==(len(partition)-1)):

            plt.figure(1)
            plt.plot(r,T_analytic,'b-',label='Analytical solution')
            plt.plot(r,T,'r^',label='Numerical solution')
            plt.xlabel('r')
            plt.ylabel('T(r)')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc='lower left')

    plt.figure(2)
    plt.grid()
    plt.plot(numpy.log10(partition),numpy.log10(error),'b-')
    plt.xlabel('log (dr)')
    plt.ylabel('log (Error)')

    print('Error ratios for 2**n spacial partition')
    print(error[:-1]/error[1:])
    print(' ')
    print('Errors')
    print(error)

    plt.show()


def test_time_error():

    delta_t = (numpy.array([0.5,1,2,4,8,16,32,64,128])*1e5)
    error = numpy.zeros(len(delta_t))
    index = 0
    for i in delta_t:

        init(number_of_zones=100,delta_t=i,max_time=1e9)
        main()
        error[index] = numpy.max(numpy.abs((T_analytic - T)))
        index += 1

        if(index==4):

            plt.figure(1)
            plt.plot(r,T_analytic,'b-',label='Analytical solution')
            plt.plot(r,T,'r^',label='Numerical solution')
            plt.xlabel('r')
            plt.ylabel('T(r)')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend(loc='lower left')

    plt.figure(2)
    plt.grid()
    plt.plot(numpy.log10(delta_t),numpy.log10(error),'b-')
    plt.xlabel('log (dt)')
    plt.ylabel('log (Error)')

    print('Error ratios for dt*(2**n) time steps')
    print(error[:-1]/error[1:])
    print(' ')
    print('Errors')
    print(error)

    plt.show()

#test_time_error()
#test_partition_error()
test_starter()



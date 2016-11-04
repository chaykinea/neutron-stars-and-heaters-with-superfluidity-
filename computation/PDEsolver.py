#     Neutron star thermal evolution
# --------------------------------------
#            PDEsolver.py
# --------------------------------------
# This module provides numerical
# computation algorithms for
# solving PDE (heat equation
# for neutron star).

from __future__ import division
from other import routines
from sys import exit
from physics.physics import *
import matplotlib.pyplot as plt

def init():                                                  

    global Rho, r_b, r, rho_r,dr,dr_b

    # --------------------------------------- mesh ----------------------------------------
    # -------------------------------------------------------------------------------------
    if (rho_vary):
        Rho = numpy.logspace( numpy.log10(rho_max),
                              numpy.log10(rho_min), Nzones + 1)  # density on cells' boundaries
    else:
        Rho = numpy.logspace( numpy.log10(loaddata.star_model()[2,3]),
                              numpy.log10(loaddata.star_model()[-2,3]), Nzones + 1)
    r_b = loaddata.radii(Rho)                                    # boundaries of the cells
    r = (r_b[1:] + r_b[:-1])/2                                   # average cell's radius
    dr = r[1:] - r[:-1]                                          # distance between cells
    dr_b = r_b[1:] - r_b[:-1]                                    # distance between cells' boundaries
    rho_r = loaddata.rho(r)                                      # cells' average density

    global T, t, dt, redshift, N

    T = T_0*numpy.ones(Nzones)                                   # initial temperature of neutron star
    t = t_0                                                      # initial time
    dt = dt_0                                                    # initial time step
    redshift = relativity_sqrt(r_b[-1])                          # redshift at infinity
    N = Nzones - 1                                               # We introduce N to prevent us from writing 'Nzones-1' a lot of times in code below

    print ('Initialization of initial parameters completed \n')
    print ('Initial parameters:')
    print ('----------------------')
    print ('Nzones:    %d' % (Nzones))
    print ('t_step:    %1.1e sec' % (dt))
    print ('t_max:     %1.3e years\n' % (t_max))
    print ('NS Radius: %1.4e cm' % (r_b[-1]))
    print ('NS Mass:   %1.4e MSun' % (loaddata.mass(r_b[-1])/MSun))
    print ('rho_max:   %1.4e gm/cm-3' % (rho_r[0]))
    print ('rho_min:   %1.4e gm/cm-3' % (rho_r[-1]))
    print ('redshift:  %5.4f \n' % (redshift))
    print ('Magnetic field:  %1.2e Gauss' % (MagField ))
    print ('Surface gravity: %1.2e cm/s^2\n' %(loaddata.g_surface()*1e14))
    if(SUPERFLUIDITY):
        print ('Superfluidity is ON')
    else:
        print ('Superfluidity is OFF')
    print ('T_init:   %1.1e K' % (T_0))
    print ('T_min:    %1.1e K \n' % (T_min))

    if(regime==0):
        print ('Implicit Euler scheme is used for t < %-8.1e years'% t_iso)
        print ('If t > %-8.1e years, isothermal scheme is used.\n' % t_iso)
    else:
        print ('Implicit Euler scheme is used for numerical computation. \n')


def time_derivative_init():               

    global A, B, C_coef, kappa_coef, L_ph_coef, T_loc

    A = numpy.zeros((Nzones,3))
    B = numpy.zeros(Nzones)

    C_coef     = relativity_sqrt(r)/(4*pi*routines.sqr(r)*dr_b)
    kappa_coef = 4*pi*routines.sqr(r_b[1:-1])*numpy.exp(loaddata.Phi(r_b[1:-1]))*relativity_sqrt(r_b[1:-1])/dr
    L_ph_coef  = 4*pi*r_b[-1]*r_b[-1]*sigma*numpy.exp(2*loaddata.Phi(r_b[-1]))

    T_loc = numpy.exp(-loaddata.Phi(r_b[-1]))
    routines.Tri_diag_matrix_solver_init(Nzones)


def time_derivative_iso_regime_init():                   

    global coef_L_photon_total, coef_C_total, coef_Q_total

    coef_L_photon_total = 4*pi*r_b[-1]*r_b[-1]*sigma*numpy.exp(2*loaddata.Phi(r_b[-1]))
    coef_Q_total = 4*pi*r*r*(r_b[1:] - r_b[:-1])/relativity_sqrt(r)
    coef_C_total = 4*pi*r*r*(r_b[1:] - r_b[:-1])/relativity_sqrt(r)

def relativity_sqrt(r):                                    
    return numpy.sqrt(1 - 2*G*loaddata.mass(r)/(( c ** 2 )*r))


def T_update_1(T,T_func,dt):                            

    T_b = (T_func[1:] + T_func[:-1])/2
    C_temp = C(T_func,rho_r)
    kappa_temp = k(T_b,Rho[1:N+1])

    B = T - dt*Q(T_func,rho_r)/C_temp
    B[N] -= dt*C_coef[N]/C_temp[N] * L_ph_coef*numpy.power(loaddata.T_e(T_func[N]*T_loc),4)

    A[0,2] = dt*C_coef[0]/C_temp[0]*kappa_coef[0]*kappa_temp[0]
    A[0,0] = 0.

    A[N,2] = 0.
    A[N,0] = dt*C_coef[N]/C_temp[N]*kappa_coef[N-1]*kappa_temp[N-1]

    A[1:N,2] = dt*C_coef[1:N]/C_temp[1:N]*kappa_coef[1:N]*kappa_temp[1:N]
    A[1:N,0] = dt*C_coef[1:N]/C_temp[1:N]*kappa_coef[0:N-1]*kappa_temp[0:N-1]

    A[:,1] = (1 + A[:,2] + A[:,0])

    for i in range(1,N):
        if B[i] < 0:
            print ('dt is too large. Simulation is terminated.')
            if (t<0.1*yrtosec):
                print('current dt value = %-8.1e sec' % dt)
            else:
                print('current dt value = %-8.1e years' % (dt/yrtosec))
            return time_step_control()

    return routines.Tri_diag_matrix_solver(-A[:,0], A[:,1], -A[:,2],B,Nzones)


def dT_isothermal(T):                                  

    T_local = T*T_loc
    Q_total = numpy.sum(Q(T, rho_r)*coef_Q_total)
    C_total = numpy.sum(C(T, rho_r)*coef_C_total)
    L_photon_total = coef_L_photon_total*(loaddata.T_e(T_local) ** 4)

    return (Q_total + L_photon_total)/C_total


def solve_PDE(output_1, output_2):                          
                                                             
    global t, dt, T, time_step_counter, T_save, t_save

    Temperature_profile_data = numpy.zeros((len(time_points)+1,Nzones))  # creating output file (for temperature profiles)
    Temperature_profile_data[0,:] = rho_r                                # writing densities into the output file

    snapshot = numpy.logspace(-1, numpy.log10(t_max*yrtosec), N_output)  # it tells at what time to make a snapshot
    snapshot_counter  = 0                                                # to count the number of the made snapshots
    time_step_counter = 0                                                # ot tells how many times we have changed time step dt

    while snapshot_counter < N_output:                                   # N_output = total number of snapshots

        t += dt

        if regime == 0:
            if (t/yrtosec < t_iso):
                x = T_update_1(T,T,dt)
                T = T_update_1(T,x,dt)
            else:
                T[-1] = T[-1] - dT_isothermal(T[-1])*dt
                T = T[-1]*numpy.ones(Nzones)
        else:
            x = T_update_1(T,T,dt)
            T = T_update_1(T,x,dt)

        if t >= snapshot[snapshot_counter]:

            snapshot_counter += 1

            print('current time   = %-8.3e years' % ((t)/yrtosec))
            print('number of step = %-4d' % snapshot_counter)

            if (t<0.1*yrtosec):
                print('time step      = %-8.1e sec' % dt)
            else:
                print('time step      = %-8.1e years' % (dt/yrtosec))

            print('T_surface_prop = %-8.1e K' % (T[-1]*numpy.exp(-loaddata.Phi(r[-1]))))
            if regime!=0:
               print('T_diff[1 iter] = %-8.1e K \n' % numpy.max(numpy.abs(x-T)))

            data = numpy.vstack([loaddata.T_e(T[-1]*numpy.exp(-loaddata.Phi(r[-1])))*redshift,
                                 t/yrtosec,
                                 T[-1]*numpy.exp(-loaddata.Phi(r_b[-1])),
                                 numpy.sum(4.*pi*numpy.power(r,2)*dr_b*Q(T,rho_r)/relativity_sqrt(r)),
                                 4.*pi*numpy.power(r_b[-1],2)*sigma*numpy.power(loaddata.T_e(T[-1]*numpy.exp(-loaddata.Phi(r[-1]))),4)*redshift*redshift])

            numpy.savetxt(output_1, data.T, fmt='%1.6e')

        if time_step_counter<len(time_points):

            if t/yrtosec > time_points[time_step_counter]:               

                dt = time_steps[time_step_counter]                   
                Temperature_profile_data[time_step_counter+1,:] = T      
                time_step_counter += 1
                T_save = T
                t_save = t

        if t > t_max*yrtosec or T[-1]*numpy.exp(-loaddata.Phi(r_b[-1])) < T_min:
            print ('Simulation is terminated. One of the critical conditions is achieved [t > t max or T < T min.]\n')
            break

    numpy.savetxt(output_2, Temperature_profile_data.T, fmt='%1.6e')

def time_step_control():

    global i, dt, time_steps, time_points, t, time_step_counter
    
    time_step_counter -= 1
    time_points[time_step_counter] *= 2
    time_steps[time_step_counter] /= 2

    print ('Time step control works.\n')
    print ('Time old = %-8.1e sec' % t)
    print ('dt old   = %-8.1e sec\n' % dt)

    t = t_save
    dt = time_steps[time_step_counter-1]

    print ('Time new = %-8.1e sec' % t)
    print ('dt new   = %-8.1e sec\n' % dt)

    return T_save


# ---------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------HEAT-SOURCE-----------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------


def source_initialization(duration,power,Left,Right):

    R_L_computator(Left,Right)

    global S, D

    D = duration*yrtosec
    S = numpy.zeros(Nzones)

    for i in range(L,R):
        S[i] = power
        
        
def R_L_computator(Left,Right):

    global L, R

    L, R = 0, 0

    for i in range(0,N+1):
        if (Left >= rho_r[i]):
            L = i
            break

    for i in range(N,-1,-1):
        if (Right <= rho_r[i]):
            R = i
            break

    if L>R or L<0 or R>N or R-L==0:
        print ('Failed to properly determine left and right boundary of the source. '
               'Simulation is terminated. Please check the source input parameters.\n')
        print ('L = %i, R = %i, N = %i\n' % (L, R, N))
        exit(0)

    print ('Left source boundary (number of element) = %i, Right source boundary = %i, '
           'Number of elements in the star = %in' % (L, R, N+1))

    R += 1                                                    


def source_visualisation(source_number):

    temp = numpy.ones(N)

    plt.plot(rho_r[0:L], temp[0:L], 'bo', linewidth=3, label='no source')
    plt.plot(rho_r[L:R], temp[L:R], 'ro', linewidth=7, label='source')
    plt.plot(rho_r[R:N], temp[R:N], 'bo', linewidth=3)
    plt.xscale('log')
    plt.yscale('linear',fontsize=24)
    plt.xlabel('Density',fontsize=24)
    plt.title('Source placement in the star',fontsize=24)
    plt.legend(loc=4)
    plt.savefig('output/source_placement_' + str(source_number) + 'pdf')


def time_step_corrector(dt, t, time_point, error, error_min):

    number_of_iterations = 0

    while (t/yrtosec < time_point and (t+dt)/yrtosec > time_point):
        number_of_iterations += 1
        error -= 1
        dt = time_steps[error]

    if (number_of_iterations>0):

        print ('Time step corrector works:\n')
        print ('Number of iterations = %i' % number_of_iterations)
        print ('Error = %i, Error_min = %i\n' % (error,error_min))

    return dt, error

def T_update_source(T,T_func,dt):                                    

    T_b = (T_func[1:] + T_func[:-1])/2
    C_temp = C(T_func,rho_r)
    kappa_temp = k(T_b,Rho[1:N+1])

    B = T - dt*(Q(T,rho_r)-S*f(t))/C_temp
    B[N] -= dt*C_coef[N]/C_temp[N] * L_ph_coef*numpy.power(loaddata.T_e(T_func[N]*T_loc),4)

    A[0,2] = dt*C_coef[0]/C_temp[0]*kappa_coef[0]*kappa_temp[0]
    A[0,0] = 0.

    A[N,2] = 0.
    A[N,0] = dt*C_coef[N]/C_temp[N]*kappa_coef[N-1]*kappa_temp[N-1]

    A[1:N,2] = dt*C_coef[1:N]/C_temp[1:N]*kappa_coef[1:N]*kappa_temp[1:N]
    A[1:N,0] = dt*C_coef[1:N]/C_temp[1:N]*kappa_coef[0:N-1]*kappa_temp[0:N-1]

    A[:,1] = (1 + A[:,2] + A[:,0])

    for i in range(1,N):
        if B[i] < 0:
            print ('dt is too large. Simulation is terminated.')
            if (t<0.1*yrtosec):
                print('current dt value = %-8.1e sec' % dt)
            else:
                print('current dt value = %-8.1e years' % (dt/yrtosec))
            return time_step_control()

    return routines.Tri_diag_matrix_solver(-A[:,0], A[:,1], -A[:,2],B,Nzones)


def solve_PDE_with_source(source_number,source_number_max): 

    global t, dt, T, time_step_counter_source, time_step_counter, T_save, t_save, D,S

    Temperature_profile_data_source_on      = numpy.zeros((len(t_source_points)+1,Nzones))
    Temperature_profile_data_source_on[0,:] = rho_r

    snapshot   = numpy.logspace(-1, numpy.log10(turn_on_time*yrtosec), N_output)
    snapshot_2 = numpy.logspace(numpy.log10(turn_on_time*yrtosec+1), numpy.log10((turn_on_time+10)*yrtosec), N_output)
    snapshot_3 = numpy.logspace(numpy.log10((turn_on_time+10)*yrtosec), numpy.log10(t_source_max*yrtosec), N_output)
    snapshot   = numpy.concatenate([snapshot,snapshot_2[1:],snapshot_3[1:]])

    snapshot_counter         = 0
    time_step_counter        = 0
    dt_source_on_off_switch  = 0
    time_step_counter_source = 0

    temp_output_file =  open('output/file_' + str(source_number) + '_cooling.dat', 'wb')

    while snapshot_counter < len(snapshot):

        t += dt

        if (time_step_counter>error_min):
            if(t > (turn_on_time*yrtosec-1.2*dt) and t <= turn_on_time*yrtosec):
                dt_source_on_off_switch = 1
                dt,time_step_counter = time_step_corrector(dt, t, turn_on_time, time_step_counter, error_min)

        if t >= turn_on_time*yrtosec:
            x = T_update_source(T,T,dt)
            T = T_update_source(T,x,dt)
        else:
            x = T_update_1(T,T,dt)
            T = T_update_1(T,x,dt)

        if t >= snapshot[snapshot_counter]:

            snapshot_counter += 1

            print('current time   = %-8.3e years' % ((t)/yrtosec))
            print('number of step = %-4d' % snapshot_counter)
            print('time step      = %-8.3e years' % (dt/yrtosec))
            print('T_surface_prop = %-8.1e K' % (T[-1]*numpy.exp(-loaddata.Phi(r[-1]))))
            print('T_diff[1 iter] = %-8.1e K' % numpy.max(numpy.abs(x-T)))
            print('Source Power   = %-8.3e '  % numpy.sum(4.*pi*numpy.power(r,2)*dr_b*S*f(t)))
            print('Surface Power  = %-8.3e '  % (4.*pi*numpy.power(r_b[-1],2)*sigma*numpy.power(loaddata.T_e(T[-1]*numpy.exp(-loaddata.Phi(r[-1]))),4)))
            print('Source N       = %-4i out of %-4i\n' % (source_number, source_number_max-1))

            data = numpy.vstack([loaddata.T_e(T[-1]*numpy.exp(-loaddata.Phi(r[-1])))*redshift,
                                 t/yrtosec,
                                 T[-1]*numpy.exp(-loaddata.Phi(r_b[-1])),
                                 numpy.sum(4.*pi*numpy.power(r,2)*dr_b*numpy.exp(2*loaddata.Phi(r))*S*f(t)/relativity_sqrt(r)),
                                 4.*pi*numpy.power(r_b[-1],2)*sigma*numpy.power(loaddata.T_e(T[-1]*numpy.exp(-loaddata.Phi(r[-1]))),4)*redshift*redshift])

            numpy.savetxt(temp_output_file, data.T, fmt='%1.8e')

        if(not dt_source_on_off_switch):
            if time_step_counter<len(time_points):
                if t/yrtosec > time_points[time_step_counter]:          
                    dt = time_steps[time_step_counter]                            
                    T_save = T
                    t_save = t
                    time_step_counter += 1
        else:
            if time_step_counter_source<len(t_source_points):
                if t/yrtosec > t_source_points[time_step_counter_source]:
                    dt = t_source_steps[time_step_counter_source]
                    Temperature_profile_data_source_on[time_step_counter_source+1,:] = T
                    time_step_counter_source += 1
                    T_save = T
                    t_save = t

        if (t > t_source_max*yrtosec or T[-1]*numpy.exp(-loaddata.Phi(r_b[-1])) < T_min):
            print ('Simulation is terminated. One of the critical conditions is achieved [t > t max or T < T min.]\n')
            break

    temp_output_file.close()

    output_file_source = open('output/file_' + str(source_number) + '.dat','wb')
    numpy.savetxt(output_file_source, Temperature_profile_data_source_on.T, fmt='%1.6e')
    output_file_source.close()

def re_init():

    global T, t, dt, redshift, N

    T = T_0*numpy.ones(Nzones)
    t = t_0
    dt = dt_0


def time_step_control_source():

    global time_step_counter_source, dt, t_source_steps, t_source_points, t

    time_step_counter_source -= 1
    t_source_points[time_step_counter_source] *= 1.2
    t_source_steps[time_step_counter_source] /= 2

    print ('Time step control works.\n')
    print ('Time old = %-8.3e years' % (t/yrtosec))
    print ('dt old   = %-8.3e sec\n' % dt)

    t = t_save
    dt = t_source_steps[time_step_counter_source-1]

    print ('Time new = %-8.3e years' % (t/yrtosec))
    print ('dt new   = %-8.3e sec\n' % dt)
    print ('time parameter = %2i\n' % time_step_counter_source)

    if(time_step_counter_source<1):
        print ('Simulation is terminated. dt is too large (((')
        exit(0)

    return T_save

def f(t):

    if(t>=turn_on_time*yrtosec and t<=(turn_on_time*yrtosec+D)):
        return 1.
    else:
        return 0.

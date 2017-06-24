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
        Rho = numpy.logspace( numpy.log10(loaddata.star_model()[2,3]),
                              numpy.log10(rho_min), Nzones + 1)  # density on cells' boundaries
    else:
        #Rho = numpy.logspace( numpy.log10(loaddata.star_model()[2,3]),
        #                      numpy.log10(loaddata.star_model()[-2,3]), Nzones + 1)
        Rho = numpy.logspace(numpy.log10(loaddata.star_model()[0, 3]),
                         numpy.log10(loaddata.star_model()[-1, 3]), Nzones + 1)
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
    print ('Magnetic field:  %1.2e Gauss' % (MagField))
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

    global A, B, C_coef, kappa_coef, L_ph_coef, T_loc, H_loc

    A = numpy.zeros((Nzones,3))
    B = numpy.zeros(Nzones)

    C_coef     = relativity_sqrt(r)/(4*pi*routines.sqr(r)*dr_b)
    kappa_coef = 4*pi*routines.sqr(r_b[1:-1])*numpy.exp(loaddata.Phi(r_b[1:-1]))*relativity_sqrt(r_b[1:-1])/dr
    L_ph_coef  = 4*pi*r_b[-1]*r_b[-1]*sigma*numpy.exp(2*loaddata.Phi(r_b[-1]))

    T_loc = numpy.exp(-loaddata.Phi(r_b[-1]))
    H_loc = numpy.exp(2*loaddata.Phi(r))
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

def source_power(model_parm=model_param):

    global H_f

    r_0 = numpy.array([1.20413015e+06, 1.21000074e+06 , 1.21203074e+06, 1.21269574e+06,  1.21353249e+06, 1.21448985e+06, 1.21486413e+06, 1.21161929e+06])
    sigma_g = numpy.array([4.38291763e+03, 3.97737319e+03, 3.78693639e+03, 3.71309773e+03, 3.60523861e+03, 3.43057815e+03, 3.26194154e+03, 2.68581589e+03])
    const_1 = numpy.array([1.66225679807e16, 1.86808408312e16, 1.98349866253e16, 2.03226874478e16, 2.10820556759e16, 2.2441279112e16, 2.39298439565e16 , 3.80131342314e16])
    const_2 = numpy.array([6.39198181e+04, 6.44242849e+04, 6.47562763e+04, 6.49105274e+04, 6.51666607e+04,  6.56611002e+04, 6.62462312e+04,6.64638989e+04 ])
    model_param_arr = numpy.array([1.40, 1.50, 1.55, 1.57, 1.60, 1.65, 1.70, 1.85])

    r_0_f = scipy.interpolate.interp1d(model_param_arr, r_0)
    sigma_g_f = scipy.interpolate.interp1d(model_param_arr, sigma_g)
    const_1_f = scipy.interpolate.interp1d(model_param_arr, const_1)
    const_2_f = scipy.interpolate.interp1d(model_param_arr, const_2)
    H_f = lambda r: const_1_f(model_parm) + const_2_f(model_parm) * 1e17 * numpy.exp(-(r-r_0_f(model_parm))**2/sigma_g_f(model_parm)**2 / 2)/numpy.sqrt(2*numpy.pi*sigma_g_f(model_parm)**2 ) # M_dot = 1e-9 M_solar/yr


def source_initialization(duration,power,Left,Right):

    R_L_computator(Left,Right)

    global S, D, H_0

    #D = duration*yrtosec
    D = duration
    S = numpy.zeros(Nzones)
    H_0 = power

    for i in range(L,R):
        #  S[i] = 1e17  steps from 1 to 6
        S[i] = H_0 * H_f(r[i]) #  step 7

        
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
           'Number of elements in the star = %i' % (L, R, N+1))

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
    plt.savefig('output/source_placement_' + str(source_number) + '.pdf',format='pdf')


def time_step_corrector(dt, t, time_point, error, error_min):

    number_of_iterations = 0

    while (t/yrtosec < time_point and (t+dt)/yrtosec > time_point):
        number_of_iterations += 1
        error -= 1
        print(error)
        dt = time_steps[error]
        if error <= 1:
            break

    if (number_of_iterations>0):

        print ('Time step corrector works:\n')
        print ('Number of iterations = %i' % number_of_iterations)
        print ('Error = %i, Error_min = %i\n' % (error,error_min))

    return dt, error

def T_update_source(T,T_func,dt):                                    

    T_b = (T_func[1:] + T_func[:-1])/2
    C_temp = C(T_func,rho_r)
    kappa_temp = k(T_b,Rho[1:N+1])

    B = T - dt*(Q(T,rho_r)-S*f(t)*H_loc)/C_temp
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
            print('dt is too large. Simulation is terminated.')
            if (t<0.1*yrtosec):
                print('current dt value = %-8.1e sec' % dt)
            else:
                print('current dt value = %-8.1e years' % (dt/yrtosec))
            return time_step_control()

    return routines.Tri_diag_matrix_solver(-A[:,0], A[:,1], -A[:,2],B,Nzones)


def solve_PDE_with_source(source_number,source_number_max): 

    global t, dt, T, time_step_counter_source, time_step_counter, T_save, t_save, D,S

    Temperature_profile_data_source_on      = numpy.zeros((len(t_points_save_data)+1,Nzones))
    Neutrino_profile_data_source_on         = numpy.zeros((len(t_points_save_data)+1,Nzones))
    Flux_profile_data_source_on             = numpy.zeros((len(t_points_save_data)+1,Nzones-1))

    Temperature_profile_data_source_on[0,:] = rho_r
    Neutrino_profile_data_source_on[0,:]    = rho_r
    Flux_profile_data_source_on[0,:]        = Rho[1:-1]

    snapshot   = numpy.logspace(-1, numpy.log10(turn_on_time*yrtosec), N_output*2)
    snapshot_2 = numpy.logspace(numpy.log10(turn_on_time*yrtosec), numpy.log10((turn_on_time+10)*yrtosec), N_output*4) 
    snapshot_3 = numpy.logspace(numpy.log10((turn_on_time+10)*yrtosec), numpy.log10(t_source_max*yrtosec), N_output*2) 
    snapshot   = numpy.concatenate([snapshot,snapshot_2[1:],snapshot_3[1:]])

    snapshot_counter         = 0
    time_step_counter        = 0
    dt_source_on_off_switch  = 0
    time_step_counter_source = 0

    data_save_counter        = 0

    temp_output_file =  open('output/cooling' + name + str(source_number) + '.dat', 'wb')

    while snapshot_counter < len(snapshot):

        t += dt

        if (time_step_counter>error_min):
            if(t > (turn_on_time*yrtosec-1.2*dt) and t <= turn_on_time*yrtosec):
                dt_source_on_off_switch = 1
                dt,time_step_counter = time_step_corrector(dt, t, turn_on_time, time_step_counter, error_min)

        if t >= turn_on_time_0*yrtosec:
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
                                 T[-1],
                                 numpy.sum(4.*pi*numpy.power(r,2)*dr_b*numpy.exp(2*loaddata.Phi(r))*S*f(t)/relativity_sqrt(r)),
                                 numpy.sum(4.*pi*numpy.power(r,2)*dr_b*Q(T,rho_r)/relativity_sqrt(r)),
                                 4.*pi*numpy.power(r_b[-1],2)*sigma*numpy.power(loaddata.T_e(T[-1]*numpy.exp(-loaddata.Phi(r[-1]))),4)*redshift*redshift,
                                 numpy.sum(4.*pi*numpy.power(r,2)*dr_b*C(T,rho_r)/relativity_sqrt(r)*T)])

            numpy.savetxt(temp_output_file, data.T, fmt='%1.10e')

        if(not dt_source_on_off_switch):
            if time_step_counter<len(time_points):
                if t/yrtosec > time_points[time_step_counter]:

                    dt = time_steps[time_step_counter]
                    time_step_counter += 1

                    T_save = T
                    t_save = t

        else:
            if time_step_counter_source<len(t_source_points):
                if t/yrtosec > t_source_points[time_step_counter_source]:

                    dt = t_source_steps[time_step_counter_source]
                    time_step_counter_source += 1

                    T_save = T
                    t_save = t

        if data_save_counter<len(t_points_save_data):
            if t/yrtosec >= t_points_save_data[data_save_counter]:

                Temperature_profile_data_source_on[data_save_counter+1,:] = T
                Neutrino_profile_data_source_on[data_save_counter+1,:]    = Q(T,rho_r)*numpy.exp(-2*loaddata.Phi(r))
                Flux_profile_data_source_on[data_save_counter+1,:]        = -kappa_coef*k(((T[:-1] + T[1:])/2),Rho[1:-1])*(T[1:] - T[:-1])

                data_save_counter += 1


        if (t > t_source_max*yrtosec or T[-1]*numpy.exp(-loaddata.Phi(r_b[-1])) < T_min):
            print ('Simulation is terminated. One of the critical conditions is achieved [t > t max or T < T min.]\n')
            break

    temp_output_file.close()

    output_file_source  = open('output/temperature' + name + str(source_number) + '.dat','wb')
    output_file_source2 = open('output/neutrino' + name + str(source_number) + '.dat','wb')
    output_file_source3 = open('output/flux' + name + str(source_number) + '.dat','wb')

    numpy.savetxt(output_file_source , Temperature_profile_data_source_on.T, fmt='%1.6e')
    numpy.savetxt(output_file_source2,    Neutrino_profile_data_source_on.T, fmt='%1.6e')
    numpy.savetxt(output_file_source3,        Flux_profile_data_source_on.T, fmt='%1.6e')

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

'''
def f(t):

    if(t>=turn_on_time*yrtosec and t<=(turn_on_time*yrtosec+D)):
        return numpy.power(numpy.sin(numpy.pi*(t - turn_on_time*yrtosec)/D),2)
    else:
        return 0.
'''

def params(time,source_name):

    global name, turn_on_time_0, turn_on_time, t_source_max, t_source_points, t_source_steps, t_points_save_data


    # ENERGY RELEASES
    '''
    turn_on_time_0 = time
    turn_on_time = turn_on_time_0 + 30000
    t_source_max = 40000+5000
    #t_source_max = turn_on_time + 100
    '''

    # ACCRETION
    turn_on_time_0 = 1e3
    turn_on_time = turn_on_time_0 + 1.45
    t_source_max = turn_on_time + 5
    name = source_name

    t_source_points     = numpy.array([turn_on_time +0.01 ,turn_on_time + 0.1,turn_on_time + 1.,turn_on_time + 5., turn_on_time + 10.,
                                   turn_on_time + 50., turn_on_time + 100., turn_on_time + 500., turn_on_time + 1000.])       # in years
    t_source_steps      = [time_steps[error_min], time_steps[error_min]*10  , time_steps[error_min]*20 , time_steps[error_min]*50,
                           time_steps[error_min]*100 , time_steps[error_min]*300 , 5.e6, 5.e7, 5.e7]

    t_points_save_data  = numpy.array([turn_on_time + 0.00001, turn_on_time + 0.01, turn_on_time + 0.1, turn_on_time + 0.5, turn_on_time + 1., turn_on_time + 1.5,
                                       turn_on_time + 2.,   turn_on_time + 2.5,  turn_on_time + 3,  turn_on_time + 3.5, turn_on_time + 4., turn_on_time + 4.5, turn_on_time + 5])

# ENERGY RELEASES
'''
def f(t):

    global H_0

    if(t>=turn_on_time_0*yrtosec):
        if(t>=turn_on_time*yrtosec and t<=(turn_on_time*yrtosec+D)):
            return 1 + (H_0-1) * numpy.power(numpy.sin(numpy.pi*(t - turn_on_time*yrtosec)/D),2)
        else:
            return 1.
    else:
        return 0.0
'''

# ACCRETION
def f(t):

    if(t>=turn_on_time_0*yrtosec):
        if t>=turn_on_time*yrtosec:
            return  numpy.exp(-(t - turn_on_time*yrtosec)/(D * 24 * 60 * 60))
        else:
            return 1.0
    else:
        return 0.0
from data import loaddata
from control.manager import *
from control.constants import *
import numpy as np
from scipy import interpolate

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def init():                                                  

    global Rho, r_b, r, rho_r,dr,dr_b

    # --------------------------------------- mesh ----------------------------------------
    # -------------------------------------------------------------------------------------


    Rho = numpy.logspace( numpy.log10(loaddata.star_model()[2,3]),
                          numpy.log10(loaddata.star_model()[-2,3]), 100 * Nzones + 1)

    r_b = loaddata.radii(Rho)                                    # boundaries of the cells
    r = (r_b[1:] + r_b[:-1])/2                                   # average cell's radius
    dr = r[1:] - r[:-1]                                          # distance between cells
    dr_b = r_b[1:] - r_b[:-1]                                    # distance between cells' boundaries
    rho_r = loaddata.rho(r)                                      # cells' average density

    global T, t, dt, redshift, N

    T = T_0*numpy.ones(Nzones)                                   # initial temperature of neutron star
    t = t_0                                                      # initial time
    dt = dt_0                                                    # initial time step
    N = 100 * Nzones - 1                                         # We introduce N to prevent us from writing 'Nzones-1' a lot of times in code below


def relativity_sqrt(r):
    return numpy.sqrt(1 - 2*G*loaddata.mass(r)/(( c ** 2 )*r))

def compute_power(rho1,rho2):

    #cfg_file = open('data/config.dat', 'wb')

    power = 0
    counter = 0

    for i in range(len(r)-1,-1,-1):
        if(rho_r[i]<=rho2 and rho_r[i] >= rho1):
            power += 4.*numpy.pi*numpy.power(r[i],2)*dr_b[i]*numpy.exp(2*loaddata.Phi(r[i]))/relativity_sqrt(r[i])
            counter += 1

    rho1_array = numpy.power(10,numpy.linspace(10,14,9))
    periods      = numpy.power(10,numpy.linspace(-1,2,7))
    source_power = numpy.power(10,numpy.linspace(18,21,4))

    rho2_array = numpy.zeros(len(rho1_array))
    counter_array = numpy.zeros(len(rho1_array))
    power_array = numpy.zeros(len(rho1_array))

    for j in range(0,len(rho1_array)):

        for i in range(len(r)-1,-1,-1):

            if rho_r[i] >= rho1_array[j]:
                power_array[j] += 4.*numpy.pi*numpy.power(r[i],2)*dr_b[i]*numpy.exp(2*loaddata.Phi(r[i]))/relativity_sqrt(r[i])
                counter_array[j] += 1
            if (power_array[j]>power):
                rho2_array[j] = rho_r[i]
                break

    for source in range(len(source_power)):
        for rho in range(len(rho1_array)):
            for period in range(len(periods)):
                save = numpy.vstack([rho2_array[rho],rho1_array[rho],source_power[source],periods[period]])
                #numpy.savetxt(cfg_file, save.T, fmt='%1.6e')

    print(power_array)

#loaddata.star_model_data_init()
#init()
#compute_power(1e11,1e12)



def find_rho():

    model_file = 'data/BSK21_1.85.dat'
    #model_file = 'data/model_rho9.dat'
    _star_model = np.loadtxt(model_file, skiprows=2)

    _rho = interpolate.interp1d(np.log(_star_model[:, 1] * 1e5), np.log(_star_model[:, 3]), kind='linear')
    _radii = interpolate.interp1d(np.log(_star_model[:, 3]), np.log(_star_model[:, 1] * 1e5), kind='linear')
    _mass = interpolate.interp1d(np.log(_star_model[:, 1] * 1e5), np.log(_star_model[:, 0] * 1.98892e33), kind='linear')
    _Phi = interpolate.interp1d(np.log(_star_model[:, 1] * 1e5), _star_model[:, 4], kind='linear')

    def mass(a):  # mass(radius)
        return np.exp(_mass(np.log(a)))
    def rho(a):  # density(radius)
        return np.exp(_rho(np.log(a)))
    def radii(a):  # density(radius)
        return np.exp(_radii(np.log(a)))
    def relativity_sqrt(r):
        return np.sqrt(1 - 2*6.67259e-8 *mass(r)/(( 2.99792458e10  ** 2 )*r))
    def Phi(a):  # dimensionless gravitational potential(radius)
        return _Phi(np.log(a))

    rho2 = [2.512293e11, 1e12, 1.275997e+13, 5.063467e+13, 6.190351e12,       2.2099571e13]
    rho3 = [2.563587e11, 1e12, 1.230001e+13, 4.756411e+13, 6.01725143048e+12, 2.1110614303e+13]
    rho1 = [3.16227766e+10, 1e11, 1e12, 1e13, 3.16227766e+11, 3.16227766e+12]

    for i in range(0,6):

        rho_sample = np.logspace(np.log10(rho1[i]), np.log10(rho3[i]), 500)
        r_sample = radii(rho_sample)[::-1]
        V = 4*np.pi*np.trapz(np.power(r_sample,2)*np.exp(2*Phi(r_sample))/relativity_sqrt(r_sample), r_sample)
        print(V)

    V_ref = 1.62838159069e+17 # M = 1.85
    rho_test = np.logspace(np.log10(1.1e10),np.log10(3e14), 10000)
    r_test = radii(rho_test)

    for i in [0,1,2,3,4,5]:
        num = find_nearest(rho_test,rho1[i])
        print(num)
        V = np.zeros(len(rho_test))
        for j in range(num+5,len(rho_test)):
            r_sample = r_test[num:j][::-1]
            V[j] = 4*np.pi*np.trapz(np.power(r_sample,2)*np.exp(2*Phi(r_sample))/relativity_sqrt(r_sample), r_sample)
        num2 = find_nearest(V,V_ref)
        print(rho_test[num2])


find_rho()




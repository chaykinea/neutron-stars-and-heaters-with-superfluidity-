__author__ = 'maryhallow'


import numpy
import timeit
import time

from data              import loaddata
from physics           import tite_1e9
from physics           import tite
from computation       import PDEsolver
from control.manager   import *
from control.constants import *

import argparse


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
args = parser.parse_args()
external_param = numpy.array(args.integers)


def data_init():

    print ('Simulation of neutron star cooling starts... \n')

    loaddata.star_model_data_init()
    if(model_file):
        tite_1e9.init(log_rho_b,rho_star)
    else:
        tite.tite_init()
    loaddata.TiTe_data_init()
    loaddata.superfluid_data_init()
    loaddata.eff_mass_data_init()

    print ('All data files are loaded.\n')

def config_reader():

    global left_boundary, right_boundary, power, duration

    config = numpy.loadtxt(source_cfg)

    left_boundary  = config[:,0]
    right_boundary = config[:,1]
    power          = config[:,2]
    duration       = config[:,3]

def main():

    print ('Files for output data are created.\n')
    print ('Generating C,Q,k tables...\n')

    PDEsolver.TableCreator()
    PDEsolver.interpolation_on_tables()
    PDEsolver.init()
    PDEsolver.time_derivative_init()
    PDEsolver.time_derivative_iso_regime_init()

    if(source_trigger):

        print ('Heat source is turned on\n')

        config_reader()
        number_of_simulations = len(power)

        print ('Number of simulations = %i\n \n' % (number_of_simulations))
        print ('Sources\' properties:\n')
        print ('Number Left Right Power Duration(Left and Right source boundaries in density units)')
        print ('-------------------------')

        for i in range(external_param[0], external_param[1]):
            print ('%2i %1.2e %1.2e %1.2e %1.2e' % (i, left_boundary[i], right_boundary[i], power[i], duration[i]))
        print ('-------------------------\n')

        for i in range(10,0,-1):
            print ('Simulation begins in %d sec...' % i)
            time.sleep(1)
        print (' ')

        start = timeit.default_timer()

        for simulation_number in range(external_param[0],external_param[1]):

            if simulation_number != external_param[0]:
                PDEsolver.re_init()
            PDEsolver.source_initialization(period[simulation_number],power[simulation_number],
                                            left_boundary[simulation_number],right_boundary[simulation_number])
            if(source_visualise):
                PDEsolver.source_visualisation(simulation_number)
            PDEsolver.solve_PDE_with_source(simulation_number,number_of_simulations)

    else:

        for i in range(10,0,-1):
            print ('Simulation begins in %d sec...' % i)
            time.sleep(1)
        print ('')

        start = timeit.default_timer()
        PDEsolver.solve_PDE(output_1, output_2)

    stop = timeit.default_timer()

    print (' ')
    print ('calculation time: %3.5f sec' % (stop - start))


def generate_cfg():

    cfg_file = open('data/config.dat', 'wb')

    loaddata.star_model_data_init()

    rho_1 = 1e11
    rho_2 = 1e12

    V = lambda rho_1,rho_2: 4*numpy.pi/3 * ( numpy.power(loaddata.radii(rho_1),3) - numpy.power(loaddata.radii(rho_2),3))
    rho_f = lambda r: numpy.power(r,2)/numpy.sqrt(1 - (2*G*loaddata.mass(r)/c/c/r))
    V_prop = lambda rho_1,rho_2: 4*numpy.pi*integrate.quad(rho_f,loaddata.radii(rho_2),loaddata.radii(rho_1))[0]

    V_const = V(rho_1,rho_2) # reference volume (corresponds to rho1 = 1e12, rho2 = 1e13)

    rho_1_array  = numpy.array([1e10,5e10,1e11,5e11,1e12,5e12,1e13])
    V_array      = numpy.array([V_const/32,V_const/16,V_const/8,V_const/4,V_const/2,V_const,V_const*2])
    Power_array  = numpy.array([1e17,5e17,1e18,5e18,1e19,5e19,1e20])
    Period_array = numpy.array([1e-1,5e-1,1e0,5e0,1e1,5e1,1e2,5e2,1e3])

    for i in range(0,len(rho_1_array)): # Volume

        rho_2_array = loaddata.rho( numpy.power((4*numpy.pi/3*numpy.power(loaddata.radii(rho_1_array[i]),3) - V_array[-2])*3/numpy.pi/4,1/3) )

        for j in range(0,len(Power_array)): # Power

            for k in range(0,len(Period_array)):

                A = numpy.vstack([rho_2_array,rho_1_array[i],Power_array[j],Period_array[k]])
                numpy.savetxt(cfg_file, A.T, fmt='%1.6e')

        print(V(rho_1_array[i],rho_2_array))
        print(V_prop(rho_1_array[i],rho_2_array))
        print(rho_2_array)

#generate_cfg()
data_init()
main()

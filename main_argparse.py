__author__ = 'maryhallow'


import numpy
import timeit
import time
import numpy as np

from scipy             import integrate
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

    ns = np.array(['AWP2', 'AWP3', 'CCDK', 'CLS', 'GIPSF', 'MSH', 'SCLBL', 'SFB', 'WAP'])
    ps  = np.array(['AO', 'BCLL', 'BS', 'CCDK', 'CCYms', 'CCYps', 'EEHO', 'EEHOr', 'T'])
    nt = np.array(['AO', 'BEEHS', 'EEHO', 'EEHOr', 'SYHHP', 'T', 'TTav', 'TToa'])

    if external_param[2]!=0:
        params = numpy.hstack([external_param[2],
                               models['neutron_singlet'][ns[external_param[3]]],
                               models['neutron_triplet'][nt[external_param[4]]],
                               models['proton_singlet'][ps[external_param[5]]]])

        #print('Turn on time:', time_roots[external_param[3]],' yr\n')
        print('Superfluidity parameters for simulation:\n')
        print(ns[external_param[3]],np.round(params[1:6],2))
        print(nt[external_param[4]],np.round(params[6:11],2))
        print(ps[external_param[5]],np.round(params[11:],2))
        print( )

        PDEsolver.add_params(*params)
    else:
        PDEsolver.add_params(SF=0)

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
            print('Simulation begins in %d sec...' % i)
            time.sleep(1)
        print(' ')

        start = timeit.default_timer()

        if external_param[2]!=0:
            file_name = '_' + ns[external_param[3]] + '_'
            #file_name = '_' + ns[external_param[3]] + '_' + nt[external_param[4]] + '_' + ps[external_param[5]] + '_'
            #PDEsolver.params(time=time_roots[external_param[3]], source_name=file_name)
            PDEsolver.params(time=10000, source_name=file_name)
        else:
            file_name = '_SF0_'
            #PDEsolver.params(time=time_roots[external_param[3]], source_name=file_name)
            PDEsolver.params(time=10000, source_name=file_name)

        for simulation_number in range(external_param[0],external_param[1]):

            if simulation_number != external_param[0]:
                PDEsolver.re_init()
            PDEsolver.source_initialization(duration[simulation_number],power[simulation_number],
                                            left_boundary[simulation_number],right_boundary[simulation_number])
            if(source_visualise):
                PDEsolver.source_visualisation(simulation_number)
            PDEsolver.solve_PDE_with_source(simulation_number,external_param[1])

    else:

        #output_1, output_2 = loaddata.create_output_files()

        if external_param[2]!=0:
            output_1 = open('output/cooling_'
                            + ns[external_param[3]] + '_'
                            + nt[external_param[4]] + '_'
                            + ps[external_param[5]] + '.dat', 'wb')
            output_2 = open('output/profiles_'
                            + ns[external_param[3]] + '_'
                            + nt[external_param[4]] + '_'
                            + ps[external_param[5]] + '.dat', 'wb')
        else:
             output_1 = open('output/cooling_SF0.dat', 'wb')
             output_2 = open('output/profiles_SF0.dat', 'wb')

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

    V = lambda rho_1, rho_2: 4*numpy.pi/3 * ( numpy.power(loaddata.radii(rho_1),3) - numpy.power(loaddata.radii(rho_2),3))
    rho_f = lambda r: numpy.power(r,2)/numpy.sqrt(1 - (2*G*loaddata.mass(r)/c/c/r))
    V_prop = lambda rho_1, rho_2: 4*numpy.pi*integrate.quad(rho_f,loaddata.radii(rho_2),loaddata.radii(rho_1))[0]

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


models =  {'neutron_singlet': {'AWP2' : [  28, 0.20,   1.5,  1.7,   2.5],
                               'AWP3' : [  50, 0.20,   2.0,  1.4,   2.0],
                               'CCDK' : [ 127, 0.18,   4.5, 1.08,   1.1],
                               'CLS'  : [ 2.2, 0.18,  0.06,  1.3,  0.03],
                               'GIPSF': [ 8.8, 0.18,   0.1,  1.2,   0.6],
                               'MSH'  : [2.45, 0.18,  0.05,  1.4,   0.1],
                               'SCLBL': [ 4.1, 0.35,   1.7, 1.67,  0.06],
                               'SFB'  : [  45, 0.10,   4.5, 1.55,   2.5],
                               'WAP'  : [  69, 0.15,   3.0,  1.4,   3.0],
                               'NONE' : [ 0.0, 0.15,   3.0,  1.4,   3.0]},
           'proton_singlet' : {'AO'   : [  14, 0.15,  0.22, 1.05,   3.8],
                               'BCLL' : [1.69, 0.05,  0.07, 1.05,  0.16],
                               'BS'   : [  17,  0.0,   2.9,  0.8,  0.08],
                               'CCDK' : [ 102,  0.0,   9.0,  1.3,   1.5],
                               'CCYms': [  35,  0.0,   5.0,  1.1,   0.5],
                               'CCYps': [  34,  0.0,   5.0, 0.95,   0.3],
                               'EEHO' : [ 4.5,  0.0,  0.57,  1.2,  0.35],
                               'EEHOr': [  61,  0.0,   6.0,  1.1,   0.6],
                               'T'    : [  48, 0.15,   2.1,  1.2,   2.8],
                               'NONE' : [ 0.0, 0.15,   2.1,  1.2,   2.8]},
           'neutron_triplet':  {'AO'  : [ 4.0,  1.2,  0.45,  3.3,   5.0],
                               'BEEHS': [0.45,  1.0,  0.40,  3.2,  0.25],
                               'EEHO' : [0.48, 1.28,   0.1, 2.37,  0.02],
                               'EEHOr': [0.23,  1.2, 0.026,  1.6, 0.008],
                               'SYHHP': [ 1.0, 2.08,  0.04,  2.7, 0.013],
                               'T'    : [ 1.2, 1.55,  0.05, 2.35,  0.07],
                               'TTav' : [ 3.0,  1.1,  0.60, 2.92,   3.0],
                               'TToa' : [ 2.1,  1.1,  0.60,  3.2,   2.4],
                               'NONE' : [ 0.0,  1.1,  0.60,  3.2,   2.4]}}

time
#generate_cfg()
data_init()
main()


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
    config = numpy.loadtxt(source_cfg)
    loaddata.star_model_data_init(file='data/BSK21_1.' + str(external_param[6]) + '.dat')
    print('File name = ', 'data/BSK21_1.' + str(external_param[6]) + '.dat', '\n')
    print('File name = ', 'data/BSK21_1.' + str(external_param[6]) + '.dat', '\n')
    print('File name = ', 'data/BSK21_1.' + str(external_param[6]) + '.dat', '\n')
    loaddata.superfluid_data_init()
    loaddata.eff_mass_data_init()

    print ('All data files are loaded.\n')

def config_reader():

    global left_boundary, right_boundary, power, duration, env_model, rho_star_arr

    config = numpy.loadtxt(source_cfg)

    left_boundary  = config[:,0]
    right_boundary = config[:,1]
    power          = config[:,2]
    duration       = config[:,3]
    env_model      = config[:,4]
    rho_star_arr   = config[:,5]

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

        start = timeit.default_timer()

        if external_param[2]!=0:
            file_name = '_' + ns[external_param[3]] + '_'
            PDEsolver.params(time=10000, source_name=file_name)
        else:
            file_name = '_SF0_'
            PDEsolver.params(time=10000, source_name=file_name)

        for simulation_number in range(external_param[0],external_param[1]):
            
            if (tite_model):
                print('new')
                tite_1e9.init(rho_bound=log_rho_b, rho_star=rho_star_arr[simulation_number], comp=int(env_model[simulation_number]))
            else:
                print('old')
                tite.tite_init()
            

            print('Comp = ', int(env_model[simulation_number]), '\n')
            print('Comp = ', int(env_model[simulation_number]), '\n')
            print('Comp = ', int(env_model[simulation_number]), '\n')
            print('Rho star = ', rho_star_arr[simulation_number], '\n')
            print('Rho star = ', rho_star_arr[simulation_number], '\n')
            print('Rho star = ', rho_star_arr[simulation_number], '\n')

            loaddata.TiTe_data_init()

            if simulation_number != external_param[0]:
                PDEsolver.re_init()

            a = 1.00 + float(external_param[6])/100
            print('model param:', a, '\n')
            print('model param:', a, '\n')
            print('model param:', a, '\n')

            PDEsolver.source_power(model_parm=a)
            PDEsolver.source_initialization(duration[simulation_number],power[simulation_number],
                                            left_boundary[simulation_number],right_boundary[simulation_number])
            if(source_visualise):
                PDEsolver.source_visualisation(simulation_number)
            PDEsolver.solve_PDE_with_source(simulation_number, external_param[1])

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

        print ('')

        start = timeit.default_timer()
        PDEsolver.solve_PDE(output_1, output_2)

    stop = timeit.default_timer()

    print (' ')
    print ('calculation time: %3.5f sec' % (stop - start))


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

data_init()
main()



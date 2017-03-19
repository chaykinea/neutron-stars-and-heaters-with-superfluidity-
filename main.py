#     Neutron star thermal evolution
# --------------------------------------
#              main.py
# --------------------------------------
# This is the main module of
# the code. Use it to start
# the simulation.

import matplotlib.pyplot as plt
import timeit
import time

from data            import loaddata
from physics         import tite_1e9
from computation     import PDEsolver
from other           import plot
from control.manager import *
from physics         import tite


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

    output_1, output_2 = loaddata.create_output_files()

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

        for i in range(0,number_of_simulations):
            print ('%2i %1.2e %1.2e %1.2e %1.2e' % (i, left_boundary[i], right_boundary[i], power[i], duration[i]))
        print ('-------------------------\n')

        for i in range(10,0,-1):
            print ('Simulation begins in %d sec...' % i)
            time.sleep(1)
        print (' ')

        start = timeit.default_timer()

        for simulation_number in range(0,number_of_simulations):

            if simulation_number != 0:
                PDEsolver.re_init()
            PDEsolver.source_initialization(duration[simulation_number],power[simulation_number],
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


def show_cooling_curves():

    temp_data_1 = numpy.loadtxt(output_file_1)
    temp_data_2 = numpy.loadtxt(output_file_2)

    plt.title('Cooling curves')
    labels = []
    for i in range (0,len(time_points)):
        labels.append(str(time_points[i]))

    line_styles = ['b-','y-','g-','r-','m-','c-','b--','y--','g--', 'r--', 'm--', 'c--']
    X = numpy.outer(temp_data_2[:,0],numpy.ones(len(time_points)))
    Y = temp_data_2[:,1:]

    plt.plot(numpy.log10(temp_data_1[:,1]) ,numpy.log10(temp_data_1[:,0]),'b--')
    plt.show()

    plot.multi_plot(numpy.log10(X), numpy.log10(Y),'$\mathrm{Temperature profile}$', labels,  1,
               '$\mathrm{log}_{10} \\thinspace \\rho \\thinspace \mathrm{(g} \\thinspace \mathrm{cm^{-3})}$', '$\mathrm{log}_{10} \\thinspace \mathrm{T}^{\infty} \mathrm{(K)}$', line_styles[0:len(time_points)], xslace = 'linear', yslace = 'linear')

#config_reader()
#config_reader_periodic()

data_init()
main()
#show_cooling_curves()




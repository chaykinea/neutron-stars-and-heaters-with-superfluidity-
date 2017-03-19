#     Neutron star thermal evolution
# --------------------------------------
#            physics.py
# --------------------------------------
# This module contains microphysics
# required to create functions
# C(T,rho), Q(T,rho), kappa(T,rho)
# for numerical computation,
# where T - local temperature
# rho - density
# C - heat capacity
# Q - redshifted neutrino emissivity
# kappa - thermal conductivity
# (all per unit volume) for
# specific neutron star model

from physics           import heatcapacity
from physics           import neutrino
from physics           import thermalconductivity
from data              import loaddata
from control.constants import *
from control.manager   import *

import scipy.interpolate
import numpy

def add_params(SF=0,AT=0,A0=0,A1=0,A2=0,A3=0,BT=0,B0=0,B1=0,B2=0,B3=0,CT=0,C0=0,C1=0,C2=0,C3=0):

    global coeffsnk0, coeffsnk1, coeffsnk2, coeffsnk3, coeffsnT0, \
           coefftnk0, coefftnk1, coefftnk2, coefftnk3, coefftnT0, \
           coeffspk0, coeffspk1, coeffspk2, coeffspk3, coeffspT0, \
           SUPERFLUIDITY

    SUPERFLUIDITY = int(SF)

    coeffsnT0 = AT * 0.5669 * MeV_erg / kB
    coeffsnk0 = A0
    coeffsnk1 = numpy.sqrt(A1)
    coeffsnk2 = A2
    coeffsnk3 = numpy.sqrt(A3)

    coefftnT0 = 0. * BT * 0.1187 * MeV_erg / kB
    coefftnk0 = B0
    coefftnk1 = numpy.sqrt(B1)
    coefftnk2 = B2
    coefftnk3 = numpy.sqrt(B3)

    coeffspT0 = 0. * CT * 0.5669 * MeV_erg / kB
    coeffspk0 = C0
    coeffspk1 = numpy.sqrt(C1)
    coeffspk2 = C2
    coeffspk3 = numpy.sqrt(C3)

# critical temperatures in K of neutron triplet and singlet superfluidity 
def _Tcn(nn):

    kf=numpy.power(3.*pi*pi*nn, 1./3.)
    if ( kf<= coefftnk0 or kf>= coefftnk2 ):
        Tc=0.0
    else:
        a=(kf-coefftnk0)*(kf-coefftnk0)
        b=(coefftnk2-kf)*(coefftnk2-kf)
        Tc=coefftnT0*a/(a+coefftnk1*coefftnk1)*b/(b+coefftnk3*coefftnk3)

    Tc = Tc+1.0e-6

    return Tc


# critical temperature in K of proton singlet superfluidity in the core 
def _Tcp(rho, np):

    if (rho <= CoreCrustBound):
        Tc=0.0
    else:
        if (coeffspk2 < 0.0):
            Tc=coeffspT0
        else:
            kf=numpy.power(3.*pi*pi*np, 1./3.)
            if ( kf<= coeffspk0 or kf>= coeffspk2 ):
                Tc=0.0
            else:
                a=(kf-coeffspk0)*(kf-coeffspk0)
                b=(coeffspk2-kf)*(coeffspk2-kf)
                Tc=coeffspT0*a/(a+coeffspk1*coeffspk1)*b/(b+coeffspk3*coeffspk3)

    Tc = Tc+1.0e-6

    return Tc


def _Tcs(nn,rho):

    if (coeffsnk2 < 0.0):
        Tcs = coeffsnT0
    else:
        kf = numpy.power(3.*pi*pi*nn, 1./3.)
        if (kf <= coeffsnk0 or kf >= coeffsnk2 ):
            Tcs = 0.0
        else:
            if(rho>neutron_drip):
                a = (kf-coeffsnk0)*(kf-coeffsnk0)
                b = (coeffsnk2-kf)*(coeffsnk2-kf)
                Tcs = coeffsnT0*a/(a+coeffsnk1*coeffsnk1)*b/(b+coeffsnk3*coeffsnk3)
            else:
                Tcs = 0.0


    Tcs = Tcs + 1.0e-6

    return Tcs


# particle number densities in the crust using analytical fitting expressions by DG Yakovlev:
# Wigner-Zeitz cell model in the inner NS crust
def _CrustData(rho, nb, nn, ne):

    # rho = density in g/cc
    # nb = number density of baryons in fm^-3
    # nn = number density of free neutrons in fm^-3
    # ne = number density of electrons in fm^-3
    # Z = number of protons in a nucleus
    # A = number of nucleons per nucleus
    # Anuc = number of nucleons in a nucleus, <= A

    if (rho > CoreCrustBound):                                   # this function does not handle the NS core
                                                                 #  => do not change particle densities
        Z, A, Anuc, Vion = 0., 0., 0., 0.
        nb_, nn_, ne_ = nb, nn, ne
        return nb_, nn_, ne_, Z, A, Anuc, Vion

    n_b=rho/m_u/1.e39                                            # baryon number density in fm^-3

    if (rho<=rho_nd):                                            # densities lower than the neutron drip

        f = numpy.log(1.+n_b/5.e-9)
        Rp = 5.688 + 0.02628*f + 0.009468*f*f                    # proton core radius
        Rn = 5.788 + 0.02077*f + 0.01489*f*f                     # neutron core radius
        np_in = 0.0738 + 1.22e-4*f - 1.641e-4*f*f
        nn_in = 0.0808 + 1.688e-4*f + 9.439e-5*f*f
        nn_out = 0.

        tn = 6.
        tp = 6.
        Nin = 4.*pi*numpy.power(Rn,3.)*nn_in*(1./3.-3./(3.+tn)+3./(3.+2.*tn)-1./(3.+3.*tn))
        Z = 4.*pi*numpy.power(Rp,3.)*np_in*(1./3.-3./(3.+tp)+3./(3.+2.*tp)-1./(3.+3.*tp))
        Anuc = Z + Nin
        A = Anuc
        Rws = numpy.power(A*3./4./pi/n_b,1./3.)                  # Wigner-Seitz radius in fm

    else:                                                        # spheres after drip

        g = n_b*100.
        f = numpy.log(g)
        Rws = 31.68-8.400*f-0.2380*f*f+0.1152*f*f*f              # Wigner-Seitz radius in fm
        tn = 1./(0.2027+0.004506*g)
        Rn = 9.406+1.481*f+0.4625*f*f+0.05738*f*f*f
        dn_n = (9.761-1.322*f-0.5544*f*f-0.07624*f*f*f)/100.
        Nin = 4.*pi*numpy.power(Rn,3.)*dn_n*(1./3.-3./(3.+tn)+3./(3.+2.*tn)-1./(3.+3.*tn))
        tp = 1./(0.1558+2.225e-3*g+9.452e-4*g*g)
        Rp = 8.345+0.7767*f+0.1333*f*f+0.008707*f*f*f
        np_in = (4.040-1.097*f-0.0723*f*f+0.0225*f*f*f)/100.
        Z = 4.*pi*numpy.power(Rp,3)*np_in*(1./3.-3./(3.+tp)+3./(3.+2.*tp)-1./(3.+3.*tp))

        Nfree = n_b*4.*pi/3.*numpy.power(Rws,3.)-Z-Nin           # number of free neutrons per Wigner-Seitz cell
        nn_out = Nfree/(4.*pi*Rws*Rws*Rws/3.)                    # number density of free neutrons in fm^-3
        nn_in = nn_out+dn_n
        np_out = 0.0

        A = Z+Nin+Nfree
        Anuc = Z+Nin+Nfree*numpy.power(Rn/Rws,3.)

    if (rho>1.0e13):                                             # to achieve the smooth core-crust transition July 2001
        nn_out = n_b

    nn_ = nn_out
    ne_ = n_b*Z/A
    nb_ = n_b
    v = numpy.power(Rn/Rws,3.)

    if(v<0.):
        v = 0.
    if(v>1.):
        v = 1.

    Vion = v

    return nb_, nn_, ne_, Z, A, Anuc, Vion

def TableCreator():

    global log_rho, log_T, C_table, Q_table, kappa_table

    if (tables_compute):

        T = numpy.logspace(numpy.log10(T_table_min),numpy.log10(T_table_max),T_table_N)
        log_T = numpy.log(T)

        rho = loaddata.star_model()[:, 3]
        log_rho = numpy.log(rho)
        Phi = loaddata.star_model()[:, 4]

        nb = loaddata.star_model()[:, 5]
        ne = loaddata.star_model()[:, 6]
        nm = loaddata.star_model()[:, 7]

        np=ne+nm
        nn=nb-np

        C_table = numpy.zeros((len(rho), len(T)))
        Q_table = numpy.zeros((len(rho),len(T)))
        kappa_table = numpy.zeros((len(rho),len(T)))

        temp_var = 0.

        for i in range (0,len(rho)):

            nb_, nn_, ne_, Z, A, Anuc, Vion = _CrustData(rho[i], nb[i], nn[i],ne[i])
            Tc_nt=_Tcn(nn[i])
            Tc_p=_Tcp(rho[i],np[i])
            Tc_ns = _Tcs(nn[i],rho[i])

            for j in range (0,T_table_N):

                T_local = T[j]*numpy.exp(-Phi[i])

                C_table[i,j] = heatcapacity._C(T_local, rho[i], nn_, ne_, nm[i], Z, Anuc, Tc_nt, Tc_ns, Tc_ns)
                Q_table[i,j] = neutrino._Q(T_local, nn_, ne_, nm[i], rho[i], Phi[i], Tc_nt, Tc_ns, Tc_p)
                kappa_table[i,j] = thermalconductivity._kappa(T_local, rho[i], ne_, nm[i], nn_, Z, Tc_p)

            if (int(i/(float(len(rho))/100))) > temp_var and numpy.mod((int(i/(float(len(rho))/100))),10)==0:

                temp_var = (int(i/(float(len(rho))/100)))
                print ('Creating tables...  %3d %%' % temp_var)

        print ('Creating tables...  %3d %%\n' % 100)
        print ('Tables are created.\n')

        if (tables_write):

            numpy.savetxt('file1.dat', T, fmt='%1.6e')
            numpy.savetxt('file2.dat', rho, fmt='%1.6e')
            numpy.savetxt('file3.dat', Q_table, fmt='%1.6e')
            numpy.savetxt('file4.dat', C_table, fmt='%1.6e')
            numpy.savetxt('file5.dat', kappa_table, fmt='%1.6e')

            print ('Tables have been written out.\n')
        else:
            print ('Tables will not be written out.\n')

    else:

        T = numpy.loadtxt('file1.dat')
        log_T = numpy.log(T)

        rho = numpy.loadtxt('file2.dat')
        log_rho = numpy.log(rho)

        Q_table = numpy.loadtxt('file3.dat')
        C_table = numpy.loadtxt('file4.dat')
        kappa_table = numpy.loadtxt('file5.dat')

        print ('Tables are loaded.\n')


def interpolation_on_tables():                                   # We create C(T,rho), Q(T,rho), kappa(T,rho) functions
                                                                 # using tables which are created in the function above
    print ('Interpolation has started.\n')
    x, y = numpy.meshgrid(log_T, log_rho)
    p = numpy.vstack([x.flatten(), y.flatten()]).T

    global log_C, log_Q, log_kappa

    log_C = scipy.interpolate.LinearNDInterpolator(p, numpy.log(C_table).flatten())
    log_Q = scipy.interpolate.LinearNDInterpolator(p, numpy.log(Q_table).flatten())
    log_kappa = scipy.interpolate.LinearNDInterpolator(p, numpy.log(kappa_table).flatten())

    print ('Interpolation completed.\n')

def C(a, b):                                                     # heat capacity(temperature, density)
    return numpy.exp(log_C(numpy.log(a), numpy.log(b)))


def k(a, b):                                                     # thermal conductivity(temperature, density)
    return numpy.exp(log_kappa(numpy.log(a), numpy.log(b)))


def Q(a, b):                                                     # neutrino emissivity(temperature, density)
    return numpy.exp(log_Q(numpy.log(a), numpy.log(b)))

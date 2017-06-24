#     Neutron star thermal evolution
# --------------------------------------
#            loadtata.py
# --------------------------------------
# This module provides the whole data
# that has to be read from input files
# before simulation starts.
# Furthermore, it contains information
# about star model functions that built
# via interpolation procedure


from control.constants import *
from control.manager import *

from sys import exit
from physics import sf_gap
from scipy import interpolate

# --This function creates output files, in which cooling data will be written-
def create_output_files():
    
    output_1 = open(output_file_1, 'wb')  # Cooling data will be written to this file
    output_2 = open(output_file_2, 'wb')  # Temperature profile data will be written to this file

    return output_1, output_2



# -----------This blog contains data for specific neutron star model----------
def star_model_data_init(file = model_file):

    global _g_surface, _star_model, _rho, _radii, _mass, _Phi, _nb, _ne, _nm, _nm_smooth, _pressure

    try:
        temp = open(file, 'r')
    except(IOError):
        print('Star_model file is not found. \n')
        print('Simulation is terminated.')
        exit(0)

    _g_surface = float(temp.readline().split()[1])
    temp.close()

    _star_model = numpy.loadtxt(file, skiprows=1)

    _rho = interpolate.interp1d(numpy.log(_star_model[:, 1] * 1e5), numpy.log(_star_model[:, 3]), kind='linear')
    _pressure = interpolate.interp1d(numpy.log(_star_model[:, 1] * 1e5), numpy.log(_star_model[:, 2]), kind='linear')
    _radii = interpolate.interp1d(numpy.log(_star_model[:, 3]), numpy.log(_star_model[:, 1] * 1e5), kind='linear',fill_value='extrapolate')
    _mass = interpolate.interp1d(numpy.log(_star_model[:, 1] * 1e5), numpy.log(_star_model[:, 0] * MSun), kind='linear')
    _Phi = interpolate.interp1d(numpy.log(_star_model[:, 1] * 1e5), _star_model[:, 4], kind='linear')
    _nb = interpolate.interp1d(numpy.log(_star_model[:, 3]), _star_model[:, 5], kind='linear')
    _ne = interpolate.interp1d(numpy.log(_star_model[:, 3]), _star_model[:, 6], kind='linear')
    _nm = interpolate.interp1d(numpy.log(_star_model[:, 3]), _star_model[:, 7], kind='linear')

    rho_eff = numpy.concatenate((_star_model[:158, 3],_star_model[212:,3]))
    mu_eff = numpy.concatenate((_star_model[:158, 7],_star_model[212:,7]))
    _nm_smooth = interpolate.interp1d(numpy.log(rho_eff),mu_eff, kind='linear')



def g_surface():  # surface gravity in 10^14 cm/s^2
    return _g_surface

def star_model():
    return _star_model

def Phi(a):  # dimensionless gravitational potential(radius)
    return _Phi(numpy.log(a))

def rho(a):  # density(radius)
    return numpy.exp(_rho(numpy.log(a)))

def pressure(a):  # Pressure(radius)
    return numpy.exp(_pressure(numpy.log(a)))

def radii(a):  # radius(density)
    return numpy.exp(_radii(numpy.log(a)))

def mass(a):  # mass(radius)
    return numpy.exp(_mass(numpy.log(a)))

def nb(a):
    return _nb(a)

def ne(a):
    return _ne(a)

def nm(a):
    return _nm_smooth(a)


# ---------------reading data that provides Te(Ti) dependence-----------------
def TiTe_data_init():

    global _T_e

    try:
        TiTe = numpy.power(10 ,numpy.loadtxt(tite_file,skiprows=3))
    except(IOError):
        print('TiTe file is not found. \n')
        print('Simulation is terminated.')
        exit(0)

    _T_e = interpolate.interp1d(numpy.log(TiTe[:, 0]), numpy.log(TiTe[:, 1]), kind='linear', fill_value='extrapolate')


def T_e(a):
    return numpy.exp(_T_e(numpy.log(a)))



# ---------------reading double superfluidity data----------------------------
def superfluid_data_init():

    global _lgR
    Nsf = 35                                                     # number of lines and rows of both n & p
                                                                 # superfluidity data
    try:
        npsf = numpy.loadtxt(npsf_file)
    except(IOError):
        print('Superfluidity data file is not found. \n')
        print('Simulation is terminated.')
        exit(0)

    tauP = numpy.power(npsf[0,:],10)
    tauN = numpy.power(npsf[1,:],10)
    Ngap_, Pgap_ = [None]*Nsf, [None]*Nsf

    for i in range(0,35):
        Ngap_[i] = sf_gap._TripletGap(tauN[i])
        Pgap_[i] = sf_gap._SingletGap(tauP[i])

    lgRfactor = npsf[2:,:]
    x, y = numpy.meshgrid(Ngap_, Pgap_)
    p = numpy.vstack([x.flatten(), y.flatten()]).T
    _lgR = interpolate.LinearNDInterpolator(p, lgRfactor.flatten())


def lgR(a,b):
    return _lgR(a,b)


# --------------------reading effective mass data-----------------------------
def eff_mass_data_init():

    global _kf

    try:
        effmass = numpy.loadtxt(effmass_file, skiprows=9)
    except(IOError):
        print('Effmass file is not found. \n')
        print('Simulation is terminated.')
        exit(0)

    kf_m = effmass[:,0]  # Fermi momentum /h
    meff = effmass[:,1]  # effective nucleon mass
    _kf = interpolate.interp1d(kf_m, meff, kind='linear', fill_value='extrapolate')  # Fermi momentum in fm^-1


def kf(a):
    return _kf(a)


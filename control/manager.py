#     Neutron star thermal evolution
# --------------------------------------
#            manager.py
# --------------------------------------
# Using this module one can control
# the whole simulation process

# ----------------------------------------------------------------
# ----------------------general parameters------------------------
# ----------------------------------------------------------------
model_file     = 'data/BSK21_1.40.dat'                         # file containing data of the model of NS
model_param    =  1.40
tite_file      = 'data/tite.dat'                               # TiTe data file
npsf_file      = 'data/npsf.dat'                               # superfluidity data file
effmass_file   = 'data/effmass.dat'                            # effective mass data file

output_file_1  = 'data/output_1.dat'                           # file containing output data (cooling curves [T_e (redshifted), time, T_i])
output_file_2  = 'data/output_2.dat'                           # file containing output data (Temperature profile [rho,T] at time steps
                                                               # which are defined below)
regime         = 1                                             # if 1 - Implicit Euler scheme is used for calculation,
                                                               # 0 - Implicit Euler scheme while t< 7.e2 years
                                                               # and after isothermal scheme is used.

dM_accreted    = 1.e-6                                         # mass of accreted envelope in units of NS mass
MagField       = 0.0e12                                        # Magnetic filed at pole in Gauss
log_rho_b      = 9.0                                           # to determine the boundary for computing NS envelope log(rho_b) in g/cm3
rho_star       = 10**8                                         # [0<=rho_star<=10**log_rho_b] if 0 then Fe envelope is used, if max then C, and mixture in between
tite_model     = 0                                             # 0 for A. Y. Potekhin et al 1997 and 1 for M. V. Beznogov et al 2016 [C-Fe mixture]

T_0            = 0.569e8                                       # 1e8 for general case     # initial redshifted temperature of the star (T(r) = const at time t = 0) in K
# ----------------------------------------------------------------
# track 83 <--> 0.569e8
# track 84.3 <--> 0.5865e8

# ----------------------------------------------------------------
# -----------parameters for creating mesh and solving PDE---------
# ----------------------------------------------------------------
t_0            = 1e3*31536000                                  # initial time (in sec)
dt_0           = 0.001                                         # initial time step (in sec)

rho_vary       = 0                                             # if 0 - rho_min and rho_max are the max and min density values
                                                               # which are avaible in the NS model file. if 1 - you can
                                                               # set them on your own.
rho_min        = 1.0044775515833761e+010                       # max density value in the mesh for solving PDE (in gramms)

t_max          = 4.e6                                          # time when simlation stops (in years)
t_iso          = 7.e2                                          # for time t > t_iso NS has iswe assume, that NS has isothermal profile (in years)
T_min          = 2.e5                                          # Temperature below which simulation stops

Nzones         = 500                                           # number of zones star will be logarithmically divided into
N_output       = 2000                                          # number of data points in output file containing cooling curves

# ----------------------------------------------------------------
time_points    = [1e-7, 1.e-6, 1.e-3, 1.e-2, 1.e-1, 1.e0, 1.e1, 1.e2, 3.e2, 1.e3, 2.e3, 1.e5, 1.e6] # points in time where we change timestep and write T(rho) data into file
#time_steps     = [10,  1.e2, 4.e2, 1.e4, 3.e4, 1.e5, 3.e5, 1.e6, 5.e6, 1.e8,  1.e8, 1.e8,1.e10]  # timestep values
#time_steps     = [10,  1.e2, 5.e3, 5.e3, 3.e4, 1.e5, 1.e6, 1.e7, 1.e8,  1.e8, 1.e8,1.e10]  # timestep values
time_steps     = [10,    1.e2, 1.e3,   2.e3,  1.e4, 5.e4, 1.e5, 1.e5, 1.e5, 1.e5, 1.e8,1.e10]  # timestep values
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# -------parameters required to compute C,Q and kappa-------------
# ----------------------------------------------------------------
CoreCrustBound = 1.5e14    # Core-Crust bound                  : to specify value
neutron_drip   = 4.e11     # neutron drip density              : to specify value
SUPERFLUIDITY  = 1         # Superfuidity                      : 1-include any superfluid effects, 0-no SF at all
var_meff       = 2         # effective mass variable           : 0-fixed effective masses, 1-use variable meff, 2-fitted masses
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# ----------parameters which construct C,Q,kappa tables-----------
# ----------------------------------------------------------------
tables_compute = 1         # compute tables before the simulation starts or not (in case when you already have tables written out)
tables_write   = 0         # write C, Q, kappa tables to file  : 1 - write, 0 - do not write
T_table_N      = 150       # max number of Temperature points  :
T_table_max    = 1.e11     # for C(T,rho),Q(T,rho),kappa(T,rho):
T_table_min    = 1.e5      # tables. Max and min Temperature   :
                           # values for the table              :
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# ---------------------Source parameters--------------------------
# ----------------------------------------------------------------

import numpy

source_cfg          = 'data/config.dat'
source_trigger      = 1
source_visualise    = 0
error_min           = 2




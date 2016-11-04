#     Neutron star thermal evolution
# --------------------------------------
#            manager.py
# --------------------------------------
# Using this module one can control
# the whole simulation process

# ----------------------------------------------------------------
# ----------------------general parameters------------------------
# ----------------------------------------------------------------
model_file     = 'data/model_rho9.dat'                         # file containing data of the model of NS
tite_file      = 'data/tite2.dat'                              # TiTe data file
npsf_file      = 'data/npsf.dat'                               # superfluidity data file
effmass_file   = 'data/effmass.dat'                            # effective mass data file

output_file_1  = 'data/output_1.dat'                           # file containing output data (cooling curves [T_e (redshifted), time, T_i])
output_file_2  = 'data/output_2.dat'                           # file containing output data (Temperature profile [rho,T] at time steps
                                                               # which are defined below)
regime         = 1                                             # if 1 - Implicit Euler scheme is used for calculation,
                                                               # 0 - Implicit Euler scheme while t< 7.e2 years
                                                               # and after isothermal scheme is used.

dM_accreted    = 1.e-15                                        # mass of accreted envelope in units of NS mass
MagField       = 0.0e12                                        # Magnetic filed at pole in Gauss
log_rho_b      = 9.0                                           # to determine the boundary for computing NS envelope log(rho_b) in g/cm3
rho_star       = 0.0                                           # [0<=rho_star<=10**log_rho_b] if 0 then Fe envelope is used, if max then C, and mixture in between
tite_model     = 0                                             # 0 for A. Y. Potekhin et al 1997 and 1 for M. V. Beznogov et al 2016 [C-Fe mixture]

T_0            = 1.e10                                         # initial redshifted temperature of the star (T(r) = const at time t = 0) in K
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# -----------parameters for creating mesh and solving PDE---------
# ----------------------------------------------------------------
t_0            = 0.                                            # initial time (in sec)
dt_0           = 1                                             # initial time step (in sec)

rho_vary       = 0                                             # if 0 - rho_min and rho_max are the max and min density values
                                                               # which are avaible in the NS model file. if 1 - you can
                                                               # set them on your own.
rho_min        = 1.0044775515833761e+010                       # max density value in the mesh for solving PDE (in gramms)
rho_max        = 7.2926887254333875e+014                       # min density value in the mesh for solving PDE (in gramms)

t_max          = 4.e6                                          # time when simlation stops (in years)
t_iso          = 7.e2                                          # for time t > t_iso NS has iswe assume, that NS has isothermal profile (in years)
T_min          = 2.e5                                          # Temperature below which simulation stops

Nzones         = 350                                           # number of zones star will be logarithmically divided into
N_output       = 500                                           # number of data points in output file containing cooling curves

# ----------------------------------------------------------------
time_points    = [1.e-6, 1.e-3, 1.e-2, 1.e-1, 1.e0, 1.e1, 3.e2, 1.e3, 2.e3, 1.e5, 1.e6] # points in time where we change timestep and write T(rho) data into file
time_steps     = [1.e2, 5.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e10,  1.e10, 1.e10,1.e10]  # timestep values
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
# ----------------parameters for superfluidity--------------------
# ----------------------------------------------------------------
coeffsnT0      = 10.2e+9                                       # K singlet neutron critical temperature
coeffsnk0      = 0.0                                           # 1/fm  Tc=T0*(k-k0)^2/((k-k0)^2+k1^2)*
coeffsnk1      = 0.6                                           # 1/fm     *(k2-k)^2/((k2-k)^2+k3^2)
coeffsnk2      = 1.45                                          # 1/fm     for k0<=k<=k2
coeffsnk3      = 0.1                                           # 1/fm  or Tc=T0=const for negative k2

coefftnT0      = 6.461e-9                                      # K triplet neutron critical temperature
coefftnk0      = 1.0                                           # 1/fm  Tc=T0*(k-k0)^2/((k-k0)^2+k1^2)*
coefftnk1      = 1.961                                         # 1/fm     *(k2-k)^2/((k2-k)^2+k3^2)
coefftnk2      = 2.755                                         # 1/fm     for k0<=k<=k2
coefftnk3      = 1.3                                           # 1/fm  or Tc=T0=const for negative k2

coeffspT0      = 20.29e+9                                      # K singlet proton critical temperature
coeffspk0      = 0.0                                           # 1/fm  Tc=T0*(k-k0)^2/((k-k0)^2+k1^2)*
coeffspk1      = 1.117                                         # 1/fm     *(k2-k)^2/((k2-k)^2+k3^2)
coeffspk2      = 1.241                                         # 1/fm     for k0<=k<=k2
coeffspk3      = 0.1473                                        # 1/fm  or Tc=T0=const for negative k2
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# ---------------------Source parameters--------------------------
# ----------------------------------------------------------------

import numpy

source_cfg          = 'data/config.dat'
source_trigger      = 0
source_visualise    = 1
turn_on_time        = 4.0e4
t_source_max        = turn_on_time + 10
error_min           = 1

# ----------------------------------------------------------------
t_source_points     = numpy.array([turn_on_time +0.01 ,turn_on_time + 0.1,turn_on_time + 1.,turn_on_time + 5., turn_on_time + 10.,
                                   turn_on_time + 50., turn_on_time + 100., turn_on_time + 500., turn_on_time + 1000.])       # in years
t_source_steps      = [time_steps[error_min], time_steps[error_min]*10  , time_steps[error_min]*20 , time_steps[error_min]*50,
                       time_steps[error_min]*100 , time_steps[error_min]*150 , time_steps[error_min]*200 , time_steps[error_min]*200 , time_steps[error_min]*200]
# ----------------------------------------------------------------


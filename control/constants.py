#     Neutron star thermal evolution
# --------------------------------------
#            constants.py
# --------------------------------------
# This module provides physics constants
# for numerical computation

# ---------------general physics constants---------------
# -------------------------------------------------------
c = 2.99792458e10                                                # speed of light (cm/s) [exact]
h = 1.05457266e-27                                               # Planck's constant hbar (g*cm^2/s)
G = 6.67259e-8                                                   # grav. constant (cm^3/g/s^2)
e = 4.8032068e-10                                                # electron charge (CGS unit)
kB = 1.380658e-16                                                # Boltzmann constant (erg/K)
sigma = 5.67051e-5                                               # Stefan-Boltzmann const (erg/cm^2/s/K^4)
alpha = (1./137.0359895)                                         # the fine structure constant
pi = 3.141592653589793238                                        # number pi
Ln10 = 2.30258509299405                                          # natural logarithm of 10
s_in_yr = 3.1556925e7                                            # seconds in tropical year
yrtosec = 31536000                                               # from years to seconds factor
# -------------------------------------------------------


# ---------------nuclear physics constants---------------
# -------------------------------------------------------
fm = 1.0e-13                                                     # Fermi (cm)
fm3 = 1.0e-39                                                    # fmi^3 (cm^3)
m_u = 1.6605402e-24                                              # atomic mass unit
Mn0 = 1.6749286e-24                                              # neutron mass (g)
Mp0 = 1.6726231e-24                                              # proton  mass (g)
MH = 1.67353e-24                                                 # hydrogen atom mass (g)
Me = 9.1093897e-28                                               # electron mass (g)
MeV_erg = 1.60217733e-6                                          # transfer factor from MeV to erg
sigma_Th = 6.6524616e-25                                         # Thompson scattering cross-section(cm^2)
Mmu_MeV = 105.658389                                             # muon mass (MeV)
Mlambda_MeV = 1115.683                                           # Lambda-hyperon mass (MeV)
Msigmam_MeV = 1197.449                                           # Sigma_minus-hyperon mass (MeV)
Mmu = (Mmu_MeV * MeV_erg / c / c)                                # Mmu (g)
Mlambda = (Mlambda_MeV * MeV_erg / c / c)                        # Mlambda (g)
Msigmam = (Msigmam_MeV * MeV_erg / c / c)                        # Msigmam (g)
Mn = (0.7 * Mn0)                                                 # neutron effective mass in NS (g)
Mp = (0.7 * Mp0)                                                 # proton effective mass in NS (g)
rho0 = 2.8                                                       # standard nucl. mass density(10^14 g/cc)
n0 = 0.16                                                        # standard nuclear number density (fm^-3)
rho_nd = 4.3e11                                                  # neutron drip density in g/cc
# -------------------------------------------------------


# ---------------Astronomical Constants------------------
# -------------------------------------------------------

MSun = 1.98892e33                                                # solar mass (g)
LSun = 3.846e33                                                  # solar luminosity (erg/s)
RSun = 6.9599e10                                                 # solar radius (cm)
Rearth = 6.37814e8                                               # Earth radius (cm)
Mearth = 5.9737e27                                               # Earth mass (g)
pc = 3.0856775807e18                                             # parsec (cm)
au = 1.4959787066e13                                             # astronomical unit
# -------------------------------------------------------


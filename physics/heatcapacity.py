#     Neutron star thermal evolution
# --------------------------------------
#            heatcapacity.py
# --------------------------------------
# This module provides specific heat capacity
# per unit volume
# _C (erg/K/cm^3)
# for degenetate Fermi-gas (npe)

from other import routines
from control.manager import *
from control.constants import *
import numpy
from physics import sf_gap

def RCp(gap):                                                    # proton superfluidity reduction factor
                                                                 # energy gap for singlet superfluidity

    a = 0.501
    b = 0.4186
    g = 1.007
    d = 1.456

    if(gap<=0.0):                                                # no superfluidity
        R = 1.0
    elif(gap>200.0):
        R = 0.0                                                  # superstrong superfluidity
    elif(gap<=0.01):
        R = 2.426*gap/0.01                                       # smoothing
    else:
        R = numpy.power(b+numpy.sqrt(g*g+(a*gap)*(a*gap)), 2.5)*numpy.exp(d-numpy.sqrt(d*d+gap*gap))

    return R

def RCn(gap):                                                    # neutron superfluidity reduction factor
                                                                 # energy gap for singlet superfluidity

    b = 0.6893
    d = 0.2824
    p1 = 1.934
    p2 = 0.79

    if(gap<=0.0):                                                # no superfluidity
        R = 1.0
    elif(gap>200.0):
        R = 0.0                                                  # superstrong superfluidity
    elif(gap<=0.01):
        R = 2.188*gap/0.01                                       # smoothing
    else:
        R = numpy.power(b+numpy.sqrt(p2*p2+(d*gap)*(d*gap)), 2.0)*numpy.exp(p1-numpy.sqrt(p1*p1+gap*gap))

    return R

# relativistic electron specific heat capacity (erg K^-1 cm^-3) 
def _Ce(T, ne):

    # T -temperatire in K
    # ne - electron number density in fm^-3

    x = routines.pF(ne)/Me/c                                     # pF - Fermi momentum
    Me_eff = Me*numpy.sqrt(1.+x*x)                               # relativistic effective mass

    Ce = Me_eff*routines.pF(ne)*(kB/h)*(kB/h)/3./h*T

    return Ce

def _Cm(T, ne, nm):

    # T - temperatire in K
    # ne - electron number density in fm^-3

    x = routines.pF(ne)/Me/c                                     # pF - Fermi momentum
    Me_eff = Me*numpy.sqrt(1.+x*x)                               # relativistic effective mass

    if(nm <=  1.e-8):
        Cm = 0.
    else:
        Cm = Me_eff*routines.pF(nm)*(kB/h)*(kB/h)/3./h*T
    return Cm

# CORE baryon (neutron + proton) specific heat capacity (erg K^-1 cm^-3)
def _Cb(T, nn, np, Tc_nt, Tc_ns, Tc_p):

    # T - temperature in K
    # nn & np - neutron & proton number densities in fm^-3
    # Tc_n & Tc_p  -  superfluidity critical temperatures for n & p in K
    # Cn, Cp,                          neutron and proton contributions
    # mn_eff, mp_eff,
    # nbb - Feb 2015: inserted
    # Ngap_t, Ngap_s, Pgap, Rn, Rp - superfluidity factors
    # Mn & Mp - neutron & proton effective masses

    if(SUPERFLUIDITY):

        if(T<Tc_nt):
            Ngap_t = sf_gap._TripletGap(T/Tc_nt)                 # neutron triplet energy gap
        else: Ngap_t = 0.

        if(T<Tc_ns):
            Ngap_s = sf_gap._SingletGap(T/Tc_ns)                 # neutron singlet energy gap
        else: Ngap_s = 0.

        if(Ngap_t>Ngap_s):
            Rn = RCn(Ngap_t)
        else:
            Rn = RCp(Ngap_s)
        if(T<Tc_p):
            Pgap = sf_gap._SingletGap(T/Tc_p)                    # proton singlet energy gap
        else:
            Pgap = 0.

        Rp = RCp(Pgap)

    else:

        Rn = Rp = 1.

    if(var_meff == 1):
        mn_eff = routines._Meff(nn)
        mp_eff = routines._Meff(np)
    elif(var_meff == 2):
        nbb = nn+np
        mn_eff = 1.369/(1+3.133*nbb+3.570*nbb*nbb-0.9797*nbb*nbb*nbb)
        mp_eff = 0.9937/(1+2.217*nbb+0.8850*nbb*nbb)
    else:
        mn_eff = Mn/Mn0
        mp_eff = Mp/Mp0

    Cn = (mn_eff * Mn0 * routines.pF(nn) * (kB/h)*(kB/h)/3./h *
          T * Rn)                                                # neutron heat capacity
    Cp = (mp_eff * Mp0 * routines.pF(np) * (kB/h)*(kB/h)/3./h *
          T * Rp)                                                # proton  heat capacity

    return Cn + Cp

# squared ion plasma frequency in s^-2
def _wpi2( ni, Z, A ):

    # ni = ion number density in cm^-3
    # Z = number of protons in nuclei
    # A = number of nucleons in nuclei

    w = 4.*pi*e*e*ni/m_u*Z*Z/A                                   # m_u - atomic mass unit
    return w

def _Ccrust_fit2( t ):
    # t = T/T_pi

    cut = 150.
    par1 = 1.392180E+0002
    par2 = 7.776458E-0002
    par3 = 7.007812E+0001
    par4 = 1.924043E-0001
    par5 = 4.466056E+0003
    par6 = 2.164839E-0005
    par7 = -6.206786E-0003
    par8 = 4.662765E-0001
    par9 = 2.046169E-0001
    par10 = 4.777220E-0006

    v = par4/t

    if(v>cut):
        ee = 0.0
    else:
        ee = numpy.exp(-v)

    v1 = par2/t

    if(v1>cut):
        ee1 = 0.0
    else: ee1 = numpy.exp(-v1)
    v2 = par8/t

    if(v2>cut):
        ee2 = 0.0
    else:
        ee2 = numpy.exp(-v2)

    v3 = par9/t

    if(v3>cut):
        ee3 = 0.0
    else:
        ee3 = numpy.exp(-v3)

    t2 = t*t
    t3 = t2*t
    t4 = t2*t2
    AA = 1.+par1*t2+par5*t4
    C = 1.+2.296*t*ee

    D = ((12.*par3*t3)/AA
    -(7.*par3*t3)*(2.*par1*t2+4.*par5*t4)/AA/AA
    -(par3*t3)*(4.*par1*t2+16.*par5*t4)/AA/AA
    +2.*(par3*t3)*routines.sqr(2.*par1*t2+4.*par5*t4)/AA/AA/AA)

    f = (D+2.296*ee*(2.*t+2.*par4+par4*par4/t)/C
    -routines.sqr(2.296*ee*(t+par4)/C)
    +par6*ee1/t*(-2.*v1+v1*v1)
    +par7*ee2/t/t*(2.-4.*v2+v2*v2)
    +par10*ee3/t/t/t*(6.-6.*v3+v3*v3))

    return( f )


# ion specific heat capacity in the crust   (by D Baiko & DG Yakovlev)
def _Cion(T, ni, Z, A):

    # T = temperature in K
    # ni = ion number density in cm^-3
    # Z = number of protons in nuclei
    # A = number of nucleons in nuclei

    wpi2 = _wpi2(ni, Z, A)
    tt = kB*T/(h*numpy.sqrt(wpi2))                               # temperature in units of ion plasma temperature

    f = _Ccrust_fit2(tt)

    Cion = 3.*kB*ni*f

    return Cion

# free neutron specific heat capacity (for the crust)
def _Cn_crust( T, nn, Tc_nt, Tc_ns ):

    # T - temperature in K
    # nn - free neutron number density in fm^-3
    # Tc_nt - critical temperature of neutron triplet sf
    # Tc_ns - critical temperature of neutron singlet sf

    if(SUPERFLUIDITY):

        if(T<Tc_nt):
            Ngap_t = sf_gap._TripletGap(T/Tc_nt)                 # neutron triplet energy gap
        else:
            Ngap_t = 0.

        if(T<Tc_ns):
            Ngap_s = sf_gap._SingletGap(T/Tc_ns)                 # neutron singlet energy gap
        else:
            Ngap_s = 0.

        if(Ngap_t>Ngap_s):
            R = RCn(Ngap_t)
        else:
            R = RCp(Ngap_s)
    else:
        R = 1.

    mn_eff = Mn/Mn0

    Cn_crust = mn_eff*Mn0*routines.pF(nn)*routines.sqr(kB/h)/3./h*T*R

    return Cn_crust

# Heat capacity (total per unit volume)
def _C(T, rho, nn, ne, nm, Z, Anuc, Tc_nt, Tc_ns, Tc_p):

    # T = temperature in K
    # rho = mass density in g/cm^3
    # nn = number density of free neutrons in fm^-3
    # ne = electron number density in fm^-3                     
    # nm = muon number density in fm^-3
    # Z = number of protons in a nucleus
    # Anuc = number of nucleons in a nucleus, <=  A
    # Tc_nt = crit. temp. of neutron triplet superfluidity in K
    # Tc_ns = crit. temp. of neutron singlet superfluidity in K
    # Tc_p = crit. temp. of proton singlet superfluidity in K
    # Ce = _Ce(T, ne, &dCe) - electron contribution
    # Cm = _Cm(T, ne, nm, &dCm) -  muon  contribution

    Ce = _Ce(T, ne)
    Cm = _Cm(T, ne, nm)

    if(rho > CoreCrustBound):                                    # CORE

        np = ne+nm
        Cb = _Cb(T, nn, np, Tc_nt, Tc_ns, Tc_p)                  # baryon contribution
        Cv = Ce+Cm+Cb

    else:                                                        # CRUST

        ni = ne/Z/fm3
        Cion = _Cion(T, ni, Z, Anuc)                             # ion contribution
        Cn_crust = _Cn_crust(T, nn, Tc_nt, Tc_ns)                # free neutron contribution

        Cv = Ce+Cion+Cn_crust

    return Cv


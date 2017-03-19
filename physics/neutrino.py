#     Neutron star thermal evolution
# --------------------------------------
#            neutrino.py
# --------------------------------------
# This module provides specific neutrino
# emissivity per unit volume
# via URCA process (Direct and Modified)
# _Q (erg/s/cm^3)

import numpy

from physics           import sf_gap
from other             import routines
from control.manager   import *
from control.constants import *
from data              import loaddata


Nnu1 = 2.                                                        # number of neutrino species other than the electron
sin2TW = 0.2319                                                  # sine squared of the Weinberg angle
Cv = (0.5+2.*sin2TW)
Ca2 = 0.75
Cv2 = (routines.sqr(Cv)+Nnu1*routines.sqr(1.-(Cv)))
Tr = 5.9302e9                                                    # electron relativistic temperature
Qc = 1.0226e23                                                   # Compton neutrino loss unit
Nnu = 3.                                                         # number of neutrino types


def fun_1(t):

    p1 = 7.662048E+00
    p2 = 1.919637E+00

    F = 2.*t/routines.sqr(2.*pi)*numpy.sqrt((2.*pi*t+p1*t*t+p2*t*t*t)/(1.+p2/4.*t))*numpy.exp(-1./t)

    return F

def fun0(t):

    p1 = 2.361340E+01
    p2 = 3.210783E+01

    t2 = t*t
    t3 = t2*t
    t4 = t3*t

    u = 2.*pi*t + p1*t2 + p2*t3 + 16.*t4
    F = 2.*t/routines.sqr(2.*pi)*numpy.sqrt(u)*numpy.exp(-1./t)

    return F

def fun1(t):

    p1 = 4.243553E+01
    p2 = 1.407638E+02
    p3 = 2.651583E+02
    p4 = 2.878642E+02

    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    t6 = t5*t

    u = 2.*pi*t + p1*t2 + p2*t3 + p3*t4 + p4*t5 + 144.*t6
    F = 2.*t/routines.sqr(2.*pi)*numpy.sqrt(u)*numpy.exp(-1./t)

    return F

def fun2(t):

    p1 = 6.133312E+01
    p2 = 3.219052E+02
    p3 = 1.153307E+03
    p4 = 2.624424E+03
    p5 = 4.468395E+03
    p6 = 4.600400E+03

    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    t6 = t5*t
    t7 = t6*t
    t8 = t7*t

    u = 2.*pi*t + p1*t2 + p2*t3 + p3*t4 + p4*t5 + p5*t6 + p6*t7 + 144.*16.*t8
    F = 2.*t/routines.sqr(2.*pi)*numpy.sqrt(u)*numpy.exp(-1./t)

    return F


def _Qpair(T_K, x):                                              # Relativistic version, cold plasma

    t = T_K/Tr
    mu = numpy.sqrt(1.+x*x)

    if(mu/t > 500.):
        F = 0.

    else:

        f_1 = fun_1(t)                                           # positron functions
        f0 = fun0(t)
        f1 = fun1(t)
        f2 = fun2(t)
        u_1 = 0.5/routines.sqr(pi)*(x*mu-numpy.log(x+mu))        # electron functions
        u0 = x*x*x/3./routines.sqr(pi)
        u1 = (mu*x*(x*x+mu*mu)-numpy.log(x+mu))/8./routines.sqr(pi)
        u2 = x*x*x*(1.+0.6*x*x)/3./routines.sqr(pi)
        u = ((Cv2+Ca2)*(8.*(u1*f2+u2*f1)+7.*(u1*f0+u0*f1)-
                      2.*(u2*f_1+u_1*f2)+5.*(u_1*f0+u0*f_1))+
                      9.*(Cv2-Ca2)*(f0*(u1+u_1)+u0*(f_1+f1)))

        F = Qc*numpy.exp(-mu/t)*u/(36.*pi)

    return F

# Neutrino emission due to Cooper pairing of neutrons and protons
def _Qcooper(T, nn, np, Tc_nt, Tc_ns, Tc_p):

    # T = temperature in K
    # rho = matter density in g cm^-3
    # nn = neutron number density in fm^-3
    # np = proton number density in fm^-3
    # Tc_nt = crit. T of neutron triplet superfluidity in K
    # Tc_ns = crit. T of neutron singlet superfluidity in K
    # Tc_p = crit. T of proton superfluidity in K
    # n0 = 0.16 fm^-3 -- nuclear number density

    T7 = numpy.power(T/1.e9,7)

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

    v1 = v2 = R_recomb = 0.                                      # neutron Cooper's neutrions

    if(T<Tc_ns):
        v1 = sf_gap._SingletGap(T/Tc_ns)                         # neutron singlet gap
    v1=0

    if(T<Tc_nt):
        v2 = sf_gap._TripletGap(T/Tc_nt)                         # neutron triplet gap

    if(v1>v2):                                                   # recombination in singlet state
                                                                 # rms error = 0.0045, max error = 0.0096 at v1 = 4

        if(v1>0. and v1<30.):

            R_recomb = (v1*v1*(0.602 + 0.5942*v1*v1 + 0.2880*v1*v1*v1*v1)*
                        numpy.sqrt(0.5547+numpy.sqrt(routines.sqr(0.4453) + routines.sqr(0.1063)*v1*v1))*
                        numpy.exp(-numpy.sqrt(4*v1*v1 + routines.sqr(2.245)) + 2.245))

    else:                                                        # recombination in triplet state
                                                                 # rms error = 0.0102, max error = 0.0338 at v2 = 2

        if(v2 > 0. and v2 < 30.):
            R_recomb = (v2*v2*(0.602*2.+3.733*v2*v2+0.3191*v2*v2*v2*v2)*
                        routines.sqr(0.7591+numpy.sqrt(routines.sqr(0.2409)+0.3145*v2*v2))
                        /(1.+0.3511*v2*v2)*numpy.exp(-numpy.sqrt(4.*v2*v2+routines.sqr(0.4616))+0.4616))
        R_recomb *=  4.17                                        # to account for axial-vector term

    Qn = 3.*1.17e21*T7*0.353*numpy.power(nn/n0,1./3.)*mn_eff*R_recomb

    if(T<Tc_p):                                                  # proton Cooper's neutrions
        v =  sf_gap._SingletGap(T/Tc_p)                          # proton singlet gap

    else:
        v = 0.

    if(v>0. and v<100.):                                         # pp recombination
                                                                 # rms error = 0.0045, max error = 0.0096 at v = 4
        R_recomb = (v*v*(0.602 + 0.5942*v*v + 0.2880*v*v*v*v)*
                    numpy.sqrt(0.5547+numpy.sqrt(routines.sqr(0.4453) + routines.sqr(0.1063)*v*v))*
                    numpy.exp(-numpy.sqrt(4*v*v + routines.sqr(2.245)) + 2.245))
    else:
        R_recomb = 0.

    Qp = (3.*1.17e21*T7*0.353*numpy.power(np/n0,1./3.)*mp_eff*R_recomb*
          (routines.sqr(0.08) + routines.sqr(1.26)*(11./42. + routines.sqr(mp_eff))*
          routines.sqr(0.353/mp_eff)*numpy.power(np/n0,2./3.)))

    return Qn*0.8+Qp

def _Qbremss_eZ(T, rho):

    # T = temperature of matter in K
    # rho = mass density in g/cc

    tau = numpy.log10(T*1.e-8)
    r = numpy.log10(rho*1.e-12)

    lgQ = (11.204 + 7.304*tau + 0.2976*r - 0.37*tau*tau + 0.188*tau*r - 0.103*r*r +
           0.0547*tau*tau*r - 6.77*numpy.log10(1. + 0.228*rho/2.8e14))

    Q = numpy.power(10,lgQ)

    return Q


def _Qcrust(T, rho, ne, nn, Tc_nt, Tc_ns):

    Qplasma = Qbremss = Qpair = Qnn = 0.

    x = routines.pF(ne)/(Me*c)                                   # electron relativistic parameter

    Qbremss = _Qbremss_eZ(T, rho)
    Qplasma = _Qplasma(T, x)
    Qpair = _Qpair(T, x)

    if((T<Tc_nt or T<Tc_ns) and SUPERFLUIDITY):

        np = Tc_p = 0.
        Qcooper = _Qcooper(T, nn, np, Tc_nt, Tc_ns, Tc_p)

    else:
        Qcooper = 0.

    Q = Qbremss + Qnn + Qpair + Qplasma + Qcooper

    return Q

# Neutrino emission due to plasmon decay.
def _Qplasma(T, x):                                              # fit by D. G. Yakovlev, 1999

    tt = T/Tr

    ne = x*x*x/(3.*pi*pi)*numpy.power((Me*c/h),3.)
    wpe = numpy.sqrt(4.*pi*e*e*ne/Me/numpy.sqrt(1.+x*x))
    fp = h*wpe/(kB*T)

    Ipl = numpy.power(tt,9.)*(16.23*numpy.power(fp,6.)+4.604*numpy.power(fp,7.5))*numpy.exp(-fp)

    Q = Qc/(96.*numpy.power(pi,4.)*alpha)*Cv2*Ipl

    return Q

def _Rn(Tgap):                                                   # neutron superfluidity reduction factor for DURCA
                                                                 # Tgap - triplet-state energy gap
    
    a = 0.2546206
    b = 0.128407
    d = 2.701395

    if(Tgap == 0.0):
        R = 1.0                                                  # no superfluidity
    elif(Tgap>200.0):
        R = 0.0                                                  # superstrong superfluidity
    else:
        R = (numpy.power(a + numpy.sqrt(routines.sqr(1.0 - a) + routines.sqr(b*Tgap)), 5.0)*
             numpy.exp(d-numpy.sqrt(d*d + Tgap*Tgap)))

    return R

def _Rp(Sgap):                                                   # proton superfluidity reduction factor for DURCA
                                                                 # Sgap - singlet-state energy gap

    a = 0.2311517
    b = 0.1438319
    d = 3.427262

    if(Sgap == 0.0):
        R = 1.0                                                  # no superfluidity
    elif(Sgap>200.0):
        R = 0.0                                                  # superstrong superfluidity
    else:
        R = (numpy.power(a+numpy.sqrt(routines.sqr(1.0-a)+routines.sqr(b*Sgap)), 5.5)*
             numpy.exp(d-numpy.sqrt(d*d+Sgap*Sgap)))
    return R

def _Rnp(Ngap, Pgap):

    if((Ngap>200.0) or (Pgap>200.0)):
        return(0.0)

    lgR = loaddata.lgR(Ngap, Pgap)
    R = numpy.power(10,lgR)

    return R

# reduction of neutron-branch of MURCA by singlet superfluidity 
def _Rmodn_s(gap):

    if(gap == 0.0):
        R = 1.0                                                  # no superfluidity
    elif(gap>200.0):
        R = 0.0                                                  # superstrong superfluidity
    else :

        a = 0.1477+numpy.sqrt(routines.sqr(0.8523)+routines.sqr(0.1175*gap))
        b = 0.1477+numpy.sqrt(routines.sqr(0.8523)+routines.sqr(0.1297*gap))
        R = 0.5*(numpy.power(a,7.5)+numpy.power(b,5.5))*numpy.exp(3.437-numpy.sqrt(routines.sqr(3.437)+gap*gap))

    return R

# reduction of proton-branch of MURCA by singlet superfluidity 
def _Rmodp_s(gap):

    if(gap == 0.0):
        R = 1.0                                                  # no superfluidity
    elif(gap>200.0):
        R = 0.0                                                  # superstrong superfluidity
    else:

        a = 0.2414+numpy.sqrt(routines.sqr(0.7586)+routines.sqr(0.1318*gap))
        R = numpy.power(a, 7.0)*numpy.exp(5.339-numpy.sqrt(routines.sqr(5.339)+4.0*gap*gap))

    return R

# reduction of proton-branch of MURCA by triplet superfluidity 
def _Rmodp_t(gap):

    if(gap == 0.0):
        R = 1.0                                                  # no superfluidity
    elif(gap>200.0):
        R = 0.0                                                  # superstrong superfluidity
    else:

        a = 0.1612+numpy.sqrt(routines.sqr(0.8388)+routines.sqr(0.1117*gap))
        b = 0.1612+numpy.sqrt(routines.sqr(0.8388)+routines.sqr(0.1274*gap))
        R = 0.5*(numpy.power(a,7.0)+numpy.power(b,5.0))*numpy.exp(2.398-numpy.sqrt(routines.sqr(2.398)+gap*gap))

    return R

# reduction of neutron-branch of MURCA by triplet superfluidity
# asymptote for T << Tc (tau = T/Tc << 1)
def _Rmodn_t_asy(tau):

    if(tau>1):
        R = 1.0                                                  # no superfluidity
    elif(tau<1.0e-3):
        R = 0.0                                                  # superstrong superfluidity
    else:
        R = 1.56e-4*numpy.power(tau,-6.0)*numpy.exp(-2.376/tau)  # asymptote!
    return R

# reduction of p-MURCA by combined neutron and proton superfluidity 
def _Rmodp(Ngap, Pgap):

    if((Ngap>200.0)or(Pgap>200.0)):
        return(0.0)

    Rn = _Rn(Ngap)
    Rnp = _Rnp(Ngap, 2.0*Pgap)

    if(Rn>1.0e-10):
        R = Rnp*_Rmodp_t(Ngap)/Rn
    else:
        R = Rnp*0.0023*Ngap*Ngap

    return R

# reduction of n-MURCA by combined neutron and proton superfluidity 
def _Rmodn(Ngap, Pgap):

    if((Ngap>200.0) or (Pgap>200.0)):
        return(0.0)

    Rp = _Rp(Pgap)
    Rnp = _Rnp(2.0*Ngap, Pgap)

    if(Rp>1.0e-10):
        R = Rnp*_Rmodn_s(Pgap)/Rp
    else:
        R = Rnp*0.0023*Pgap*Pgap

    return R

# reduction of np-bremsstrahlung by singlet superfluidity 
def _Rbremss_np_s(gap):

    if(gap<= 0.0):
        R = 1.0                                                  # no superfluidity
    elif(gap>100.0):
        R = 0.0                                                  # superstrong superfluidity
    else:

        a = 0.9982+numpy.sqrt(routines.sqr(0.0018)+routines.sqr(0.3815*gap))
        b = 0.3949+numpy.sqrt(routines.sqr(0.6051)+routines.sqr(0.2666*gap))
        R = 1.0/2.732*( a*numpy.exp(1.306-numpy.sqrt(routines.sqr(1.306)+gap*gap)) +
                      1.732*numpy.power(b, 7.0)*numpy.exp(3.303-numpy.sqrt(routines.sqr(3.303)+4.0*gap*gap)) )

    return R

# reduction of pp-bremsstrahlung by singlet superfluidity 
def _Rbremss_pp_s(gap):

    if(gap<= 0.0):
        R = 1.0                                                  # no superfluidity
    elif(gap>100.0):
        R = 0.0                                                  # superstrong superfluidity
    else:
        a = 0.1747+numpy.sqrt(routines.sqr(0.8253)+routines.sqr(0.07933*gap))
        b = 0.7333+numpy.sqrt(routines.sqr(0.2667)+routines.sqr(0.1678*gap))
        R = 0.5*( a*a*numpy.exp(4.228-numpy.sqrt(routines.sqr(4.228)+4.0*gap*gap)) +
                numpy.power(b, 7.5)*numpy.exp(7.762-numpy.sqrt(routines.sqr(7.762)+9.0*gap*gap)) )

    return R

# reduction of np-bremssstrahlung by neutron and proton superfluidity 
# *** UNDER CONSTRUCTION *** 
def _Rbremss_np(Ngap, Pgap):

    if((Ngap>200.0) or (Pgap>200.0)):
        return(0.0)
    Rp = _Rp(Pgap)
    Rnp = _Rnp(Ngap, Pgap)

    if(Rp>1.0e-30):
        R = Rnp*_Rbremss_np_s(Pgap)/Rp
    else:
        R = 0.0

    return R

def Threshold( fB1, fB2, fL ):

    cond1 = ( fB2 + fL >=  fB1 )
    cond2 = ( fB1 >=  numpy.fabs(fB2-fL) )
    cond3 = ( (fB1>0.) and (fB2>0.) and (fL>0.))
    cond = (cond1 and cond2 and cond3)

    if(cond1):
        return 1
    else:
        return 0

def _Qurca(T, nn, ne, nm, Tc_nt, Tc_ns, Tc_p):

    # T = temperature in K
    # nn = number density of free neutrons in fm^-3
    # ne = electron number density in fm^-3
    # nm = muon number density in fm^-3
    # Tc_nt = crit. temp. of neutron triplet superfluidity in K
    # Tc_ns = crit. temp. of neutron singlet superfluidity in K
    # Tc_p = critical temperature of proton superfluidity in K
    # a_n, a_p, a_nn, a_np, a_pp - Born corrections
    # b_n, b_p, b_nn, b_np, b_pp - non-Born corrections
    # Ngap, Pgap, Ngap1 - neutron & proton SF gaps
    # Rdir, Rmodn, Rmodp - reduction factors for URCAs
    # Rbremss_nn, Rbremss_np, Rbremss_pp - reduction factors for Bremss

    if(T<Tc_nt):
        Ngap = sf_gap._TripletGap(T/Tc_nt)                       # neutron triplet energy gap

    else:
        Ngap = 0.

    #  Approximate treatment of singlet n sf in the core!
    if(T<Tc_ns):
        Ngap1 = sf_gap._SingletGap(T/Tc_ns)                      # neutron singlet energy gap
    else:
        Ngap1 = 0.

    if (Ngap1>Ngap):
        Ngap = Ngap1

    if(T<Tc_p):
        Pgap = sf_gap._SingletGap(T/Tc_p)      # proton singlet  energy gap
    else:
        Pgap = 0.

    Rdir = Rmodn = Rmodp = Rbremss_nn = Rbremss_np = Rbremss_pp = 1.

    if(SUPERFLUIDITY):

        if((T<Tc_nt) and (T>= Tc_p)):                            # neutron superfluidity

            Rdir = _Rn(Ngap)                                     # reduction of  DURCA  by neutrons
            Rmodn = _Rmodn(Ngap, 0.)                             # reduction of n-MURCA by neutrons
            # asymptote is _Rmodn_t_asy(T/Tc_n)
            Rmodp = _Rmodp_t(Ngap)                               # reduction of p-MURCA by neutrons
            Rbremss_nn = _Rbremss_pp_s(Ngap)                     # *** UNDER CONSTRUCTION ***
            Rbremss_np = _Rbremss_np(Ngap, 0.)                   # *** UNDER CONSTRUCTION ***
            Rbremss_pp = 1.                                      # no reduction of pp-bremss

        elif((T>= Tc_nt) and (T<Tc_p)):                          # proton superfluidity

            Rdir = _Rp(Pgap)                                     # reduction of  DURCA  by protons
            Rmodn = _Rmodn_s(Pgap)                               # reduction of n-MURCA by protons
            Rmodp = _Rmodp_s(Pgap)                               # reduction of p-MURCA by protons
            Rbremss_nn = 1.                                      # no reduction of nn-bremss
            Rbremss_np = _Rbremss_np_s(Pgap)                     # reduction of np-bremss by protons
            Rbremss_pp = _Rbremss_pp_s(Pgap)                     # reduction of pp-bremss by protons

        elif((T<Tc_nt) and (T<Tc_p)):                            # neutron and proton superfluidity

            Rdir = _Rnp(Ngap, Pgap)                              # reduction of  DURCA  by the two
            Rmodn = _Rmodn(Ngap, Pgap)                           # reduction of n-MURCA by the two
            Rmodp = _Rmodp(Ngap, Pgap)                           # reduction of p-MURCA by the two
            Rbremss_nn = _Rbremss_pp_s(Ngap)                     # *** UNDER CONSTRUCTION ***
            Rbremss_np = _Rbremss_np(Ngap, Pgap)                 # *** UNDER CONSTRUCTION ***
            Rbremss_pp = _Rbremss_pp_s(Pgap)                     # *** UNDER CONSTRUCTION ***

    T6 = numpy.power(T*1.e-9, 6.)
    T8 = numpy.power(T*1.e-9, 8.)

    np = ne+nm
    onethird = 1./3.
    nn13 = numpy.power(nn/n0, onethird)

    ne13 = numpy.power(ne/n0,onethird)
    nm13 = numpy.power(nm/n0,onethird)
    np13 = numpy.power((ne+nm)/n0,onethird)

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

    # Direct URCA process with electrons:
    #  n -> p + e + nu_e~,  p + e -> n + nu_e
    if(Threshold(nn13, np13, ne13)):
        Qdurca_e = 4.0e27*(mn_eff)*(mp_eff)*ne13*T6              # DURCA emissivity
        #  Qdurca_e = 1.0e23*T6    # to simulate pi-K-condensation

    else:
        Qdurca_e = 0.

    Qdurca_e *= Rdir                                             # superfluidity reduction factor

    #  *****  Direct URCA process with muons:   *****
    #  n -> p + mu + nu_mu~,  p + mu -> n + nu_mu
    if(Threshold(nn13, np13, nm13)):
        Qdurca_m = 4.0e27*(mn_eff)*(mp_eff)*ne13*T6              # DURCA emissivity
        #  here we take into account that me_eff = mm_eff = chem. potential_e/c^2
        #  Qdurca_m = 1.0e25*T6  # to simulate pi-condensation

    else:
        Qdurca_m = 0.

    Qdurca_m *= Rdir                                             # superfluidity reduction factor

    Qdurca = Qdurca_e + Qdurca_m

    #  Modified URCA process with electrons:
    #  n + n -> n + p + e + nu_e~,  n + p + e -> n + n + nu_e
    #  p + n -> p + p + e + nu_e~,  p + p + e -> p + n + nu_e

    a_n = a_p = 1.76-0.63/(nn13*nn13)
    b_n = b_p = 0.68

    Qmurca_n_e = 8.55e21*numpy.power(mn_eff, 3.)*(mp_eff)*np13*T8*a_n*b_n

    if(nn13 >= 3.*np13 - ne13):

        p_factor1 = routines.sqr(ne13+3.*np13-nn13)/(8.*ne13*np13)
        Qmurca_p_e  = 8.53e21*numpy.power(mp_eff, 3.)*(mn_eff)*np13*T8*a_p*b_p*p_factor1

    else:

        p_factor2 = 1./2.*(3. - nn13/np13)
        Qmurca_p_e  = 8.53e21*numpy.power(mp_eff, 3.)*(mn_eff)*np13*T8*a_p*b_p*p_factor2

    Qmurca_n_e*= Rmodn                                           # superfluidity reduction factor
    Qmurca_p_e*= Rmodp                                           # superfluidity reduction factor

    #  if(ne < nn/64.) Qmurca_p = 0. - momentum conservation condition

    if(nn13 > 3.*np13 + ne13):
        Qmurca_p_e  = 0.

    Qmurca_e = Qmurca_n_e + Qmurca_p_e        # total electron MURCA emissivity

    #  Modified URCA process with muons:
    #  n + n -> n + p + mu + nu_m~,  n + p + mu -> n + n + nu_m
    #  p + n -> p + p + mu + nu_m~,  p + p + mu -> p + n + nu_m

    Corrector = nm13/ne13
    Qmurca_n_m = 8.55e21*numpy.power(mn_eff, 3.)*(mp_eff)*np13*T8*a_n*b_n*Corrector

    if(nn13 >=  3.*np13 - nm13):
        Qmurca_p_m = (8.53e21*numpy.power(mp_eff, 3.)*(mn_eff)*np13*T8*a_p*
                      b_p*routines.sqr(nm13+3.*np13-nn13)/(8.*ne13*np13))

    else:

        p_factor2 = 1./2.*(3. - nn13/np13)
        Qmurca_p_m  = 8.53e21*numpy.power(mp_eff, 3.)*(mn_eff)*np13*T8*a_p*b_p*p_factor2*Corrector

    Qmurca_n_m *= Rmodn                                           # superfluidity reduction factor
    Qmurca_p_m *= Rmodp                                           # superfluidity reduction factor

    if(nn13 > 3.*np13 + nm13):
        Qmurca_p_m  = 0.

    Qmurca_m = Qmurca_n_m + Qmurca_p_m                            # total electron MURCA emissivity
    Qmurca = Qmurca_e + Qmurca_m

    #  Bremsstrahlung:
    #  n + n -> n + n + nu + nu~
    #  n + p -> n + p + nu + nu~
    #  p + p -> p + p + nu + nu~

    a_nn = 0.59
    a_np = 1.06
    a_pp = 0.11
    b_nn = 0.56
    b_np = 0.66
    b_pp = 0.7

    QBnn = 7.4e19*numpy.power(mn_eff, 4.)*nn13*T8*a_nn*b_nn*Nnu
    QBpp = 7.4e19*numpy.power(mp_eff, 4.)*np13*T8*a_pp*b_pp*Nnu
    QBnp = 1.5e20*routines.sqr(mn_eff)*routines.sqr(mp_eff)*np13*T8*a_np*b_np*Nnu

    QBnn*= Rbremss_nn                                            # superfluidity reduction factor
    QBpp*= Rbremss_pp                                            # superfluidity reduction factor
    QBnp*= Rbremss_np                                            # superfluidity reduction factor

    Qbremss = QBnn + QBpp + QBnp                                 # total Bremsstrahlung emissivity

    Q = Qdurca + Qmurca + Qbremss

    return Q

# Neutrino emissivity (total, redshifted, per unit volume)
def _Q(T, nn, ne, nm, rho, phi, Tc_nt, Tc_ns, Tc_p):

    # T = temperature in K
    # nn = number density of free neutrons in fm^-3
    # ne = electron number density in fm^-3
    # nm = muon number density in fm^-3
    # Vion = fraction of volume occupied by ions in the crust
    # rho =  mass density in g/cm^3
    # phi = dimensionless gravitational potential
    # Tc_nt = crit. temp. of neutron triplet superfluidity in K
    # Tc_nt = crit. temp. of neutron singlet superfluidity in K
    # Tc_p = crit. temp. of proton singlet superfluidity in K

    CoreCrustBound2 = 0.9999667*CoreCrustBound
    redshift = numpy.exp(2.*phi)                                 # gravitational redshift factor
    np = ne + nm

    if(rho < CoreCrustBound):                                    # CRUST

        Qcrust = _Qcrust(T, rho, ne, nn, Tc_nt, Tc_ns)
        Qcrust *= redshift

    if(rho > CoreCrustBound2):                                   # CORE

        Qurca = _Qurca(T, nn, ne, nm, Tc_nt, Tc_ns, Tc_p)

        if((T<Tc_nt or T<Tc_ns or T<Tc_p) and SUPERFLUIDITY):    # if superfluidity
            Qcooper = _Qcooper(T, nn, np, Tc_nt, Tc_ns, Tc_p)

        else:
            Qcooper = 0.

        Qcore = (Qurca + Qcooper)*redshift

    # Smoothing the core-crust boundary jump

    if(rho <=  CoreCrustBound2):
        Q = Qcrust
    elif (rho >=  CoreCrustBound):
        Q = Qcore
    else:
        Q = ((rho-CoreCrustBound2)/(CoreCrustBound-CoreCrustBound2)*Qcore +
             (CoreCrustBound-rho)/(CoreCrustBound-CoreCrustBound2)*Qcrust)

    return Q


#     Neutron star thermal evolution
# --------------------------------------
#          thermalconductivity.py
# --------------------------------------
# This module provides thermal
# conductivity per unit volume
# _kappa(erg/cm^3/K/s)
# for degenetate Fermi-gas (npe)

from other import routines
from physics import sf_gap
import numpy
from control.manager import *
from control.constants import *

def coulog2(PF0, gamma, Zion, CMI, xnuc, xnuct):

    # Version 24.02.00
    # Fit to the Coulomb logarithms CLel (for electric conductivity)
    # and CLth (for thermal conductivity) of dense Coulomb plasmas
    # Input: PF0 = p_F / mc - relativity (density) parameter,
    # gamma - Coulomb coupling parameter of ions,
    # Zion - mean charge of the ion,
    # CMI - mean atomic mass number
    # Dimensional quantities are in the relativistic units (m_e = \hbar = c = 1)

    Uminus1 = 2.8
    Uminus2 = 12.973
    AUM = 1822.8
    BOHR = 137.036

    # Uminus1,Uminus2 - dimensionless frequency moments of phonons
    # AUM - atomic mass unit divided by the electron mass
    # BOHR - radius of the first Bohr orbit in the rel.units

    #   ----------------------   Preliminaries   ---------------------------

    DENS = PF0**3/3./pi**2                                       # number density of electrons
    DENSI = DENS/Zion                                            # number density of ions (rel.)
    SPHERION = (.75/pi/DENSI)**.3333333                          # Ion sphere radius
    EF0 = numpy.sqrt(1.+PF0**2)                                  # Fermi energy
    VF0 = PF0/EF0                                                # Fermi velocity
    PM2 = (2.*PF0)**2                                            # squared max.momentum transfer
    Q2icl = 3.*gamma/SPHERION**2                                 # squared Debye momentum
    Q2e = 4./pi/BOHR*PF0*EF0                                     # k_{TF}^2 (strong degen.)
    TRP = Zion/gamma*numpy.sqrt(CMI*AUM*SPHERION/3./BOHR)        # = T/T_p
    BORNCOR = VF0*Zion*pi/BOHR                                   # first non-Born correction

    #   ---------------------------------------------------------------

    XS = ((Q2icl*(1.+.06*gamma)*numpy.exp(- numpy.sqrt(gamma))+Q2e)/PM2*
       numpy.exp(-BORNCOR))                                      # eff.screening param.
    XW = Uminus2*PM2/Q2icl*(1.+.3333*BORNCOR)                    # eff. Debye-Waller param.

    # Modification WITH FINITE SIZES OF NUCLEI; xnuc = r_{nuc}/a_i

    XW1 = 14.7327*xnuc**2                                        # = 4(9\pi/4)^{2/3} x_{nucl}^2  = coeff.at q^2
    XW1 = XW1*(1.+.3333*BORNCOR)*(1.+Zion/13.* numpy.sqrt(xnuc))
    CL = coulan2(XS,XW,VF0,XW1)
    A0 = 1.683* numpy.sqrt(PF0/CMI/Zion)                         # zero-vibr.param.(Baiko&Yakovlev95)
    VIBRCOR =  numpy.exp(-A0/4.*Uminus1* numpy.exp(-9.1*TRP))    # corr.for zero-vibr.
    T0 = .19/Zion**.16667                                        # critical T/T_p parameter
    G0 = TRP/ numpy.sqrt(TRP**2+T0**2)*(1.+(Zion/125.)**2)       # mod.10.01.99
    CLel = CL*G0*VIBRCOR                                         # 1st FIT
    G2 = TRP/ numpy.sqrt(.0081+TRP**2)**3
    THtoEL = 1.+G2/G0*(1.+BORNCOR*VF0**3)*.0105*(1.-1./Zion)*(1.e0+xnuct**2* numpy.sqrt(2.e0*Zion))
    TRU = TRP*3.*VF0*BOHR/Zion**.3333333                         # T/Tu
    if(TRU<20.):                                                 # correction for dying-out umklapp processes

        CLhigh = CLel
        EU = numpy.exp(-1.e0/TRU)
        CLlowK = 50.* numpy.sqrt(PF0/CMI)/Zion*TRP**3            # low-T lim.for kappa
        CLlowS = CLlowK/VF0/BOHR/.75*TRP**2                      # low-T lim.for sigma
        CLel = CLhigh*EU+CLlowS*(1.-EU)
        THtoEL = (CLhigh*THtoEL*EU+CLlowK*(1.-EU))/CLel

    CLth = CLel*THtoEL

    return CLel, CLth

def coulan2(XS, XW0, V, XW1):

    # Analytic expression for Coulomb logarithm - Version 23.05.00
    # XS = (q_s/2p)^2, where p - momentum, q_s - eff.scr.momentum
    # XW0 = u_{-2} (2p/\hbar q_D)^2, where u_{-2} = 13, q_D^2 = 3\gamma/a_i^2
    # V = p/(mc)
    # XW1 = s1*(2p/\hbar)^2, s1 \approx r_{nuc}^2

    EPS = 1.e-2
    EPS1 = 1.e-3
    EULER = 0.5772156649

    if(XS<0. or XW0<0. or V<0. or XW1<0.):
        return 0
    
    for I in range(0,2):
        
        if(I== 0):
            
            XW = XW0+XW1
            B = XS*XW
            
        else:                                                    # to do the 2nd term
            XW = XW1
            B = XS*XW
        if(I == 0 or key == 2):                                  # 23.05.00: for key = 2 re-check
            
            # Check applicability of asymptotes:
            if(XW<EPS):
                key = 1
                pass
            else:
                if(XW > 1./EPS and B > 1./EPS):
                    key = 2
                elif(XS<EPS1 and B<EPS1/(1.+XW)):
                    key = 3
                else:
                    key = 4
                    
        EA = numpy.exp(-XW)
        E1 = 1.-EA
        
        if(key != 1):
            E2 = (XW-E1)/XW

        if(key == 1):

            CL0 = numpy.log((XS+1.)/XS)
            CL1 = .5*XW*(2.-1./(XS+1.)-2.*XS*CL0)
            CL2 = .5*XW*(1.5-3.*XS-1./(XS+1.)+3.*XS**2*CL0)

        elif(key == 2):

            CL0 = numpy.log(1.e0+1./XS)
            CL1 = (CL0-1.e0/(1.+XS))/2.
            CL2 = (2.*XS+1.)/(2.*XS+2.)-XS*CL0

        elif(key == 3):
            CL1 = .5*(EA*expint(XW,0)+numpy.log(XW)+EULER)
            CL2 = .5*E2

        elif(key == 4):

            CL0 = numpy.log((XS+1.)/XS)
            EL = expint(B,0)-expint(B+XW,0)*EA
            CL1 = .5*(CL0+XS/(XS+1.)*E1-(1.+B)*EL)
            CL2 = .5*(E2-XS*XS/(1.+XS)*E1-2.*XS*CL0+XS*(2.+B)*EL)

        else:
            return 0
            
        if(I == 0):                                              # 1st term calculated
            coulan2 = CL1-V**2*CL2

            if(XW1<EPS1):
                return 0                                         # don't calculate the 2nd term
        else:

            coulan2 = coulan2-(CL1-V**2*CL2)                     # 2nd term calculated
            
    return coulan2

def expint(XI,L):                                                # expint = e^XI E_{L+1}(XI)

    gamma = .5772156649e0
    Nrep = 21

    if(XI > 1.):

        CL = L
        CI = Nrep
        C = 0.

        for I in range (Nrep,0,-1):

            C = CI/(XI+C)
            C = (CL+CI)/(1.+C)
            CI = CI-1.

        Q0 = 1./(XI+C)

    else:
        
        PSI = -gamma

        for K in range(1,L+1):
            PSI = PSI+1./K

        Q0 = 0.
        CMX = 1.
        CL = L
        CM = -1.

        for M in range(0,Nrep+1):
            
            CM = CM+1.
            
            if(M != 0):
                CMX = -CMX*XI/CM
            if(M != L) :
                DQ = CMX/(CM-CL)
            else:
                DQ = CMX*(numpy.log(XI+1.e-20)-PSI)
                
            Q0 = Q0-DQ
            
        Q0 = numpy.exp(XI)*Q0
        
    return Q0

def soyam(t,x):

    # This factor corrects the volume for nucl.shape acc.to Oyamatsu:
    # nn_in = nn_out+dn_n*(1-(r/Rn)^t); x = min(1,Rws/Rn); same for protons.

    return x**3-9.*x**(3.+t)/(3.+t)+9.*x**(3.+2.*t)/(3.+2.*t)-x**(3.+3.*t)/(1.+t)

def oyaform(bard,index):

    # This block has been copied from "conrt.pas" (D.G.Yakovlev)
    # and converted into Fortran. It realizes the SMOOTH COMPOSITION model
    # Version 24.02.99
    # Input: bard (baryon density in fm^{-3}), index (of phase)
    # Output: Z (total number of protons inside the nucleus)
    # Anuc (num.of barions within the nucleus), A (within the cell)
    # tp (smoothness param.of proton core)
    # xnuc (effective proton-core radius divided by the WS cell radius)
    # xnuct (a second proton-core parameter for use in a quantum crystal)
    # Internal variables:
    # Nin = A-Z - total number of neutrons inside the nucleus
    # Nfree - total number of free neutrons (incl.ones penetrating nuclei)
    
    if(index == 30) : # {densities lower than the neutron drip}
        
        f = numpy.log(1.e0+bard/5.0e-9)
        Rp = 5.688+0.02628*f+0.009468*f*f                        # max.proton core radius
        Rn = 5.788+0.02077*f+0.01489*f*f                         # max.neutron core radius
        np_in = 0.0738+1.22e-4*f-1.641e-4*f*f                    # centr.num.dens.of protons
        nn_in = 0.0808+1.688e-4*f+9.439e-5*f*f                   # same for neutrons
        nn_out = 0.0
        tp = 6.e0
        tn = tp
        Nin = pi/.75*Rn**3*nn_in*soyam(tn,1.e0)
        Z = pi/.75*Rp**3*np_in*soyam(tp,1.e0)
        Anuc = Z+Nin                                             # {nucleons within a nucleus}
        A = Anuc
        Rws = (A*.75/pi/bard)**.333333

        if(Rws<Rn):
            print ('Too large Rn for outer envelope! oyaform has returned zero')
            return 0

        aa = (A/bard)**.333333                                   # {cube size}
        
    elif(index == 3):                                            # {spheres after drip}
        
        g = bard*100.
        f = numpy.log(g)
        Rws = 31.68-8.400*f-0.2380*f*f+0.1152*f**3
        tn = 1./(0.2027+0.004506*g)                              # param.of shape
        Rn = 9.406+1.481*f+0.4625*f*f+0.05738*f**3
        dn_n = (9.761-1.322*f-0.5544*f*f-0.07624*f**3)/100.      # n-height
        Nin = pi/.75*Rn**3*dn_n*soyam(tn,numpy.min([1.e0,Rws/Rn]))
        tp = 1./(0.1558+2.225e-3*g+9.452e-4*g*g)
        Rp = 8.345+0.7767*f+0.1333*f*f+0.008707*f**3
        np_in = (4.040-1.097*f-0.0723*f*f+0.0225*f**3)/100.
        Z = pi/.75*Rp**3*np_in*soyam(tp,numpy.min([1.e0,Rws/Rp]))
        Nfree = bard*pi/.75*Rws**3-Z-Nin                         # free neutrons outside+under nuc.
        nn_out = Nfree/(pi/.75*Rws**3)                           # free neutron density
        nn_in = nn_out+dn_n                                      # max.n-density
        A = Z+Nfree+Nin                                          # total num.of barions in the cell (A')
        Anuc = Z+Nin+Nfree*(Rn/Rws)**3                           # number of barions within Rn
        
        if(Rn > Rws):
            Anuc = A
            
        aa = (A/bard)**.333333 # {cube size}

    else:
        return

    Rp0eff = (Z/pi*.75/np_in)**.333333
    Rp2eff = Rp* numpy.sqrt((1.-15./(5.+tp)+15./(5.+2.*tp)-5./(5.+3.*tp))/soyam(tp,1.e0))
    Rp1eff = Rp*(1.-12./(4.+tp)+12./(4.+2.*tp)-4./(4.+3.*tp))/ soyam(tp,1.e0)
    Rp3eff = Rp*((1.-18./(6.+tp)+18./(6.+2.*tp)-6./(6.+3.*tp))/soyam(tp,1.e0))**.333333
    xnuc = Rp2eff/Rws
    xnuct = xnuc*tp/(.6+tp)

    return Z, Anuc, A, xnuc, xnuct

def coulogmain(rho,T):

    # The demonstrative MAIN program
    # Input: Zion ("Z") - ion charge number,
    # CMI ("A") - atomic mass number,
    # CMI1 ("A'") - number of nucleons per Wigner-Seitz cell,
    # rho - density in g/cm^3
    # xnuc0, tp - nuclear form parameters
    # Output: CLel, CLth
    
    # AUM - atomic mass unit divided by the electron mass
    # AUD - relativistic unit of density in g/cm^3
    
    Uminus2 = 12.973
    AUM = 1822.8
    AUD = 15819.7
    UNISIG = 7.763e20
    UNIKAP = 2.778e15
    BOHR = 137.036
    TEMP = T/5.93e9
    
    if(rho > rho_nd):
        index = 3
    else:
        index = 30

    bard = rho/1.66054e15                                        # number density of barions in fm^{-3}
    
    Zion,CMI,CMI1,xnuc,xnuct = oyaform(bard,index)

    DENSI = rho/(AUD*AUM*CMI1)                                   # number density of ions (rel.)
    SPHERION = (.75/pi/DENSI)**.3333333                          # Ion sphere radius
    DENS = DENSI*Zion                                            # number density of electrons
    PF0 = (3.*pi**2*DENS)**.3333333                              # Fermi momentum ( = x_r)
    VF0 = PF0/ numpy.sqrt(1.+PF0**2)                             # Fermi velocity
    gamma = Zion**2/BOHR/TEMP/SPHERION                           # Ion coupling parameter

    CLel, CLth = coulog2(PF0, gamma, Zion, CMI, xnuc, xnuct)
    
    return CLth

# electron thermal conductivity in the crust
def _kappa_crust(T, rho, ne, Z):
 
    # T = temperature in K
    # rho = mass density in g/cm^3
    # ne = electron number density in fm^-3
    # Z = number of protons in nuclei
    # A = number of nucleons in nuclei (Anuc)
    
    x = routines.pF(ne)/Me/c                                     # parameter of electron relativism
    Me_eff = Me*numpy.sqrt(1.+x*x)                               # effective electron mass

    try:
        Lth = coulogmain(rho,T)
    except ValueError or ZeroDivisionError:
        return 0.

    nu = 4.*Z/(3.*pi*h)*routines.sqr(alpha)*Me_eff*c*c*Lth

    if (nu == 0.):
        return 1

    kappa = pi*pi*kB*kB*T*(ne/fm3)/(3.*Me_eff*nu)

    return kappa 


# Core thermal conductivity (electrons)
def _kappa_core_e (T, ne, nm, nn, Tc_p):                         # thermal conductivity in the core.

    # Electron and muon contribution
    # Transverse screening is included
    # T = temperature in K
    # ne = electron number density in fm^-3
    # nm = muon number density in fm^-3
    # Tc_p = critical temperature of proton sf
    
    T8 = T*1.e-8                                                 # normalized temperature

    np = ne+nm
    pp = fm*routines.pF(np)/h
    pe = fm*routines.pF(ne)/h
    pm = fm*routines.pF(nm)/h                                    # Fermi momenta in fm^-1

    # effective masses in Mp0

    me_eff = numpy.sqrt(0.295e-6+0.04415*pe*pe)
    mm_eff = numpy.sqrt(0.01267+0.04415*pm*pm)

    if(var_meff == 1):
        mp_eff =  routines._Meff(np)
    elif(var_meff == 2):
        nbb = nn+np
        mp_eff = 0.9937/(1+2.217*nbb+0.8850*nbb*nbb)
    else:
        mp_eff = Mp/Mp0

    ql = 0.21032*numpy.sqrt(pe*me_eff+pm*mm_eff+pp*mp_eff)       # screening momenta
    qt = 0.09640*numpy.sqrt(pe*pe+pm*pm+pp*pp)
                                      
    if(SUPERFLUIDITY and T<Tc_p):                                # Superfluid suppression factors

        y = sf_gap._SingletGap(T/Tc_p)
        r = (pe*pe+pm*pm)/pp/pp
        Rp = ((0.7694 + numpy.sqrt(routines.sqr(0.2306)+routines.sqr(0.07207*y))
              + (27.0*y*y+0.1476*numpy.power(y,4.))*numpy.exp(-numpy.sqrt(routines.sqr(4.273)+y*y))
              + 0.5051*(numpy.exp(4.273-numpy.sqrt(routines.sqr(4.273)+y*y))-1.))*
              numpy.exp(1.187-numpy.sqrt(routines.sqr(1.187)+y*y)))
        Rt = ((0.6-0.16*r)*numpy.exp(-0.14*y*y)+(1-0.6+0.16*r)/
             numpy.sqrt(1+y*y*routines.sqr(1.37*(1-0.6+0.16*r)/r)))
        Rtl = numpy.power(r+1.,1./3.)/numpy.power(routines.sqr(r+1.)-0.7568*y+routines.sqr(0.5061*y),1./6.)

    else:

        Rp = 1.
        Rt = 1.
        Rtl = 1.

    # Collision frequencies and relaxation times
    # Effective masses excluded in denominators 

    nuet = 0.14675 *T8*pe*Rt                                     # transverse collision frequencies

    nueec = 0.2950e-3*numpy.power(T8,5/3)*me_eff*me_eff*pe/ql/ql*numpy.power(qt,-2/3)*Rtl
    # crossed ee

    nueel = 0.2591e-3 *T8*T8*me_eff*me_eff*me_eff*me_eff/pe/ql/ql/ql                      # long ee
    nuepl = 0.2591e-3 *T8*T8*me_eff*me_eff*mp_eff*mp_eff/pe/ql/ql/ql*Rp                   # long ep

    if(nm != 0):

        nueml = 0.2591e-3 *T8*T8*me_eff*me_eff*mm_eff*mm_eff/pe/ql/ql/ql
        # long emu

        numt = 0.14675 *T8*pm*Rt                                                        # transv mu
        nummc = 0.2950e-3*numpy.power(T8,5/3)*mm_eff*mm_eff*pm/ql/ql*numpy.power(qt,-2/3)*Rtl
        # crossed mumu
        numml = 0.2591e-3 *T8*T8*mm_eff*mm_eff*mm_eff*mm_eff/pm/ql/ql/ql
        # long mumu
        numpl = 0.2591e-3 *T8*T8*mm_eff*mm_eff*mp_eff*mp_eff/pm/ql/ql/ql*Rp
        # long mup
        numel = 0.2591e-3 *T8*T8*me_eff*me_eff*mm_eff*mm_eff/pm/ql/ql/ql
        # long mue

        numec = 0.2950e-3*numpy.power(T8,5/3)*me_eff*mm_eff*pe*pe/pm/ql/ql*numpy.power(qt,-2/3)*Rtl
        # crossed me'
        nuemc = 0.2950e-3*numpy.power(T8,5/3)*mm_eff*me_eff*pm*pm/pe/ql/ql*numpy.power(qt,-2/3)*Rtl
        # crossed em'

        nue = nuet+nueec+nueel+nuepl+nueml                       # total collision frequency
        num = numt+nummc+numml+numpl+numel

        taue = (num-nuemc)/(num*nue-numec*nuemc)
        taum = (nue-numec)/(num*nue-numec*nuemc)

    else:

        nue = nuet+nueec+nueel+nuepl
        taue = 1./nue
        taum = 0.

    kappa = 0.3747e23*T8*(ne*taue+nm*taum)

    return kappa

# Core thermal conductivity (neutrons)
def _kappa_core_n(T, nn, np):                                    # thermal conductivity in the core

    # Neutron contribution:
    # Baiko, Haensel & Yakovlev, unpublished
    # Superfluidity effects are not taken into account
    # T = temperature in K
    # nn = neutron number density in fm^-3
    # np = proton number density in fm^-3

    corr = 1.2                                                   # correction to the relaxation time approximation

    T8 = T*1.e-8

    kFn = routines.pF(nn)/h*fm
    kFp = routines.pF(np)/h*fm

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

    # "in-vacuum" integrals
    try:
        Sn0 = (7.88/routines.sqr(kFn)*(1.-0.2241*kFn+0.2006*routines.sqr(kFn))/(1.-0.1742*kFn))
        Sp1 = (0.8007*kFp/routines.sqr(kFn)*(1.+31.28*kFp-0.0004285*routines.sqr(kFp)+26.85*kFn
                                   +0.08012*routines.sqr(kFn))/(1.-0.5898*kFn+0.2368*routines.sqr(kFn)+0.5838*routines.sqr(kFp)
                                                       +0.884*kFn*kFp))
        Sp2 = (0.383*numpy.power(kFp,4.)/numpy.power(kFn,5.5)*(1.+102.*kFp+53.91*kFn)/
        (1.-0.7087*kFn+0.2537*routines.sqr(kFn)+9.404*routines.sqr(kFp)-1.589*kFn*kFp))

        # "in-medium" corrections
        u = kFn-1.556
        Kn = (1./routines.sqr(mn_eff)*(0.4891+1.111*u*u-0.2283*u*u*u+0.01589*kFp
                             -0.02099*routines.sqr(kFp)+0.2773*u*kFp))
        u = kFn-2.126
        Kp1 = (1./routines.sqr(mp_eff)*(0.04377+1.1*u*u+0.118*u*u*u+0.1626*kFp
                              +0.3871*u*kFp-0.299*numpy.power(u,4.)))
        u = kFn - 2.116
        Kp2 = (1./routines.sqr(mp_eff)*(0.0001313+1.248*u*u+0.2403*u*u*u+0.3257*kFp
                              +0.5536*u*kFp-0.3237*numpy.power(u,4.)+0.09786*u*u*kFp))

        nu_nn = (3.48e15*numpy.power(mn_eff,3.)*routines.sqr(T8)*Sn0*Kn)
        nu_np = (3.48e15*(mn_eff)*routines.sqr(mp_eff)*routines.sqr(T8)*(Sp1*Kp1+Sp2*Kp2))
        tau = corr/(nu_nn+nu_np)

        kappa = pi*pi*kB*kB*T*(nn/fm3)/(3.*mn_eff*Mn0)*tau

        return kappa

    except ZeroDivisionError:
        return 0.

# Thermal conductivity (total per unit volume)
def _kappa(T, rho, ne, nm, nn, Z, Tc_p):

    # T = temperature in K
    # rho = mass density in g/cm^3
    # ne = electron number density in fm^-3
    # nm = muon number density in fm^-3
    # nn = neutron number density in fm^-3
    # Z = number of protons in a nucleus
    # A = number of nucleons per nucleus
    # Tc_nt = crit. temp. of neutron triplet sf
    # Tc_ns = crit. temp. of neutron singlet sf
    # Tc_p = critical temperature of proton sf

    np = ne+nm

    if(rho > CoreCrustBound):                                    # CORE

        kappa = _kappa_core_e(T, ne, nm, nn, Tc_p) + _kappa_core_n(T, nn, np)

    else:                                                        # CRUST

        kappa = _kappa_crust(T, rho, ne, Z)

    return kappa 


#------------------------------------



def coulogmain_for_grad(rho,T):

    if(rho > rho_nd):
        index = 3
    else:
        index = 30

    bard = rho/1.66054e15

    Zion,CMI,CMI1,xnuc,xnuct = oyaform(bard,index)


    return Zion,CMI
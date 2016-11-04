
# DECIMAL LOG of the Mean flux ratio (by the Simpson's algorithm)
#                                                       Version 21.08.02
#
# Modified version; August 2002
import numpy
from control.constants import *
from control.manager import *
from data import loaddata

def tite_init():

    Ms=loaddata.star_model()[-1, 0]*MSun                         # Neutron Star Mass in g
    Rs=loaddata.star_model()[-1,1]*1.e5                          # Neutron Star Radius in cm
    rgr=2.953*Ms/MSun/Rs*1.0e5                                   # r_g/R

    titspalex(rgr,MagField ,loaddata.g_surface(),dM_accreted)

    print('TiTe file is created.\n')

def titspalex(rgr,Bpole,g14,DeltaM_M):

    MAXTb = 36
    NM0 = 8
    EPS = 1.e-6

    tite = open('data/tite.dat', 'w')
    B12=Bpole/1.0e12
    tite.write('It works\n')
    tite.write('%6.3f %7.5f\n' %(rgr, B12))
    tite.write('------------------\n')

    for IT in range(1,MAXTb+1):

        Tb9=10.**(4.9+0.2*(IT-1)-9.0)

        # amu = cos polar angle on the surface
        # g14 is surface gravity in units of 10^14 cm/s^2

        eta=g14*g14*DeltaM_M
        Ti=Tb9*1.e9

        fr=-3.0/rgr**3*(numpy.log(1.0-rgr)+rgr+rgr**2/2.0) # GR corr of B_r
        pr=numpy.sqrt(1.0-rgr)*(3.0/(1.0-rgr)-2.0*fr)      # GR corr of B_t

        NM=NM0

        while True:

            NM=NM*2
            H=1.e0/NM
            S=0.0    # future integral

            NBLINK=2
            NSIGN=2

            for I in range(0,NM+1):
                amu=H*I # "mu" (=cos[polar angle on the surface])

                B12=Bpole/1.0e12*numpy.sqrt(amu**2+(pr/2./fr)**2*(1.-amu**2))# local B

                zeta=Tb9-numpy.sqrt(7.*Tb9*numpy.sqrt(g14))/1.0e3
                Ts6Fe4=g14*((7.*zeta)**2.25+(zeta/3.)**1.25)
                Ts6a4=(g14*(18.1*Tb9)**2.42*(0.447+0.075*numpy.log10(Ti)/(1.0+(6.2*Tb9)**4))+3.2*Tb9**1.67*Ts6Fe4)/(1+3.2*Tb9**1.67)

                R4Fe=1.0  # /* magnetic field correction */
                R4a =1.0

                if (B12>1.0e-4):

                    frt=2.0*fr*amu
                    cos2=frt**2/(frt**2+pr**2*(1.0-amu**2)) # angle btw B and normal
                    sin2=1.0-cos2

                    #  For iron envelope

                    a1=1.76e-4
                    a2=0.038
                    a3=1.5
                    a4=0.0132
                    a5=0.620
                    a6=0.318
                    a7=2.3e-9
                    a8=3.0
                    a9=0.160
                    a10=21.0
                    a11=4.7e5
                    b1=159.0
                    b2=270.0
                    b3=172.0
                    b4=110.0
                    b5=0.363
                    b6=0.181
                    b7=0.50
                    b8=0.619

                    chi1=(a1+a2*Tb9**a3)/(Tb9**2+a4*Tb9**a5)
                    chi1=1.0+chi1*B12**a6/(1.+a7*B12/Tb9**a8)**a9
                    cc=1.+1.0/(3.7+(a10+a11/B12**1.5)*Tb9**2)
                    chi1=chi1/cc

                    be=1.0/(1.0+b5*Tb9**b6)
                    chi2=numpy.sqrt(1.0+b1*B12/(1.0+b2*Tb9**b7))
                    chi2=chi2/(1.0+b3*B12/(1.0+b4*Tb9**b8))**be

                    al=4.0+numpy.sqrt(chi2/chi1)
                    chi=(chi1**al*cos2+chi2**al*sin2)**(1.0/al)
                    R4Fe=chi**4

                    #  Accreted envelope

                    a1=4.50e-3
                    a2=0.055
                    a3=2.0
                    a4=0.0595
                    a5=0.328
                    a6=0.237
                    a7=6.8e-7
                    a8=2.0
                    a9=0.113
                    a10=163.0
                    a11=3.4e5
                    b1=172.0
                    b2=155.0
                    b3=383.0
                    b4=94.0
                    b5=0.383
                    b6=0.367
                    b7=2.28
                    b8=1.69

                    chi1=(a1+a2*Tb9**a3)/(Tb9**2+a4*Tb9**a5)
                    chi1=1.+chi1*B12**a6/(1.0+a7*B12/Tb9**a8)**a9
                    cc=1.+1.0/(3.7+(a10+a11/B12**1.5)*Tb9**2)
                    chi1=chi1/cc

                    be=1.0/(1.0+b5*Tb9**b6)
                    chi2=numpy.sqrt(1.0+b1*B12/(1.0+b2*Tb9**b7))
                    chi2=chi2/(1.0+b3*B12/(1.0+b4*Tb9**b8))**be

                    al=(2.0+chi2/chi1)**2
                    chi=(chi1**al*cos2+chi2**al*sin2)**(1.0/al)
                    R4a=chi**4

                if(DeltaM_M<1.e-30):
                    Ts6_4=Ts6Fe4*R4Fe  # nonaccreted envelope
                elif(eta>1.0e-6):
                    Ts6_4=Ts6a4*R4a    # fully accreted envelope
                else:
                    ak=-numpy.log10(1.0e6*eta)
                    ga=1./(1.+3.8*(0.1*ak)**9)/(1.+0.171*ak**3.5*Tb9)

                    Ts6_4=Ts6a4*R4a*ga+Ts6Fe4*R4Fe*(1.0-ga)

                NSTAT=NBLINK   # Simpson statweight 1 4 2 4 2 4 1

                if (I==0 or I==NM):
                    NSTAT=1

                S=S+NSTAT*Ts6_4

                NBLINK=NBLINK+NSIGN
                NSIGN=-NSIGN

            S=S*H/3.
            SM=S     #for comparison

            if(numpy.abs(S-SM)<=EPS*S): # end of integration
                break

        Ts6_4=S  # surface-average flux
        Te=1.e6*Ts6_4**0.25  #/* surface temp. in K, nonredshifted */
        F1=numpy.log10(Te) # LOG of the magnetic Te

        tite.write('%6.3f %7.5f \n' %(numpy.log10(Tb9)+9.0, F1))



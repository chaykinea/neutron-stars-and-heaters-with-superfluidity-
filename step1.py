__author__ = 'maryhallow'


import matplotlib.pylab as plt
from physics import heatcapacity
from physics import neutrino
from physics import thermalconductivity
from physics import physics
from data    import loaddata
from matplotlib.ticker import AutoMinorLocator
import numpy as npp


# report 1

def test_f():

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    T_0 = [1,0]
    lw  = [1.5,3]

    loaddata.superfluid_data_init()
    loaddata.star_model_data_init()

    T_Q = npp.array([1e7,5e7,2.3e8,1e9,5e9,1e10])
    T_kappa      = npp.array([1e7,1e8,1e9,1e10])

    rho = loaddata.star_model()[:, 3]
    log_rho = npp.log10(rho)
    Phi = loaddata.star_model()[:, 4]

    nb = loaddata.star_model()[:, 5]
    ne = loaddata.star_model()[:, 6]
    nm = loaddata.nm(npp.log(rho))

    np=ne+nm
    nn=nb-np

    C_table = npp.zeros((len(rho), len(T_kappa),2))
    Q_table = npp.zeros((len(rho),len(T_Q),2))
    kappa_table = npp.zeros((len(rho),len(T_kappa),2))

    for k in range(0,2):

        for i in range (0,len(rho)):

            nb_, nn_, ne_, Z, A, Anuc, Vion = physics._CrustData(rho[i], nb[i], nn[i],ne[i])

            Tc_nt = physics._Tcn(nn[i])*T_0[k]
            Tc_p  = physics._Tcp(rho[i],np[i])*T_0[k]
            Tc_ns = physics._Tcs(nn[i],rho[i])*T_0[k]

            for j in range(len(T_Q)):
                T_local = T_Q[j]*npp.exp(-Phi[i])
                Q_table[i,j,k] = neutrino._Q(T_local, nn_, ne_, nm[i], rho[i], Phi[i], Tc_nt, Tc_ns, Tc_p)*npp.exp(-2*Phi[i])

            for j in range(len(T_kappa)):
                T_local = T_kappa[j]*npp.exp(-Phi[i])
                C_table[i,j,k] = heatcapacity._C(T_local, rho[i], nn_, ne_, nm[i], Z, Anuc, Tc_nt, Tc_ns, Tc_ns)
                kappa_table[i,j,k] = thermalconductivity._kappa(T_local, rho[i], ne_, nm[i], nn_, Z, Tc_p)

    linestyles = npp.array([['r-','b-','g-','y-','m-','c-','k-'],
                            ['r--','b--','g--','y--','m--','c--','k--']])
    labels_Q = npp.array(['T = 10$^{7}$ K','T = 5 $\\times$ 10$^{7}$ K',
                          'T = 2.3 $\\times$ 10$^{8}$ K','T = 10$^{9}$ K','T = 5 $\\times$ 10$^{9}$ K',
                          'T = 10$^{10}$ K'])
    labels_kappa = npp.array(['T = 10$^{7}$ K','T = 10$^{8}$ K',
                          'T = 10$^{9}$ K','T = 10$^{10}$ K'])

    pf1, a1 = plt.subplots()

    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    a1.xaxis.set_minor_locator(x_minor_locator)
    a1.yaxis.set_minor_locator(y_minor_locator)

    plt.title('Neutrino emissivity',fontsize=22)
    plt.ylim(2,32)
    plt.xlim(9,npp.log10(7.2946929931640625e+014))
    plt.xlabel('log $\\rho$ [g cm$^{-3}$]',fontsize=22)
    plt.ylabel('log $\\rm Q$ [erg s$^{-1}$ cm$^{-3}$]',fontsize=22)
    plt.yticks([5,10,15,20,25,30],fontsize=22)
    plt.xticks([9,10,11,12,13,14],fontsize=22)
    plt.xlim(9,14.7)

    for i in range (len(T_Q)-1,-1,-1):

        plt.plot(log_rho,npp.log10(Q_table[:,i,0]),linestyles[0,i],lw=lw[0],label=labels_Q[i])
        plt.plot(log_rho,npp.log10(Q_table[:,i,1]),linestyles[1,i],lw=lw[1])


    plt.figure(1)
    plt.legend(loc='upper left',fontsize=16)
    plt.savefig('neutrino.pdf',format='pdf')


    pf2, a2 = plt.subplots()

    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    a2.xaxis.set_minor_locator(x_minor_locator)
    a2.yaxis.set_minor_locator(y_minor_locator)

    plt.title('Heat capacity',fontsize=22)
    plt.ylim(15,22)
    plt.xlim(9,npp.log10(7.2946929931640625e+014))
    plt.xlabel('log $\\rho$ [g cm$^{-3}$]',fontsize=22)
    plt.ylabel('log $\\rm C$ [erg K$^{-1}$ cm$^{-3}$]',fontsize=22)
    plt.yticks([15,16,17,18,19,20,21,22],fontsize=22)
    plt.xticks([9,10,11,12,13,14],fontsize=22)
    plt.xlim(9,14.7)

    for i in range (len(T_kappa)-1,-1,-1):

        plt.plot(log_rho,npp.log10(C_table[:,i,0]),linestyles[0,i],lw=lw[0],label=labels_kappa[i])
        plt.plot(log_rho,npp.log10(C_table[:,i,1]),linestyles[1,i],lw=lw[1])

    plt.figure(2)
    plt.legend(loc='upper left',fontsize=16)
    plt.savefig('heat_cap.pdf',format='pdf')

    pf3, a3 = plt.subplots()

    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    a3.xaxis.set_minor_locator(x_minor_locator)
    a3.yaxis.set_minor_locator(y_minor_locator)

    plt.title('Thermal conductivity',fontsize=22)
    plt.ylim(17,25)
    plt.xlim(9,npp.log10(7.2946929931640625e+014))
    plt.xlabel('log $\\rho$ [g cm$^{-3}$]',fontsize=22)
    plt.ylabel('log $\kappa$ [erg K$^{-1}$ cm$^{-1}$ s$^{-1}$]',fontsize=22)
    plt.yticks([17,18,19,20,21,22,23,24,25],fontsize=22)
    plt.xticks([9,10,11,12,13,14],fontsize=22)
    plt.xlim(9,14.7)

    for i in range (len(T_kappa)-1,-1,-1):

        plt.plot(log_rho,npp.log10(kappa_table[:,i,0]),linestyles[0,i],lw=lw[0],label=labels_kappa[i])
        plt.plot(log_rho,npp.log10(kappa_table[:,i,1]),linestyles[1,i],lw=lw[1])


    plt.figure(3)
    plt.legend(loc='upper left',fontsize=16)
    plt.savefig('kappa.pdf',format='pdf')

    plt.show()


def T_crit():

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    figg, axx = plt.subplots()
    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    axx.xaxis.set_minor_locator(x_minor_locator)
    axx.yaxis.set_minor_locator(y_minor_locator)
    plt.title('Critical temperature',fontsize=22)

    loaddata.superfluid_data_init()
    loaddata.star_model_data_init()

    rho = loaddata.star_model()[:, 3]
    log_rho = npp.log10(rho)

    nb = loaddata.star_model()[:, 5]
    ne = loaddata.star_model()[:, 6]
    nm = loaddata.nm(npp.log(rho))

    np=ne+nm
    nn=nb-np

    Tc_ns = npp.zeros_like(rho)
    Tc_nt = npp.zeros_like(rho)
    Tc_p  = npp.zeros_like(rho)

    for i in range(len(rho)):

        Tc_nt[i] = npp.log10(physics._Tcn(nn[i]))
        Tc_p[i]  = npp.log10(physics._Tcp(rho[i],np[i]))
        Tc_ns[i] = npp.log10(physics._Tcs(nn[i],rho[i]))

    plt.xlim(11,15.2)
    plt.ylim(7.4,10.2)

    plt.xlabel('log $\\rho$ [g cm$^-3$]',fontsize=22)
    plt.ylabel('log $T_{c}$ [K]',fontsize=22)
    plt.yticks([7.5,8,8.5,9,9.5,10],fontsize=22)
    plt.xticks([11,12,13,14,15],fontsize=22)

    plt.plot(log_rho,Tc_ns,lw=2,label='1ns')
    plt.plot(log_rho,Tc_nt,lw=2,label='1nt')
    plt.plot(log_rho,Tc_p ,lw=2,label='1p')

    plt.legend(loc='upper left',fontsize=20)
    #plt.savefig('Tcrit.pdf',format='pdf')

    plt.show()


def cooling_curve():

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    line_styles = ['b-','r--','y--','k--']
    dashes = npp.array([[14,4],[20,5],[5,9]])
    line_thickness = [4,3.5,3,2.5]

    cooling_1 = npp.loadtxt('output/cooling_curve.dat')
    cooling_2 = npp.loadtxt('output/cooling_curve_ps.dat')
    cooling_3 = npp.loadtxt('output/cooling_curve_ns_ps.dat')
    cooling_4 = npp.loadtxt('output/cooling_curve_ns_nt_ps.dat')

    fig, ax = plt.subplots()

    plt.xlim(1,6.7)
    plt.ylim(5.4,6.6)

    plt.title('Cooling curves',fontsize=22)
    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(4)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    ax.xaxis.set_minor_locator(x_minor_locator)
    ax.yaxis.set_minor_locator(y_minor_locator)

    plt.xlabel('log t [yr]',fontsize=22)
    plt.ylabel('log $\\rm T^{\infty}_{\\rm s}$ [K]',fontsize=22)
    plt.yticks([5.4,5.6,5.8,6,6.2,6.4,6.6],fontsize=22)
    plt.xticks([1,2,3,4,5,6],fontsize=22)

    plt.text(1.4,5.7, '$M \\thinspace = \\thinspace1.4\\thinspace  M_{\odot}$',fontsize=22)
    plt.text(1.4,5.8 , '$R \\thinspace = \\thinspace12.6\\thinspace  \mathrm{km}$',fontsize=22)
    plt.text(1.4,5.6 , '$\\rho_{\\rm b} = 10^{9} \\thinspace \\rm g \\thinspace  \\rm cm^{-3}$',fontsize=22)

    plt.plot(npp.log10(cooling_1[:,1]),npp.log10(cooling_1[:,0]),line_styles[0], linewidth=line_thickness[0],label='No SF')
    plt.plot(npp.log10(cooling_2[:,1]),npp.log10(cooling_2[:,0]),line_styles[1], linewidth=line_thickness[1], dashes = (dashes[0,0],dashes[0,1],4,dashes[0,1]),label='1p')
    plt.plot(npp.log10(cooling_3[:,1]),npp.log10(cooling_3[:,0]),line_styles[2], linewidth=line_thickness[2], dashes = (dashes[1,0],dashes[1,1]),label='1p + 1ns')
    plt.plot(npp.log10(cooling_4[:,1]),npp.log10(cooling_4[:,0]),line_styles[3], linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]),label='1p + 1ns + 1nt')

    plt.legend(loc='upper right',fontsize=18)
    plt.savefig('cooling_curves.pdf',format='pdf')
    plt.show()


def show_source_curves(a):

    placement = npp.array([[9,6,3,0],[18,15,12,0]])

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    line_styles = ['b-','r--','y--','k--']
    dashes = npp.array([[14,4],[20,5],[5,9]])
    line_thickness = [4,3,2,1]

    figgg, axxx = plt.subplots()

    plt.title('Constant source [cooling curves]',fontsize=22)
    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(4)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    axxx.xaxis.set_minor_locator(x_minor_locator)
    axxx.yaxis.set_minor_locator(y_minor_locator)

    plt.xlabel('log t [yr]',fontsize=22)
    plt.ylabel('log $\\rm T^{\infty}_{\\rm s}$ [K]',fontsize=22)
    plt.yticks([6.1,6.2,6.3,6.4,6.5],fontsize=22)
    plt.xticks([-2,-1,0,1,2,3],fontsize=22)

    plt.xlim(-2,2.7)
    plt.ylim(6.1,6.5)

    for i in placement[a,:]:

        temp_1 = npp.loadtxt('../data/1/cooling curves_profiles/no_sf/file_' + str(i) + '_cooling.dat')
        temp_2 = npp.loadtxt('../data/1/cooling curves_profiles/1p/file_' + str(i) + '_cooling.dat')
        temp_3 = npp.loadtxt('../data/1/cooling curves_profiles/1p_1ns/file_' + str(i) + '_cooling.dat')

        temp_1[:,1] -= 1e3
        temp_2[:,1] -= 1e3
        temp_3[:,1] -= 1e3

        if i!=0:
            plt.plot(npp.log10(temp_1[:,1]), npp.log10(temp_1[:,0]), line_styles[0], linewidth=line_thickness[0])
            plt.plot(npp.log10(temp_2[:,1]), npp.log10(temp_2[:,0]), line_styles[1], linewidth=line_thickness[1], dashes = (dashes[0,0],dashes[0,1],4,dashes[0,1]))
            plt.plot(npp.log10(temp_3[:,1]), npp.log10(temp_3[:,0]), line_styles[2], linewidth=line_thickness[2], dashes = (dashes[1,0],dashes[1,1]))
        else:
            plt.plot(npp.log10(temp_1[:,1]), npp.log10(temp_1[:,0]), line_styles[0], linewidth=line_thickness[0], label='no SF')
            plt.plot(npp.log10(temp_2[:,1]), npp.log10(temp_2[:,0]), line_styles[1], linewidth=line_thickness[1], dashes = (dashes[0,0],dashes[0,1],4,dashes[0,1]), label='1p')
            plt.plot(npp.log10(temp_3[:,1]), npp.log10(temp_3[:,0]), line_styles[2], linewidth=line_thickness[2], dashes = (dashes[1,0],dashes[1,1]), label='1p + 1ns')

    if(a==0):
        plt.text(-0.92,6.203 , '3', fontsize=20,color='b')
        plt.text(-1.01,6.243 , '3', fontsize=20,color='y')
        plt.text(-0.63,6.38 , '3', fontsize=20,color='r')

        plt.text(-0.569,6.20 , '2', fontsize=20,color='b')
        plt.text(-0.56,6.23 , '2', fontsize=20,color='y')
        plt.text(-0.35,6.35 , '2', fontsize=20,color='r')

        plt.text(0.9,6.185 , '1', fontsize=20,color='b')
        plt.text(-0.1,6.24 , '1', fontsize=20,color='y')
        plt.text(1.5,6.3 , '1', fontsize=20,color='r')

        plt.text(1,6.127 , '0', fontsize=20,color='b')
        plt.text(0.2,6.21 , '0', fontsize=20,color='y')
        plt.text(1.55,6.27 , '0', fontsize=20,color='r')

        plt.text(-1.5,6.45 ,'$\\rho_{1} = 10^{11} \\thinspace \\rm g \\thinspace cm^{-3}$', fontsize=23)
    else:
        plt.text(-0.055,6.203 , '3', fontsize=20,color='b')
        plt.text(-0.13,6.243 , '3', fontsize=20,color='y')
        plt.text(0.27,6.37 , '3', fontsize=20,color='r')

        plt.text(0.35,6.20 , '2', fontsize=20,color='b')
        plt.text(0.35,6.245 , '2', fontsize=20,color='y')
        plt.text(0.8,6.35 , '2', fontsize=20,color='r')

        plt.text(1.1,6.145 , '1', fontsize=20,color='b')
        plt.text(0.75,6.235 , '1', fontsize=20,color='y')
        plt.text(1.25,6.275 , '1', fontsize=20,color='r')

        plt.text(1.25,6.126 , '0', fontsize=20,color='b')
        plt.text(0.9,6.21 , '0', fontsize=20,color='y')
        plt.text(1.55,6.248 , '0', fontsize=20,color='r')

        plt.text(-1.5,6.45 ,'$\\rho_{1} = 10^{12} \\thinspace \\rm g \\thinspace cm^{-3}$', fontsize=23)

    plt.legend(loc='upper right',fontsize=16)
    plt.savefig('source_curves_' + str(a) + '.pdf',format='pdf')


def show_source_profiles(a):

    placement = npp.array([[9,6,3,0],[18,15,12,0]])

    line_styles = ['b-','r--','y--','k--']
    dashes = npp.array([[14,4],[20,5],[5,9]])
    line_thickness = [4,3,2,2.5]

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    prf, ax1 = plt.subplots()

    plt.title('Constant source [temperature profiles]',fontsize=22)
    plt.ylim(8.2,9.3)
    plt.xlim(9,npp.log10(7.2946929931640625e+014))
    plt.xlabel('log $\\rho$ [g cm$^{-3}$]',fontsize=22)
    plt.ylabel('log $\\rm T^{\infty}$ [K]',fontsize=22)
    plt.yticks([8.2,8.4,8.6,8.8,9.0,9.2,9.4],fontsize=22)
    plt.xticks([9,10,11,12,13,14],fontsize=22)

    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(4)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    ax1.xaxis.set_minor_locator(x_minor_locator)
    ax1.yaxis.set_minor_locator(y_minor_locator)

    for i in placement[a,:]:

        temp_1 = npp.loadtxt('../data/1/cooling curves_profiles/no_sf/file_' + str(i) + '.dat')
        temp_2 = npp.loadtxt('../data/1/cooling curves_profiles/1p/file_' + str(i) + '.dat')
        temp_3 = npp.loadtxt('../data/1/cooling curves_profiles/1p_1ns/file_' + str(i) + '.dat')

        if i!=0:
            plt.plot(npp.log10(temp_1[:,0]), npp.log10(temp_1[:,-1]), line_styles[0], linewidth=line_thickness[0])
            plt.plot(npp.log10(temp_2[:,0]), npp.log10(temp_2[:,-1]), line_styles[1], linewidth=line_thickness[1],  dashes = (dashes[0,0],dashes[0,1],4,dashes[0,1]))
            plt.plot(npp.log10(temp_3[:,0]), npp.log10(temp_3[:,-1]), line_styles[2], linewidth=line_thickness[2],  dashes = (dashes[1,0],dashes[1,1]))

        else:
            plt.plot(npp.log10(temp_1[:,0]), npp.log10(temp_1[:,-1]), line_styles[0], linewidth=line_thickness[0], label='no SF')
            plt.plot(npp.log10(temp_2[:,0]), npp.log10(temp_2[:,-1]), line_styles[1], linewidth=line_thickness[1], dashes = (dashes[0,0],dashes[0,1],4,dashes[0,1]), label='1p')
            plt.plot(npp.log10(temp_3[:,0]), npp.log10(temp_3[:,-1]), line_styles[2], linewidth=line_thickness[2], dashes = (dashes[1,0],dashes[1,1]), label='1p + 1ns')

    if(a==0):
        plt.axvline(x=11,     linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')
        plt.axvline(x=12, linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')

        plt.text(13.15,8.573 , '0', fontsize=16,color='r')
        plt.text(12.2,8.617 , '1', fontsize=17,color='r')
        plt.text(13.6,8.65  , '2', fontsize=20,color='r')
        plt.text(13.4,8.8 , '3', fontsize=20,color='r')

        plt.text(10,8.53 , '0', fontsize=20,color='y')
        plt.text(10,8.61 , '1', fontsize=20,color='y')
        plt.text(11.65,8.78 , '2', fontsize=20,color='y')
        plt.text(12.17,8.85 , '3', fontsize=20,color='y')

        plt.text(13.8,8.27  , '0', fontsize=20,color='b')
        plt.text(13.5,8.37 , '1', fontsize=20,color='b')
        plt.text(12.6,8.775 , '2', fontsize=20,color='b')
        plt.text(12.47,8.94 , '3', fontsize=20,color='b')
    else:
        plt.axvline(x=12,     linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')
        plt.axvline(x=13.104, linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')

        plt.text(13.5,8.573 , '0', fontsize=15,color='r')
        plt.text(13.2,8.645 , '1', fontsize=20,color='r')
        plt.text(13.4,8.76  , '2', fontsize=20,color='r')
        plt.text(13.7,8.89 , '3', fontsize=20,color='r')

        plt.text(10,8.47 , '0', fontsize=20,color='y')
        plt.text(10,8.558 , '1', fontsize=13,color='y')
        plt.text(10,8.7 , '2', fontsize=20,color='y')
        plt.text(10,8.78 , '3', fontsize=20,color='y')

        plt.text(13.8,8.27  , '0', fontsize=20,color='b')
        plt.text(13.5,8.37 , '1', fontsize=20,color='b')
        plt.text(13.4,8.64 , '2', fontsize=20,color='b')
        plt.text(13.2,8.94 , '3', fontsize=20,color='b')


    plt.legend(loc='upper right',fontsize=16)
    plt.savefig('temperature_profiles_' + str(a) + '.pdf',format='pdf')


def show_source_neutrino_profiles(a):

    placement = npp.array([[9,6,3,0],[18,15,12,0]])

    source = npp.array([1e-19,1e-19,1,1,1,1e-19,1e-19])
    rho_source = npp.array([[10,10.999,11,11.5,12,12.001,13],[10,11.999,12,12.5,13.104,13.105,14]])

    H_c    = npp.array([0,
                        1e17,5e17,
                        1e18,5e18,
                        1e19,5e19,
                        1e20,5e20,
                        1e21])

    line_styles = ['b-','r--','y--','k--']
    dashes = npp.array([[14,4],[20,5],[5,9]])
    line_thickness = [4,3,2,2.5]

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    prf, ax1 = plt.subplots()

    plt.title('Constant source [neutrino profiles]',fontsize=22)
    plt.ylim(12,22)
    plt.xlim(9,npp.log10(7.2946929931640625e+014))
    plt.xlabel('log $\\rho$ [g cm$^{-3}$]',fontsize=22)
    plt.ylabel('log $\\rm Q_{\\nu}$ [erg cm$^{-3}$ s$^{-1}$]',fontsize=22)
    plt.yticks([12,14,16,18,20,22],fontsize=22)
    plt.xticks([9,10,11,12,13,14],fontsize=22)

    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(4)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    ax1.xaxis.set_minor_locator(x_minor_locator)
    ax1.yaxis.set_minor_locator(y_minor_locator)

    for i in placement[a,:]:

        temp_1 = npp.loadtxt('../data/1/cooling curves_profiles/no_sf/file_' + str(i) + '_neutrino.dat')
        temp_2 = npp.loadtxt('../data/1/cooling curves_profiles/1p/file_' + str(i) + '_neutrino.dat')
        temp_3 = npp.loadtxt('../data/1/cooling curves_profiles/1p_1ns/file_' + str(i) + '_neutrino.dat')

        if i!=0:
            plt.plot(rho_source[a,:],npp.log10(H_c[i-a*9]*source),line_styles[3], linewidth=line_thickness[3],  dashes = (dashes[2,0],dashes[2,1]))

            plt.plot(npp.log10(temp_1[:,0]), npp.log10(temp_1[:,-1]), line_styles[0], linewidth=line_thickness[0])
            plt.plot(npp.log10(temp_2[:,0]), npp.log10(temp_2[:,-1]), line_styles[1], linewidth=line_thickness[1],  dashes = (dashes[0,0],dashes[0,1],4,dashes[0,1]))
            plt.plot(npp.log10(temp_3[:,0]), npp.log10(temp_3[:,-1]), line_styles[2], linewidth=line_thickness[2],  dashes = (dashes[1,0],dashes[1,1]))
        else:
            plt.plot(npp.log10(temp_1[:,0]), npp.log10(temp_1[:,-1]), line_styles[0], linewidth=line_thickness[0], label='no SF')
            plt.plot(npp.log10(temp_2[:,0]), npp.log10(temp_2[:,-1]), line_styles[1], linewidth=line_thickness[1], dashes = (dashes[0,0],dashes[0,1],4,dashes[0,1]), label='1p')
            plt.plot(npp.log10(temp_3[:,0]), npp.log10(temp_3[:,-1]), line_styles[2], linewidth=line_thickness[2], dashes = (dashes[1,0],dashes[1,1]), label='1p + 1ns')


    if(a==0):

        plt.text(12.17,15.7 , '0', fontsize=20,color='r')
        plt.text(12.4,16.8 , '1', fontsize=20,color='r')
        plt.text(12.17,18.1 , '2', fontsize=20,color='r')
        plt.text(12.75,18.95 , '3', fontsize=20,color='r')

        plt.text(10.2,14.65 , '0', fontsize=20,color='y')
        plt.text(10.5,15.45 , '1', fontsize=20,color='y')
        plt.text(11.2,18.8 , '2', fontsize=20,color='y')
        plt.text(11.55,20.4 , '3', fontsize=20,color='y')

        plt.text(13,14 , '0', fontsize=20,color='b')
        plt.text(12.8,14.5 , '1', fontsize=20,color='b')
        plt.text(13.6,15.5 , '2', fontsize=20,color='b')
        plt.text(13.55,16.73 , '3', fontsize=20,color='b')

        plt.text(11.4,18.1 , '1', fontsize=20,color='k')
        plt.text(11.4,19.8 , '2', fontsize=20,color='k')
        plt.text(11.4,21.1 , '3', fontsize=20,color='k')

    else:

        plt.text(12.17,15.7 , '0', fontsize=20,color='r')
        plt.text(12.4,16.75 , '1', fontsize=20,color='r')
        plt.text(12.92,19.1 , '2', fontsize=20,color='r')
        plt.text(13.4,19.4 , '3', fontsize=20,color='r')

        plt.text(9.7,14.1 , '0', fontsize=20,color='y')
        plt.text(9.7,14.9 , '1', fontsize=18,color='y')
        plt.text(9.7,16.2 , '2', fontsize=20,color='y')
        plt.text(9.7,16.85 , '3', fontsize=20,color='y')

        plt.text(12.9,14 , '0', fontsize=20,color='b')
        plt.text(12.8,14.5 , '1', fontsize=20,color='b')
        plt.text(13.55,17.3 , '2', fontsize=20,color='b')
        plt.text(13.8,17.62 , '3', fontsize=20,color='b')

        plt.text(12.4,18.1 , '1', fontsize=20,color='k')
        plt.text(12.4,19.8 , '2', fontsize=20,color='k')
        plt.text(12.4,21.1 , '3', fontsize=20,color='k')


    plt.legend(loc='upper left',fontsize=16)
    plt.savefig('neutrino_profiles_' + str(a) + '.pdf',format='pdf')


#test_f()
T_crit()
#cooling_curve()

#show_source_curves(1)
#show_source_profiles(0)
#show_source_neutrino_profiles(0)
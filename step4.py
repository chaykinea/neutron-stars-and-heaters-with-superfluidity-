from control.constants import *
from control.manager   import *
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import matplotlib.cm as cm


neutron_s_model_names = np.array(['AWP2', 'AWP3', 'CCDK', 'CLS', 'GIPSF', 'MSH', 'SCLBL', 'SFB', 'WAP'])
neutron_s_model_names2 = np.array(['AWP2', 'GIPSF', 'SCLBL'])
proton_s_model_names  = np.array(['AO', 'BCLL', 'BS', 'CCDK', 'CCYms', 'CCYps', 'EEHO', 'EEHOr', 'T'])
neutron_t_model_names = np.array(['AO', 'BEEHS', 'EEHO', 'EEHOr', 'SYHHP', 'T', 'TTav', 'TToa'])

dashes = np.array([[12,4],[22,4],[7,4], [10,10], [12,4],[5,12], [2,2]])
line_thickness = np.array([4.5, 4.1, 3.3, 1.7, 2.9, 2.5, 2.1, 3.7, 2])
#colors = np.array(['r', 'purple', 'b', 'm', 'c', 'chocolate', 'darkgreen', 'orange', 'k'])
colors = np.array(['r', 'purple', 'b', 'm', 'c', 'chocolate', 'darkgreen', 'orange', 'k'])
colors = np.array(['red', 'green','darkblue'])

source_rho_1 = np.array([1e11, 1e12, 1e13])
source_rho_2 = np.array([1.000022e12, 1.275997e13, 5.063467e13])

# model 1
ns = np.array(['AWP2', 'AWP3', 'CCDK', 'CLS', 'GIPSF', 'MSH', 'SCLBL', 'SFB', 'WAP'])
nt = np.array(['AO'])
ps = np.array(['CCDK'])
time_roots = np.array([817.0204, 972.3984, 991.4243, 969.2274, 978.7404, 962.8855, 801.1655, 928.0047, 972.3984])
time_root_no_sf = 19506.77

# time_roots = np.array([244186.4, 259977.9,  267302.8,  261785.3, 263624.5, 259977.9,  245898.7, 254587.2, 259977.9])
# time_root_no_sf = 325.3953


# model 2
# ns = np.array(['GIPSF'])
# nt = np.array(['AO', 'BEEHS', 'EEHO', 'EEHOr', 'SYHHP', 'T', 'TTav', 'TToa'])
# ps = np.array(['CCDK'])
# time_roots = np.array([978.7404, 18444.49, 3747.005, 288548.5, 156477.2, 1080.212, 1926.863, 1273.641])
# time_root_no_sf = 19506.77

# time_roots = np.array([263624.5, 172458.8, 477.6023, 404.6698, 633233.8, 275.1297, 122674.5, 227792.4])
# time_root_no_sf = 325.3953

# model 3
# ns = np.array(['GIPSF'])
# nt = np.array(['AO'])
# ps = np.array(['AO', 'BCLL', 'BS', 'CCDK', 'CCYms', 'CCYps', 'EEHO', 'EEHOr', 'T'])
# time_roots = np.array([880.44, 896.2949, 886.782, 978.7404, 908.9788, 889.953, 928.0047,  921.6628, 893.1239])
# time_root_no_sf = 19506.77

# time_root_no_sf = 325.3953
# time_roots = np.array([378477.3, 378477.3, 378477.3, 263624.5, 378477.3, 378477.3, 378477.3, 378477.3, 378477.3])



models = { 'neutron_singlet': {'AWP2' : [  28, 0.20,   1.5,  1.7,   2.5],
                               'AWP3' : [  50, 0.20,   2.0,  1.4,   2.0],
                               'CCDK' : [ 127, 0.18,   4.5, 1.08,   1.1],
                               'CLS'  : [ 2.2, 0.18,  0.06,  1.3,  0.03],
                               'GIPSF': [ 8.8, 0.18,   0.1,  1.2,   0.6],
                               'MSH'  : [2.45, 0.18,  0.05,  1.4,   0.1],
                               'SCLBL': [ 4.1, 0.35,   1.7, 1.67,  0.06],
                               'SFB'  : [  45, 0.10,   4.5, 1.55,   2.5],
                               'WAP'  : [  69, 0.15,   3.0,  1.4,   3.0]},
           'proton_singlet' : {'AO'   : [  14, 0.15,  0.22, 1.05,   3.8],
                               'BCLL' : [1.69, 0.05,  0.07, 1.05,  0.16],
                               'BS'   : [  17,  0.0,   2.9,  0.8,  0.08],
                               'CCDK' : [ 102,  0.0,   9.0,  1.3,   1.5],
                               'CCYms': [  35,  0.0,   5.0,  1.1,   0.5],
                               'CCYps': [  34,  0.0,   5.0, 0.95,   0.3],
                               'EEHO' : [ 4.5,  0.0,  0.57,  1.2,  0.35],
                               'EEHOr': [  61,  0.0,   6.0,  1.1,   0.6],
                               'T'    : [  48, 0.15,   2.1,  1.2,   2.8]},
           'neutron_triplet':  {'AO'  : [ 4.0,  1.2,  0.45,  3.3,   5.0],
                               'BEEHS': [0.45,  1.0,  0.40,  3.2,  0.25],
                               'EEHO' : [0.48, 1.28,   0.1, 2.37,  0.02],
                               'EEHOr': [0.23,  1.2, 0.026,  1.6, 0.008],
                               'SYHHP': [ 1.0, 2.08,  0.04,  2.7, 0.013],
                               'T'    : [ 1.2, 1.55,  0.05, 2.35,  0.07],
                               'TTav' : [ 3.0,  1.1,  0.60, 2.92,   3.0],
                               'TToa' : [ 2.1,  1.1,  0.60,  3.2,   2.4]}}


def _Tcn(kf):

    if ( kf<= coefftnk0 or kf>= coefftnk2 ):
        Tc=0.0
    else:
        a=(kf-coefftnk0)*(kf-coefftnk0)
        b=(coefftnk2-kf)*(coefftnk2-kf)
        Tc=coefftnT0*a/(a+coefftnk1)*b/(b+coefftnk3)

    return Tc


# critical temperature in K of proton singlet superfluidity in the core
def _Tcp(kf):

    if (coeffspk2 < 0.0):
        Tc=coeffspT0
    else:
        if kf<= coeffspk0 or kf>= coeffspk2:
            Tc=0.0
        else:
            a=(kf-coeffspk0)*(kf-coeffspk0)
            b=(coeffspk2-kf)*(coeffspk2-kf)
            Tc=coeffspT0*a/(a+coeffspk1)*b/(b+coeffspk3)

    return Tc


def _Tcs(kf):

    if (coeffsnk2 < 0.0):
        Tcs = coeffsnT0
    else:
        if (kf <= coeffsnk0 or kf >= coeffsnk2 ):
            Tcs = 0.0
        else:
            a = (kf-coeffsnk0)*(kf-coeffsnk0)
            b = (coeffsnk2-kf)*(coeffsnk2-kf)
            Tcs = coeffsnT0*a/(a+coeffsnk1)*b/(b+coeffsnk3)

    return Tcs

def pF(nd):

    if nd > 0.:
        p=h*numpy.power( 3.*pi*pi*nd, 1./3. )/fm
    else:
        p=0.
    return p


def plot_style(xticks=5,yticks=5):

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    fig,ax = plt.subplots()
    x_minor_locator = AutoMinorLocator(xticks)
    y_minor_locator = AutoMinorLocator(yticks)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    ax.xaxis.set_minor_locator(x_minor_locator)
    ax.yaxis.set_minor_locator(y_minor_locator)


def T_crit():

    global coeffsnk0, coeffsnk1, coeffsnk2, coeffsnk3, coeffsnT0, \
        coefftnk0, coefftnk1, coefftnk2, coefftnk3, coefftnT0,\
        coeffspk0, coeffspk1, coeffspk2, coeffspk3, coeffspT0

    k = np.linspace(0,3.7,400)

    Tc_ns = np.zeros_like(k)
    Tc_nt = np.zeros_like(k)
    Tc_p  = np.zeros_like(k)

    plot_style()
    plt.title('$\\rm T_{c} \\thinspace neutron \\thinspace  ^{1}S_{0}$', fontsize=22)

    idx = 0
    for model in neutron_s_model_names2:

        coeffsnT0 = models['neutron_singlet'][model][0] * 0.5669 * MeV_erg / kB * 1e-9
        coeffsnk0 = models['neutron_singlet'][model][1]
        coeffsnk1 = models['neutron_singlet'][model][2]
        coeffsnk2 = models['neutron_singlet'][model][3]
        coeffsnk3 = models['neutron_singlet'][model][4]

        print(colors[idx])
        for i in range(len(k)):

            Tc_ns[i] = _Tcs(k[i])

        if idx>2:
            plt.plot(k,Tc_ns,color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx-2,0],dashes[idx-2,1]),label=model)
        else:
            plt.plot(k,Tc_ns,color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1],4,dashes[idx,1]),label=model)
        idx += 1

    plt.legend(loc='upper right',fontsize=15)

    plt.xlabel('$\\rm k_{F} \\thinspace (fm^{-1})$',fontsize=22)
    plt.ylabel('$\\rm T_{c} \\thinspace (10^{9} K)$',fontsize=22)
    plt.yticks([0,2,4,6,8,10,12,14],fontsize=22)
    plt.xticks([0,0.5,1,1.5,2,2.5,3,3.5],fontsize=22)
    plt.xlim(0,1.7)
    plt.ylim(0,14)
    plt.savefig('Tcrit_ns.pdf',format='pdf')
    '''
    plot_style()
    idx = 0
    plt.title('$\\rm T_{c} \\thinspace neutron \\thinspace  ^{3}P_{2}$', fontsize=22)

    for model in neutron_t_model_names:

        coefftnT0 = models['neutron_triplet'][model][0] * 0.1187 * MeV_erg / kB * 1e-9
        coefftnk0 = models['neutron_triplet'][model][1]
        coefftnk1 = models['neutron_triplet'][model][2]
        coefftnk2 = models['neutron_triplet'][model][3]
        coefftnk3 = models['neutron_triplet'][model][4]

        for i in range(len(k)):

            Tc_nt[i] = _Tcn(k[i])

        if idx>2:
            plt.plot(k,Tc_nt,color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx-2,0],dashes[idx-2,1]),label=model)
        else:
            plt.plot(k,Tc_nt,color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1],4,dashes[idx,1]),label=model)
        idx += 1

    plt.legend(loc='upper right',fontsize=15)

    plt.xlabel('$\\rm k_{F} \\thinspace (fm^{-1})$',fontsize=22)
    plt.ylabel('$\\rm T_{c} \\thinspace (10^{9} K)$',fontsize=22)
    plt.yticks([0.2,0.4,0.6,0.8,1.0],fontsize=22)
    plt.xticks([0,0.5,1,1.5,2,2.5,3,3.5],fontsize=22)
    plt.xlim(0,3.4)
    plt.ylim(0,1)
    plt.savefig('Tcrit_nt.pdf',format='pdf')

    plot_style()
    idx = 0
    plt.title('$\\rm T_{c} \\thinspace proton \\thinspace  ^{1}S_{0}$', fontsize=22)

    for model in proton_s_model_names:

        coeffspT0 = models['proton_singlet'][model][0] * 0.5669 * MeV_erg / kB * 1e-9
        coeffspk0 = models['proton_singlet'][model][1]
        coeffspk1 = models['proton_singlet'][model][2]
        coeffspk2 = models['proton_singlet'][model][3]
        coeffspk3 = models['proton_singlet'][model][4]

        for i in range(len(k)):

            Tc_p[i]  = _Tcp(k[i])

        if idx>2:
            plt.plot(k,Tc_p,color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx-2,0],dashes[idx-2,1]),label=model)
        else:
            plt.plot(k,Tc_p,color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1],4,dashes[idx,1]),label=model)
        idx += 1

    plt.legend(loc='upper right',fontsize=15)

    plt.xlabel('$\\rm k_{F} \\thinspace (fm^{-1})$',fontsize=22)
    plt.ylabel('$\\rm T_{c} \\thinspace (10^{9} K)$',fontsize=22)
    plt.yticks([0,2,4,6,8,10,12,14],fontsize=22)
    plt.xticks([0,0.5,1,1.5,2,2.5,3,3.5],fontsize=22)
    plt.xlim(0,1.5)
    plt.ylim(0,7)
    plt.savefig('Tcrit_ps.pdf',format='pdf')
    '''

T_crit()

def show_cooling_curves():

    plot_style(yticks=5)

    plt.axhline(y=5.4, linewidth=3, dashes = (dashes[2,0],dashes[2,1]), color='k')
    #plt.text(2.3, 5.51, ' When log $ \\rm T^{\infty}_{s}$ crosses',fontsize=17)
    #plt.text(2.3, 5.47, ' the line heat source',fontsize=17)
    #plt.text(2.3, 5.43, ' starts to work',fontsize=17)

    data = np.loadtxt('output/cooling_SF0.dat')
    plt.plot(np.log10(data[:,1]),np.log10(data[:,0]), color='y', lw=3,label='No SF')

    for i in [0]:
        for j in [0]:
            for k in [0,1,2,3,4,5,6,7,8]:

                idx = i + j + k
                data = np.loadtxt('output/cooling_' + ns[i] + '_' + nt[j] + '_' + ps[k] + '.dat')
                if(idx>2):
                    plt.plot(np.log10(data[:,1]),np.log10(data[:,0]),color=colors[idx],
                             linewidth=line_thickness[idx], dashes = (dashes[idx-2,0],dashes[idx-2,1]),label=ps[idx])
                else:
                    plt.plot(np.log10(data[:,1]),np.log10(data[:,0]),color=colors[idx],
                             linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1],4,dashes[idx,1]),label=ps[idx])

    plt.xlabel('log t [yr]',fontsize=22)
    plt.ylabel('log $\\rm T^{\infty}_{\\rm s}$ [K]',fontsize=22)
    plt.yticks([5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6],fontsize=22)
    plt.xticks([0,1,2,3,4,5,6],fontsize=22)
    plt.xlim(0.8,6.0)
    plt.ylim(5.3,6.44)
    plt.legend(loc='upper right',fontsize=17)
    plt.savefig('profile_SF_model_3_hm_cooling.pdf', format='pdf')

    plt.show()


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


def find_roots(value):

    data = np.loadtxt('output/cooling_SF0.dat')
    root = find_nearest(np.log10(data[:,0]),value)
    print(np.log10(data[root,0]),data[root,1])

    for k in [0]:
        for j in [0]:
            for i in [0,1,2,3,4,5,6,7,8]:

                data = np.loadtxt('output/cooling_' + ns[i] + '_' + nt[j] + '_' + ps[k] + '.dat')
                root = find_nearest(np.log10(data[:,0]),value)
                print(np.log10(data[root,0]),data[root,1])

def plot_1(b,c):

    plot_style(yticks=5,xticks=2)

    #yy = numpy.array([30.5,31,31.5,32,32.5,33,33.5,34.,34.5,35,35.5,36,36.5,37])
    #xx = numpy.array([0,10,20,30,40,50])
    yy = np.array([30,31,32,33,34,35,36,37])
    xx = np.array([0,1,2,3,4,5,6,7,8,9,10])

    plt.ylabel('$\\rm log \\thinspace L^{\infty} [erg \\thinspace s^{-1}]$',fontsize=32)
    plt.xlabel('$\\rm t \\thinspace [yr]$',fontsize=32)

    density_label = numpy.array(['$10^{11}$','$10^{12}$','$10^{13}$'])
    #period_label = numpy.array(['$10^{0}$','$10^{1}$'])
    period_label = numpy.array(['$10^{0}$'])

    plt.yticks(yy, fontsize=24)
    plt.xticks(xx, fontsize=24)

    data = np.loadtxt('output/_SF0_' + str(b*len(period_label)+c) + '_cooling.dat')
    data[:,1] -= time_root_no_sf
    plt.plot((data[:,1]),numpy.log10(data[:,5]),'y-',lw=3,label='NO SF')

    for i in [0]:
        for j in [0]:
            for k in [0,1,2,3,4,5,6,7,8]:

                data = np.loadtxt('output/_' + ns[i] + '_' + nt[j] + '_' + ps[k] +
                                  '_' + str(b*len(period_label)+c) + '_cooling.dat')
                idx = i + j + k
                data[:,1] -= time_roots[idx]
                if(idx>2):
                    plt.plot(data[:,1], numpy.log10(data[:,5]),color=colors[idx],
                             linewidth=line_thickness[idx], dashes = (dashes[idx-2,0],dashes[idx-2,1]),label=ps[idx])
                else:
                    plt.plot(data[:,1], numpy.log10(data[:,5]),color=colors[idx],
                             linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1],4,dashes[idx,1]),label=ps[idx])

    plt.ylim(y_min,y_max)
    plt.xlim(x_min,x_max)

    plt.text((x_max-x_min)/15, y_min + (y_max - y_min)*0.93 , '$ \\rm H_{0}  =  $' + '$10^{20}$' + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=22)
    plt.text((x_max-x_min)/15, y_min + (y_max - y_min)*0.88 , '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=22)
    plt.text((x_max-x_min)/15, y_min + (y_max - y_min)*0.83 , '$ \\rm \\Delta t  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=22)

    print(density_label[b])
    print(period_label[c])

    plt.legend(loc='upper right', fontsize = 16)
    plt.savefig('temperature_SF_model_3_hm_' + str(b) + '_' + str(c) + '.pdf', format='pdf')
    plt.show()


def temperature(b,c,time_index):

    plot_style()

    yy = numpy.array([0,1,2,3,4,5,6,7,8,9,10,11])
    xx = numpy.array([9,10,11,12,13,14,15])

    plt.yticks(yy, fontsize=24)
    plt.xticks(xx, fontsize=24)

    plt.ylabel('$\\rm T^{\infty}_{r}/10^{8} $',fontsize=32)
    plt.xlabel('$\\rm log \\thinspace \\rho_{1} [g \\thinspace cm^{-3}]$',fontsize=32)

    density_label = numpy.array(['$10^{11}$','$10^{12}$','$10^{13}$'])
    period_label = numpy.array(['$10^{0}$'])
    time_label  = numpy.array(['$0.01$','$0.1$','$0.5$','$1.0$','$2.0$','$3.0$',
                               '$4.0$','$5.0$','$6.0$','$7.0$','$8.0$','$9.0$','$10.0$'])

    plt.axvline(x=numpy.log10(source_rho_1[b]), linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')
    plt.axvline(x=numpy.log10(source_rho_2[b]), linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')

    data = np.loadtxt('output/_SF0_' + str(b*len(period_label)+c) + '_temperature.dat')
    data[:,1] -= time_root_no_sf
    plt.plot((numpy.log10(data[:,0])),data[:,time_index+1]/1e8,'y-',lw=3,label='NO SF')

    for i in [0]:
        for j in [0]:
            for k in [0,1,2,3,4,5,6,7,8]:

                data = np.loadtxt('output/_' + ns[i] + '_' + nt[j] + '_' + ps[k] +
                                  '_' + str(b*len(period_label)+c) + '_temperature.dat')
                idx = i + j + k
                data[:,1] -= time_roots[idx]
                if(idx>2):
                    plt.plot((numpy.log10(data[:,0])),data[:,time_index+1]/1e8,color=colors[idx],
                             linewidth=line_thickness[idx], dashes = (dashes[idx-2,0],dashes[idx-2,1]),label=ps[idx])
                else:
                    plt.plot((numpy.log10(data[:,0])),data[:,time_index+1]/1e8,color=colors[idx],
                             linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1],4,dashes[idx,1]),label=ps[idx])

    plt.text(10.1, y_min + (y_max-y_min)*0.93 , '$ \\rm H_{0}  =  $' + '$10^{20}$' + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=17)
    plt.text(10.1, y_min + (y_max-y_min)*0.86 , '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=17)
    plt.text(10.1, y_min + (y_max-y_min)*0.79 , '$ \\rm \\Delta t  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=17)

    plt.ylim(y_min,y_max)
    plt.xlim(x_min,x_max)

    print(density_label[b])
    print(period_label[c])
    print(time_label[time_index])

    plt.legend(loc='upper right',fontsize=15)
    plt.savefig('profile_SF_model_3_hm_' + str(b) + '_' + str(c) + '.pdf', format='pdf')
    plt.show()


def flux(b,c,time_index):

    plot_style()

    yy = numpy.array([-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,0,500,1000,1500,2000,2500])
    xx = numpy.array([9,10,11,12,13,14,15])

    plt.yticks(yy, fontsize=24)
    plt.xticks(xx, fontsize=24)

    plt.ylabel('$\\rm L^{\infty}_{r}/L_{\odot}$',fontsize=32)
    plt.xlabel('$\\rm log \\thinspace \\rho_{1} [g \\thinspace cm^{-3}]$',fontsize=32)

    density_label = numpy.array(['$10^{11}$','$10^{12}$','$10^{13}$'])
    period_label = numpy.array(['$10^{0}$','$10^{1}$'])
    time_label  = numpy.array(['$0.01$','$0.1$','$0.5$','$1.0$','$2.0$','$3.0$',
                               '$4.0$','$5.0$','$6.0$','$7.0$','$8.0$','$9.0$','$10.0$'])

    plt.axvline(x=numpy.log10(source_rho_1[b]), linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')
    plt.axvline(x=numpy.log10(source_rho_2[b]), linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')

    data = np.loadtxt('output/_SF0_' + str(b*len(period_label)+c) + '_flux.dat')
    data[:,1] -= time_root_no_sf
    plt.plot((numpy.log10(data[:,0])),data[:,time_index+1]/LSun,'y-',lw=3,label='NO SF')

    for i in [0]:
        for j in [0]:
            for k in [0,1,2,3,4,5,6,7,8]:
                print(k)
                data = np.loadtxt('output/_' + ns[i] + '_' + nt[j] + '_' + ps[k] +
                                  '_' + str(b*len(period_label)+c) + '_flux.dat')
                idx = i + j + k
                data[:,1] -= time_roots[idx]
                if(idx>2):
                    plt.plot((numpy.log10(data[:,0])),data[:,time_index+1]/LSun,color=colors[idx],
                             linewidth=line_thickness[idx], dashes = (dashes[idx-2,0],dashes[idx-2,1]),label=ps[idx])
                else:
                    plt.plot((numpy.log10(data[:,0])),data[:,time_index+1]/LSun,color=colors[idx],
                             linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1],4,dashes[idx,1]),label=ps[idx])

    plt.text(10.1, y_min + (y_max-y_min)*0.93 , '$ \\rm H_{0}  =  $' + '$10^{20}$' + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=17)
    plt.text(10.1, y_min + (y_max-y_min)*0.86 , '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=17)
    plt.text(10.1, y_min + (y_max-y_min)*0.79 , '$ \\rm \\Delta t  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=17)

    plt.ylim(y_min,y_max)
    plt.xlim(x_min,x_max)

    print(density_label[b])
    print(period_label[c])
    print(time_label[time_index])

    plt.legend(loc='upper right',fontsize=15)
    #plt.savefig('flux_SF_model_3_lm_' + str(b) + '_' + str(c) + '.pdf', format='pdf')
    plt.show()


def max_min_plot():

    plot_style()
    levels1 = numpy.array([1.1, 1.5, 2.2])
    x = numpy.linspace(-1, 2, 7)
    y = numpy.linspace(10, 14, 9)

    periods = numpy.power(10, numpy.linspace(-1, 2, 7))

    x = numpy.linspace(-1, 2, 7)
    y = numpy.linspace(10, 14, 9)
    '''
    z_no = numpy.zeros((len(x),len(y)))
    z_sf = numpy.zeros((len(x),len(y)))

    for j in range(0,len(y)):
        for i in range(0,len(x)):

            temp1 = numpy.loadtxt('output/_SF0_' + str(len(periods)*j + i) + '_cooling.dat')
            temp2 = numpy.loadtxt('output/_GIPSF_AO_CCDK_' + str(len(periods)*j + i) + '_cooling.dat')

            L_max_no   = 2000 + numpy.argmax(temp1[2000:,5])
            L_max_sf_1 = 2000 + numpy.argmax(temp2[2000:,5])

            z_sf[i,j] = temp2[L_max_sf_1,5]/temp2[2000,5]
            z_no[i,j] = temp1[L_max_no,5]  /temp1[2000,5]

            print(i, j)

    numpy.savetxt('contour_max_min_nosf.dat',z_no)
    numpy.savetxt('contour_max_min_sf.dat', z_sf)
    '''
    z_no = numpy.loadtxt('contour_max_min_nosf.dat')
    z_sf = numpy.loadtxt('contour_max_min_sf.dat')

    ratio1 = np.zeros_like(z_no)

    for j in range(0,len(y)):
        for i in range(0,len(x)):

            ratio1[i,j] = z_sf[i,j]/z_no[i,j]

    X, Y = numpy.meshgrid(x, y)

    plt.figure(1)
    im = plt.imshow(numpy.rot90(ratio1)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1, 2, 10, 14),vmax=2.5,vmin=1)
    CS1 = plt.contour(X, Y, numpy.transpose(ratio1),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS1, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,14,9),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)

    #plt.ylabel('$\mathrm{log} \\thinspace \\rho_{1} \\thinspace [\mathrm{g \\thinspace   cm^{-3}}]$',fontsize=32)
    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)

    CBI2 = plt.colorbar(im, orientation='vertical')
    CBI2.ax.tick_params(labelsize=24)
    CBI2.set_label('$\\rm L^{\infty}_{max} (SF) / L^{\infty}_{max} (non-SF)$', labelpad=5, y=0.45,fontsize=22)

    plt.savefig('contour_ratios.pdf',format='pdf')


def max_plot():

    plot_style()
    levels1 = numpy.array([1.1, 2, 10, 25, 35, 50, 75, 100])
    x = numpy.linspace(-1, 2, 7)
    y = numpy.linspace(10, 14, 9)

    z_no = numpy.loadtxt('contour_max_min_nosf.dat')
    z_sf = numpy.loadtxt('contour_max_min_sf.dat')

    X, Y = numpy.meshgrid(x, y)

    plt.figure(1)
    plt.title('No SF',fontsize=22)
    im = plt.imshow(numpy.rot90(z_no)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,2,10,14),vmax=100,vmin=1)
    CS1 = plt.contour(X, Y, numpy.transpose(z_no),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS1, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,14,9),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)

    plt.ylabel('$\mathrm{log} \\thinspace \\rho_{1} \\thinspace [\mathrm{g \\thinspace   cm^{-3}}]$',fontsize=32)
    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)

    plt.savefig('contour_max_min_nosf.pdf',format='pdf')

    plt.figure(2)
    plt.title('ns: GIPSF, nt: AO, ps: CCDK',fontsize=22)
    im = plt.imshow(numpy.rot90(z_sf)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,2,10,14),vmax=100,vmin=1)
    CS1 = plt.contour(X, Y, numpy.transpose(z_sf),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS1, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,14,9),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)

    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)

    CBI2 = plt.colorbar(im, orientation='vertical')
    CBI2.ax.tick_params(labelsize=24)
    CBI2.set_label('$\\rm  L^{\infty}_{max} / L^{\infty}_{min} $', labelpad=5, y=0.45,fontsize=22)

    plt.savefig('contour_max_min_sf.pdf',format='pdf')

    plt.show()


def delay_plot():

    plot_style()

    periods = numpy.power(10,numpy.linspace(-1,2,7))

    levels1 = numpy.array([0.1,0.3,1,2,3,5,10,17])

    x = numpy.linspace(-1,2,7)
    y = numpy.linspace(10,13.5,8)

    '''
    z_no = numpy.zeros((len(x),len(y)))
    z_sf = numpy.zeros((len(x),len(y)))


    for j in range(0,len(y)):
        for i in range(0,len(x)):

            temp1 = numpy.loadtxt('output/_SF0_' + str(len(periods)*j + i) + '_cooling.dat')
            temp2 = numpy.loadtxt('output/_GIPSF_AO_CCDK_' + str(len(periods)*j + i) + '_cooling.dat')

            L_max_no = 2000 + numpy.argmax(temp1[2000:,5])
            L_max_sf = 2000 + numpy.argmax(temp2[2000:,5])
            H_max = 2000 + numpy.argmax(temp1[2000:,3])

            time_diff1 = temp1[L_max_no,1] - temp1[H_max,1]
            time_diff2 = temp2[L_max_sf,1] - temp2[H_max,1]

            z_no[i,j] = time_diff1
            z_sf[i,j] = time_diff2

            print(i,j)

    numpy.savetxt('contour_time_delay_nosf.dat',z_no)
    numpy.savetxt('contour_time_delay_sf.dat',z_sf)
    '''
    z_no = numpy.loadtxt('contour_time_delay_nosf.dat')
    z_sf = numpy.loadtxt('contour_time_delay_sf.dat')

    X, Y = numpy.meshgrid(x, y)

    plt.figure(1)
    im_no = plt.imshow(numpy.rot90(z_no)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,2,10,13.5),vmin=0,vmax=15)
    CS1 = plt.contour(X, Y, numpy.transpose(z_no),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS1, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,13.5,8),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)

    plt.ylabel('$\mathrm{log} \\thinspace \\rho_{1} \\thinspace [\mathrm{g \\thinspace   cm^{-3}}]$',fontsize=32)
    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)
    plt.savefig('contour_time_delay_nosf.pdf',format='pdf')

    levels1 = numpy.array([0.1,0.3,1,2,3,5,17])
    plt.figure(2)
    im_sf = plt.imshow(numpy.rot90(z_sf)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,2,10,13.5),vmin=0,vmax=15)
    CS2 = plt.contour(X, Y, numpy.transpose(z_sf),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS2, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,13.5,8),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)
    CBI2 = plt.colorbar(im_sf, orientation='vertical',ticks=numpy.array([15,10,5,0]))
    CBI2.ax.tick_params(labelsize=24)
    CBI2.set_label('$ \mathrm{t_{retarded}}$ [yr]', labelpad=0, y=0.45,fontsize=24)

    #plt.ylabel('$\mathrm{log} \\thinspace \\rho_{1} \\thinspace [\mathrm{g \\thinspace   cm^{-3}}]$',fontsize=32)
    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)
    plt.savefig('contour_time_delay_sf.pdf',format='pdf')

    plt.show()



def delay_plot_dif():

    plot_style()

    periods = numpy.power(10,numpy.linspace(-1,2,7))

    x = numpy.linspace(-1,2,7)
    y = numpy.linspace(10,13.5,8)

    '''
    z_no = numpy.zeros((len(x),len(y)))
    z_sf = numpy.zeros((len(x),len(y)))


    for j in range(0,len(y)):
        for i in range(0,len(x)):

            temp1 = numpy.loadtxt('output/_SF0_' + str(len(periods)*j + i) + '_cooling.dat')
            temp2 = numpy.loadtxt('output/_GIPSF_AO_CCDK_' + str(len(periods)*j + i) + '_cooling.dat')

            L_max_no = 2000 + numpy.argmax(temp1[2000:,5])
            L_max_sf = 2000 + numpy.argmax(temp2[2000:,5])
            H_max = 2000 + numpy.argmax(temp1[2000:,3])

            time_diff1 = temp1[L_max_no,1] - temp1[H_max,1]
            time_diff2 = temp2[L_max_sf,1] - temp2[H_max,1]

            z_no[i,j] = time_diff1
            z_sf[i,j] = time_diff2

            print(i,j)

    numpy.savetxt('contour_time_delay_nosf.dat',z_no)
    numpy.savetxt('contour_time_delay_sf.dat',z_sf)
    '''
    z_no = numpy.loadtxt('contour_time_delay_nosf.dat')
    z_sf = numpy.loadtxt('contour_time_delay_sf.dat')
    difference = z_no-z_sf

    X, Y = numpy.meshgrid(x, y)

    levels1 = numpy.array([1,4,7,10, 15])
    plt.figure(1)
    im_sf = plt.imshow(numpy.rot90(difference)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,2,10,13.5),vmin=0,vmax=20)
    CS2 = plt.contour(X, Y, numpy.transpose(difference),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS2, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,13.5,8),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)
    CBI2 = plt.colorbar(im_sf, orientation='vertical',ticks=numpy.array([20,15,10,5,0]))
    CBI2.ax.tick_params(labelsize=24)
    CBI2.set_label('t$_{\\rm nosf}$ - t$_{\\rm sf}$ [yr]', labelpad=0, y=0.45,fontsize=28)

    #plt.ylabel('$\mathrm{log} \\thinspace \\rho_{1} \\thinspace [\mathrm{g \\thinspace   cm^{-3}}]$',fontsize=32)
    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)
    plt.savefig('contour_time_delay_diff.pdf',format='pdf')

    plt.show()


def halfwidth_plot():

    plot_style()

    periods = numpy.power(10,numpy.linspace(-1,2,7))

    levels1 = numpy.array([1.7,2.5,3.1,3.7])
    x = numpy.linspace(-1,1.5,6)
    y = numpy.linspace(10,13.5,8)
    z = numpy.zeros((len(x),len(y)))
    '''
    for j in range(0,len(y)):
        for i in range(0,len(x)):

            temp1 = numpy.loadtxt('output/_SF0_' + str(len(periods)*j + i) + '_cooling.dat')
            temp2 = numpy.loadtxt('output/_GIPSF_AO_CCDK_' + str(len(periods)*j + i) + '_cooling.dat')

            L_max_no = 2000 + numpy.argmax(temp1[2000:,5])
            L_max_sf = 2000 + numpy.argmax(temp2[2000:,5])

            A = 2000

            value1 = (temp1[A,5]+temp1[L_max_no,5])/2
            value2 = (temp2[A,5]+temp2[L_max_sf,5])/2

            left1 = A + find_nearest(temp1[A:L_max_no,5] ,value1)
            right1 = L_max_no + find_nearest(temp1[L_max_no:,5] ,value1)

            left2 = A + find_nearest(temp2[A:L_max_sf,5] ,value2)
            right2 = L_max_sf + find_nearest(temp2[L_max_sf:,5] ,value2)

            z[i,j] = (temp1[right1,1] - temp1[left1,1])/(temp2[right2,1] - temp2[left2,1])

            print(i,j)

    numpy.savetxt('contour_halfwidth_ratios.dat',z)
    '''
    z = numpy.loadtxt('contour_halfwidth_ratios.dat')

    X, Y = numpy.meshgrid(x, y)

    plt.figure(1)
    im = plt.imshow(numpy.rot90(z)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,1.5,10,13.5),vmin=1,vmax=4)
    CS1 = plt.contour(X, Y, numpy.transpose(z),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS1, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,1.5,6),fontsize=24)
    plt.yticks(numpy.linspace(10,13.5,8),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)

    plt.ylabel('$\mathrm{log} \\thinspace \\rho_{1} \\thinspace [\mathrm{g \\thinspace   cm^{-3}}]$',fontsize=32)
    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)

    CBI2 = plt.colorbar(im, orientation='vertical',ticks=numpy.array([1,2,3,4,5]))
    CBI2.ax.tick_params(labelsize=24)
    CBI2.set_label('$\\rm \\frac{halfwidth (non-SF)}{halfwidht (SF)}$', labelpad=5, y=0.45,fontsize=30)

    plt.savefig('contour_halfwidth_ratios.pdf',format='pdf')

    plt.show()


def visual(a,b,c):

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams.update({'figure.autolayout': True})

    f, (axarr) = plt.subplots(1, sharex=True)

    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    axarr.xaxis.set_minor_locator(x_minor_locator)
    axarr.yaxis.set_minor_locator(y_minor_locator)

    yy = numpy.array([32,33,34,35,36,37,38,39,40,42,44,46])
    #xx = numpy.array([0,50,100,150])
    xx = numpy.array([0,5,10,15,20,25,30,35,40,45,50])
    #xx = numpy.array([0,10,20,30,40,50,60,70,80,90,100])

    plt.yticks(yy,fontsize=24)
    plt.xticks(xx,fontsize=24)

    axarr.set_ylim(y_min,y_max)
    #axarr.set_xlim(-10,110)
    axarr.set_xlim(x_min,x_max)

    axarr.set_ylabel('$\\rm log \\thinspace L^{\infty} [erg \\thinspace s^{-1}]$',fontsize=32)
    axarr.set_xlabel('$\\rm t \\thinspace [yr]$',fontsize=32)


    rho1_array = numpy.power(10,numpy.linspace(10,14,9))
    periods      = numpy.power(10,numpy.linspace(-1,2,7))
    source_power = numpy.power(10,numpy.linspace(18,21,4))

    source_label  = numpy.array(['$10^{18}$','$10^{19}$','$10^{20}$','$10^{21}$'])
    density_label = numpy.array(['$10^{10}$','$10^{10.5}$','$10^{11}$','$10^{11.5}$','$10^{12}$','$10^{12.5}$','$10^{13}$','$10^{13.5}$','$10^{14}$'])
    period_label = numpy.array(['$10^{-1}$','$10^{-0.5}$','$10^{0}$','$10^{0.5}$','$10^{1}$','$10^{1.5}$','$10^{2}$'])

    for i in numpy.array([c]):

        cooling_sf = numpy.loadtxt('output/_GIPSF_AO_CCDK_' + str(len(periods)*b + i) + '_cooling.dat')
        cooling_sf[:,1] -= 978.7404

        cooling_no = numpy.loadtxt('output/_SF0_' + str(len(periods)*b + i) + '_cooling.dat')
        cooling_no[:,1] -= 19506.77

        print(len(periods)*len(rho1_array)*a + len(periods)*b + i)

        plt.plot((cooling_no[:,1]),numpy.log10(cooling_no[:,3]),'y-',lw=3,label='Heater')
        plt.plot((cooling_no[:,1]),numpy.log10(cooling_no[:,5]),'g-',lw=3,label='Surface (non-superfluid)')
        plt.plot((cooling_sf[:,1]),numpy.log10(cooling_sf[:,5]),'r--',lw=3.5,label='Surface (superfluid)')

        L_max_sf = 2000 + numpy.argmax(cooling_sf[2000:,5])
        print(L_max_sf)
        L_max_no = 2000 + numpy.argmax(cooling_no[2000:,5])
        print(L_max_no)
        H_max = 2000 + numpy.argmax(cooling_no[2000:,3])
        print(H_max)

        plt.plot(cooling_sf[L_max_sf,1], numpy.log10(cooling_sf[L_max_sf,5]),  'ko', markersize=12)
        plt.plot(cooling_no[L_max_no,1], numpy.log10(cooling_no[L_max_no,5]),  'ko', markersize=12)
        plt.plot(cooling_no[H_max,1], numpy.log10(cooling_no[H_max,3]),  'ko', markersize=12)

        A = 000

        value1 = (cooling_no[A,5]+cooling_no[L_max_no,5])/2
        value2 = (cooling_sf[A,5]+cooling_sf[L_max_sf,5])/2

        left1 = A + find_nearest(cooling_no[A:L_max_no,5] ,value1)
        right1 = L_max_no + find_nearest(cooling_no[L_max_no:,5] ,value1)

        left2 = A + find_nearest(cooling_sf[A:L_max_sf,5] ,value2)
        right2 = L_max_sf + find_nearest(cooling_sf[L_max_sf:,5] ,value2)

        plt.plot(cooling_sf[left2,1], numpy.log10(cooling_sf[left2,5]),  'co', markersize=12)
        plt.plot(cooling_no[left1,1], numpy.log10(cooling_no[left1,5]),  'mo', markersize=12)

        plt.plot(cooling_sf[right2,1], numpy.log10(cooling_sf[right2,5]),  'co', markersize=12)
        plt.plot(cooling_no[right1,1], numpy.log10(cooling_no[right1,5]),  'mo', markersize=12)

    plt.text((x_max-x_min)/2.2, y_max-2 , '$ \\rm H_{0}  =  $' + source_label[a] + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=22)
    plt.text((x_max-x_min)/2.2, y_max-2.5 , '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=22)
    plt.text((x_max-x_min)/2.2, y_max-3 , '$ \\rm P  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=22)

    plt.legend(loc='upper right',fontsize=22)
    #plt.savefig('profile_' + str(a) + '_' + str(b) + '_' + str(c) + '.pdf', format='pdf')
    plt.show()


x_max=150
x_min=-1
y_min=32.5
y_max=36
#visual(2,7,0)
#np.log10(7.2946929931640625e+014)
x_max=15
x_min=10
y_min=0
y_max=2

#max_min_plot()
#max_plot()

#y_min = -2000
#y_max = +2000
#x_max=8
#x_min=-1
#y_min=30.5
#y_max=34

#find_roots(5.9736530)
#plot_1(2,0)
#flux(2,0,2)
#temperature(2,0,2)
#show_cooling_curves()
#check_models()
#T_crit()
#halfwidth_plot()
#delay_plot()
#delay_plot_dif()

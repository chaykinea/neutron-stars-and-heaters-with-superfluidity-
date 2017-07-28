import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from data    import loaddata
from scipy import interpolate
from control.constants import *

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

names = np.array(['SF0','AWP2', 'GIPSF', 'SCLBL'])
labels = np.array(['No SF','AWP2', 'GIPSF', 'SCLBL'])
colors = np.array(['black','darkblue','red','darkorange'])
shape = np.array(['s','^','o','d'])
line_thickness = np.array([4.8, 2.7,2.2,3.2])
dashes = np.array([[2,1e-15],[4,8],[20,3],[9,6]])
order = np.array([1,4,2,3])


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

def plot_style(xticks=5,yticks=5):
    
    global ax
    
    plt.rc('text', usetex=True)
    #plt.rcParams['mathtext.fontset'] = 'cm'
    #plt.rcParams['mathtext.rm'] = 'serif'
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
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
    ax.tick_params(axis='both', which='both', pad=8, left='on', right='on',top='on',bottom='on')

    plt.rcParams['lines.linewidth'] = 1.0
    plt.rcParams['lines.dashed_pattern'] = [6, 6] 
    plt.rcParams['lines.dashdot_pattern'] = [3, 5, 1, 5]
    plt.rcParams['lines.dotted_pattern'] = [1, 3]
    plt.rcParams['lines.scale_dashes'] = False
    plt.rcParams['errorbar.capsize'] = 6

def shift(Pow=5, change=0):

    rho_values = np.array([3.16227766e+10, 1e11, 3.16227766e+11, 1e12, 3.16227766e+12, 1e13])
    rho_sample = np.linspace(10.5,13, 200)
    rho_change = np.array([0,1,-1,2,-2,3])
    y_sample = np.zeros((3,200,3))
    Pow_idx = 0

    plot_style()

    for Pow in [3,5,4]:

        array_hm = np.zeros((3, 6))

        for rho in range(0, 6):

            if rho==2 or rho==4:
                num = + 2*Pow + int((rho-2)/2)
                data = np.loadtxt('output_hm2/cooling_' + names[0] + '_' + str(num) + '.dat')
                ref = data[4000 + np.argmax(data[4000:, -2]),-2]

                for i in range(1, 4):

                    data = np.loadtxt('output_hm2/cooling_' + names[i] + '_' + str(num) + '.dat')
                    max_arg = 4000 + np.argmax(data[4000:, -2])
                    print('output_hm/cooling_' + names[i] + '_' + str(num) + '.dat')
                    array_hm[int(i-1), rho] += data[max_arg, -2]/ref
                    #plt.scatter(np.log10(rho_values[rho]), data[max_arg, -2]/ref, marker=shape[0], s=170,facecolor=colors[i], edgecolor='black',linewidth=2.8, zorder=3)

            else:
                num = 8*Pow + 2*rho_change[rho]
                data = np.loadtxt('output_hm/cooling_' + names[0] + '_' + str(num) + '.dat')
                ref = data[4000 + np.argmax(data[4000:, -2]),-2]

                for i in range(1, 4):

                    data = np.loadtxt('output_hm/cooling_' + names[i] + '_' + str(num) + '.dat')
                    max_arg = 4000 + np.argmax(data[4000:, -2])
                    print('output_hm/cooling_' + names[i] + '_' + str(num) + '.dat')
                    array_hm[int(i-1), rho] += data[max_arg, -2]/ref
                    #plt.scatter(np.log10(rho_values[rho]), data[max_arg, -2]/ref, marker=shape[0], s=170,facecolor=colors[i], edgecolor='black',linewidth=2.8, zorder=3)

        for i in range(0, 3):
            hm = interpolate.interp1d(np.log10(rho_values),array_hm[i,:], kind='cubic')
            y_sample[i,:,Pow_idx] = hm(rho_sample)

        Pow_idx += 1
    
    idx = 1
    plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),label=labels[idx])
    idx += 1
    plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),label=labels[idx])
    idx += 1
    plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),label=labels[idx])

    plt.xticks([10.5,11,11.5,12,12.5,13],fontsize=20)
    plt.xlim(11.0,13)

    zorder = np.array([3,2,2,2])
    for i in range(0,3):
        idx = i+1
        #plt.fill_between(rho_sample, y_sample[i,:,0], y_sample[i,:,1], facecolor=colors[i+1], interpolate=True, alpha=0.5)
        plt.plot(rho_sample, y_sample[i,:,2], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]))
        #plt.plot(rho_sample, y_sample[i,:,0],lw=3,color=colors[i+1],zorder=zorder[i+1])
        plt.plot(rho_sample, y_sample[i,:,1], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]))

    plt.yticks([1,1.5,2,2.5,3,3.5,4],fontsize=20)
    plt.ylim(0.9,4)

    plt.xlabel('$\\rm log \\thinspace$$\\rho_{1}   \\thinspace \\thinspace \\thinspace \\rm [g \\thinspace cm^{-3}]$',fontsize=24)
    plt.ylabel('$L^{\infty, \\rm sf}_{\\rm max}/L^{\infty}_{\\rm max}$',fontsize=24)
    plt.text(11.05,3.1,'$ M = 1.85 \\thinspace \\rm M \odot $', fontsize=23)

    plt.text(11.96, 3.7, '$ L_{\\rm h0} = 10^{3} L_{\\rm hc}$', fontsize=21,color=colors[2])
    plt.text(12.3, 3.13, '$ L_{\\rm h0} = 5 \\times 10^{2} L_{\\rm hc}$', fontsize=21,color=colors[2])
    plt.text(12.24, 2.6, '$ L_{\\rm h0} = 10^{3} L_{\\rm hc}$', fontsize=21,color=colors[1])
    plt.text(12.28, 1.87, '$ L_{\\rm h0} = 5 \\times 10^{2} L_{\\rm hc}$', fontsize=21,color=colors[1])
    plt.text(12.45, 1.38, '$ L_{\\rm h0} = 10^{3} L_{\\rm hc}$', fontsize=21,color=colors[3])
    plt.text(12.36, 0.96, '$ L_{\\rm h0} = 5 \\times 10^{2} L_{\\rm hc}$', fontsize=20,color=colors[3])
    #plt.text(11.75, 2.0, '$ H_{\\rm 0} = 2 \\times 10^{2} L_{\\rm hc}$', fontsize=21)

    plt.legend(loc='upper left', fontsize=20,frameon=False)
    plt.savefig('ratios_density_m1.eps', format='eps')
    plt.show()

shift()

def shift2(Pow=5, change=0):

    rho_values = np.array([3.16227766e+10, 1e11, 3.16227766e+11, 1e12, 3.16227766e+12, 1e13])
    rho_sample = np.linspace(10.5,13, 200)
    rho_change = np.array([0,1,-1,2,-2,3])
    y_sample = np.zeros((3,200,3))
    Pow_idx = 0

    plot_style()

    for Pow in [3,5,4]:

        array_hm = np.zeros((3, 6))

        for rho in range(0, 6):

            if rho==2 or rho==4:
                num = + 2*Pow + int((rho-2)/2)
                data = np.loadtxt('output_lm2/cooling_' + names[0] + '_' + str(num) + '.dat')
                ref = data[4000 + np.argmax(data[4000:, -2]),-2]

                for i in range(1, 4):

                    data = np.loadtxt('output_lm2/cooling_' + names[i] + '_' + str(num) + '.dat')
                    max_arg = 4000 + np.argmax(data[4000:, -2])
                    print('output_lm/cooling_' + names[i] + '_' + str(num) + '.dat')
                    array_hm[int(i-1), rho] += data[max_arg, -2]/ref
                    #plt.scatter(np.log10(rho_values[rho]), data[max_arg, -2]/ref, marker=shape[0], s=170,facecolor=colors[i], edgecolor='black',linewidth=2.8, zorder=3)

            else:
                num = 8*Pow + 2*rho_change[rho]
                data = np.loadtxt('output_lm/cooling_' + names[0] + '_' + str(num) + '.dat')
                ref = data[4000 + np.argmax(data[4000:, -2]),-2]

                for i in range(1, 4):

                    data = np.loadtxt('output_lm/cooling_' + names[i] + '_' + str(num) + '.dat')
                    max_arg = 4000 + np.argmax(data[4000:, -2])
                    print('output_lm/cooling_' + names[i] + '_' + str(num) + '.dat')
                    array_hm[int(i-1), rho] += data[max_arg, -2]/ref
                    #plt.scatter(np.log10(rho_values[rho]), data[max_arg, -2]/ref, marker=shape[0], s=170,facecolor=colors[i], edgecolor='black',linewidth=2.8, zorder=3)

        for i in range(0, 3):
            hm = interpolate.interp1d(np.log10(rho_values),array_hm[i,:], kind='cubic')
            y_sample[i,:,Pow_idx] = hm(rho_sample)

        Pow_idx += 1


    plt.xticks([10.5,11,11.5,12,12.5,13],fontsize=20)
    plt.xlim(11.0,13)

    zorder = np.array([3,2,2,2])
    for i in range(0,3):
        idx = i+1
        #plt.fill_between(rho_sample, y_sample[i,:,0], y_sample[i,:,1], facecolor=colors[i+1], interpolate=True, alpha=0.5)
        plt.plot(rho_sample, y_sample[i,:,2], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]))
        #plt.plot(rho_sample, y_sample[i,:,0],lw=3,color=colors[i+1],zorder=zorder[i+1])
        plt.plot(rho_sample, y_sample[i,:,1], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]))

    plt.yticks([1,1.5,2,2.5,3,3.5,4],fontsize=20)
    plt.ylim(0.9,2.5)

    plt.xlabel('$\\rm log \\thinspace$$\\rho_{1}   \\thinspace \\thinspace \\thinspace \\rm [g \\thinspace cm^{-3}]$',fontsize=24)
    #plt.ylabel('$L^{\infty,\\rm sf}_{\\rm max}/L^{\infty}_{\\rm max}$',fontsize=24)
    plt.text(11.1,2.37,'$ M = 1.40 \\thinspace \\rm M \odot$', fontsize=23)

    #plt.text(11.65, 1.55, '$ H_{\\rm 0} = 2 \\times 10^{2} H_{\\rm c}$', fontsize=21)
    plt.text(11.702, 1.98, '$ L_{\\rm h0} = 5 \\times 10^{2} L_{\\rm hc}$', fontsize=21, color=colors[2])
    plt.text(11.85,2.41, '$ L_{\\rm h0} = 10^{3} L_{\\rm hc}$', fontsize=21,color=colors[2])
    plt.text(11.74, 1.66, '$ L_{\\rm h0} =  10^{3} L_{\\rm hc}$', fontsize=21, color=colors[1])
    plt.text(11.93,1.34, '$ L_{\\rm h0} = 5 \\times 10^{2} L_{\\rm hc}$', fontsize=21,color=colors[1])
    plt.text(12.45, 1.13, '$ L_{\\rm h0} = 10^{3} L_{\\rm hc}$', fontsize=21,color=colors[3])
    plt.text(12.4, 0.96, '$ L_{\\rm h0} = 5 \\times 10^{2} L_{\\rm hc}$', fontsize=20,color=colors[3])

    #plt.legend(loc='upper left', fontsize=20)
    plt.savefig('ratios_density_m2.eps', format='eps')
    plt.show()

shift2(Pow=4,change=1)

def _Tcs(nn):

    kf = np.power(3.*pi*pi*nn, 1./3.)
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


def T_crit():

    global coeffsnk0, coeffsnk1, coeffsnk2, coeffsnk3, coeffsnT0

    loaddata.superfluid_data_init()
    loaddata.star_model_data_init()

    rho = loaddata.star_model()[:, 3]
    log_rho = np.log10(rho)

    nb = loaddata.star_model()[:, 5]
    ne = loaddata.star_model()[:, 6]
    nm = loaddata.nm(np.log(rho))
    Phi = loaddata.star_model()[:, 4]

    npp = ne+nm
    nn = nb-npp

    plot_style()

    output = np.zeros((len(nn),4))
    output[:,0] = log_rho

    #plt.title('$\\rm T_{c} \\thinspace neutron \\thinspace  ^{1}S_{0}$', fontsize=22)

    for idx,model in zip(range(1,4),['AWP2', 'GIPSF', 'SCLBL']):

        coeffsnT0 = models['neutron_singlet'][model][0] * 0.5669 * MeV_erg / kB * 1e-9
        coeffsnk0 = models['neutron_singlet'][model][1]
        coeffsnk1 = models['neutron_singlet'][model][2]
        coeffsnk2 = models['neutron_singlet'][model][3]
        coeffsnk3 = models['neutron_singlet'][model][4]

        print(colors[idx])

        Tc_ns = np.zeros(len(nn))

        for i in range(len(nn)):

            Tc_ns[i] = _Tcs(nn[i])*np.exp(Phi[i])

        output[:, idx] = Tc_ns
        plt.plot(log_rho,Tc_ns,color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),label=labels[idx],zorder=order[idx])


    plt.legend(loc='upper left',fontsize=23,frameon=False)

    plt.xlabel('$\\rm log \\thinspace$$\\rho \\thinspace \\thinspace \\thinspace \\rm [g \\thinspace cm^{-3}]$',fontsize=24)
    plt.ylabel('$T_{\\rm c}$  $/(\\rm 10^{9} \\thinspace K)$',fontsize=22)
    plt.axvline(np.log10(1.5e14),dashes=(32,6),lw=2,color='black')
    #plt.axvline(np.log10(4.e11),lw=2,color='black')
    plt.yticks([0,2,4,6,8,10,12,14],fontsize=22)
    plt.xticks([10,11,12,13,14,15],fontsize=22)
    plt.xlim(10.9,14.5)
    plt.ylim(0,14)

    np.save('Tcrit_data.npy', output)
    #plt.show()
    plt.savefig('Tcrit_ns.eps',format='eps')

T_crit()

def heat_cap():

    plot_style()

    T = np.loadtxt('tables_hm/file1.dat')
    log_T = np.log(T)
    rho = np.loadtxt('tables_hm/file2.dat')
    log_rho = np.log(rho)

    C = np.loadtxt('tables_hm/file4.dat')
    C_AWP2 = np.loadtxt('tables_hm/file40.dat')
    C_GIPSF = np.loadtxt('tables_hm/file44.dat')
    C_SCLBL = np.loadtxt('tables_hm/file46.dat')
    C_e = np.loadtxt('tables_hm/file4e.dat')
    C_eb = np.loadtxt('tables_hm/file4eb.dat')

    x, y = np.meshgrid(log_T, log_rho)
    p = np.vstack([x.flatten(), y.flatten()]).T

    log_C = interpolate.LinearNDInterpolator(p, np.log(C).flatten())
    log_C2 = interpolate.LinearNDInterpolator(p, np.log(C_AWP2).flatten())
    log_C3 = interpolate.LinearNDInterpolator(p, np.log(C_GIPSF).flatten())
    log_C4 = interpolate.LinearNDInterpolator(p, np.log(C_SCLBL).flatten())
    log_C5 = interpolate.LinearNDInterpolator(p, np.log(C_e).flatten())
    log_C6 = interpolate.LinearNDInterpolator(p, np.log(C_eb).flatten())

    heat_capacity = lambda a,b: np.exp(log_C(np.log(a), np.log(b)))
    heat_capacity_AWP2 = lambda a,b: np.exp(log_C2(np.log(a), np.log(b)))
    heat_capacity_GIPSF = lambda a,b: np.exp(log_C3(np.log(a), np.log(b)))
    heat_capacity_SCLBL = lambda a,b: np.exp(log_C4(np.log(a), np.log(b)))
    heat_capacity_e = lambda a,b: np.exp(log_C5(np.log(a), np.log(b)))
    heat_capacity_eb = lambda a,b: np.exp(log_C6(np.log(a), np.log(b)))

    rho = np.logspace(10,np.log10(1e15),300)

    for T in [1e7,1e9]:
        idx = 0
        plt.plot(np.log10(rho),np.log10(heat_capacity(T,rho)), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]))
        plt.plot(np.log10(rho),np.log10(heat_capacity_eb(T,rho)),dashes=[3,3],color='black',lw=2.5,zorder=order[idx])
        idx += 1
        plt.plot(np.log10(rho),np.log10(heat_capacity_AWP2(T,rho)), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
        idx += 1
        plt.plot(np.log10(rho),np.log10(heat_capacity_GIPSF(T,rho)), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
        idx += 1
        plt.plot(np.log10(rho),np.log10(heat_capacity_SCLBL(T,rho)), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
        #plt.plot(np.log10(rho),np.log10(heat_capacity_e(T,rho)),dashes=[3,3],color='black',lw=2.5)
        
    
    idx = 0
    plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),label=labels[idx])
    idx += 1
    plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),label=labels[idx])
    idx += 1
    plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),label=labels[idx])
    idx += 1
    plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),label=labels[idx])
    plt.plot([-10,-20], [-10,-20], lw=2.5, color='black',dashes=[3,3], label='$\\rm C_{e} + C_{ion}$')

    plt.yticks(np.array([15,16,17,18,19,20,21]),fontsize=20)
    plt.xticks(np.array([10,11,12,13,14,15]),fontsize=20)
    plt.xlim(10,15)
    plt.ylim(14.9,21)
    plt.legend(loc='upper left',fontsize=19, frameon=False)

    plt.text(10.3,15.6,'$T = 10^{7} \\thinspace \\rm K$',fontsize=21)
    #plt.text(10.6,16.6,'$T = 10^{8} \\thinspace \\rm K$',fontsize=21)
    plt.text(10.5,18,'$T = 10^{9} \\thinspace \\rm K$',fontsize=21)

    plt.xlabel('$\\rm log \\thinspace$$\\rho   \\thinspace \\thinspace \\thinspace \\rm [g \\thinspace cm^{-3}]$',fontsize=24)
    plt.ylabel('$\\rm log \\thinspace$$ C_{V}   \\thinspace \\thinspace \\thinspace \\thinspace \\rm [erg \\thinspace K^{-1} cm^{-3}]$',fontsize=24)
    
    plt.savefig('heatcap.eps', format='eps')

heat_cap()








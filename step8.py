import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from scipy import integrate
from scipy import stats
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator

names = np.array(['SF0','AWP2', 'GIPSF', 'SCLBL'])
labels = np.array(['No SF','AWP2', 'GIPSF', 'SCLBL'])
colors = np.array(['black','darkblue','red','darkorange'])
shape = np.array(['s','^','o','d'])
line_thickness = np.array([4.8, 2.7,2.2,3.2])
dashes = np.array([[2,1e-15],[4,8],[20,3],[9,6]])
order = np.array([1,4,2,3])

order2 = np.array([1,2,3,4,5])
colors2 = np.array(['green','orange','red','blue','yellow'])
dashes2 = np.array([[2,1e-15],[7,1],[3,2],[1.5,3]])
line_thickness2 = np.array([3.5,3,2.5,2.,1.5])

colors3 = np.array(['green','orange','red','blue','yellow','brown','gray','black'])

def fixlogax(ax, a='x'):
    if a == 'x':
        labels = [item.get_text() for item in ax.get_xticklabels()]
        positions = ax.get_xticks()
        # print positions
        # print labels
        for i in range(len(positions)):
            labels[i] = '$10^{'+str(int(np.log10(positions[i])))+'}$'
        if np.size(np.where(positions == 1)) > 0:
            labels[np.where(positions == 1)[0][0]] = '$1$'
        if np.size(np.where(positions == 10)) > 0:
            labels[np.where(positions == 10)[0][0]] = '$10$'
        if np.size(np.where(positions == 0.1)) > 0:
            labels[np.where(positions == 0.1)[0][0]] = '$0.1$'
        # print positions
        # print labels
        ax.set_xticklabels(labels)
    if a == 'y':
        labels = [item.get_text() for item in ax.get_yticklabels()]
        positions = ax.get_yticks()
        # print positions
        # print labels
        for i in range(len(positions)):
            labels[i] = '$10^{'+str(int(np.log10(positions[i])))+'}$'
        if np.size(np.where(positions == 1)) > 0:
            labels[np.where(positions == 1)[0][0]] = '$1$'
        if np.size(np.where(positions == 10)) > 0:
            labels[np.where(positions == 10)[0][0]] = '$10$'
        if np.size(np.where(positions == 0.1)) > 0:
            labels[np.where(positions == 0.1)[0][0]] = '$0.1$'
        # print positions
        # print labels
        ax.set_yticklabels(labels)

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

MSun = 1.98892e33
LSun = 3.846e33
c = 2.99792458e10
G = 6.67259e-8
from_ev_to_erg = 1.6021766208e-12
from_yr_to_sec = 31536000

T_IGR = np.array([99.7, 91.5, 89.2, 84.8, 88.5, 84.6, 82.8])
t_IGR = np.array([55609, 55680.5, 55810.5, 56060, 56187.5, 56228, 56340]) - 55556
err_IGR = np.array([1.6, 1.5, 1.5, 1.5, 1.9, 2.0, 1.2])

def plot1():
    plot_style()
    k_b = 8.617330350e-5

    data1 = np.loadtxt('output/cooling_GIPSF_' + str(0+28) + '.dat')
    outdata = np.zeros([len(data1[::5,0]),15])
    outdata[:,0] = (data1[::5, 1]-1.001000e3)*365

    for idx in [6,5,4,3,2,1,0]:
 
        data1 = np.loadtxt('output/cooling_GIPSF_' + str(idx+28) + '.dat')
        data2 = np.loadtxt('output/cooling_SF0_' + str(idx+28) + '.dat')

        plt.plot((data1[:, 1]-1.0010e3)*365, data1[:, 0]*k_b, color=colors3[idx])
        plt.plot((data2[:, 1]-1.0010e3)*365, data2[:, 0]*k_b, '--',color=colors3[idx])

        outdata[:,idx*2+1] = data1[::5, 0]*k_b
        outdata[:,idx*2+2] = data2[::5, 0]*k_b

    np.savetxt('relaxation_curves.dat',outdata,fmt='%1.5e')

    plt.plot([10,10],[20,20],color='black',lw=2,label='GIPSF model')
    plt.plot([10,10],[20,20],'--',color='black',lw=2,label='No SF')

    plt.xscale('log')
    plt.xticks([1,10,100,1000,10000], fontsize=20)

    plt.text(1.50,155,'$\dot{M} = 10^{-8} \\thinspace \\rm M \\odot / yr$',fontsize=22)
    plt.text(1.50,151,'$\widetilde{T}_{0} = 10^{8} \\thinspace \\rm K$',fontsize=22)
    
    plt.yticks([50,60,70,80,90,100,110,120,130,140,150,160], fontsize=20)
    plt.ylim(105, 160)
    plt.xlabel('$\\rm Days \\thinspace  \\thinspace  since  \\thinspace  \\thinspace  '
               'onset \\thinspace  \\thinspace  of  \\thinspace  \\thinspace quiescence$',fontsize=22)
    plt.ylabel('$kT^{\infty}_{\\rm s}$ $\\rm eV$',fontsize=22)
    plt.legend(loc='upper right',fontsize=19,scatterpoints=1,frameon=False)
    fixlogax(plt.gca(), a='x')
    plt.xlim(1,10000)
    plt.savefig('relaxation_curves_H_Delta_M_acc_10-6.eps',format='eps')

    plt.show()

#plot1()


def tite_plot():

    old_fe = np.loadtxt('data/tite3.dat',skiprows=3)
    old_acc = np.loadtxt('data/tite.dat',skiprows=3)
    new_he = np.loadtxt('data/tite2.dat',skiprows=3)

    plot_style()
    
    plt.plot(old_fe[:,0], old_fe[:,1],lw=3,color='red',label='Potekhin 1997, Fe only ($\Delta M_{\\rm acc} = 0$)')
    plt.plot(old_acc[:,0], old_acc[:,1],lw=3,color='green',label='Potekhin 1997, $\Delta M_{\\rm acc}/ M_{\\rm star} = 10^{-6}$')
    plt.plot(new_he[:,0], new_he[:,1],lw=3,linestyle='--',color='blue',label='Beznogov 2016, He only')
    plt.xticks([5,6,7,8,9,10,11], fontsize=20)
    plt.yticks([4,5,6,7,8], fontsize=20)
    plt.ylim(5,7)
    plt.xlim(6,9)
    plt.ylabel('$\\rm log \\thinspace $$T_{\\rm s} \\thinspace \\thinspace \\thinspace \\rm [K]$',fontsize=24)
    plt.xlabel('$\\rm log \\thinspace $$T_{\\rm b} \\thinspace \\thinspace \\thinspace \\rm [K]$',fontsize=24)
    plt.legend(loc='upper left',fontsize=19,scatterpoints=1,frameon=False)
    plt.savefig('tstb_models.pdf',format='pdf')
     
    plt.show()

#tite_plot()


def plot2(): 
    plot_style()
    k_b = 8.617330350e-5
    config = np.loadtxt('data/config.dat')
    data1 = np.loadtxt('output/cooling_GIPSF_' + str(0+35) + '.dat')
    outdata = np.zeros([len(data1[::5,0]),15])
    outdata[:,0] = (data1[::5, 1]-1.000200e3)*365

    for idx in [6,5,4,3,2,1,0]:
 
        data1 = np.loadtxt('output/cooling_GIPSF_' + str(idx+35) + '.dat')
        data2 = np.loadtxt('output/cooling_SF0_' + str(idx+35) + '.dat')

        plt.plot((data1[:, 1]-1.000200e3)*365, data1[:, 0]*k_b, color=colors3[idx],linewidth='3')
        plt.plot((data2[:, 1]-1.000200e3)*365, data2[:, 0]*k_b, '--',color=colors3[idx],linewidth='3')

        outdata[:,idx*2+1] = data1[::5, 0]*k_b
        outdata[:,idx*2+2] = data2[::5, 0]*k_b

    plt.plot([10,10],[20,20],color='black',lw=2,label='GIPSF model')
    plt.plot([10,10],[20,20],'--',color='black',lw=2,label='No SF')
    np.savetxt('IGR.dat',outdata,fmt='%1.5e')
    data11 = np.loadtxt('output/cooling_SF0_' + str(0+43) + '.dat')
    data111 = np.loadtxt('output/cooling_SF0_' + str(0+44) + '.dat')
    plt.plot((data11[:, 1]-1.000200e3)*365, data11[:, 0]*k_b, '--',color='k',linewidth='3')
    plt.plot((data111[:, 1]-1.000200e3)*365, data111[:, 0]*k_b, '--',color='k',linewidth='3')
    plt.scatter(t_IGR, T_IGR, s=100, color='black', marker='^',label='IGR J17480--2446', zorder=6)
    plt.errorbar(x=t_IGR, y=T_IGR, yerr=err_IGR, color='black', fmt=' ',zorder=6)

    plt.xscale('log')
    plt.xticks([1,10,100,1000,10000], fontsize=20)
  
    plt.text(1.50,80,'$\dot{M} = 3.5\\times 10^{-9} \\thinspace \\rm M \\odot / yr$',fontsize=22)
    plt.text(1.50,73,'$kT_{\\rm s 0}^{\infty} = 69 \\thinspace \\thinspace \\rm eV$',fontsize=22)
    plt.text(1.50,69,'$\Delta t = 12.5 \\thinspace \\thinspace \\rm yr$',fontsize=22)
    
    plt.yticks(np.array([90,100,110,120,130,140,150,160,170,180])-30, fontsize=20)
    plt.ylim(70, 110)
    plt.xlabel('$\\rm Days \\thinspace  \\thinspace  since  \\thinspace  \\thinspace  '
               'onset \\thinspace  \\thinspace  of  \\thinspace  \\thinspace quiescence$',fontsize=22)
    plt.ylabel('$kT^{\infty}_{\\rm s}$ $\\rm eV$',fontsize=22)
    plt.legend(loc='upper right',fontsize=19,scatterpoints=1,frameon=False)
    fixlogax(plt.gca(), a='x')
    plt.xlim(1,100000)
    #plt.savefig('IGR.pdf',format='pdf')
    plt.show()

plot2()


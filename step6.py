import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from scipy import interpolate

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

names = np.array(['SF0','AWP2', 'GIPSF', 'SCLBL'])
labels = np.array(['No SF','AWP2', 'GIPSF', 'SCLBL'])
colors = np.array(['black','red','green','darkblue'])
shape = np.array(['s','^','o','d'])


def plot_style(xticks=5,yticks=5):

    plt.rcParams.update({'figure.autolayout': True})
    #plt.tight_layout()
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

def show_curves_1(Pow=4, dt=0, timeshift=4e4):

    plot_style()

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams.update({'figure.autolayout': True})
    zorders = np.array([2,1,1,1])
    zorders2 = np.array([4,3,3,3])

    for rho in range(1,2):
        num = dt + 8*Pow + 2*rho
        for i in range(0,4):
            data = np.loadtxt('output_lm/cooling_' + names[i] + '_' + str(num) + '.dat')
            data2 = np.loadtxt('output_hm/cooling_' + names[i] + '_' + str(num) + '.dat')

            max_arg = 4000 + np.argmax(data[4000:, -2])
            max_arg2 = 4000 + np.argmax(data2[4000:, -2])

            power_ref = data[4000,3]
            power_ref2 = data2[4000,3]

            plt.plot(data[:,1]-timeshift, np.log10(data[:,-2]/power_ref), dashes=(17,4,4,4), color=colors[i], lw=3.5, zorder=zorders[i])
            plt.scatter(data[max_arg,1]-timeshift, np.log10(data[max_arg, -2]/power_ref), marker=shape[0], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])
            plt.plot(data2[:,1]-timeshift, np.log10(data2[:,-2]/power_ref2),'--', color=colors[i], lw=3, zorder=zorders[i])
            plt.scatter(data2[max_arg2,1]-timeshift, np.log10(data2[max_arg2, -2]/power_ref2), marker=shape[1], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])

            if i==0:

                max_arg = 4000 + np.argmax(data[4000:, 3])
                plt.plot(data[:,1]-timeshift, np.log10(data[:,3]/power_ref), color=colors[i], lw=3.5)
                plt.scatter(data[max_arg,1]-timeshift, np.log10(data[max_arg, 3]/power_ref), marker=shape[2], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])

                plt.plot([-10,-20], [-10,-20],'-', lw=3, color=colors[0], label=labels[0])
                plt.plot([-10,-20], [-10,-20], '-', lw=3, color=colors[1], label=labels[1])
                plt.plot([-10,-20], [-10,-20], '-', lw=3, color=colors[2], label=labels[2])
                plt.plot([-10,-20], [-10,-20], '-', lw=3, color=colors[3], label=labels[3])

        plt.legend(loc='upper right', fontsize=20)

        plt.xticks(np.linspace(0,35,8),fontsize=22)
        plt.yticks(np.linspace(-4,3,8),fontsize=22)
        plt.ylim(-3.4, 3.1)
        plt.xlim(-1,21)
        plt.ylabel('$\\rm log $$L^{\infty}_{\\rm{s,h}}/L_{\\rm{c}}^{\infty}$', fontsize=26)
        plt.xlabel('$t \\thinspace \\rm{yr}$', fontsize=26)

        plt.text(6,0.2,'$\\rm HEATER$',fontsize=22)
        #plt.text(6,1.5,'$\\rm H_{0} = 5\\times 10^{2} H_{c} $',fontsize=22)
        plt.text(6.5,2.6,'$\\rm \\rho_{1} = 10^{11} \\thinspace g \\thinspace cm^{-3}$',fontsize=22)
        plt.text(8,-1.3,'$\\rm SURFACE, \\thinspace M = 1.40 M_{\odot}$',fontsize=22)
        plt.text(5,-2.8,'$\\rm SURFACE, \\thinspace M = 1.85 M_{\odot}$',fontsize=22)

        plt.savefig('fig1.pdf',format='pdf')
        plt.show()

def show_curves_2(Pow=4, dt=0, timeshift=4e4):

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams.update({'figure.autolayout': True})
    zorders = np.array([2,1,1,1])
    zorders2 = np.array([4,3,3,3])

    for rho in range(2,3):
        num = dt + 8*Pow + 2*rho
        for i in range(0,4):
            data = np.loadtxt('output_lm/cooling_' + names[i] + '_' + str(num) + '.dat')
            data2 = np.loadtxt('output_hm/cooling_' + names[i] + '_' + str(num) + '.dat')

            max_arg = 4000 + np.argmax(data[4000:, -2])
            max_arg2 = 4000 + np.argmax(data2[4000:, -2])

            power_ref = data[4000,3]
            power_ref2 = data2[4000,3]

            ax2 = plt.subplot(3,1,2)
            x_minor_locator = AutoMinorLocator(5)
            y_minor_locator = AutoMinorLocator(5)
            plt.tick_params(which='both', width=1.7)
            plt.tick_params(which='major', length=9)
            plt.tick_params(which='minor', length=5)
            plt.xlim(-1,21)
            plt.yticks(np.linspace(-1,-0.5,2),fontsize=22)
            plt.ylim(-1.1,-0.4)
            ax2.xaxis.set_minor_locator(x_minor_locator)
            ax2.yaxis.set_minor_locator(y_minor_locator)

            plt.plot(data[:,1]-timeshift, np.log10(data[:,-2]/power_ref), dashes=(17,4,4,4), color=colors[i], lw=3.5, zorder=zorders[i])
            plt.scatter(data[max_arg,1]-timeshift, np.log10(data[max_arg, -2]/power_ref), marker=shape[0], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])

            ax3 = plt.subplot(3,1,3)
            x_minor_locator = AutoMinorLocator(5)
            y_minor_locator = AutoMinorLocator(5)
            plt.tick_params(which='both', width=1.7)
            plt.tick_params(which='major', length=9)
            plt.tick_params(which='minor', length=5)
            plt.xlim(-1,21)
            plt.yticks(np.linspace(-3.,-1,3),fontsize=22)
            plt.ylim(-3.3,-0.5)
            ax3.xaxis.set_minor_locator(x_minor_locator)
            ax3.yaxis.set_minor_locator(y_minor_locator)
            plt.plot(data2[:,1]-timeshift, np.log10(data2[:,-2]/power_ref2),'--', color=colors[i], lw=3, zorder=zorders[i])
            plt.scatter(data2[max_arg2,1]-timeshift, np.log10(data2[max_arg2, -2]/power_ref2), marker=shape[1], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])

            if i==0:

                ax1 = plt.subplot(3,1,1)
                x_minor_locator = AutoMinorLocator(5)
                y_minor_locator = AutoMinorLocator(5)
                plt.tick_params(which='both', width=1.7)
                plt.tick_params(which='major', length=9)
                plt.tick_params(which='minor', length=5)
                plt.yticks(np.linspace(0,3,4),fontsize=22)
                plt.ylim(-0.1,3.1)
                plt.xlim(-1,21)
                ax1.xaxis.set_minor_locator(x_minor_locator)
                ax1.yaxis.set_minor_locator(y_minor_locator)
                max_arg = 4000 + np.argmax(data[4000:, 3])
                plt.plot(data[:,1]-timeshift, np.log10(data[:,3]/power_ref), color=colors[i], lw=3.5)
                plt.scatter(data[max_arg,1]-timeshift, np.log10(data[max_arg, 3]/power_ref), marker=shape[2], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])

                plt.plot([-10,-20], [-10,-20],'-', lw=3, color=colors[0], label=labels[0])
                plt.plot([-10,-20], [-10,-20], '-', lw=3, color=colors[1], label=labels[1])
                plt.plot([-10,-20], [-10,-20], '-', lw=3, color=colors[2], label=labels[2])
                plt.plot([-10,-20], [-10,-20], '-', lw=3, color=colors[3], label=labels[3])

        plt.legend(loc='upper right', fontsize=20)

        plt.xticks(np.linspace(-0,35,8),fontsize=22)
        plt.xlim(-1,21)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.ylabel('$\\rm log $$L^{\infty}_{\\rm{s}}/L_{\\rm{c}}^{\infty}$', fontsize=26)
        plt.xlabel('$t \\thinspace \\rm{yr}$', fontsize=26)

        plt.subplot(3,1,1)
        plt.text(2,0.4,'$\\rm HEATER$',fontsize=22)
        #plt.text(13,2.4,'$\\rm H_{0} = 5\\times 10^{2} H_{c} $',fontsize=22)
        plt.text(13,2.4,'$\\rm \\rho_{1} = 10^{12} \\thinspace g \\thinspace cm^{-3}$',fontsize=22)
        plt.ylabel('$\\rm log $$L^{\infty}_{\\rm{h}}/L_{\\rm{c}}^{\infty}$', fontsize=26)
        plt.subplot(3,1,2)
        plt.text(9,-0.55,'$\\rm SURFACE, \\thinspace M = 1.40 M_{\odot}$',fontsize=22)
        plt.ylabel('$\\rm log $$L^{\infty}_{\\rm{s}}/L_{\\rm{c}}^{\infty}$', fontsize=26)
        plt.subplot(3,1,3)
        plt.text(9,-1.1,'$\\rm SURFACE, \\thinspace M = 1.85 M_{\odot}$',fontsize=22)

        plt.savefig('fig2.pdf',format='pdf')
        plt.show()

#show_curves_1()
#show_curves_2()

def temperature_profile_profile(rho=2, Pow=4, dt=0):

    time_points = np.array([0.01, 0.1, 0.5, 1., 2., 3., 4., 5, 6., 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22])
    Tcrit = np.load('Tcrit_data.npy')

    plot_style()

    num = dt + 8*Pow + 2*rho
    numbers = np.array([1, 2, 3, 6, 18])
    arg_num = np.array([0, 3, 4])
    dashes = np.array([[20, 4, 20, 4],[3, 3, 3, 3],[10,1e-15,10,1e-15],[10,3,2,3],[12,3,12,3]])
    thickness = np.array([5.5, 3, 3.5, 2.5, 2])
    idx = 0

    ax1 = plt.subplot(2,1,1)
    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    ax1.xaxis.set_minor_locator(x_minor_locator)
    ax1.yaxis.set_minor_locator(y_minor_locator)
    plt.yticks(np.array([8,8.2,8.4,8.6,8.8,9]),fontsize=20)
    plt.xticks(np.array([9,10,11,12,13,14,15]),fontsize=20)
    if(rho==1):
        plt.ylabel('$\\rm log \\thinspace$$ T^{\infty} \\thinspace \\rm K$',fontsize=24)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.xlim(9,15)

    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,1]*1e9), color=colors[1], lw=1.5)
    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,2]*1e9), color=colors[2], lw=1.5)
    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,3]*1e9), color=colors[3], lw=1.5)

    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,1]*1e9), facecolor=colors[1], interpolate=True, alpha=0.25)
    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,2]*1e9), facecolor=colors[2], interpolate=True, alpha=0.25)
    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,3]*1e9), facecolor=colors[3], interpolate=True, alpha=0.25)

    if rho==2:
        plt.ylim(8.17,8.73)
        plt.text(9.2,8.23,'$\\rho_{1} = 10^{12} \\thinspace \\rm g \\thinspace cm^{-3}$', fontsize=23)
    elif rho==3:
        plt.ylim(8.17,8.63)
        plt.text(9.2,8.22,'$\\rho_{1} = 10^{13} \\thinspace \\rm g \\thinspace cm^{-3}$', fontsize=23)
    else:
        plt.ylim(8.17,8.83)
        plt.text(9.1,8.74,'$\\rho_{1} = 10^{11} \\thinspace \\rm g \\thinspace cm^{-3}$', fontsize=23)
        plt.text(9.1,8.67,'$\\rm M = 1.40 \\thinspace M_{\\odot}$',fontsize=23)

    for time in numbers[:3]:
        for i in range(3,-1,-1):
            if time==1 and dt==0:
                data = np.loadtxt('temperature_lm/temperature_' + names[i] + '_' + str(num+1) + '.dat')
                plt.plot(np.log10(data[:,0]), np.log10(data[:,time]),  dashes = dashes[idx,:], color=colors[i], lw=thickness[idx])
            else:
                data = np.loadtxt('temperature_lm/temperature_' + names[i] + '_' + str(num) + '.dat')
                plt.plot(np.log10(data[:,0]), np.log10(data[:,time]), dashes = dashes[idx,:], color=colors[i], lw=thickness[idx])
        idx += 1

    if rho==1:
        plt.text(12.45, 8.5, '$T^{\\infty} < T_{\\rm cr}^{\\infty}$', fontsize=23)
    elif rho==2:
        plt.text(13.11, 8.5, '$T^{\\infty} < T_{\\rm cr}^{\\infty}$', fontsize=23)
        plt.plot([-10,-20], [-10,-20], lw=thickness[4], color='black', dashes=dashes[4,:], label='$t = 15.0 \\rm \\thinspace yr$')
        plt.plot([-10,-20], [-10,-20], lw=thickness[3], color='black', dashes=dashes[3,:], label='$t = 3.0 \\rm \\thinspace yr$')
        plt.plot([-10,-20], [-10,-20], lw=thickness[2], color='black', dashes=dashes[2,:], label='$t = 0.5 \\rm \\thinspace yr$')
        plt.plot([-10,-20], [-10,-20], lw=thickness[1], color='black', dashes=dashes[1,:], label='$t = 0.1 \\rm \\thinspace yr$')
        plt.plot([-10,-20], [-10,-20], lw=thickness[0], color='black', dashes=dashes[0,:], label='$t = 0.0 \\rm \\thinspace yr$')
        plt.legend(loc='upper left',fontsize=20)
    else:
        plt.text(13.11, 8.5, '$T^{\\infty} < T_{\\rm cr}^{\\infty}$', fontsize=23)

        plt.plot([-10,-20], [-10,-20], lw=thickness[4], color='black', dashes=dashes[4,:], label='$t = 15.0 \\rm \\thinspace yr$')
        plt.plot([-10,-20], [-10,-20], lw=thickness[3], color='black', dashes=dashes[3,:], label='$t = 3.0 \\rm \\thinspace yr$')
        plt.plot([-10,-20], [-10,-20], lw=thickness[2], color='black', dashes=dashes[2,:], label='$t = 0.5 \\rm \\thinspace yr$')
        plt.plot([-10,-20], [-10,-20], lw=thickness[1], color='black', dashes=dashes[1,:], label='$t = 0.1 \\rm \\thinspace yr$')
        plt.plot([-10,-20], [-10,-20], lw=thickness[0], color='black', dashes=dashes[0,:], label='$t = 0.0 \\rm \\thinspace yr$')
        plt.legend(loc='upper left',fontsize=20)

    idx = 0
    ax2 = plt.subplot(2,1,2)
    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    ax2.xaxis.set_minor_locator(x_minor_locator)
    ax2.yaxis.set_minor_locator(y_minor_locator)
    if(rho==2 or rho==3):
        plt.yticks(np.array([8,8.1,8.2,8.3,8.4,8.5,8.6]),fontsize=20)
    else:
        plt.yticks(np.array([8,8.2,8.4,8.6,8.8,9]),fontsize=20)
        plt.ylabel('$\\rm log \\thinspace$$ T^{\infty} \\thinspace \\rm K$',fontsize=24)
    plt.xlabel('$\\rm log \\thinspace$$\\rho \\thinspace \\rm g \\thinspace cm^{-3}$',fontsize=24)
    plt.xticks(np.array([9,10,11,12,13,14,15]),fontsize=20)
    plt.xlim(9, 15)

    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,1]*1e9), color=colors[1], lw=1.5)
    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,2]*1e9), color=colors[2], lw=1.5)
    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,3]*1e9), color=colors[3], lw=1.5)

    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,1]*1e9), facecolor=colors[1], interpolate=True, alpha=0.25)
    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,2]*1e9), facecolor=colors[2], interpolate=True, alpha=0.25)
    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,3]*1e9), facecolor=colors[3], interpolate=True, alpha=0.25)

    if rho==2:
        plt.ylim(8.17,8.53)
    elif rho==3:
        plt.ylim(8.17,8.35)
    else:
        plt.ylim(8.17,8.68)

    for time in numbers[arg_num]:
        for i in range(3,-1,-1):
            if time==1 and dt==0:
                data = np.loadtxt('temperature_lm/temperature_' + names[i] + '_' + str(num+1) + '.dat')
                plt.plot(np.log10(data[:,0]), np.log10(data[:,time]),  dashes = dashes[arg_num[idx],:], color=colors[i], lw=thickness[arg_num[idx]])
            else:
                data = np.loadtxt('temperature_lm/temperature_' + names[i] + '_' + str(num) + '.dat')
                plt.plot(np.log10(data[:,0]), np.log10(data[:,time]), dashes = dashes[arg_num[idx],:], color=colors[i], lw=thickness[arg_num[idx]])
        idx += 1

    if rho==1:
        plt.text(12.45, 8.5, '$T^{\\infty} < T_{\\rm cr}^{\\infty}$', fontsize=23)
    elif rho==3:
        plt.text(12.5, 8.31, '$T^{\\infty} < T_{\\rm cr}^{\\infty}$', fontsize=23)
    else:
        plt.text(13, 8.4, '$T^{\\infty} < T_{\\rm cr}^{\\infty}$', fontsize=23)

    plt.savefig('profile_2.pdf', format='pdf')
    plt.show()

#temperature_profile_profile(rho=2)
#temperature_profile_profile()


def timedelay(Pow=5, dt=0):

    rho_values = np.array([3.16227766e+10, 1e11, 3.16227766e+11, 1e12, 3.16227766e+12, 1e13])
    rho_change = np.array([0,1,-1,2,-2,3])
    plot_style()
    Pow_idx = 0

    y_sample = np.zeros((4, 200, 3))
    rho_sample = np.linspace(10.5,13, 200)

    for Pow in [0, 5, 4]:

        array_lm = np.zeros((4, 6))

        for rho in range(0, 6):

            if rho == 2 or rho == 4:

                num = + 2*Pow + int((rho-2)/2)

                for i in range(0, 4):

                    data = np.loadtxt('output_lm2/cooling_' + names[i] + '_' + str(num) + '.dat')
                    max_arg = 4000 + np.argmax(data[4000:, -2])
                    max_arg3 = 4000 + np.argmax(data[4000:, 3])

                    print('output_lm2/cooling_' + names[i] + '_' + str(num) + '.dat')

                    array_lm[i, rho] += data[max_arg, 1]-data[max_arg3, 1]
                    #plt.scatter(np.log10(rho_values[rho]), array_lm[i, rho], marker=shape[0], s=170,facecolor=colors[i], edgecolor='black',linewidth=2.8, zorder=5)

            else:

                num = dt + 8*Pow + 2*rho_change[rho]

                for i in range(0, 4):

                    data = np.loadtxt('output_lm/cooling_' + names[i] + '_' + str(num) + '.dat')
                    max_arg = 4000 + np.argmax(data[4000:, -2])
                    max_arg3 = 4000 + np.argmax(data[4000:, 3])

                    print('output_lm/cooling_' + names[i] + '_' + str(num) + '.dat')

                    array_lm[i, rho] += data[max_arg, 1]-data[max_arg3, 1]
                    #plt.scatter(np.log10(rho_values[rho]), array_lm[i, rho], marker=shape[0], s=170,facecolor=colors[i], edgecolor='black',linewidth=2.8, zorder=5)

        for i in range(0, 4):

            lm = interpolate.interp1d(np.log10(rho_values),array_lm[i,:], kind='cubic')
            y_sample[i,:,Pow_idx] = lm(rho_sample)

        Pow_idx += 1

    alpha = np.array([0.5, 0.5, 0.4, 0.5])
    zorder = np.array([4,2,2,3])

    for i in range(0,4):
        plt.fill_between(rho_sample, y_sample[i,:,0], y_sample[i,:,1], facecolor=colors[i], interpolate=True, alpha=alpha[i])
        plt.plot(rho_sample, y_sample[i,:,2],lw=3,color=colors[i],zorder=zorder[i])

    plt.xticks([10.5,11,11.5,12,12.5,13],fontsize=20)
    plt.xlim(11.4,13)

    plt.yticks([0,2,4,6,8,10],fontsize=20)
    plt.ylim(1,10)

    plt.xlabel('$\\rm log \\thinspace$$\\rho_{1} \\thinspace \\rm g \\thinspace cm^{-3}$',fontsize=24)
    plt.ylabel('$\Delta t_{\\rm r}$',fontsize=24)

    #plt.text(11.45,6.8,'$ H_{\\rm{0}}/H_{\\rm{c}} = 10^{2}$', fontsize=23)
    plt.text(11.45,9,'$ M = 1.4 M_{\odot}$', fontsize=25)

    plt.legend(loc='upper left', fontsize=20)
    #plt.savefig('timedelay.pdf', format='pdf')
    plt.show()

def kappa():

    plot_style()
    rho_values = np.array([1e10,10**10.5,1e11,10**11.5,1e12,10**12.5,1e13])
    A = np.zeros((2,len(rho_values)))
    idx = 0

    T = np.loadtxt('tables_hm/file1.dat')
    log_T = np.log(T)
    rho = np.loadtxt('tables_hm/file2.dat')
    log_rho = np.log(rho)
    x, y = np.meshgrid(log_T, log_rho)
    p = np.vstack([x.flatten(), y.flatten()]).T

    C = np.loadtxt('tables_hm/file4.dat')

    log_C = interpolate.LinearNDInterpolator(p, np.log(C).flatten())
    heat_capacity = lambda a,b: np.exp(log_C(np.log(a), np.log(b)))
    rho = np.logspace(10,np.log10(1e15),300)
    for T in [1e8,1e9]:
        plt.plot(np.log10(rho),np.log10(heat_capacity(T,rho)), dashes=[20,8],color='black',lw=3,zorder=3)
        A[idx,:] = heat_capacity(T*np.ones(len(rho_values)),rho_values)
        idx += 1

    #print(A[0,:])
    #print(A[0,:])
    print(A[1,:]/A[0,:])
    idx = 0

    C = np.loadtxt('tables_hm/file44.dat')

    log_C = interpolate.LinearNDInterpolator(p, np.log(C).flatten())
    heat_capacity = lambda a,b: np.exp(log_C(np.log(a), np.log(b)))
    rho = np.logspace(10,np.log10(1e15),300)
    for T in [1e8,1e9]:
        plt.plot(np.log10(rho),np.log10(heat_capacity(T,rho)), dashes=[20,8],color='black',lw=3,zorder=3)
        A[idx,:] = heat_capacity(T*np.ones(len(rho_values)),rho_values)
        idx += 1

    #print(A[0,:])
    #print(A[0,:])
    print(A[1,:]/A[0,:])
    plt.yticks(np.array([15,16,17,18,19,20,21]),fontsize=20)
    plt.xticks(np.array([10,11,12,13,14,15]),fontsize=20)
    plt.xlim(10,15)
    #plt.ylim(14.9,27)

    plt.xlabel('$\\rm log \\thinspace$$\\rho \\thinspace \\rm g \\thinspace cm^{-3}$',fontsize=24)
    plt.ylabel('$\\rm log \\thinspace$$ C_{V}  \\thinspace \\thinspace \\rm erg \\thinspace K^{-1} cm^{-3}$',fontsize=24)
    plt.show()

#kappa()
#flux()
#timedelay()


A = np.loadtxt('output/cooling_AWP2_0.dat')
B = np.loadtxt('output/cooling_SF0_0.dat')

plt.plot(A[:,1]-2e3,np.log10(A[:,5]))
plt.plot(B[:,1]-2e3,np.log10(B[:,5]))
plt.show()



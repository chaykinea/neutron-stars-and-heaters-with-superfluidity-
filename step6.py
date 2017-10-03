import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from scipy import interpolate

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

names = np.array(['SF0','AWP2', 'GIPSF', 'SCLBL'])
labels = np.array(['No SF','AWP2', 'GIPSF', 'SCLBL'])
colors = np.array(['black','darkblue','red','darkorange'])
colors2 = np.array(['black','deepskyblue','salmon','#ffcc99'])
shape = np.array(['s','^','o','d'])
line_thickness = np.array([4.8, 2.7,2.2,3.2])
dashes = np.array([[2,1e-15],[4,8],[20,3],[9,6]])
order = np.array([1,4,2,3])

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

def show_curves_1(Pow=4, dt=0, timeshift=4e4):

    plot_style()

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams.update({'figure.autolayout': True})
    zorders = np.array([2,1,1,1]) + 5
    zorders2 = np.array([4,3,3,3]) + 5

    for rho in range(1,2):
        num = dt + 8*Pow + 2*rho
        for i in range(0,4):
            idx = i
            data = np.loadtxt('output_lm/cooling_' + names[i] + '_' + str(num) + '.dat')
            data2 = np.loadtxt('output_hm/cooling_' + names[i] + '_' + str(num) + '.dat')

            max_arg = 4000 + np.argmax(data[4000:, -2])
            max_arg2 = 4000 + np.argmax(data2[4000:, -2])

            power_ref = data[4000,3]
            power_ref2 = data2[4000,3]

            plt.plot(data[:,1]-timeshift, np.log10(data[:,-2]/power_ref), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
            plt.scatter(data[max_arg,1]-timeshift, np.log10(data[max_arg, -2]/power_ref), marker=shape[0], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])
            plt.plot(data2[:,1]-timeshift, np.log10(data2[:,-2]/power_ref2), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
            plt.scatter(data2[max_arg2,1]-timeshift, np.log10(data2[max_arg2, -2]/power_ref2), marker=shape[1], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])

            if i==0:

                max_arg = 4000 + np.argmax(data[4000:, 3])
                plt.plot(data[:,1]-timeshift, np.log10(data[:,3]/power_ref), color=colors[i], lw=3.5)
                plt.scatter(data[max_arg,1]-timeshift, np.log10(data[max_arg, 3]/power_ref), marker=shape[2], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])

                plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])
                idx += 1
                plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])
                idx += 1
                plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])
                idx += 1
                plt.plot([-10,-20], [-10,-20], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])

        plt.legend(loc='upper right', fontsize=21, frameon=False)

        plt.xticks(np.linspace(0,35,8),fontsize=29)
        plt.yticks(np.linspace(-4,3,8),fontsize=29)
        plt.ylim(-3.4, 3.1)
        plt.xlim(-1,21)
        plt.ylabel('$\\rm log$ $L^{\infty}_{\\rm{s,h}}/L_{\\rm{h0}}^{\infty}$', fontsize=29)
        plt.xlabel('$t, \\thinspace \\thinspace \\rm{yr}$', fontsize=29)

        plt.text(6,0.2,'$\\rm HEATER$',fontsize=22)
        #plt.text(6,1.5,'$\\rm H_{0} = 5\\times 10^{2} H_{c} $',fontsize=22)
        plt.text(5.8,2.5,'$\\rm \\rho_{1} = 10^{11} \\thinspace \\thinspace g \\thinspace cm^{-3}$',fontsize=26)
        plt.text(8,-1.3,'$\\rm SURFACE,$ $M = 1.40 \\thinspace \\rm M \odot$',fontsize=22)
        plt.text(5,-2.8,'$\\rm SURFACE,$ $M = 1.85 \\thinspace \\rm M \odot$',fontsize=22)

        plt.savefig('fig1.eps',format='eps')
        plt.show()

def show_curves_2(Pow=4, dt=0, timeshift=4e4):
    plot_style(xticks=5,yticks=5)

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams.update({'figure.autolayout': True})
    zorders = np.array([2,1,1,1])
    zorders2 = np.array([4,3,3,3])

    for rho in range(2,3):
        num = dt + 8*Pow + 2*rho
        for i in range(0,4):

            idx = i
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
            plt.yticks(np.linspace(-1,-0.5,2),fontsize=29)
            plt.ylim(-1.1,-0.4)
            ax2.xaxis.set_minor_locator(x_minor_locator)
            ax2.yaxis.set_minor_locator(y_minor_locator)
            ax2.tick_params(axis='both', which='both', pad=8, left='on', right='on',top='on',bottom='on')

            plt.plot(data[:,1]-timeshift, np.log10(data[:,-2]/power_ref), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
            plt.scatter(data[max_arg,1]-timeshift, np.log10(data[max_arg, -2]/power_ref), marker=shape[0], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=order[idx]+1)

            ax3 = plt.subplot(3,1,3)
            x_minor_locator = AutoMinorLocator(5)
            y_minor_locator = AutoMinorLocator(5)
            plt.tick_params(which='both', width=1.7)
            plt.tick_params(which='major', length=9)
            plt.tick_params(which='minor', length=5)

            plt.xlim(-1,21)
            plt.yticks(np.linspace(-3.,-1,3),fontsize=29)
            plt.ylim(-3.3,-0.5)

            ax3.xaxis.set_minor_locator(x_minor_locator)
            ax3.yaxis.set_minor_locator(y_minor_locator)
            ax3.tick_params(axis='both', which='major', pad=8)
            ax3.tick_params(axis='both', which='both', pad=8, left='on', right='on',top='on',bottom='on')
            plt.plot(data2[:,1]-timeshift, np.log10(data2[:,-2]/power_ref2), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
            plt.scatter(data2[max_arg2,1]-timeshift, np.log10(data2[max_arg2, -2]/power_ref2), marker=shape[1], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=order[idx]+1)

            if i==0:

                ax1 = plt.subplot(3,1,1)

                x_minor_locator = AutoMinorLocator(5)
                y_minor_locator = AutoMinorLocator(5)
                plt.tick_params(which='both', width=1.7)
                plt.tick_params(which='major', length=9)
                plt.tick_params(which='minor', length=5)
                plt.yticks(np.linspace(0,3,4),fontsize=29)
                plt.ylim(-0.1,3.1)
                plt.xlim(-1,21)

                ax1.xaxis.set_minor_locator(x_minor_locator)
                ax1.yaxis.set_minor_locator(y_minor_locator)
                ax1.tick_params(axis='both', which='major', pad=8)
                ax1.tick_params(axis='both', which='both', pad=8, left='on', right='on',top='on',bottom='on')
                max_arg = 4000 + np.argmax(data[4000:, 3])

                plt.plot(data[:,1]-timeshift, np.log10(data[:,3]/power_ref), color=colors[i], lw=3.5)
                plt.scatter(data[max_arg,1]-timeshift, np.log10(data[max_arg, 3]/power_ref), marker=shape[2], s=170,
                        facecolor=colors[i], edgecolor='grey',linewidth=1.5, zorder=zorders2[i])

        plt.xticks(np.linspace(-0,35,8),fontsize=29)
        plt.xlim(-1,21)

        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax1.get_xticklabels(), visible=False)

        plt.ylabel('$\\rm log $ $L^{\infty}_{\\rm{s}}/L_{\\rm{h0}}^{\infty}$', fontsize=29, labelpad=24)
        plt.xlabel('$t, \\thinspace \\thinspace \\rm{yr}$', fontsize=29)

        plt.subplot(3,1,1)

        plt.text(2,0.4,'$\\rm HEATER$',fontsize=22)
        plt.text(12,2.16,'$\\rm \\rho_{1} = 10^{12} \\thinspace \\thinspace g \\thinspace cm^{-3}$',fontsize=26)
        plt.ylabel('$\\rm log $ $L^{\infty}_{\\rm{h}}/L_{\\rm{h0}}^{\infty}$', fontsize=29, labelpad=49)

        plt.subplot(3,1,2)
        plt.text(8,-0.55,'$\\rm SURFACE,$ $M = 1.40 \\thinspace \\rm M \odot$',fontsize=22)
        plt.ylabel('$\\rm log $ $L^{\infty}_{\\rm{s}}/L_{\\rm{h0}}^{\infty}$', fontsize=29)

        plt.subplot(3,1,3)
        plt.text(8,-1.1,'$\\rm SURFACE,$ $M = 1.85 \\thinspace \\rm M \odot$',fontsize=22)

        plt.savefig('fig2.eps',format='eps')
        plt.show()

#show_curves_1()
show_curves_2()

def temperature_profile_profile(rho=2, Pow=4, dt=0):

    time_points = np.array([0.01, 0.1, 0.5, 1., 2., 3., 4., 5, 6., 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22])
    Tcrit = np.load('Tcrit_data.npy')

    plot_style()

    num = dt + 8*Pow + 2*rho
    plt.rc('text', usetex=True)
    plt.rcParams.update({'figure.autolayout': True})

    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    ax1 = plt.subplot(2,1,1)

    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)

    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)

    ax1.xaxis.set_minor_locator(x_minor_locator)
    ax1.yaxis.set_minor_locator(y_minor_locator)
    ax1.tick_params(axis='both', which='both', pad=8, left='on', right='on',top='on',bottom='on')

    plt.yticks(np.array([8,8.2,8.4,8.6,8.8,9]),fontsize=29)
    plt.xticks(np.array([9,10,11,12,13,14,15]),fontsize=29)

    if(rho==1):
        plt.ylabel('$\\rm log \\thinspace$$ \\widetilde{T} \\thinspace \\thinspace \\thinspace \\rm [K]$',fontsize=29)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.xlim(9,np.log10(1.5e14))

    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,2]*1e9), facecolor=colors2[2], interpolate=True, alpha=1)
    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,1]*1e9), facecolor=colors2[1], interpolate=True, alpha=1)
    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,3]*1e9), facecolor=colors2[3], interpolate=True, alpha=1)

    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,1]*1e9), color=colors[1], lw=2)
    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,2]*1e9), color=colors[2], lw=2)
    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,3]*1e9), color=colors[3], lw=2)

    if rho==2:
        plt.ylim(8.17,8.73)
        plt.text(9.1,8.32,'$\\rho_{1} = 10^{12} \\thinspace \\thinspace \\rm g \\thinspace cm^{-3}$', fontsize=24)
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
        plt.legend(loc='upper left',fontsize=21, frameon=False)
    else:
        plt.ylim(8.17,8.83)
        plt.text(9.1,8.71,'$\\rho_{1} = 10^{11} \\thinspace \\thinspace \\rm g \\thinspace cm^{-3}$', fontsize=24)
        plt.text(9.1,8.63,'$M = 1.40 \\thinspace \\rm M \\odot$',fontsize=24)

    for time in [1,3]:
        for i in range(3,-1,-1):
            idx = i
            if time==1 and dt==0:
                data = np.loadtxt('temperature_lm/temperature_' + names[i] + '_' + str(num+1) + '.dat')
                plt.plot(np.log10(data[:,0]), np.log10(data[:,time]), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
            else:
                data = np.loadtxt('temperature_lm/temperature_' + names[i] + '_' + str(num) + '.dat')
                plt.plot(np.log10(data[:,0]), np.log10(data[:,time]), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])

    if rho==1:
        plt.text(13, 8.52, '$T < T_{\\rm cr}$', fontsize=25)
        plt.text(9.5,8.25, '$t=0.0 \\thinspace \\thinspace \\rm yr$',fontsize=23)
        plt.text(10.45,8.42, '$t=0.5 \\thinspace \\thinspace \\rm yr$',fontsize=23)

        plt.text(11.82,8.67,'GIPSF',fontsize=17,rotation=86)
        plt.text(12.23,8.68,'AWP2',fontsize=17,rotation=84)
        plt.text(13.03,8.75,'SCLBL',fontsize=17,rotation=83)

        plt.text(13.72,8.42,'GIPSF',fontsize=17,rotation=-90)
        #plt.text(14.18,8.47,'AWP2, SCLBL',fontsize=14,rotation=-89)
        
    elif rho==2:
        plt.text(13.1, 8.49, '$T < T_{\\rm cr}$', fontsize=25)
        plt.text(12.1,8.23, '$t=0.0 \\thinspace \\thinspace \\rm yr$',fontsize=23)
        plt.text(11.07,8.52, '$t=0.5 \\thinspace \\thinspace \\rm yr$',fontsize=23)
        
        plt.text(11.82,8.67,'GIPSF',fontsize=17,rotation=86)
        plt.text(12.23,8.69,'AWP2',fontsize=17,rotation=84)
        plt.text(13.00,8.69,'SCLBL',fontsize=17,rotation=83)

        plt.text(13.74,8.42,'GIPSF',fontsize=17,rotation=-90)
        #plt.text(14.20,8.43,'AWP2, SCLBL',fontsize=14,rotation=-89)
    idx = 0

    ax2 = plt.subplot(2,1,2)

    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)

    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)

    ax2.xaxis.set_minor_locator(x_minor_locator)
    ax2.yaxis.set_minor_locator(y_minor_locator)
    ax2.tick_params(axis='both', which='both', pad=8, left='on', right='on',top='on',bottom='on')

    if rho==2:
        plt.yticks(np.array([8,8.1,8.2,8.3,8.4,8.5,8.6]),fontsize=29)
    else:
        plt.yticks(np.array([8,8.2,8.4,8.6,8.8,9]),fontsize=29)
        plt.ylabel('$\\rm log \\thinspace$$ \\widetilde{T} \\thinspace \\thinspace \\thinspace \\rm [K]$',fontsize=29)
    plt.xlabel('$\\rm log \\thinspace$$\\rho \\thinspace \\thinspace \\thinspace \\rm [g \\thinspace cm^{-3}]$',fontsize=29)
    plt.xticks(np.array([9,10,11,12,13,14,15]),fontsize=29)
    plt.xlim(9, np.log10(1.5e14))

    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,2]*1e9), facecolor=colors2[2], interpolate=True, alpha=1)
    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,1]*1e9), facecolor=colors2[1], interpolate=True, alpha=1)
    plt.fill_between(Tcrit[:,0], np.ones(len(Tcrit[:,0])), np.log10(Tcrit[:,3]*1e9), facecolor=colors2[3], interpolate=True, alpha=1)

    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,1]*1e9), color=colors[1], lw=2)
    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,2]*1e9), color=colors[2], lw=2)
    plt.plot(Tcrit[:,0], np.log10(Tcrit[:,3]*1e9), color=colors[3], lw=2)

    if rho==2:
        plt.ylim(8.17,8.53)
    else:
        plt.ylim(8.17,8.68)

    for time in [6,18]:
        for i in range(3,-1,-1):
            idx = i
            if time==1 and dt==0:
                data = np.loadtxt('temperature_lm/temperature_' + names[i] + '_' + str(num+1) + '.dat')
                plt.plot(np.log10(data[:,0]), np.log10(data[:,time]), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
            else:
                data = np.loadtxt('temperature_lm/temperature_' + names[i] + '_' + str(num) + '.dat')
                plt.plot(np.log10(data[:,0]), np.log10(data[:,time]), color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
        idx += 1

    if rho==1:
        plt.text(13, 8.47, '$T < T_{\\rm cr}$', fontsize=25)
        plt.text(9.66,8.35, '$t=15.0 \\thinspace \\thinspace \\rm yr$',fontsize=23)
        plt.text(9.5,8.54, '$t=3.0 \\thinspace \\thinspace \\rm yr$',fontsize=23)
      
    else:
        plt.text(13, 8.40, '$T < T_{\\rm cr}$', fontsize=25)
        plt.text(9.5,8.4, '$t=3.0 \\thinspace \\thinspace \\rm yr$',fontsize=23)
        plt.text(11.71,8.29, '$t=15.0 \\thinspace \\thinspace \\rm yr$',fontsize=23)
    if rho==1:
        plt.savefig('profile_1.eps', format='eps')
        #plt.savefig('profile_1.jpeg', format='jpeg')
    else:
        plt.savefig('profile_2.eps', format='eps')
        #plt.savefig('profile_2.jpeg', format='jpeg')
    plt.show()

#temperature_profile_profile(rho=2)
#temperature_profile_profile(rho=1)


def timedelay(Pow=5, dt=0):

    rho_values = np.array([3.16227766e+10, 1e11, 3.16227766e+11, 1e12, 3.16227766e+12, 1e13])
    # 0 -- 10, 1 -- 50, 2 -- 100, 3 -- 200, 4 -- 500, 5 -- 1000
    rho_change = np.array([0,1,-1,2,-2,3])
    plot_style()

    y_sample = np.zeros((4, 200))
    rho_sample = np.linspace(10.5,13, 200)

    for Pow in [4]:

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
            y_sample[i,:] = lm(rho_sample)


    alpha = np.array([0.5, 0.5, 0.4, 0.5])
    zorder = np.array([4,2,2,3])

    for idx in range(4):
        plt.plot(rho_sample, y_sample[idx,:], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])

    plt.xticks([10.5,11,11.5,12,12.5,13],fontsize=29)
    plt.xlim(11.4,13)

    plt.yticks([0,2,4,6,8,10],fontsize=29)
    plt.ylim(1,10)

    plt.xlabel('$\\rm log \\thinspace$$\\rho_{1} \\thinspace \\thinspace  \\thinspace \\rm [g \\thinspace cm^{-3}]$',fontsize=29)
    plt.ylabel('$\Delta t_{\\rm r},\\thinspace  \\thinspace \\rm yr$',fontsize=29)
    plt.text(11.44,6.6,'$M = 1.40 \\thinspace \\rm M \\odot$',fontsize=24)

    plt.legend(loc='upper left', fontsize=23,frameon=False)
    plt.savefig('timedelay.eps', format='eps')
    plt.show()

#timedelay(Pow=5, dt=0)



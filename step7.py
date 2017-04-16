import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from scipy import integrate
from scipy import stats
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator

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


MSun = 1.98892e33
LSun = 3.846e33
c = 2.99792458e10
G = 6.67259e-8
from_ev_to_erg = 1.6021766208e-12
from_yr_to_sec = 31536000
Mn0 = 1.6749286e-24
Mp0 = 1.6726231e-24
Mb = (Mn0 + Mp0) / 2

M_dot = 1e-9 * MSun / from_yr_to_sec

data = np.loadtxt('heaters_data.dat')
rho_1 = data[:, 0]
rho_2 = data[:, 0] * (1 + data[:, 2])
y = data[:, 1]


_star_model = np.loadtxt('data/model_rho9.dat', skiprows=2)
_nb = interpolate.interp1d(np.log10(_star_model[:, 3]), _star_model[:, 5], kind='linear')
_ne = interpolate.interp1d(np.log10(_star_model[:, 3]), _star_model[:, 6], kind='linear')
_nm = interpolate.interp1d(np.log10(_star_model[:, 3]), _star_model[:, 7], kind='linear')
_rho = interpolate.interp1d(np.log10(_star_model[:, 1] * 1e5), np.log10(_star_model[:, 3]), kind='linear')
_pressure = interpolate.interp1d(np.log10(_star_model[:, 1] * 1e5), np.log10(_star_model[:, 2]), kind='linear')
_radii = interpolate.interp1d(np.log10(_star_model[:, 3]), np.log10(_star_model[:, 1] * 1e5), kind='linear')
_mass = interpolate.interp1d(np.log10(_star_model[:, 1] * 1e5), np.log10(_star_model[:, 0] * MSun), kind='linear')
_Phi = interpolate.interp1d(np.log10(_star_model[:, 1] * 1e5), _star_model[:, 4], kind='linear')

def radii(a):  # radius(density)
    return np.power(10, _radii(np.log10(a)))

def mass(a):  # mass(radius)
    return np.power(10, _mass(np.log10(a)))

def relativity_sqrt(r):
    return np.sqrt(1 - 2*G*mass(r)/(( c ** 2 )*r))

def rho(a):  # density(radius)
    return np.power(10, _rho(np.log10(a)))

H_guess = 1e17

def f2(x,a,b,c,d):
    r_0 = radii(1e9)
    output = np.zeros_like(x)
    for i in range(len(output)):
        r = np.linspace(r_0, x[i], 200)
        output[i] = H_guess * a * np.abs(np.trapz((stats.norm.pdf(x=r, loc=b, scale=c)+d)*dV(r)*np.exp(2*_Phi(np.log10(r))), r))
    return output

dV = lambda r: 4*np.pi*r*r/relativity_sqrt(r)
V = lambda r1,r2: integrate.quad(dV,r1,r2)[0]

def power():

    r_1 = radii(rho_1)
    r_2 = radii(rho_2)

    Q_total = (y * 1e3 * from_ev_to_erg) * M_dot / Mb
    Q_total_cum = np.cumsum(Q_total)

    param, mtrx = optimize.curve_fit(f2, r_1, Q_total_cum, p0=[6.39198181e+04,   1.20413015e+06,   4.38291763e+03,   2.60061730e-06])
    print(param)

    #param, mtrx = optimize.curve_fit(f2, r_1, Q_total_cum, p0=[  6.64642971e+04,   1.21161931e+06,   2.68584250e+03,   5.71917036e-06])
    #print(param)

    rho_sample = np.linspace(9,14,1000)
    print(param[3] * H_guess * param[0]/1e16)

    r_0 = 1.20413015e+06 # cm
    r_sample = radii(np.power(10,rho_sample))
    sigma_g = 4.38291763e+03 # cm
    H_f = lambda r: 1.66225679807e16 + 6.39198181e+04 * 1e17 * np.exp(-(r-r_0)**2/sigma_g**2 / 2)/np.sqrt(2*np.pi*sigma_g**2 )

    #r_0 = 1.21161929e+06 # cm
    #r_sample = radii(np.power(10,rho_sample))
    #sigma_g = 2.68581589e+03  # cm
    #H_f = lambda r: 3.80131342314e16 +  6.64638989e+04 * 1e17 * np.exp(-(r-r_0)**2/sigma_g**2 / 2)/np.sqrt(2*np.pi*sigma_g**2 )

    plot_style()
    plt.plot(rho_sample,H_f(r_sample)/1e17,ls='-',lw=3,c='r')
    plt.xlabel('$\\rm log \\thinspace \\rho $ $\\rm g \\thinspace cm^{-3}$',fontsize=22)
    plt.ylabel('$\\rm Q_{h}/10^{17}$ $\\rm erg \\thinspace s^{-1} \\thinspace cm^{-3}$',fontsize=22)
    plt.xticks([9,10,11,12,13,14],fontsize=20)
    plt.yticks([0,1,2,3,4,5,6],fontsize=20)
    plt.savefig('fig2.pdf',format='pdf')
    plt.show()

    check = np.zeros_like(rho_sample)
    for i in range(len(check)):
        check[i] = np.abs(np.trapz(H_f(r_sample[:i+1])*np.exp(2*_Phi(np.log10(r_sample[:i+1]))) * dV(r_sample[:i+1]), r_sample[:i+1]))

    plot_style()
    plt.step(np.log10(data[:,0]), Q_total_cum/1e34,lw=3)
    plt.plot(rho_sample, check/1e34,lw=3,color='red')
    plt.xlabel('$\\rm log \\thinspace \\rho $ $\\rm g \\thinspace cm^{-3}$',fontsize=22)
    plt.ylabel('$\\rm Q_{< \\rho}/10^{34}$ $\\rm erg \\thinspace s^{-1}$',fontsize=22)
    plt.xticks([9,10,11,12,13,14],fontsize=20)
    plt.yticks([0,2,4,6,8,10,12],fontsize=20)
    plt.savefig('fig1.pdf',format='pdf')
    plt.show()


def plot1():
    plot_style()

    t_MAXI = np.array([5.5, 16.1, 23.3, 31.9, 51.1, 85.4, 104.5, 134.9 ,  150.9,  291.8, 497.1])
    T_MAXI_model_II = np.array([308,  307,  298,  365,  276, 320, 252.6, 243.6 ,  241,  208,   184.3])
    err_MAXI = np.array([11, 4, 3, 7, 3, 5, 1, 1, 2, 2, 1.1])

    T_XTE = np.array([164.2, 159.5, 156.8, 150.0,  129.1,  159.3,  136.0,  126.3,  125.4,  129.6,  124.0,  123.9, 124.1])
    t_XTE = np.array([ 2.77, 10.63, 16.06, 49.31, 174.15, 225.54, 298.12, 430.89, 539.90, 592.50, 652.44, 705.20, 795.45])
    err_XTE = np.array([3.6,   2.5,   1.3,   1.2,    4.7,    2.0,    2.0,    3.1,    1.5,    2.2,    2.2,    2.0, 1.7])

    k_b = 8.617330350e-5
    config = np.loadtxt('data/config.dat')
    for i,col in zip(range(0,2),['green','orange']):
        print(config[i,:])
        #data2 = np.loadtxt('output/cooling_SF0_' + str(i) + '.dat')
        data = np.loadtxt('output/cooling_GIPSF_' + str(i) + '.dat')
        plt.plot((data[:, 1]-4e3)*365, data[:, 0]*k_b, lw=3, color=col, ls='-')
        #plt.plot((data2[:, 1]-4e3)*365, data2[:, 0]*k_b, lw=3, color='red', ls='--')

    plt.scatter(t_MAXI, T_MAXI_model_II, s=100, color='magenta', marker='o',label='MAXI J0556–332')
    plt.errorbar(x=t_MAXI, y=T_MAXI_model_II, yerr=err_MAXI, color='magenta', fmt=' ')

    plt.scatter(t_XTE, T_XTE, s=100, color='black', marker='^',label='XTE J1701–462')
    plt.errorbar(x=t_XTE, y=T_XTE, yerr=err_XTE, color='black', fmt=' ')

    plt.plot([1,2], [1,2], lw=3, color='black', ls='-', label='cooling code')

    plt.legend(loc='upper right',fontsize=20,scatterpoints=1)
    plt.xscale('log')
    plt.xticks([1,10,100,1000],fontsize=20)
    plt.yticks([100,150,200,250,300,350],fontsize=20)
    plt.ylim(100, 350)
    plt.xlabel('$\\rm Days \\thinspace  \\thinspace  since  \\thinspace  \\thinspace  '
               'end \\thinspace  \\thinspace  of  \\thinspace  \\thinspace outburst$',fontsize=22)
    plt.ylabel('$\\rm kT^{\infty}_{s}$ $\\rm eV$',fontsize=22)

    plt.xlim(1, 1e3)
    plt.savefig('fig3.pdf',format='pdf')

    plt.show()


def plot2():

    plt.rcParams.update({'figure.autolayout': True})
    #plt.tight_layout()
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    ax = plt.subplot(2,1,1)
    x_minor_locator = AutoMinorLocator(5)
    y_minor_locator = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    ax.xaxis.set_minor_locator(x_minor_locator)
    ax.yaxis.set_minor_locator(y_minor_locator)

    for i,col in zip(range(0,2),['green','orange']):
        data = np.loadtxt('output/cooling_GIPSF_' + str(i) + '.dat')
        plt.plot((data[:, 1]-4e3)*365, np.log10(data[:, 3]), lw=3, color=col, ls='-')
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.ylim(30,37)
    plt.yticks([30,31,32,33,34,35,36,37],fontsize=18)
    plt.ylabel('$\\rm log \\thinspace L^{\infty}_{h}$ $\\rm erg \\thinspace s^{-1}$', fontsize=22)
    plt.xlim(-500,2000)

    axx = plt.subplot(2,1,2)
    x_minor_locator2 = AutoMinorLocator(5)
    y_minor_locator2 = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    axx.xaxis.set_minor_locator(x_minor_locator2)
    axx.yaxis.set_minor_locator(y_minor_locator2)

    t_MAXI = np.array([5.5, 16.1, 23.3, 31.9, 51.1, 85.4, 104.5, 134.9 ,  150.9,  291.8, 497.1])
    T_MAXI_model_II = np.array([308,  307,  298,  365,  276, 320, 252.6, 243.6 ,  241,  208,   184.3])
    err_MAXI = np.array([11, 4, 3, 7, 3, 5, 1, 1, 2, 2, 1.1])

    T_XTE = np.array([164.2, 159.5, 156.8, 150.0,  129.1,  159.3,  136.0,  126.3,  125.4,  129.6,  124.0,  123.9, 124.1])
    t_XTE = np.array([ 2.77, 10.63, 16.06, 49.31, 174.15, 225.54, 298.12, 430.89, 539.90, 592.50, 652.44, 705.20, 795.45])
    err_XTE = np.array([3.6,   2.5,   1.3,   1.2,    4.7,    2.0,    2.0,    3.1,    1.5,    2.2,    2.2,    2.0, 1.7])

    k_b = 8.617330350e-5
    for i,col in zip(range(0,2),['green','orange']):
        data = np.loadtxt('output/cooling_GIPSF_' + str(i) + '.dat')
        plt.plot((data[:, 1]-4e3)*365, data[:, 0]*k_b, lw=3, color=col, ls='-')

    plt.scatter(t_MAXI, T_MAXI_model_II, s=100, color='magenta', marker='o',label='MAXI J0556–332')
    plt.errorbar(x=t_MAXI, y=T_MAXI_model_II, yerr=err_MAXI, color='magenta', fmt=' ')

    plt.scatter(t_XTE, T_XTE, s=100, color='black', marker='^',label='XTE J1701–462')
    plt.errorbar(x=t_XTE, y=T_XTE, yerr=err_XTE, color='black', fmt=' ')

    plt.xlim(-500,2000)
    plt.xticks([-500,0,500,1000,1500,2000],fontsize=18)
    plt.yticks([100,150,200,250,300,350],fontsize=18)
    plt.ylim(100, 350)
    plt.xlabel('$\\rm Days \\thinspace  \\thinspace  since  \\thinspace  \\thinspace  '
               'end \\thinspace  \\thinspace  of  \\thinspace  \\thinspace outburst$',fontsize=22)
    plt.ylabel('$\\rm kT^{\infty}_{s}$ $\\rm eV$',fontsize=22)
    plt.savefig('fig4.pdf',format='pdf')
    plt.show()

plot1()
#plot2()

def plot3():

    plot_style()
    data = np.loadtxt('heaters_data.dat')
    rho_1 = data[:, 0]
    y = data[:, 1]
    y = np.cumsum(y)
    plt.xlabel('$\\rm log \\thinspace \\rho $ $\\rm g \\thinspace cm^{-3}$',fontsize=22)
    plt.ylabel('$\\rm Integrated \\thinspace \\thinspace heat \\thinspace \\thinspace per \\thinspace \\thinspace barion$ $\\rm MeV$',fontsize=22)
    plt.step(np.log10(rho_1), y/1e3,lw=3)
    plt.yticks([0,0.5,1.0,1.5,2.0],fontsize=18)
    plt.xticks([9,10,11,12,13,14],fontsize=18)
    plt.savefig('fig5.pdf',format='pdf')
    plt.show()

# plot3()


def plot4():

    plot_style()

    for i,col,lb in zip(range(0,2),['green','orange'],['MAXI J0556–332','XTE J1701–462']):
        data = np.loadtxt('output/cooling_GIPSF_' + str(i) + '.dat')
        plt.plot((data[:, 1]-4e3)*365, np.log10(data[:, 5]), lw=3, color=col, ls='-',label=lb)

    plt.xlim(-500,2000)
    plt.xticks([-500,0,500,1000,1500,2000],fontsize=18)
    plt.yticks([33.5,34,34.5,35,35.5,36],fontsize=18)
    plt.ylim(33.5,36)
    plt.xlabel('$\\rm Days \\thinspace  \\thinspace  since  \\thinspace  \\thinspace  '
               'end \\thinspace  \\thinspace  of  \\thinspace  \\thinspace outburst$',fontsize=22)
    plt.ylabel('$\\rm log \\thinspace L_{s}^{\infty} \\thinspace erg \\thinspace s^{-1}$',fontsize=22)
    plt.legend(loc='upper right',fontsize=20)
    plt.savefig('fig6.pdf',format='pdf')
    plt.show()

#plot4()


def f2(x,a,b,c,d):
    r_0 = radii(1e9)
    output = np.zeros_like(x)
    for i in range(len(output)):
        r = np.linspace(r_0, x[i], 200)
        output[i] = H_guess * a * np.abs(np.trapz((stats.norm.pdf(x=r, loc=b, scale=c)+d)*dV(r)*np.exp(2*_Phi(np.log10(r))), r))
    return output
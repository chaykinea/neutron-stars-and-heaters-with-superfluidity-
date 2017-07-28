import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from scipy import integrate
from scipy import stats
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc

names = np.array(['SF0','AWP2', 'GIPSF', 'SCLBL'])
labels = np.array(['No SF','AWP2', 'GIPSF', 'SCLBL'])
colors = np.array(['black','darkblue','red','darkorange'])
shape = np.array(['s','^','o','d'])
line_thickness = np.array([4.8, 2.7,2.2,3.2])
dashes = np.array([[2,1e-15],[4,8],[20,3],[9,6]])
order = np.array([1,4,2,3])

order2 = np.array([1,2,3,4])
colors2 = np.array(['green','orange','red','blue'])
dashes2 = np.array([[2,1e-15],[7,1],[3,2],[1.5,3]])
dashes2 = np.array([[2,1e-15],[7,1e-15],[3,1e-15],[1.5,1e-15]])
line_thickness2 = np.array([6,4.7,3.4,2.1])

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
Mn0 = 1.6749286e-24
Mp0 = 1.6726231e-24
Mb = (Mn0 + Mp0) / 2

M_dot = 1e-9 * MSun / from_yr_to_sec # IMPORTANT

data = np.loadtxt('heaters_data.dat')
rho_1 = data[:, 0]
rho_2 = data[:, 0] * (1 + data[:, 2])
y = data[:, 1]

_star_model = np.loadtxt('data/BSK21_1.40.dat', skiprows=2)
model_param = 1.40 # IMPORTANT
_nb = interpolate.interp1d(np.log10(_star_model[:, 3]), _star_model[:, 5], kind='linear', fill_value='extrapolate')
_ne = interpolate.interp1d(np.log10(_star_model[:, 3]), _star_model[:, 6], kind='linear', fill_value='extrapolate')
_nm = interpolate.interp1d(np.log10(_star_model[:, 3]), _star_model[:, 7], kind='linear', fill_value='extrapolate')
_rho = interpolate.interp1d(np.log10(_star_model[:, 1] * 1e5), np.log10(_star_model[:, 3]), kind='linear', fill_value='extrapolate')
_pressure = interpolate.interp1d(np.log10(_star_model[:, 1] * 1e5), np.log10(_star_model[:, 2]), kind='linear', fill_value='extrapolate')
_radii = interpolate.interp1d(np.log10(_star_model[:, 3]), np.log10(_star_model[:, 1] * 1e5), kind='linear', fill_value='extrapolate')
_mass = interpolate.interp1d(np.log10(_star_model[:, 1] * 1e5), np.log10(_star_model[:, 0] * MSun), kind='linear', fill_value='extrapolate')
_Phi = interpolate.interp1d(np.log10(_star_model[:, 1] * 1e5), _star_model[:, 4], kind='linear', fill_value='extrapolate')

t_MAXI = np.array([5.5, 16.1, 23.3, 31.9, 51.1, 85.4, 104.5, 134.9 ,  150.9,  291.8, 497.1])
T_MAXI_model_II = np.array([308,  307,  298,  365,  276, 320, 252.6, 243.6 ,  241,  208,   184.3])
err_MAXI = np.array([11, 4, 3, 7, 3, 5, 1, 1, 2, 2, 1.1])

T_XTE = np.array([164.2, 159.5, 156.8, 150.0,  129.1,  159.3,  136.0,  126.3,  125.4,  129.6,  124.0,  123.9, 124.1])
t_XTE = np.array([ 2.77, 10.63, 16.06, 49.31, 174.15, 225.54, 298.12, 430.89, 539.90, 592.50, 652.44, 705.20, 795.45])
err_XTE = np.array([3.6,   2.5,   1.3,   1.2,    4.7,    2.0,    2.0,    3.1,    1.5,    2.2,    2.2,    2.0, 1.7])

T_KS = np.array([103, 88, 76, 72, 70, 67, 70, 63.35])
t_KS = np.array([51995.1, 52165.7, 52681.6, 52859.5, 53430.5, 53500.4, 53525.4, 54978.86]) - 51930.5
err_KS = np.array([2.5, 2, 3, 4, 3 , 7.5, 3, 3])

T_IGR = np.array([99.7, 91.5, 89.2, 84.8, 88.5, 84.6, 82.8])
t_IGR = np.array([55609, 55680.5, 55810.5, 56060, 56187.5, 56228, 56340]) - 55556
err_IGR = np.array([1.6, 1.5, 1.5, 1.5, 1.9, 2.0, 1.2])

t_EXO = np.array([54755.5, 54776, 54886, 54908, 54992, 55013, 55306, 55364, 55489, 55745, 56505]) - 54714
T_EXO = np.array([129.1, 126.1, 122.6, 120.0, 117.8, 115.5, 116.8, 116.2, 115.4, 117.6, 109.9])
err_EXO = np.array([2.3, 2.2, 2.6, 2.0, 2.5, 2.2, 2.5, 1.9, 2.2, 2.2, 2.0])

t_MXB = np.array([52197.8, 52563.2, 52712.2, 52768.9, 53560.0, 53576.7, 54586.37]) - 52159.5
T_MXB = np.array([121, 85, 77 ,73, 58, 54, 56.16])
err_MXB = np.array([2, 2, 1, 1.5, 3.5, 4.5, 1.5])

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
dV2 = lambda r: 4*np.pi*r*r
V = lambda r1,r2: integrate.quad(dV,r1,r2)[0]


def power():

    r_1 = radii(rho_1)
    r_2 = radii(rho_2)

    Q_total = (y * 1e3 * from_ev_to_erg) * M_dot / Mb
    Q_total_cum = np.cumsum(Q_total)
    Q_total_cum = np.append(Q_total_cum,Q_total_cum[-1])
    data_plot = np.append(data[:,0],10**15)

    rho_sample = np.linspace(9,14,1000)
    r_sample = radii(np.power(10, rho_sample))
    

    #param, mtrx = optimize.curve_fit(f2, r_1, Q_total_cum, p0=[6.49105274e+04,   1.21269574e+06,   3.71309773e+03,   3.13087696e-06])  # M 1.4
    #print(param)
    #print(param[3] * H_guess * param[0]/1e16)

    r_0 = np.array(
        [1.20413015e+06, 1.21000074e+06, 1.21203074e+06, 1.21269574e+06, 1.21353249e+06, 1.21448985e+06, 1.21486413e+06,
         1.21161929e+06])
    sigma_g = np.array(
        [4.38291763e+03, 3.97737319e+03, 3.78693639e+03, 3.71309773e+03, 3.60523861e+03, 3.43057815e+03, 3.26194154e+03,
         2.68581589e+03])
    const_1 = np.array(
        [1.66225679807e16, 1.86808408312e16, 1.98349866253e16, 2.03226874478e16, 2.10820556759e16, 2.2441279112e16,
         2.39298439565e16, 3.80131342314e16])
    const_2 = np.array(
        [6.39198181e+04, 6.44242849e+04, 6.47562763e+04, 6.49105274e+04, 6.51666607e+04, 6.56611002e+04, 6.62462312e+04,
         6.64638989e+04])
    model_param_arr = np.array([1.40, 1.50, 1.55, 1.57, 1.60, 1.65, 1.70, 1.85])

    r_0_f = interpolate.interp1d(model_param_arr, r_0)
    sigma_g_f = interpolate.interp1d(model_param_arr, sigma_g)
    const_1_f = interpolate.interp1d(model_param_arr, const_1)
    const_2_f = interpolate.interp1d(model_param_arr, const_2)

    H_f = lambda r: const_1_f(model_param) + const_2_f(model_param) * 1e17 * np.exp(
        -(r - r_0_f(model_param)) ** 2 / sigma_g_f(model_param) ** 2 / 2) / np.sqrt(
        2 * np.pi * sigma_g_f(model_param) ** 2)  # M_dot = 1e-9 M_solar/yr

    #plot_style()
    #plt.plot(rho_sample,H_f(r_sample)/1e17,ls='-',lw=3,c='r')
    #plt.xlabel('$\\rm log \\thinspace \\rho $ $\\rm g \\thinspace cm^{-3}$',fontsize=22)
    #plt.ylabel('$\\rm Q_{h}/10^{17}$ $\\rm erg \\thinspace s^{-1} \\thinspace cm^{-3}$',fontsize=22)
    #plt.xticks([9,10,11,12,13,14],fontsize=20)
    #plt.yticks([0,1,2,3,4,5,6],fontsize=20)
    #plt.savefig('fig2.pdf',format='pdf')
    #plt.show()

    check = np.zeros_like(rho_sample)
    for i in range(len(check)):
        check[i] = np.abs(np.trapz(H_f(r_sample[:i+1])*np.exp(2*_Phi(np.log10(r_sample[:i+1]))) * dV(r_sample[:i+1]), r_sample[:i+1]))
    
    #aa = np.argmin(np.abs(np.power(10,rho_sample)-3e12))
    #bb = np.argmin(np.abs(np.power(10,rho_sample)-8e11)) 
    #print(10**rho_sample[aa])
    #print(10**rho_sample[bb])
    #aa_int = np.abs(np.trapz(H_f(r_sample[:aa+1])*np.exp(2*_Phi(np.log10(r_sample[:aa+1]))) * dV(r_sample[:aa+1]), r_sample[:aa+1]))
    #bb_int = np.abs(np.trapz(H_f(r_sample[:bb+1])*np.exp(2*_Phi(np.log10(r_sample[:bb+1]))) * dV(r_sample[:bb+1]), r_sample[:bb+1]))
    #print(aa_int)
    #print(bb_int)
    #print(3e1*aa_int)
    #print(6e1*bb_int)
    
    #plot_style()
    #plt.step(np.log10(data[:,0]), Q_total_cum/1e34,lw=3)
    #plt.plot(rho_sample, check/1e34,lw=3,color='red')
    #plt.xlabel('$\\rm log \\thinspace \\rho $ $\\rm g \\thinspace cm^{-3}$',fontsize=24)
    #plt.ylabel('$\\rm Q_{< \\rho}/10^{34}$ $\\rm erg \\thinspace s^{-1}$',fontsize=24)
    #plt.xticks([9,10,11,12,13,14],fontsize=20)
    #plt.yticks([0,2,4,6,8,10,12],fontsize=20)
    #plt.xlim(9,14)
    #plt.ylim(0,12)
    #plt.savefig('fig1.pdf',format='pdf')
    #plt.show()

    num_mass = np.zeros(4)
    plot_style()

    for rho_max,idx,lb in zip([1e13,10**12.5,1e12,1e11],[3,2,1],['$\\rho_{\\rm acc} = 10^{13} \\rm \\thinspace \\thinspace g \\thinspace cm^{-3}$',
        '$\\rho_{\\rm acc} = 10^{12.5} \\rm \\thinspace \\thinspace g \\thinspace cm^{-3}$',
        '$\\rho_{\\rm acc} = 10^{12} \\rm \\thinspace \\thinspace g \\thinspace cm^{-3}$',
        '$\\rho_{\\rm acc} = 10^{11} \\rm \\thinspace \\thinspace g \\thinspace cm^{-3}$']):
        num = np.argmin(np.abs(np.power(10,rho_sample)-rho_max))
        print(rho_sample[num])
        num_int_max = np.abs(np.trapz(H_f(r_sample[:num+1])*np.exp(2*_Phi(np.log10(r_sample[:num+1]))) * dV(r_sample[:num+1]), r_sample[:num+1]))
        num_mass[idx] = np.abs(np.trapz(10**rho_sample[:num+1] * dV2(r_sample[:num+1]), r_sample[:num+1])) / MSun
        print(num_mass[idx])
        num_int = np.ones(len(rho_sample)) * num_int_max
        for int_idx in range(0,num+1):
            num_int[int_idx] = np.abs(np.trapz(H_f(r_sample[:int_idx+1])*np.exp(2*_Phi(np.log10(r_sample[:int_idx+1]))) * dV(r_sample[:int_idx+1]), r_sample[:int_idx+1]))

        plt.plot(rho_sample, num_int/1e34,color=colors2[idx],
                     linewidth=line_thickness2[idx], dashes = (dashes2[idx,0],dashes2[idx,1]),zorder=order2[idx],label=lb)
    
    #plt.text(11.8,0.92,'$\Delta M_{\\rm acc} = 7.1 \\times 10^{-6} \\thinspace \\rm M \\odot$',fontsize=20)
    plt.text(12.1,2.65,'$\Delta M_{\\rm acc} = 6.2 \\times 10^{-5} \\thinspace \\rm M \\odot$',fontsize=18)
    plt.text(12.7,7.85,'$1.7 \\times 10^{-4} \\thinspace \\rm M \\odot$',fontsize=18)
    plt.text(12.885,10.35,'$7.3 \\times 10^{-4} \\thinspace \\rm M \\odot$',fontsize=18)
    plt.text(9.3,8.2,'$\dot{M} = 10^{-8} \\thinspace \\rm M \\odot / yr$',fontsize=22)

    plt.xlabel('$\\rm log \\thinspace$$\\rho \\thinspace \\thinspace \\thinspace \\rm [g \\thinspace cm^{-3}]$',fontsize=24)
    plt.ylabel('$L_{\\rm h}^{\infty}$  $/ (\\rm 10^{35} \\thinspace erg \\thinspace s^{-1})$',fontsize=24)
    plt.step(np.log10(data_plot), Q_total_cum/1e34,linestyle='--',lw=3,color='black',zorder=1,label='$\\rm Haensel \\thinspace \\thinspace and \\thinspace \\thinspace Zdunik \\thinspace \\thinspace (2008)$')
    output = np.vstack([data_plot,Q_total_cum/1e34])
    np.savetxt('curve_from_HZ_2008.dat',output.T,fmt='%1.5e')
    plt.xticks([9,10,11,12,13,14],fontsize=20)
    plt.yticks([0,2,4,6,8,10,12],fontsize=20)
    plt.xlim(9,14)
    plt.ylim(0,12.5)
    plt.legend(loc='upper left',fontsize=20,scatterpoints=1,frameon=False)
    plt.savefig('fig112.eps',format='eps')
    plt.show()

power()

def power2():

    r_1 = radii(rho_1)
    r_2 = radii(rho_2)

    Q_total = (y * 1e3 * from_ev_to_erg) * M_dot / Mb
    Q_total_cum = np.cumsum(Q_total)

    rho_sample = np.linspace(9,14,1000)
    r_sample = radii(np.power(10, rho_sample))

    #param, mtrx = optimize.curve_fit(f2, r_1, Q_total_cum, p0=[6.49105274e+04,   1.21269574e+06,   3.71309773e+03,   3.13087696e-06])  # M 1.4
    #print(param)
    #print(param[3] * H_guess * param[0]/1e16)

    r_0 = np.array(
        [1.20413015e+06, 1.21000074e+06, 1.21203074e+06, 1.21269574e+06, 1.21353249e+06, 1.21448985e+06, 1.21486413e+06,
         1.21161929e+06])
    sigma_g = np.array(
        [4.38291763e+03, 3.97737319e+03, 3.78693639e+03, 3.71309773e+03, 3.60523861e+03, 3.43057815e+03, 3.26194154e+03,
         2.68581589e+03])
    const_1 = np.array(
        [1.66225679807e16, 1.86808408312e16, 1.98349866253e16, 2.03226874478e16, 2.10820556759e16, 2.2441279112e16,
         2.39298439565e16, 3.80131342314e16])
    const_2 = np.array(
        [6.39198181e+04, 6.44242849e+04, 6.47562763e+04, 6.49105274e+04, 6.51666607e+04, 6.56611002e+04, 6.62462312e+04,
         6.64638989e+04])
    model_param_arr = np.array([1.40, 1.50, 1.55, 1.57, 1.60, 1.65, 1.70, 1.85])

    r_0_f = interpolate.interp1d(model_param_arr, r_0)
    sigma_g_f = interpolate.interp1d(model_param_arr, sigma_g)
    const_1_f = interpolate.interp1d(model_param_arr, const_1)
    const_2_f = interpolate.interp1d(model_param_arr, const_2)

    H_f = lambda r: const_1_f(model_param) + const_2_f(model_param) * 1e17 * np.exp(
        -(r - r_0_f(model_param)) ** 2 / sigma_g_f(model_param) ** 2 / 2) / np.sqrt(
        2 * np.pi * sigma_g_f(model_param) ** 2)  # M_dot = 1e-9 M_solar/yr

    check = np.zeros_like(rho_sample)
    for i in range(len(check)):
        check[i] = np.abs(np.trapz(H_f(r_sample[:i+1])*np.exp(2*_Phi(np.log10(r_sample[:i+1]))) * dV(r_sample[:i+1]), r_sample[:i+1]))
    
    num_mass = np.zeros(4)
    plot_style()

    for rho_max,idx,lb in zip([1e13,10**12.5,1e12,1e11],[3,2,1],['$\\rho_{\\rm acc} = 10^{13} \\rm \\thinspace \\thinspace g \\thinspace cm^{-3}$',
'$\\rho_{\\rm acc} = 10^{12.5} \\rm \\thinspace \\thinspace g \\thinspace cm^{-3}$',
'$\\rho_{\\rm acc} = 10^{12} \\rm \\thinspace \\thinspace g \\thinspace cm^{-3}$',
'$\\rho_{\\rm acc} = 10^{11} \\rm \\thinspace \\thinspace g \\thinspace cm^{-3}$']):
        num = np.argmin(np.abs(np.power(10,rho_sample)-rho_max))
        print(rho_sample[num])
        num_int_max = np.abs(np.trapz(H_f(r_sample[:num+1])*np.exp(2*_Phi(np.log10(r_sample[:num+1]))) * dV(r_sample[:num+1]), r_sample[:num+1]))
        num_mass[idx] = np.abs(np.trapz(10**rho_sample[:num+1] * dV2(r_sample[:num+1]), r_sample[:num+1])) / MSun
        times = np.array([-1000, -0.0000001, 0, 365, 365.0000001, 1000]) - 365
        powers = np.ones(len(times)) * num_int_max / 1e34
        powers[:2] = 0
        powers[-2:] = 0
        plt.plot(times, powers,color=colors2[idx],
                     linewidth=line_thickness2[idx], dashes = (dashes2[idx,0],dashes2[idx,1]),zorder=order2[idx],label=lb)
    
    plt.xlabel('$\\rm Days$',fontsize=22)
    plt.ylabel('$L_{\\rm h}^{\infty}$  $(\\rm 10^{35} \\thinspace erg \\thinspace s^{-1})$',fontsize=24)
    plt.xticks([-500,-250,0,250,500,750,1000],fontsize=20)
    plt.yticks([0,2,4,6,8,10,12],fontsize=20)
    plt.xlim(-365,500)
    plt.ylim(0,12.5)
    plt.legend(loc='upper right',fontsize=20,scatterpoints=1,frameon=False)
    plt.savefig('fig1111.pdf',format='pdf')
    plt.show()

power2()

def plot1():
    plot_style()

    k_b = 8.617330350e-5
    config = np.loadtxt('data/config.dat')

    for i in [4,7]:

        print(config[i,:])
        data1 = np.loadtxt('output/cooling_AWP2_' + str(i) + '.dat')
        data2 = np.loadtxt('output/cooling_GIPSF_' + str(i) + '.dat')
        data3 = np.loadtxth('output/cooling_SCLBL_' + str(i) + '.dat')
        data = np.loadtxt('output/cooling_SF0_' + str(i) + '.dat')

        a = 4803

        idx = 0
        plt.plot((data[:a+1, 1]-1.00145e3)*365, data[:a+1, 0]*k_b, color=colors[idx],
             linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],alpha=0.25)
        idx = 1
        plt.plot((data1[:a+1, 1]-1.00145e3)*365, data1[:a+1, 0]*k_b, color=colors[idx],
             linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],alpha=0.25)
        idx = 2
        plt.plot((data2[:a+1, 1]-1.00145e3)*365, data2[:a+1, 0]*k_b, color=colors[idx],
             linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],alpha=0.25)
        idx = 3
        plt.plot((data3[:a+1, 1]-1.00145e3)*365, data3[:a+1, 0]*k_b, color=colors[idx],
             linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],alpha=0.25)
       
        
        idx = 0
        plt.plot((data[a:, 1]-1.00145e3)*365, data[a:, 0]*k_b, color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
        idx = 1
        plt.plot((data1[a:, 1]-1.00145e3)*365, data1[a:, 0]*k_b, color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
        idx = 2
        plt.plot((data2[a:, 1]-1.00145e3)*365, data2[a:, 0]*k_b, color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
        idx = 3
        plt.plot((data3[a:, 1]-1.00145e3)*365, data3[a:, 0]*k_b, color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx])
        

    plt.axvline(365,dashes=(15,5),lw=3,color='black')
    plt.scatter(t_MAXI, T_MAXI_model_II, s=100, color='magenta', marker='o',label='MAXI J0556–332', zorder=6)
    plt.errorbar(x=t_MAXI, y=T_MAXI_model_II, yerr=err_MAXI, color='magenta', fmt=' ', zorder=6)
    
    idx = 0
    plt.plot([1,2], [1,2], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])
    idx = 1
    plt.plot([1,2], [1,2], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])
    idx = 2
    plt.plot([1,2], [1,2], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])
    idx = 3
    plt.plot([1,2], [1,2], color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx,0],dashes[idx,1]),zorder=order[idx],label=labels[idx])
    


    plt.scatter(t_XTE, T_XTE, s=100, color='blue', marker='^',label='XTE J1701–462', zorder=6)
    plt.errorbar(x=t_XTE, y=T_XTE, yerr=err_XTE, color='blue', fmt=' ',zorder=6)

    plt.legend(loc='upper right',fontsize=20,scatterpoints=1,frameon=False)
    #plt.xscale('log')
    #plt.xticks([1,10,100,1000,10000],fontsize=20)
    plt.xticks([-500,-250,0,250,500,750,1000], fontsize=20)
    plt.yticks([125,150,175,200,225,250,275,300],fontsize=20)
    plt.ylim(110, 280)
    plt.xlabel('$\\rm Days \\thinspace  \\thinspace  since  \\thinspace  \\thinspace  '
               'end \\thinspace  \\thinspace  of  \\thinspace  \\thinspace outburst$',fontsize=22)
    plt.ylabel('$\\rm k$$T^{\infty}_{\\rm s}$ $\\rm eV$',fontsize=22)

    plt.xlim(0,1000)
    plt.savefig('fig_MAXI_XTE_linear.pdf',format='pdf')

    plt.show()

#plot1()

def plot_new():
    plot_style()

    k_b = 8.617330350e-5
    config = np.loadtxt('data/config.dat')

    for i in [4]:
        print(config[i,:])
        data2 = np.loadtxt('output/temperature_SF0_' + str(i) + '.dat')
        data = np.loadtxt('output/temperature_GIPSF_' + str(i) + '.dat')
        for jj,col,lbl in zip([1, 4, 5, 9], ['red', 'orange', 'green', 'blue'], [0.0, 0.5, 1, 3]):
            plt.plot(data[:,0],     data[:, jj]/1e8, '-',color=col, lw=3, label=str(lbl) + ' years')
            plt.plot(data2[:, 0], data2[:, jj]/1e8,'--',color=col,lw=3)

    plt.plot([1,2], [1,2], lw=5, color='black', ls='-', label='SF ON')
    plt.plot([1,2], [1,2], lw=5, color='black', ls='--', label='SF OFF')

    plt.legend(loc='upper right',fontsize=18,scatterpoints=1,frameon=False)
    plt.xscale('log')
    plt.xlabel('$\\rm \\rho \\thinspace g \\thinspace cm^{3} $',fontsize=22)
    plt.ylabel('$\\rm \\tilde{T} $ $\\rm 10^{8} K$',fontsize=22)
    plt.xticks([1e9,1e10,1e11,1e12,1e13,1e14,1e15], fontsize=20)
    plt.yticks([1,2,3,4], fontsize=20)
    plt.xlim(1e9,1e15)
    plt.ylim(0.8, 2.2)
    plt.text(2e10,1.8,'XTE',fontsize=23)
    plt.savefig('XTE_temperature_profile.pdf',format='pdf')

    plt.show()

#plot_new()

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

    for i,col,lbl in zip([6,5,4],['green','red','blue'], ['MAXI model 2', 'MAXI model 1', 'XTE']):
        print(i)
        data = np.loadtxt('output/cooling_GIPSF_' + str(i) + '.dat')
        plt.plot((data[:, 1]-1.00145e3)*365, np.log10(np.abs(data[:, 3])+1e-9), lw=3, color=col, ls='-',label=lbl)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.yticks([30,31,32,33,34,35,36,37],fontsize=18)
    plt.ylim(34,37)
    plt.ylabel('$\\rm log \\thinspace L^{\infty}_{h}$ $\\rm erg \\thinspace s^{-1}$', fontsize=22)
    plt.xlim(-532,1500)
    plt.legend(loc='upper right', fontsize=18, frameon=False)

    axx = plt.subplot(2,1,2)
    x_minor_locator2 = AutoMinorLocator(5)
    y_minor_locator2 = AutoMinorLocator(5)
    plt.tick_params(which='both', width=1.7)
    plt.tick_params(which='major', length=9)
    plt.tick_params(which='minor', length=5)
    axx.xaxis.set_minor_locator(x_minor_locator2)
    axx.yaxis.set_minor_locator(y_minor_locator2)

    k_b = 8.617330350e-5
    for i,col in zip(range(4,7),['blue', 'red','green']):
        data = np.loadtxt('output/cooling_GIPSF_' + str(i) + '.dat')
        data2 = np.loadtxt('output/cooling_SF0_' + str(i) + '.dat')
        plt.plot((data[:, 1]-1.00145e3)*365, np.log10(data[:, 5]), lw=3, color=col, ls='-')
        plt.plot((data2[:, 1] -1.00145e3) * 365, np.log10(data2[:, 5]), lw=3, color=col, ls='--')

    #plt.scatter(t_KS, T_KS, s=100, color='red', marker='s', label='KS 1731-260')
    #plt.errorbar(x=t_KS, y=T_KS, yerr=err_KS, color='red', fmt=' ')

    #plt.scatter(t_MXB, T_MXB, s=100, color='darkgreen', marker='D', label='MXB 1659−29')
    #plt.errorbar(x=t_MXB, y=T_MXB, yerr=err_MXB, color='darkgreen', fmt=' ')

    plt.xticks([-500,0,500,1000,1500,2000],fontsize=18)
    plt.xlim(-532,1500)
    #plt.yticks([50,75,100,125,150],fontsize=18)
    #plt.ylim(50, 150)
    plt.xlabel('$\\rm Days \\thinspace  \\thinspace  since  \\thinspace  \\thinspace  '
               'end \\thinspace  \\thinspace  of  \\thinspace  \\thinspace outburst$',fontsize=22)
    #plt.ylabel('$\\rm kT^{\infty}_{s}$ $\\rm eV$',fontsize=22)
    plt.yticks([32, 32.5, 33, 33.5, 34, 34.5, 35, 35.5, 36], fontsize=18)
    plt.ylim(33.6, 35)
    plt.ylabel('$\\rm log \\thinspace L_{s}^{\infty} \\thinspace erg \\thinspace s^{-1}$', fontsize=22)
    plt.savefig('fig4_new.pdf',format='pdf')
    plt.show()

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

    for i,col,lb in zip(range(2,4),['blue','orange'],['KS 1731-260','MXB 1659−29']):
        data = np.loadtxt('output/cooling_GIPSF_' + str(i) + '.dat')
        data2 = np.loadtxt('output/cooling_SF0_' + str(i) + '.dat')
        plt.plot((data[:, 1]-4e3)*365, np.log10(data[:, 5]), lw=3, color=col, ls='-',label=lb)
        plt.plot((data2[:, 1] - 4e3) * 365, np.log10(data2[:, 5]), lw=3, color=col, ls='--')

    plt.xlim(-500,2000)
    plt.xticks([-500,0,500,1000,1500,2000],fontsize=18)
    plt.yticks([32,32.5,33,33.5,34,34.5,35,35.5,36],fontsize=18)
    plt.ylim(32.3,34.2)
    plt.xlabel('$\\rm Days \\thinspace  \\thinspace  since  \\thinspace  \\thinspace  '
               'end \\thinspace  \\thinspace  of  \\thinspace  \\thinspace outburst$',fontsize=22)
    plt.ylabel('$\\rm log \\thinspace L_{s}^{\infty} \\thinspace erg \\thinspace s^{-1}$',fontsize=22)
    plt.legend(loc='upper right',fontsize=20)
    plt.savefig('fig6_new.pdf',format='pdf')
    plt.show()

#plot4()





def f2(x,a,b,c,d):
    r_0 = radii(1e9)
    output = np.zeros_like(x)
    for i in range(len(output)):
        r = np.linspace(r_0, x[i], 200)
        output[i] = H_guess * a * np.abs(np.trapz((stats.norm.pdf(x=r, loc=b, scale=c)+d)*dV(r)*np.exp(2*_Phi(np.log10(r))), r))
    return output




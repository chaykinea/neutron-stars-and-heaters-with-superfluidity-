__author__ = 'maryhallow'


import numpy
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
from scipy import interpolate
import matplotlib.cm as cm


# report 2

def find_nearest(array,value):
    idx = (numpy.abs(array-value)).argmin()
    return idx


def cooling_curve():

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['figure.figsize'] = 8, 7.5

    line_styles = ['b-','r--','y--','k--']
    dashes = numpy.array([[14,4],[20,5],[5,9]])
    line_thickness = [4,3.5,3,2.5]

    cooling_1 = numpy.loadtxt('output/cooling_curve_nosf.dat')
    cooling_2 = numpy.loadtxt('output/cooling_curve_sf2.dat')

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

    plt.plot(numpy.log10(cooling_1[:,1]),numpy.log10(cooling_1[:,0]),line_styles[0], linewidth=line_thickness[0],label='No SF')
    plt.plot(numpy.log10(cooling_2[:,1]),numpy.log10(cooling_2[:,0]),line_styles[1], linewidth=line_thickness[1], dashes = (dashes[0,0],dashes[0,1],4,dashes[0,1]),label='SF')

    plt.legend(loc='upper right',fontsize=18)
    plt.savefig('cooling_curves_2.pdf',format='pdf')
    plt.show()


def interpolation():

    cooling_1 = numpy.loadtxt('output/cooling_curve_nosf.dat')
    #cooling_2 = numpy.loadtxt('output/cooling_curve_sf2.dat')

    nosf_cooling = interpolate.interp1d(numpy.log10(cooling_1[:, 1]), numpy.log10(cooling_1[:, 4]), kind='linear')
    #sf_cooling   = interpolate.interp1d(numpy.log10(cooling_2[:, 1]), numpy.log10(cooling_2[:, 4]), kind='linear')

    #nosf_cooling_2 = interpolate.interp1d(numpy.log10(cooling_1[:, 0]), numpy.log10(cooling_1[:, 1]), kind='linear')
    #sf_cooling_2   = interpolate.interp1d(numpy.log10(cooling_2[:, 0]), numpy.log10(cooling_2[:, 1]), kind='linear')

    #nosf_cooling_3 = interpolate.interp1d(numpy.log10(cooling_1[:, 1]), numpy.log10(cooling_1[:, 0]), kind='linear')
    sf_cooling_3   = interpolate.interp1d(numpy.log10(cooling_2[:, 1]), numpy.log10(cooling_2[:, 0]), kind='linear')
    '''
    time = numpy.linspace(4.51,4.55,10000)

    for i in range(len(time)):
        a = sf_cooling(time[i])-nosf_cooling(time[i])
        b = sf_cooling(time[i+1])-nosf_cooling(time[i+1])
        if(a>=0 and b<=0):
            time_root = (time[i]+time[i+1])/2
            print(time_root)
            print(numpy.power(10,time_root))
            print(nosf_cooling_3(time[i+1]))
            break

    time_1 = nosf_cooling_2(6.)
    time_2  = sf_cooling_2(6.)
    '''
    #print(time_1,time_2)
    #print(numpy.power(10,time_1),numpy.power(10,time_2))



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

        cooling_sf = numpy.loadtxt('sf/file_' + str(len(periods)*len(rho1_array)*a + len(periods)*b + i) +  '_cooling.dat')
        cooling_sf[:,1] -= 35370.228363

        cooling_no = numpy.loadtxt('output/file_' + str(len(periods)*len(rho1_array)*a + len(periods)*b + i) +  '_cooling.dat')
        cooling_no[:,1] -= 35370.228363

        print(len(periods)*len(rho1_array)*a + len(periods)*b + i)

        plt.plot((cooling_no[:,1]),numpy.log10(cooling_no[:,3]),'y-',lw=3,label='Heater')
        plt.plot((cooling_no[:,1]),numpy.log10(cooling_no[:,5]),'g-',lw=3,label='Surface (non-superfluid)')
        plt.plot((cooling_sf[:,1]),numpy.log10(cooling_sf[:,5]),'r--',lw=3.5,label='Surface (superfluid)')


        L_max_sf = 1000 + numpy.argmax(cooling_sf[1000:,5])
        print(L_max_sf)
        L_max_no = 1000 + numpy.argmax(cooling_no[1000:,5])
        print(L_max_no)
        H_max = 1000 + numpy.argmax(cooling_no[1000:,3])
        print(H_max)

        plt.plot(cooling_sf[L_max_sf,1], numpy.log10(cooling_sf[L_max_sf,5]),  'ko', markersize=12)
        plt.plot(cooling_no[L_max_no,1], numpy.log10(cooling_no[L_max_no,5]),  'ko', markersize=12)
        plt.plot(cooling_no[H_max,1], numpy.log10(cooling_no[H_max,3]),  'ko', markersize=12)

        A = 1000

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


def delay_plot():

    plt.rcParams['axes.linewidth'] = 2.5
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 8
    plt.rcParams.update({'figure.autolayout': True})

    rho1_array = numpy.power(10,numpy.linspace(10,14,9))
    periods      = numpy.power(10,numpy.linspace(-1,2,7))
    source_power = numpy.power(10,numpy.linspace(18,21,4))

    levels1 = numpy.array([0.1,0.3,1,2,3,5,10,17])

    x = numpy.linspace(-1,2,7)
    y = numpy.linspace(10,13.5,8)
    z_no = numpy.zeros((len(x),len(y)))
    z_sf = numpy.zeros((len(x),len(y)))

    '''
    for j in range(0,len(y)):
        for i in range(0,len(x)):
            temp1 = numpy.loadtxt('nosf/file_' + str(len(periods)*len(rho1_array)*2 + len(periods)*j + i) + '_cooling.dat')
            temp2 = numpy.loadtxt('sf/file_' + str(len(periods)*len(rho1_array)*2 + len(periods)*j + i) + '_cooling.dat')

            L_max_no = 1000 + numpy.argmax(temp1[1000:,5])
            L_max_sf = 1000 + numpy.argmax(temp2[1000:,5])
            H_max = 1000 + numpy.argmax(temp1[1000:,3])

            time_diff1 = temp1[L_max_no,1] - temp1[H_max,1]
            time_diff2 = temp2[L_max_sf,1] - temp2[H_max,1]

            z_no[i,j] = time_diff1
            z_sf[i,j] = time_diff2

            print(i,j)


    numpy.savetxt('contour_time_delay_1e20_nosf.dat',z_no)
    numpy.savetxt('contour_time_delay_1e20_sf.dat',z_sf)

    '''
    z_no = numpy.loadtxt('contour_time_delay_1e20_nosf.dat')
    z_sf = numpy.loadtxt('contour_time_delay_1e20_sf.dat')

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
    plt.savefig('contour_time_delay_1e20_nosf.pdf',format='pdf')

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
    plt.savefig('contour_time_delay_1e20_sf.pdf',format='pdf')

    plt.show()


def halfwidth_plot():

    plt.rcParams['axes.linewidth'] = 2.5
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 8
    plt.rcParams.update({'figure.autolayout': True})

    rho1_array = numpy.power(10,numpy.linspace(10,14,9))
    periods      = numpy.power(10,numpy.linspace(-1,2,7))
    source_power = numpy.power(10,numpy.linspace(18,21,4))

    levels1 = numpy.array([1.5,2.5,3.5,4.5])
    x = numpy.linspace(-1,1.5,6)
    y = numpy.linspace(10,13.5,8)
    z = numpy.zeros((len(x),len(y)))

    '''
    for j in range(0,len(y)):
        for i in range(0,len(x)):
            temp1 = numpy.loadtxt('nosf/file_' + str(len(periods)*len(rho1_array)*2 + len(periods)*j + i) + '_cooling.dat')
            temp2 = numpy.loadtxt('sf/file_' + str(len(periods)*len(rho1_array)*2 + len(periods)*j + i) + '_cooling.dat')

            L_max_no = 1000 + numpy.argmax(temp1[1000:,5])
            L_max_sf = 1000 + numpy.argmax(temp2[1000:,5])

            A = 1000

            value1 = (temp1[A,5]+temp1[L_max_no,5])/2
            value2 = (temp2[A,5]+temp2[L_max_sf,5])/2

            left1 = A + find_nearest(temp1[A:L_max_no,5] ,value1)
            right1 = L_max_no + find_nearest(temp1[L_max_no:,5] ,value1)

            left2 = A + find_nearest(temp2[A:L_max_sf,5] ,value2)
            right2 = L_max_sf + find_nearest(temp2[L_max_sf:,5] ,value2)

            z[i,j] = (temp1[right1,1] - temp1[left1,1])/(temp2[right2,1] - temp2[left2,1])

            print(i,j)

    numpy.savetxt('contour_halfwidth_ratios_1e20.dat',z)
    '''

    z = numpy.loadtxt('contour_halfwidth_ratios_1e20.dat')

    X, Y = numpy.meshgrid(x, y)

    plt.figure(1)
    im = plt.imshow(numpy.rot90(z)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,1.5,10,13.5),vmin=1,vmax=5)
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

    plt.savefig('contour_halfwidth_ratios_1e20_sf.pdf',format='pdf')

    plt.show()



def max_plot():

    plt.rcParams['axes.linewidth'] = 2.5
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 8
    plt.rcParams.update({'figure.autolayout': True})

    rho1_array = numpy.power(10,numpy.linspace(10,14,9))
    periods      = numpy.power(10,numpy.linspace(-1,2,7))
    source_power = numpy.power(10,numpy.linspace(18,21,4))

    levels1 = numpy.array([1.3,2,3])
    x = numpy.linspace(-1,2,7)
    y = numpy.linspace(10,13.5,8)
    z = numpy.zeros((len(x),len(y)))

    '''
    for j in range(0,len(y)):
        for i in range(0,len(x)):
            temp1 = numpy.loadtxt('nosf/file_' + str(len(periods)*len(rho1_array)*2 + len(periods)*j + i) + '_cooling.dat')
            temp2 = numpy.loadtxt('sf/file_' + str(len(periods)*len(rho1_array)*2 + len(periods)*j + i) + '_cooling.dat')

            L_max_no = 1000 + numpy.argmax(temp1[1000:,5])
            L_max_sf = 1000 + numpy.argmax(temp2[1000:,5])

            z[i,j] = (temp2[L_max_sf,5]/temp2[1000,5])/(temp1[L_max_no,5]/temp1[1000,5])

            print(i,j)
    '''

    # numpy.savetxt('contour_max_ratios_1e20.dat',z)

    z = numpy.loadtxt('contour_max_ratios_1e20.dat')

    X, Y = numpy.meshgrid(x, y)

    plt.figure(1)
    im = plt.imshow(numpy.rot90(z)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,2,10,13.5),vmin=1,vmax=3.5)
    CS1 = plt.contour(X, Y, numpy.transpose(z),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS1, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,13.5,8),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)

    plt.ylabel('$\mathrm{log} \\thinspace \\rho_{1} \\thinspace [\mathrm{g \\thinspace   cm^{-3}}]$',fontsize=32)
    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)

    CBI2 = plt.colorbar(im, orientation='vertical',ticks=[1,1.5,2,2.5,3,3.5])
    CBI2.ax.tick_params(labelsize=24)
    CBI2.set_label('$\\rm L^{\infty}_{max} (SF) / L^{\infty}_{max} (non-SF)$', labelpad=5, y=0.45,fontsize=22)

    plt.text(-0.6, 10.8 ,'$ \\rm everything\\thinspace below $',fontsize=30,color='blue')
    plt.text(-0.6, 10.5 ,'$ \\rm \\thinspace [log \\thinspace \\rho_{1} = 11.0] \\approx 1$',fontsize=30,color='blue')

    plt.savefig('contour_max_ratios_1e20.pdf',format='pdf')
    plt.show()



def max_min_plot():

    plt.rcParams['axes.linewidth'] = 2.5
    plt.rcParams['figure.figsize'] = 8, 7.5
    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 8
    plt.rcParams.update({'figure.autolayout': True})

    rho1_array = numpy.power(10,numpy.linspace(10,14,9))
    periods      = numpy.power(10,numpy.linspace(-1,2,7))
    source_power = numpy.power(10,numpy.linspace(18,21,4))

    levels1 = numpy.array([2,10,25,35,50,75,175])
    x = numpy.linspace(-1,2,7)
    y = numpy.linspace(10,13.5,8)
    z_no = numpy.zeros((len(x),len(y)))
    z_sf = numpy.zeros((len(x),len(y)))

    '''
    for j in range(0,len(y)):
        for i in range(0,len(x)):
            temp1 = numpy.loadtxt('nosf/file_' + str(len(periods)*len(rho1_array)*2 + len(periods)*j + i) + '_cooling.dat')
            temp2 = numpy.loadtxt('sf/file_' + str(len(periods)*len(rho1_array)*2 + len(periods)*j + i) + '_cooling.dat')

            L_max_no = 1000 + numpy.argmax(temp1[1000:,5])
            L_max_sf = 1000 + numpy.argmax(temp2[1000:,5])

            z_sf[i,j] = temp2[L_max_sf,5]/temp2[1000,5]
            z_no[i,j] = temp1[L_max_no,5]/temp1[1000,5]

            print(i,j)

    numpy.savetxt('contour_max_min_1e20_nosf.dat',z_no)
    numpy.savetxt('contour_max_min_1e20_sf.dat',z_sf)
    '''

    z_no = numpy.loadtxt('contour_max_min_1e20_nosf.dat')
    z_sf = numpy.loadtxt('contour_max_min_1e20_sf.dat')

    X, Y = numpy.meshgrid(x, y)

    plt.figure(1)
    im = plt.imshow(numpy.rot90(z_no)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,2,10,13.5),vmax=100,vmin=1)
    CS1 = plt.contour(X, Y, numpy.transpose(z_no),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS1, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,13.5,8),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)

    plt.ylabel('$\mathrm{log} \\thinspace \\rho_{1} \\thinspace [\mathrm{g \\thinspace   cm^{-3}}]$',fontsize=32)
    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)

    plt.savefig('contour_max_min_1e20_nosf.pdf',format='pdf')

    plt.figure(2)
    im = plt.imshow(numpy.rot90(z_sf)[::-1,:], interpolation='bicubic', origin='lower',cmap=cm.hot,extent=(-1,2,10,13.5),vmax=100,vmin=1)
    CS2 = plt.contour(X, Y, numpy.transpose(z_sf),colors='blue',levels=levels1, linewidths=2.0)
    plt.clabel(CS2, fontsize=22,fmt='%-3.1f')
    plt.xticks(numpy.linspace(-1,2,7),fontsize=24)
    plt.yticks(numpy.linspace(10,13.5,8),fontsize=24)
    plt.tick_params(width=2,color='blue',length=7)

    plt.xlabel('$\mathrm{log}   \\thinspace \\Delta t \\thinspace \\thinspace [\mathrm{yr}]$',fontsize=32)

    CBI2 = plt.colorbar(im, orientation='vertical')
    CBI2.ax.tick_params(labelsize=24)
    CBI2.set_label('$\\rm  L^{\infty}_{max} / L^{\infty}_{min} $', labelpad=5, y=0.45,fontsize=22)

    plt.savefig('contour_max_min_1e20_sf.pdf',format='pdf')


    plt.show()



x_max=71
x_min=-1
y_min=33
y_max=38
#visual(2,4,2)
#halfwidth_plot()
#max_plot()
#max_min_plot()
#delay_plot()
#interpolation()
#cooling_curve()
import numpy as np
cooling_1 = numpy.loadtxt('output/cooling_curve_nosf.dat')
function = interpolate.interp1d(np.log10(cooling_1[:,0]),cooling_1[:,1])
for i in [6.2,6.1,6.0,5.9,5.8]:
    print(function(i))

'''
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['figure.figsize'] = 8, 7.5
plt.rcParams.update({'figure.autolayout': True})

f = plt.figure()
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['figure.figsize'] = 8, 7.5
plt.rcParams.update({'figure.autolayout': True})
plt.tick_params(which='both', width=1.7)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=5)
axarr = f.add_subplot(1,1,1)
yy = numpy.array([0,0.5,1,1.5,2,2.5,3.0,3.5])
xx = numpy.array([5.8,5.9,6.0,6.1,6.2])
x_minor_locator = AutoMinorLocator(2)
y_minor_locator = AutoMinorLocator(5)
plt.tick_params(which='both', width=1.7)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=5)
axarr.xaxis.set_minor_locator(x_minor_locator)
axarr.yaxis.set_minor_locator(y_minor_locator)
plt.yticks(yy,fontsize=24)
plt.xticks(xx,fontsize=24)
plt.ylabel('$\\rm log \\thinspace L^{\infty}_{max}/L^{\infty}_{c}$',fontsize=32)
plt.xlabel('$\\rm log \\thinspace T^{\infty}_{s} \\thinspace K$',fontsize=32)

time_roots = np.array([183.0929720889007,1759.0346861045732,19974.966903050416,137916.19678829238,429396.4912770355])
numbers = np.array([2733, 3473, 3999, 3999, 3999])
temperature = np.array([6.2,6.1,6.0,5.9,5.8])
plt.ylim(0,2.5)
labels = np.array(['$\\rm H_{0} = 5 \\times 10^{17}\\thinspace  erg \\thinspace s^{-1} cm^{-3}$',
                   '$\\rm H_{0} = 5 \\times 10^{18}\\thinspace  erg \\thinspace s^{-1} cm^{-3}$',
                   '$\\rm H_{0} = 5 \\times 10^{19}\\thinspace  erg \\thinspace s^{-1} cm^{-3}$',
                   '$\\rm H_{0} = 5 \\times 10^{20}\\thinspace  erg \\thinspace s^{-1} cm^{-3}$'])
colors = ['blueviolet','lawngreen','lightcoral','darkgreen']
thick = [5,4.3,3.6,2.9]
for j in range(3,-1,-1):
    maxx = []
    for i in range(0,5):
        temp = numpy.loadtxt('output/_SF0_' + str(i) + '_' + str(j) + '_cooling.dat')
        temp[:,1] -= time_roots[i]
        maxx.append(np.max(temp[numbers[i]:,5])/temp[numbers[i],5])
    function = interpolate.interp1d(temperature[::-1],maxx[::-1],kind='cubic')
    sample = np.linspace(5.8,6.2,500)
    values = function(sample)
    plt.plot(sample,np.log10(values),color=colors[j],label=labels[j],lw=thick[j])


plt.legend(loc='upper right',fontsize=16)
plt.savefig('fig.pdf',format='pdf')
plt.show()
'''

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['figure.figsize'] = 8, 7.5
plt.rcParams.update({'figure.autolayout': True})

f = plt.figure()
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['figure.figsize'] = 8, 7.5
plt.rcParams.update({'figure.autolayout': True})
plt.tick_params(which='both', width=1.7)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=5)
axarr = f.add_subplot(1,1,1)
yy = numpy.array([0,0.5,1,1.5,2,2.5,3.0,3.5])
xx = numpy.array([0,2,4,6,8,10,12,14,16,18,20])
x_minor_locator = AutoMinorLocator(2)
y_minor_locator = AutoMinorLocator(5)
plt.tick_params(which='both', width=1.7)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=5)
axarr.xaxis.set_minor_locator(x_minor_locator)
axarr.yaxis.set_minor_locator(y_minor_locator)
plt.yticks(yy,fontsize=24)
plt.xticks(xx,fontsize=24)
plt.ylabel('$\\rm log \\thinspace L^{\infty}_{s}/L^{\infty}_{c}$',fontsize=32)
plt.xlabel('$\\rm t \\thinspace yr$',fontsize=32)
dashes = np.array([[12,4],[22,4],[7,4], [2,2], [12,4],[5,12], [2,2]])
line_thickness = np.array([4.5, 4.1, 3.3, 1.7, 3.5, 3.5, 3.1, 3,5, 2])
time_roots = np.array([183.0929720889007,1759.0346861045732,19974.966903050416,137916.19678829238,429396.4912770355])
numbers = np.array([2733, 3473, 3999, 3999, 3999])
temperature = np.array([6.2,6.1,6.0,5.9,5.8])
plt.ylim(0,2.5)
plt.xlim(-1,11)
labels = np.array(['$\\rm T_{s}^{\infty} = 10^{6.2}\\thinspace K$',
                   '$\\rm T_{s}^{\infty} = 10^{6.1}\\thinspace K$',
                   '$\\rm T_{s}^{\infty} = 10^{6.0}\\thinspace K$',
                   '$\\rm T_{s}^{\infty} = 10^{5.9}\\thinspace K$',
                   '$\\rm T_{s}^{\infty} = 10^{5.8}\\thinspace K$'])

colors = ['blueviolet','lightcoral','lawngreen', 'darkgreen','orange']
thick = [5.2,4.5,3.8,3.1,2.4]
for j in range(3,2,-1):
    for i in range(4,-1,-1):
        temp = numpy.loadtxt('output/_SF0_' + str(i) + '_' + str(j) + '_cooling.dat')
        temp[:,1] -= time_roots[i]
        idx = i
        if(idx>2 and idx<4):
            plt.plot(temp[:,1],np.log10(temp[:,5]/temp[numbers[i],5]),color=colors[idx],
                     linewidth=line_thickness[idx], dashes = (dashes[idx-2,0],dashes[idx-2,1]),label=labels[idx])
        elif(idx<=2):
            plt.plot(temp[:,1],np.log10(temp[:,5]/temp[numbers[i],5]),color=colors[idx],
                     linewidth=line_thickness[idx+2], dashes = (dashes[idx+2,0],dashes[idx+2,1],4,dashes[idx+2,1]),label=labels[idx])
        else:
            plt.plot(temp[:,1],np.log10(temp[:,5]/temp[numbers[i],5]), '-',color=colors[idx],
                 linewidth='3',label=labels[idx])


plt.text(5.2,1.2,'$\\rm \\Delta t = 10 \\thinspace days$',fontsize=20)
plt.text(5.2,1.0,'$\\rm \\rho_{1} = 10^{10} \\thinspace g cm^{-3}$',fontsize=20)
plt.text(5.2,0.8,'$\\rm H_{0} = 5\\times 10^{20} \\thinspace erg \\thinspace s^{-1} cm^{-3}$',fontsize=20)
plt.legend(loc='upper right',fontsize=18)
plt.savefig('fig8.pdf',format='pdf')
plt.show()


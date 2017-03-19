__author__ = 'maryhallow'

import numpy
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
from control.manager import *
from control.constants import *
from scipy import interpolate
import matplotlib.cm as cm
from data import loaddata


def find_nearest(array,value):
    idx = (numpy.abs(array-value)).argmin()
    return idx


def visual(b,c):

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

    density_label = numpy.array(['$10^{11}$','$10^{12}$','$10^{13}$'])
    period_label = numpy.array(['$10^{0}$','$10^{1}$'])


    cooling_sf = numpy.loadtxt('sf/file_' + str(b*len(period_label)+c) +  '_cooling.dat')
    cooling_sf[:,1] -= 35370.228363

    cooling_no = numpy.loadtxt('no_sf/file_' + str(b*len(period_label)+c) +  '_cooling.dat')
    cooling_no[:,1] -= 35370.228363

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



    plt.text((x_max-x_min)/2.2, y_max-2 , '$ \\rm H_{0}  =  $' + '$10^{20}$' + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=22)
    plt.text((x_max-x_min)/2.2, y_max-2.5 , '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=22)
    plt.text((x_max-x_min)/2.2, y_max-3 , '$ \\rm P  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=22)

    plt.legend(loc='upper right',fontsize=22)
    #plt.savefig('profile_' + str(a) + '_' + str(b) + '_' + str(c) + '.pdf', format='pdf')
    plt.show()


x_max=51
x_min=-1
y_min=33
y_max=38

#visual(2,0)


def flux(b,c,time_index,time_index2):

    dashes = numpy.array([[14,4],[20,5],[5,9]])
    line_thickness = [4,3,2,2.5]

    loaddata.superfluid_data_init()
    loaddata.star_model_data_init()

    print(t_points_save_data[time_index] - turn_on_time )

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

    yy = numpy.array([-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,0,500,1000,1500,2000,2500])
    xx = numpy.array([9,10,11,12,13,14,15])

    plt.yticks(yy,fontsize=24)
    plt.xticks(xx,fontsize=24)

    axarr.set_ylim(y_min,y_max)
    axarr.set_xlim(x_min,x_max)

    axarr.set_ylabel('$\\rm L^{\infty}_{r}/L_{\odot}$',fontsize=32)
    axarr.set_xlabel('$\\rm log \\thinspace \\rho_{1} [g \\thinspace cm^{-3}]$',fontsize=32)

    density_label = numpy.array(['$10^{11}$','$10^{12}$','$10^{13}$'])
    period_label = numpy.array(['$10^{0}$','$10^{1}$'])
    time_label  = numpy.array(['$0.01$','$0.1$','$0.5$','$1.0$','$2.0$','$3.0$',
                               '$4.0$','$5.0$','$6.0$','$7.0$','$8.0$','$9.0$','$10.0$'])

    flux_sf = numpy.loadtxt('sf/file_' + str(b*len(period_label)+c) +  '_flux.dat')
    flux_no = numpy.loadtxt('no_sf/file_' + str(b*len(period_label)+c) +  '_flux.dat')

    plt.axvline(x=13,     linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')
    plt.axvline(x=numpy.log10(5.063467e13),  linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')

    plt.plot((numpy.log10(flux_no[:,0])),flux_no[:,time_index+1]/LSun,'g-',lw=2.5,label='non-sf, $ \\rm t-t_{0}  =  $' + time_label[time_index] + '$\\rm \\thinspace \\thinspace yr $')
    plt.plot((numpy.log10(flux_sf[:,0])),flux_sf[:,time_index+1]/LSun,'r--',lw=3.5,label='sf, $ \\rm t-t_{0}  =  $' + time_label[time_index] + '$\\rm \\thinspace \\thinspace yr $')

    #plt.plot((numpy.log10(flux_no[:,0])),flux_no[:,time_index2+1]/LSun,'y-',lw=2.5,label='non-sf, $ \\rm t-t_{0}  =  $' + time_label[time_index2] + '$\\rm \\thinspace \\thinspace yr $')
    #plt.plot((numpy.log10(flux_sf[:,0])),flux_sf[:,time_index2+1]/LSun,'b--',lw=3.5,label='sf, $ \\rm t-t_{0}  =  $' + time_label[time_index2] + '$\\rm \\thinspace  \\thinspace yr $')

    #plt.text(12.3, -700 , '$ \\rm H_{0}  =  $' + '$10^{20}$' + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=18)
    #plt.text(12.3, -900 , '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=18)
    #plt.text(12.3, -1100 , '$ \\rm P  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=18)

    plt.text(9.3, -600 , '$ \\rm H_{0}  =  $' + '$10^{20}$' + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=18)
    plt.text(9.3, -1100 , '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=18)
    plt.text(9.3, -1600 , '$ \\rm P  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=18)

    print(density_label[b])
    print(period_label[c])

    plt.legend(loc='upper right',fontsize=16)
    plt.savefig('flux_' + str(b) + '_' + str(c) + '.pdf', format='pdf')
    plt.show()


def temperature(b,c,time_index,time_index2):

    dashes = numpy.array([[14,4],[20,5],[5,9]])
    line_thickness = [4,3,2,2.5]

    loaddata.superfluid_data_init()
    loaddata.star_model_data_init()

    print(t_points_save_data[time_index] - turn_on_time )

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

    yy = numpy.array([1,2,3,4,5,6,7,8,9,10,11])
    xx = numpy.array([9,10,11,12,13,14,15])

    plt.yticks(yy,fontsize=24)
    plt.xticks(xx,fontsize=24)

    axarr.set_ylim(y_min,y_max)
    axarr.set_xlim(x_min,x_max)

    axarr.set_ylabel('$\\rm T^{\infty}_{r}/10^{8} $',fontsize=32)
    axarr.set_xlabel('$\\rm log \\thinspace \\rho_{1} [g \\thinspace cm^{-3}]$',fontsize=32)

    density_label = numpy.array(['$10^{11}$','$10^{12}$','$10^{13}$'])
    period_label = numpy.array(['$10^{0}$','$10^{1}$'])
    time_label  = numpy.array(['$0.01$','$0.1$','$0.5$','$1.0$','$2.0$','$3.0$',
                               '$4.0$','$5.0$','$6.0$','$7.0$','$8.0$','$9.0$','$10.0$'])

    T_sf = numpy.loadtxt('sf/file_' + str(b*len(period_label)+c) +  '.dat')
    T_no = numpy.loadtxt('no_sf/file_' + str(b*len(period_label)+c) +  '.dat')

    plt.axvline(x=11,     linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')
    plt.axvline(x=12,  linewidth=line_thickness[3], dashes = (dashes[2,0],dashes[2,1]), color='k')

    plt.plot((numpy.log10(T_no[:,0])),T_no[:,time_index+1]/1e8,'g-',lw=2.5,label='non-sf, $ \\rm t-t_{0}  =  $' + time_label[time_index] + '$\\rm \\thinspace \\thinspace yr $')
    plt.plot((numpy.log10(T_sf[:,0])),T_sf[:,time_index+1]/1e8,'r--',lw=3.5,label='sf, $ \\rm t-t_{0}  =  $' + time_label[time_index] + '$\\rm \\thinspace \\thinspace yr $')

    #plt.plot((numpy.log10(flux_no[:,0])),flux_no[:,time_index2+1]/LSun,'y-',lw=2.5,label='non-sf, $ \\rm t-t_{0}  =  $' + time_label[time_index2] + '$\\rm \\thinspace \\thinspace yr $')
    #plt.plot((numpy.log10(flux_sf[:,0])),flux_sf[:,time_index2+1]/LSun,'b--',lw=3.5,label='sf, $ \\rm t-t_{0}  =  $' + time_label[time_index2] + '$\\rm \\thinspace  \\thinspace yr $')

    #plt.text(12.3, -700 , '$ \\rm H_{0}  =  $' + '$10^{20}$' + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=18)
    #plt.text(12.3, -900 , '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=18)
    #plt.text(12.3, -1100 , '$ \\rm P  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=18)

    plt.text(12.3, 7.8 , '$ \\rm H_{0}  =  $' + '$10^{20}$' + '$\\rm \\thinspace  erg \\thinspace  s^{-1} \\thinspace  cm^{-3}$',fontsize=19)
    plt.text(12.3, 7, '$ \\rm \\rho_{1}  =  $' + density_label[b] + '$\\rm \\thinspace  g \\thinspace  cm^{-3}$',fontsize=19)
    plt.text(12.3, 6.2 , '$ \\rm P  =  $' + period_label[c] + '$\\rm \\thinspace  yr $',fontsize=19)

    print(density_label[b])
    print(period_label[c])

    plt.legend(loc='upper right',fontsize=18)
    plt.savefig('temperature_' + str(b) + '_' + str(c) + '.pdf', format='pdf')
    plt.show()

y_min = 1
y_max = 11
x_min = 9
x_max = 15

temperature(0,1,7,12)
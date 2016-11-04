#     Neutron star thermal evolution
# --------------------------------------
#               sf_gap.py
# --------------------------------------
# This module provides superfluidity
# energy gaps for singlet  S  and triplet
# P (m=0) states, where 'm' is the
# projection of total moment of the system
# on the quantization axis.
# tau = T/Tc, where T - temperature,
# Tc - critical temperature of superfluidity.

import numpy

def _TripletGap(tau):

    if(tau >= 1.0):
        gap = 0.0
    elif(tau < 1.0e-4):
        gap = 1000.0
    else:
        gap = numpy.sqrt(1.0 - tau) * (0.7893 + 1.188 / tau)

    return gap

def _SingletGap(tau):

    if(tau >= 1.0):
        gap = 0.0
    elif(tau < 1.0e-4):
        gap = 1000.0
    else:
        gap = numpy.sqrt(1.0 - tau) * (1.456 - 0.157 / numpy.sqrt(tau) + 1.764 / tau)

    return gap



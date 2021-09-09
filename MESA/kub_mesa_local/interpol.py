#!/usr/bin/env python3 
import numpy as np
import sys
import consts

from interpol_beta import *

INTERPOL(float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]))

print('mdot HBA x1e12='+str(testCF((4.2/365)**(2/3),float(sys.argv[2]),float(sys.argv[4]),BETA0(float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[3])))/1e12))
print('lambda='+str(BETA0(float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[3]))))
# to compare:

#print(INTERPOL(float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))/testCF(consts.Rsun/(2.*(float(sys.argv[3])/5780.)**2.) / consts.au2cm,float(sys.argv[2]),float(sys.argv[4]),BETA0(float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[3]))))

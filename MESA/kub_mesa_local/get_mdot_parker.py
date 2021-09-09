#!/usr/bin/env python3

from scipy.special import lambertw
import numpy as np
import consts

def get_mdot_parker(radius,pressure,csound,f_rrbondi): 
	
	mdot = 4.*np.pi*radius**2.*pressure*csound**-1.*np.sqrt(-lambertw(-f_rrbondi)).real

	with open('temp_mdot_parker.txt', 'w') as f:
		f.writelines(str(mdot))

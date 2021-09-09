#!/usr/bin/env python3
# Tu+ 2015

import numpy as np
import consts

def lum_evo(period):
	
	omega0 = 2.*np.pi/(period*24.*60.*60.)
	omega_sun = 2.*np.pi/(27.*24.*60.*60.)
	tsat = 2.9*(omega0/omega_sun)**1.14
	b=(0.35*np.log10(omega0/omega_sun) - 0.98)**-1.
	a= 10.**(27.2)*(4570.)**-b
	
	a_euv=10.**(4.8)*10.**(0.86*np.log10(a))
	b_euv = b*0.86
	
	print(10**(-3.13)*consts.Lsun/(4.*np.pi*consts.au2cm**2.))
	print(10.**(4.8)*10.**(0.86*np.log10(10**(-3.13)*consts.Lsun))/(4.*np.pi*consts.au2cm**2.))
	
	# convert erg s^-1 to erg s^-1 cm^-2
	return [a*tsat**b/(4.*np.pi*consts.au2cm**2.),b,a_euv*tsat**b_euv/(4.*np.pi*consts.au2cm**2.),b_euv,tsat]

#!/usr/bin/env python3

pathset='/Applications/kub_mesa_local/LOGS/' #cr_local/hwchen_larogers_2016/LOGS/'

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import glob
import consts
from numpy import diff
from scipy.special import lambertw

def ad_den(r_,grad_ad_,R_p,cs_,Mcore_,rho_, gamma_): 
	return rho_ * (1. + grad_ad_ * (consts.G * Mcore_ / (cs_**2. * R_p)) * (R_p / r_ - 1.))**(1. / (gamma_ - 1.))

colnms = ['temperature', 'density', 'pressure', 'tau', 'csound', 'radius', 'sch_stable', 'ledoux_stable', 'energy', 'logE', 'luminosity', 'cv', 'cp', 'dm', 'mass','opacity','pressure_scale_height','entropy','empty', 'empty']

dist=0.1 #(/365.)**(2./3.)
totalmass=15.625
coremass=12.5
envfraction=0.2

all_files=glob.glob(os.path.join(pathset, "**"))
logfiles=[]
for fl in all_files:
	if 'profile' in np.array(fl.split('/'))[-1] and '.data' in np.array(fl.split('/'))[-1] and 'massloss' in np.array(fl.split('/'))[-1]:
		logfiles.append(fl)
	#if '_' in np.array(fl.split('/'))[-1]:
	#	if '.data' in np.array(fl.split('/'))[-1] and 'TH' not in np.array(fl.split('/'))[-1]:
	#		#if str(dist)+'_'+str(totalmass)+'_'+str(envfraction)+'_8.25' in np.array(fl.split('/'))[-1]:
	#		logfiles.append(fl)

print(logfiles)
def find_nearest(array, value):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return array[idx]

fd = np.ones((len(logfiles),19))

flc=0
for fl in logfiles:
	print(fl)
	lines_ = []
	with open(fl, 'r') as temp_f:
		lines = temp_f.readlines()
		cntr=0
		for l in lines:
			lines_.append(l)
			if cntr==2:
				strs=' '.join(l.split())
				age=strs.split(' ')[4]
				mdot=float(strs.split(' ')[21])*consts.Msun/(60*60*24*365) # cgs
			cntr+=1
	print('age [yr] = '+str(age))
	print('mdot [1e13 cgs]='+str(mdot/1e13))
	lines_=np.asarray(lines_)[5:]
	#print(lines_[0])

	data=np.ones((len(lines_),20))
	#print(lines_[1:2])
	linecount=0
	for linee in lines_[1:]:
		vcount=0
		#print(linee)
		for el in linee.split(' '):
			if len(el.split('E')) > 1.:
				if len(el.split('E')[1])>1:
					if el!='' and vcount<=23. and el[-1]!='+':
						if '\n' in el:
							el = el[:-2]
						#print(el)
						data[linecount,vcount]=el
						vcount+=1
		linecount+=1
	
	dat=pd.DataFrame(data)

	dat.columns=colnms
	dat=dat.drop(['empty'],axis=1)
	#print(np.asarray(dat['radius'])[0:-1]*consts.Rsun/consts.Rearth)

	massfrac = (np.asarray(dat['mass'])[0:-1]*(consts.Msun/consts.Mearth) - coremass)/(totalmass*envfraction)
	csound =  np.asarray(dat['csound'])[0:-1]
	temperature = np.asarray(dat['temperature'])[0:-1]
	gamma_ad = np.asarray(dat['cp'])[0:-1]/np.asarray(dat['cv'])[0:-1]
	sch_stable = np.asarray(dat['sch_stable'])[0:-1]
	pressure = np.asarray(dat['pressure'])[0:-1]
	rad_grad = np.diff(np.log(temperature))/np.diff(np.log(pressure)) # length n-1
	density = np.asarray(dat['density'])[0:-1]
	radius=np.asarray(dat['radius'])[0:-1]*consts.Rsun/consts.Rearth
	pressure_scale_height = np.asarray(dat['pressure_scale_height'])[0:-1]*consts.Rsun/consts.Rearth
	opacity = np.asarray(dat['opacity'])[0:-1]
	luminosity = np.asarray(dat['luminosity'])[0:-1]*consts.Lsun

	#plt.plot(np.asarray(dat['radius'])[0:-1]*consts.Rsun/consts.Rearth, (np.asarray(dat['mass'])[0:-1]*(consts.Msun/consts.Mearth) - coremass)/(totalmass*envfraction),label=age)
	#plt.plot(np.asarray(dat['radius'])[0:-1]*consts.Rsun/consts.Rearth, np.asarray(dat['cp'])[0:-1]/np.asarray(dat['cv'])[0:-1],label=age)
	#plt.plot(radius, (gamma_ad-1)/gamma_ad,label=age)
	plt.plot(radius, temperature,label=age)
	#plt.plot(radius, gamma_ad,label=age)
	#plt.plot(radius, massfrac,label=age)
	#plt.plot(radius, density,label=age)
	#plt.plot(radius[0:-1],rad_grad)
	#plt.plot(radius, sch_stable,label=age)
	#plt.plot(radius, csound,label=age)
	#plt.plot(radius,luminosity,label=age)
	
	r_rcb = radius[np.where(sch_stable>=1)]
	r_rcb_ = r_rcb[np.where(r_rcb>3)][-1]
	cs_rcb = csound[np.where(r_rcb==r_rcb_)]
	gamma_rcb = gamma_ad[np.where(r_rcb==r_rcb_)]
	grad_ad_rcb = (gamma_rcb-1)/gamma_rcb
	rho_rcb = density[np.where(r_rcb==r_rcb_)]
	print('rcb='+str(r_rcb_/consts.Rearth))
	print('cs, km/s='+str(csound[0]/1e5))
	Rbondi=consts.G*totalmass*consts.Mearth / (2. * csound[0]**2.)
	r_rb=radius[0]/(Rbondi/consts.Rearth)
	print('R/Rbondi='+str(r_rb))
	print('H='+str(pressure_scale_height[0]))
	Mach_=(csound[0] * opacity[0] * np.abs(mdot))/(4.*np.pi*consts.G*totalmass*consts.Mearth)
	print('Mach #='+str((csound[0] * opacity[0] * np.abs(mdot))/(4.*np.pi*consts.G*totalmass*consts.Mearth)))
	print('Mach # at photosphere, analytic='+str((radius[0]*consts.Rearth / Rbondi)**-2. * np.exp(-2. * (radius[0]*consts.Rearth / Rbondi)**-1.)))
	print('Mach # at photosphere, analytic='+str(np.sqrt(-lambertw(-1. * r_rb**-4 * np.exp(3. - 4. / r_rb)))))
	print('mdot, predicted = '+str(4.*np.pi*consts.G*totalmass*consts.Mearth*(np.sqrt(-lambertw(-1. * r_rb**-4 * np.exp(3. - 4. / r_rb))))/(opacity[0]*csound[0])))
	print('T(surface)='+str(temperature[0]))
	print('T_eq='+str((5780. * np.sqrt(0.5 * radius[0]*consts.Rearth * (dist*consts.au2cm)**-1.))))
	print('Lambda='+str(consts.G*totalmass*consts.Mearth*consts.mH/(consts.k*temperature[0]*radius[0]*consts.Rearth)))
	print('R roche='+str(dist*consts.au2cm * (totalmass*consts.Mearth/(totalmass*consts.Mearth+consts.Msun))**(1./3.)/consts.Rearth))
	
#plt.xlim(2,8)
plt.legend()
#plt.xscale('log')
plt.yscale('log')
plt.show()

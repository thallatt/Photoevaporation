#!/usr/bin/env python3

import numpy as np
from scipy import interpolate
import scipy.interpolate as sinter
import consts
import matplotlib.pyplot as plt
import pandas as pd
import math
from get_powerlaw import make_pdf
from mesa_data import *
from lum_evolution import lum_evo
import get_Protdistrb
import numpy.random as nr
from lognormal_cdf import lognormal_cdf
from powerlaw_cdf import powerlaw_cdf
from rayleigh_cdf import rayleigh_cdf

#### ANALYTIC ####

eta=0.15

def Lum_alt(Lx): return 10.**(4.8 + 0.86 * np.log10(Lx))

def Lum_tu(t_,p_):
	factor = 4.*np.pi*(consts.au2cm**2.)
	xray,xray_index,euv,euv_index,t_sat = lum_evo(p_)[0],lum_evo(p_)[1],lum_evo(p_)[2],lum_evo(p_)[3],lum_evo(p_)[4]
	# convert from flux to luminosity
	if t_/1.e6 < t_sat:
		L_ = Lum_alt(10**(30.46)) + 10**(30.46)
	else:
		L_  = xray*factor*((t_/1.e6)/t_sat)**xray_index + euv*factor*((t_/1.e6)/t_sat)**euv_index
	return L_

#print("Tu Lsateuv flux at au = "+str(Lum_alt(10**(30.46))/(4.*np.pi*(consts.au2cm**2.))))
#print("Tu Lsatxray flux at au = "+str(10**(30.46)/(4.*np.pi*(consts.au2cm**2.))))
#print("Tu Lsat total at au = "+str((10**(30.46)+Lum_alt(10**(30.46)))/(4.*np.pi*(consts.au2cm**2.))))
#print("Tu 90th euv flux at 226myr at au = "+str(1.2e36*((226)**-2.15)/(4.*np.pi*(consts.au2cm**2.))))
#print("Tu 90th xray flux at 226myr at au = "+str(2.3e36*((226)**-2.50)/(4.*np.pi*(consts.au2cm**2.))))
#print("Tu 90th euv flux at 600myr at au = "+str(1.2e36*((600)**-2.15)/(4.*np.pi*(consts.au2cm**2.))))
#print(str(3685.1546645010812*(600/226)**-2.15))
#print("Tu 90th xray flux at 600myr at au = "+str(2.3e36*((600)**-2.50)/(4.*np.pi*(consts.au2cm**2.))))
#print(str(3685.1546645010812*(600/226)**-2.15))
#print("Tu 50th euv flux at 23myr at au = "+str(4.8e32*((23)**-1.22)/(4.*np.pi*(consts.au2cm**2.))))
#print("Tu 50th xray flux at 23myr at au = "+str(2.6e32*((23)**-1.42)/(4.*np.pi*(consts.au2cm**2.))))
#print("Tu 50th euv  flux at 600myr at au = "+str(4.8e32*((600)**-1.22)/(4.*np.pi*(consts.au2cm**2.))))
#print("Tu 50th xray flux at 600myr at au = "+str(2.6e32*((23)**-1.42)/(4.*np.pi*(consts.au2cm**2.))))
#print(str(3702.890719955262*(600/23)**-1.22))
#print("Ribas flux at 23myr at au = "+str(Lum_ribas(23e6)/(4.*np.pi*(consts.au2cm**2.))))
#print("Ribas flux at 600myr at au = "+str(Lum_ribas(600e6)/(4.*np.pi*(consts.au2cm**2.))))

def mcore(GCR_,period_, radius_):

	a_=(period_/365.)**(2./3.)*consts.au2cm
	gyr_ = 1.e9*consts.s2yr
	#return interpolate.interp1d(dps,mcrit_1day)(period_)*consts.Mearth
	return np.sqrt(GCR_**-1. * (1.+GCR_)**-1. * gyr_ * eta * Lum_tu(1.3e8,1.)*(radius_*consts.Rearth)**3. * (4. * a_**2. * consts.G)**-1.)

def mcore__(corefn,period_star):

	#return interpolate.interp1d([1.,3.,5.,10.,20.],[interpolate.interp1d(dps,mcrit_1day)(period_),interpolate.interp1d(dps,mcrit_3day)(period_),interpolate.interp1d(dps,mcrit_5day)(period_),interpolate.interp1d(dps,mcrit_10day)(period_),interpolate.interp1d(dps,mcrit_20day)(period_)],fill_value='extrapolate')(period_star)
	return corefn(period_star)

scale=1.

# analytic power law 

#residuals_ = {}
#pwrlaw_forplot = {}
#for counter, numleft_ in enumerate(np.asarray(pss['y'])/scale):
#	per = dps[counter]
#	numevap_ = 1. - numleft_
#	for alpha__ in np.linspace(0,2,10):
#		if mcore(gcrs__(per),per,rad_(per)) / m0 > 1.:
#			resid_ = np.abs(numevap_ - ((mcore(gcrs__(per),per,rad_(per)) / m0)**(1.-alpha__) - 1.) / ((50./4.)**(1.-alpha__) - 1.))
#			#print('alpha='+str(alpha__)+', period='+str(per)+', numevap= '+str(numevap_)+', num less than mcrit='+str(((mcore(gcrs__(per),per,rad_(per)) / m0)**(1.-alpha__) - 1.) / ((40./4.)**(1.-alpha__) - 1.))+', mcrit='+str(mcore(gcrs__(per),per,rad_(per))/consts.Mearth))
#			if alpha__ in list(residuals_):
#				ls_ = residuals_[alpha__][0]
#				residuals_[alpha__] = np.array([np.sum([ls_,resid_])])
#			else:
#				residuals_[alpha__]=np.array([resid_])
#			if alpha__ in list(pwrlaw_forplot):
#				cu = pwrlaw_forplot[alpha__]
#				pwrlaw_forplot[alpha__] = np.append(cu, ((mcore(gcrs__(per),per,rad_(per)) / m0)**(1.-alpha__) - 1.) / ((50./4.)**(1.-alpha__) - 1.))
#			else:
#				pwrlaw_forplot[alpha__] = np.array([((mcore(gcrs__(per),per,rad_(per)) / m0)**(1.-alpha__) - 1.) / ((50./4.)**(1.-alpha__) - 1.)])
#		else:
#			resid_ = np.abs(numevap_ - ((1.0)**(1.-alpha__) - 1.) / ((50./4.)**(1.-alpha__) - 1.))
#			if alpha__ in list(residuals_):
#				ls_ = residuals_[alpha__][0]
#				residuals_[alpha__] = np.array([np.sum([ls_,resid_])])
#			else:
#				residuals_[alpha__]=np.array([resid_])
#			if alpha__ in list(pwrlaw_forplot):
#				cu = pwrlaw_forplot[alpha__]
#				pwrlaw_forplot[alpha__] = np.append(cu, ((1.0)**(1.-alpha__) - 1.) / ((50./4.)**(1.-alpha__) - 1.))
#			else:
#				pwrlaw_forplot[alpha__] = np.array([((1.0)**(1.-alpha__) - 1.) / ((50./4.)**(1.-alpha__) - 1.)])
#
#for key,dat in pwrlaw_forplot.items():
#	plt.plot(dps,(1.-dat) * (np.array(pss['y'])[-1]/scale / (1-dat[-1])),label=str(float(round(key,2))))
#	plt.errorbar(np.array(pss['x']),np.array(pss['y'])/scale,yerr=[(np.array(pss['y'])-np.array(pss['dybelow']))/scale,(np.array(pss['dyabove'])-np.array(pss['y']))/scale],marker='o',color='black')
#plt.xscale('log')
#plt.yscale('log')
#plt.legend()
#plt.show()

#### SAMPLE STELLAR ROTATION PERIODS ####

Protstar='N2362'
hiM, loM = 1.4, 0.5
if Protstar == 'N2362':
    N2362 = np.genfromtxt("/Users/Tim/Documents/Work/MSc/code/Prot_data/Irwin_ngc2362.txt", skip_header=53, delimiter='|', filling_values=np.nan)
    N2362M = N2362.transpose()[5]
    N2362Per = N2362.transpose()[3][np.where((N2362M > loM)*(N2362M <= hiM) > 0)[0]]
    invFProt = sinter.interp1d(np.linspace(0, 1, len(N2362Per)), np.sort(N2362Per))
elif Protstar == 'ONC':
    #ONC = np.genfromtxt("../Prot_data/herbstONCrPer.txt", skip_header=38, delimiter='|', filling_values=np.nan)
    #oncM = ONC.transpose()[3]
    #oncPer = ONC.transpose()[0][np.where((oncM > loM)*(oncM <= hiM) > 0)[0]]
    #oncPer = oncPer[np.where(~np.isnan(oncPer))[0]]
    #invFProt = sinter.interp1d(np.linspace(0, 1, len(oncPer)), np.sort(oncPer))
    invFProt = get_Protdistrb.getSProtFunc('ONC-comb', loM, hiM)
elif Protstar == 'N2547':
    N2547 = np.genfromtxt("/Users/Tim/Documents/Work/MSc/code/Prot_data/Irwin_ngc2547.txt", skip_header=57, delimiter='|', filling_values=np.nan)
    N2547M = N2547.transpose()[7]
    N2547Per = N2547.transpose()[5][np.where((N2547M > loM) * (N2547M <= hiM) > 0)[0]]
    invFProt = sinter.interp1d(np.linspace(0, 1, len(N2547Per)), np.sort(N2547Per))

#### EVAPORATION ####

def evaporate(dist,mean,*argv):
	planets={}
	
	if dist!='powerlaw':
		if dist=='rayleigh':
			rayl=rayleigh_cdf(mean,4.,50.)
			for counter,val in enumerate(dps):
				core_fn = interpolate.interp1d([1.,5.,10.],[interpolate.interp1d(dps,mcrit_1day_5_)(val),interpolate.interp1d(dps,mcrit_5day_5_)(val),interpolate.interp1d(dps,mcrit_10day_5_)(val)],fill_value='extrapolate')
				record_=[]
				for i in range(0,50000):
					Prot=5.83#invFProt(nr.uniform(0, 1, 1))
					mesa_crit=mcore__(core_fn,Prot)
					record_.append(mesa_crit)
				mesa_crit=np.median(mesa_crit)
				num_=rayl(mesa_crit)[0]
				# deal with bounds outside of cdf
				if num_>1.0: 
					num_=1.0
				if num_<0.:
					num_=0.
				planets[val]=1.-num_
			
			#for counter,val in enumerate(dps):
			#	sample_mesa=[]
			#	for i in range(1,1000):
			#		smp=np.random.rayleigh(mean,1)
			#		Prot=invFProt(nr.uniform(0, 1, 1))
			#		mesa_crit=mcore__(val,Prot)
			#		if smp > mesa_crit:
			#			sample_mesa.append(smp)
			#	planets[val]=sample_mesa#cm_f
		
		elif dist=='lognormal':
			lognm=lognormal_cdf(np.exp(mean),argv[0],4.,50.)
			for counter,val in enumerate(dps):
				core_fn = interpolate.interp1d([1.,3.,5.,10.,20.],[interpolate.interp1d(dps,mcrit_1day_5)(val),interpolate.interp1d(dps,mcrit_3day_5)(val),interpolate.interp1d(dps,mcrit_5day_5)(val),interpolate.interp1d(dps,mcrit_10day_5)(val),interpolate.interp1d(dps,mcrit_20day_5)(val)],fill_value='extrapolate')
				record_=[]
				for i in range(0,50000):
					Prot=5.83#invFProt(nr.uniform(0, 1, 1))
					mesa_crit=mcore__(core_fn,Prot)
					record_.append(mesa_crit)
				mesa_crit=np.median(mesa_crit)
				num_=lognm(mesa_crit)[0]
				# deal with bounds outside of cdf
				if num_>1.0: 
					num_=1.0
				if num_<0.:
					num_=0.
			#for counter,val in enumerate(dps):
			#	sample_mesa=[]
			#	for i in range(1,1000):
			#		smp=np.random.lognormal(mean,argv[0],1)
			#		Prot=invFProt(nr.uniform(0, 1, 1))
			#		mesa_crit=mcore__(val,Prot)
			#		if smp > mesa_crit:
			#			sample_mesa.append(smp)
			#	planets[val]=sample_mesa#cm_f
				planets[val]=1.-num_

		#sample_mesa=np.array(sample_mesa)
		#planets[val]=sample_mesa#cm_f
		
	else:
		#pwrlaw=powerlaw_cdf(argv[0],4.,50.)
		#for counter,val in enumerate(dps):
		#	core_fn = interpolate.interp1d([1.,3.,5.,10.,20.],[interpolate.interp1d(dps,mcrit_1day_5)(val),interpolate.interp1d(dps,mcrit_3day_5)(val),interpolate.interp1d(dps,mcrit_5day_5)(val),interpolate.interp1d(dps,mcrit_10day_5)(val),interpolate.interp1d(dps,mcrit_20day_5)(val)],fill_value='extrapolate')
		#	record_=[]
		#	for i in range(1,10000):
		#		Prot=5.83#invFProt(nr.uniform(0, 1, 1))[0]
		#		mesa_crit=mcore__(core_fn,Prot)
		#		record_.append(mesa_crit)
		#	mesa_crit=np.median(record_)
		#	print(mesa_crit)
		#	num_=pwrlaw(mesa_crit)[0]
		#	print(num_)
		#	planets[val]=1.-num_

		pwrlaw=make_pdf(argv[0], 4., 50.)
		pwrlaw_cdf=powerlaw_cdf(argv[0],4.,50.)
		omegas=[]
		print('populating power law...')
		
		for counter,val in enumerate(dps):
			core_fn = interpolate.interp1d([1.,3.,5.,10.,20.],[interpolate.interp1d(dps,mcrit_1day_5)(val),interpolate.interp1d(dps,mcrit_3day_5)(val),interpolate.interp1d(dps,mcrit_5day_5)(val),interpolate.interp1d(dps,mcrit_10day_5)(val),interpolate.interp1d(dps,mcrit_20day_5)(val)],fill_value='extrapolate')
			sample_mesa=[]
			omegas=[]
			ms=[]
			for i in range(1,10000):
				smp = pwrlaw(np.random.uniform(0, 1, 1))[0]				
				Prot=invFProt(nr.uniform(0, 1, 1))[0]
				mesa_crit=mcore__(core_fn,Prot)
				ms.append(mesa_crit)
			num_=pwrlaw_cdf(np.median(ms))[0]
			
			# deal with bounds outside of cdf
			if num_>1.0: 
				num_=1.0
			if num_<0.:
				num_=0.
			planets[val]=1.-num_

	plot_dist =False
	if plot_dist:
		if dist=='rayleigh':
			sample=np.random.rayleigh(mean,1000)
		elif dist=='lognormal':
			sample=np.random.lognormal(mean,argv[0],1000)
		else:
			sample=[]
			for i in range(1,1000):
				sample.append(pwrlaw(np.random.uniform(0, 1, 1))[0])
		x=plt.hist(sample,bins=np.logspace(np.log10(np.min(sample)),np.log10(np.max(sample)),8))
		#plt.plot(np.linspace(6,np.max(sample),15), np.max(x[0]) * (np.linspace(6,np.max(sample),15)/6.)**(-argv[0]))
		plt.xscale('log')
		plt.yscale('log')
		plt.show()

	num_pl_ar=[]
	#for id,dat in planets.items():
	#	num_pl_ar.append(len(dat))
	
	for id,dat in planets.items():
		num_pl_ar.append(np.abs(dat))

	dat = pd.DataFrame(np.array([dps,np.array(num_pl_ar)/1000 * (np.array(pss['y'])[-1] / (np.array(num_pl_ar)[-1]/1000))])).T
	dat.columns = ['p','num']

	# plot data
	plt.errorbar(dps,np.array(pss['y'])/scale,yerr=[(np.array(pss['y'])-np.array(pss['dybelow']))/scale,(np.array(pss['dyabove'])-np.array(pss['y']))/scale],marker='o',color='black')
	petys=np.array(kdat['y'])*(np.array(pss['y'])[-1]/scale/np.array(kdat['y'])[-1])
	plt.errorbar(kps,petys,yerr=[(np.array(kdat['y'])-np.array(kdat['dybelow']))*(np.array(pss['y'])[-1]/scale/np.array(kdat['y'])[-1]),(np.array(kdat['dyabove'])-np.array(kdat['y']))*(np.array(pss['y'])[-1]/scale/np.array(kdat['y'])[-1])],marker='o',color='red')
	plt.xscale('log')
	plt.yscale('log')
	if dist=='rayleigh':
		plt.loglog(dps,np.array(num_pl_ar) * (np.array(pss['y'])[-1] / np.array(num_pl_ar)[-1]),label='mcore dist = '+str(dist)+', mean = '+str(mean))
		#dat.to_csv('evaporation_'+str(dist)+'_mean='+str(mean)+'_output_HD_R_8.csv')
	elif dist=='lognormal':
		plt.loglog(dps,np.array(num_pl_ar) * (np.array(pss['y'])[-1] / np.array(num_pl_ar)[-1]),label='mcore dist = '+str(dist)+', mean = '+str(np.exp(mean))+', std = '+str(argv[0]))
		#dat.to_csv('evaporation_'+str(dist)+'_mean='+str(round(mean,1))+'_dev='+str(argv[0])+'_output_HD_R_8.csv')
	else:
		plt.loglog(dps,np.array(num_pl_ar) * (np.array(pss['y'])[-1] / np.array(num_pl_ar)[-1]),label='mcore dist = '+str(dist)+', index = '+str(argv[0]))
		#dat.to_csv('evaporation_'+str(dist)+'_index='+str(argv[0])+'_output.csv')
	#plt.show()

#evaporate('lognormal',np.log(14.0),0.5)
#evaporate('lognormal',np.log(19.0),0.5)
#evaporate('lognormal',np.log(17.0),0.5)

evaporate('lognormal',np.log(12.0),0.5)
evaporate('lognormal',np.log(10.0),0.5)
evaporate('lognormal',np.log(8.0),0.5)

#evaporate('rayleigh',15.)
#evaporate('rayleigh',12.)
#evaporate('rayleigh',10.)

#evaporate('lognormal',np.log(9.),0.65)
#evaporate('lognormal',np.log(9.5),0.65)
#evaporate('lognormal',np.log(10.),0.65)

#evaporate('lognormal',np.log(12.),0.65)
#evaporate('lognormal',np.log(12.),0.5)
#evaporate('lognormal',np.log(11.),0.5)

#evaporate('lognormal',np.log(10.0),0.7)

#evaporate('powerlaw',np.log(20.),1.6)
#evaporate('powerlaw',np.log(20.),1.8)
#evaporate('powerlaw',np.log(20.),1.4)

##
plt.text(30,0.1,'N2362; M$_{critical}$ from MESA HD')
plt.legend()
plt.show()

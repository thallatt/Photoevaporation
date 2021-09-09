#!/usr/bin/env python3
import numpy as np
from scipy import interpolate
import pandas as pd

def SMAXIS(Mss,Teq):
#INPUT: Mstar [Msun], Teq [K]
#OUTPUT: d0 [au]
    mstd = np.loadtxt('teq-sma.txt');
    mss00 = mstd[:,0];
    teq00 = mstd[:,1];
    d00 = mstd[:,2];
    #mstd=mstd1;
    Mss_g = [0.4, 0.6, 0.8, 1.0, 1.3];
    T_g = np.linspace(300, 3000, 28); #300:100:3000;

    d0_g = np.zeros( (len(Mss_g), len(Teq)));
    
    for i in range(len(Mss_g)):
	
        tmp = np.interp( Teq, T_g, d00[np.where(mss00 == Mss_g[i])]);
        d0_g[i,:] = tmp;
    #end for
    
    d0 = np.zeros( (len(Mss), 1));
    for i in range(len(Teq)):
        d0[i] = np.interp( Mss[i], Mss_g, d0_g[:,i]);
    #end for
	
    return d0

def BETA0(rp,mp,t0):
	#Re, Me, kb, mh, G0 required  #it is the gravitational parameter Lambda, we used to name it 'beta' until someone has pointed out that it's then looking as plasma parameter
	bb=G0*mp*Me/(rp*Re*(kb*t0/mh));
	return bb;

def unique_rows(a):
	a = np.ascontiguousarray(a)
	unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
	return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def ADEN(rp,mp):
	#Re, Me required
	aveg_den=3.*mp*Me/(rp*rp*rp*Re*Re*Re*3.1416*4);
	return aveg_den;
	
def ADEN_R(mp,aden0):
	#Re, Me required
	aveg_R=pow(3.*mp*Me/aden0/(Re*Re*Re*3.1416*4),1./3.);
	return aveg_R;
	
def lin_extpl(x1, x2, y1, y2, x):
	#to change!!!
	k = (y1-y2)/(x1-x2);
	b = y1-(k*x1);
	y = k*x + b;
	return y

##########################################################################################################################################

def Lborder(dau,fxuv,rpre):
	#INPUT: d0 [AU], Fxuv [erg/s/cm^2], Rpl [RE]

	zeta=-1.297796148718774+6.861843637445744;
	eta=0.884595403184073+0.009459807206476;

	beta=32.019929264625155 -16.408393523348366;
	alp1=0.422232254541188 -1;
	alp2=-1.748858849270155 +3.286179370395197;
	alp3=3.767941293231585 -2.75;

	K = (zeta + eta*np.log(dau)); 
	C = beta + alp1*np.log(fxuv) + alp2*np.log(dau) + alp3*np.log(rpre);

	ld= np.exp(C/K);
	return ld;

def testCF(dau,fxuv,rpre,ld):
	#INPUT: d0 [AU], Fxuv [erg/s/cm^2], Rpl [RE], LAMBDA
	zeta=-1.297796148718774;
	eta=0.884595403184073;

	beta=16.408393523348366 ;
	alp2=-3.286179370395197;
	alp3=2.75;

	K = (zeta + eta*np.log(dau)); 
	C = beta + np.log(fxuv) + alp2*np.log(dau) + alp3*np.log(rpre);

	mdot= np.exp(C + K*np.log(ld));
	return mdot;

def testCFJ(dau,fxuv,rpre,ld):
	#INPUT: d0 [AU], Fxuv [erg/s/cm^2], Rpl [RE], LAMBDA

	zeta=-6.861843637445744;
	eta=-0.009459807206476;

	beta=32.019929264625155;
	alp1=0.422232254541188;
	alp2=-1.748858849270155;
	alp3=3.767941293231585;

	K = (zeta + eta*np.log(dau)); 
	C = beta + alp1*np.log(fxuv) + alp2*np.log(dau) + alp3*np.log(rpre);

	mdot= np.exp(C + K*np.log(ld));
	return mdot;


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################


def INTERPOL(Mss,EUV,T_i,r_i,m_i):   #,NTEST
## INPUT: Mstar [Msun], EUV [erg/cm/s^2], Teq [K], Rpl [Re], Mpl [Me]
## import numpy as np
## this routine requires following functions (the ones in brackets are used only in the full version; now commented out):
## 		BETA0; (ROCHE); (BETA0_M); SMAXIS; (dissbar0); (lin_extpl); (Len_K)
## 
## OUTPUT: dotM, ...

	wild = 1;# 0 -- strictly prohibit extrapolation; 1 -- allow extrapolation (use only for evolution!)
	
	# set teq grid based on stellar mass
	Mssg0 = np.array([1.0]);#([0.4, 0.6, 0.8, 1.0, 1.3]);
	if(Mss>0.6):
		Tg0 = np.array([1100, 1500, 2000]);
	else:
		Tg0 = np.array([300, 700, 1100, 1500]);
	
	# set planet radius grid: 1-10 Re
	rg0 = np.linspace(1, 10, 10); 
	# set planet mass grid
	mg0 = np.array([1., 1.6, 2.1, 3.2, 4.3, 5., 6.7, 7.8, 9., 12.1, 16.2, 21.7, 29.1, 39., 45.1, 60.5, 81.0, 93.8, 108.6]);
	
	## first check of input parameters
	# check whether any inputs are beyond grid ranges
	use = ((Mss > max(Mssg0) or Mss<min(Mssg0)) or (T_i>max(Tg0) or T_i<min(Tg0)) or (r_i>max(rg0) or r_i<min(rg0)) or (m_i>max(mg0) or m_i<min(mg0)));
	# if yes, and extrapolation is turned off, shut down
	if(use and wild==0):
		print('Planet parameters are not in the grid range, 0');
		ITPL = -1;
		return ITPL;
	
	## uncomment the following if input is not set in the main file
	# read in hydro model data to interpolate between
	input0 = np.loadtxt('lhy_1Msun_1100-1500K_cut.dat');
	
	mss = input0[:,0];
	teq = input0[:,1];
	euv = input0[:,2];
	rpl = input0[:,3];
	mpl = input0[:,4];
	lhy = input0[:,5]; # = minus mdot

	# compute beta for all the hydro models
	bt = BETA0(rpl, mpl, teq); 

	drt=np.transpose(input0)
	dat = pd.DataFrame(np.array([drt[0],drt[1],drt[2],drt[3],drt[4],drt[5],bt])).T
	dat.columns = ['Mstar','Teq','FEUV','R','M','Mdot','lambda']
	#print(dat.loc[(dat['M']>20)&(dat['M']<40)&(dat['FEUV']>20000)&(dat['FEUV']<90000)])

	is__=np.where(euv==min(euv))

	#print(mss[is__])
	#print(teq[is__])
	#print(euv[is__])
	#print(rpl[is__])
	#print(mpl[is__])
	#print(lhy[is__])
	#print(bt[is__])
	## TODO: 
	
	###########[lhy,~,bt,~,rpl,mpl,mss,teq,euv,reff,tmax,nmin,vmax,rdis,ksimax]=loadallpoints();
	# compute beta for planet of interest
	beta_i = BETA0( r_i, m_i, T_i);#it is the gravitational parameter Lambda, we used to name it 'beta' until someone has pointed out that it's then looking as plasma parameter
	#d0_i=SMAXIS(Mss,T_i);
	
	## find closest M* and set Mss grid

	n_tmp, = np.where( abs(Mssg0-Mss) == 0);   # index where star mass is equal to mssg0 == array([0])

	# if mass does not match grid value, choose closest grid mass
	if(len(n_tmp) == 0):
		dmssg = np.abs(Mssg0-Mss); # difference between mssg0(=1) and star mass
		mssg = np.array(sorted(Mssg0[np.argsort(dmssg)][0:2])); # grid mass closest to user mass
	else: # otherwise, mssg=Mssg[0] == array([1.])
		mssg = Mssg0[n_tmp];

	# mass loss grid for stellar mass: np.linspace(0,0,1) == array([0.])
	lhy_ii_mss = np.linspace(0, 0, len(mssg));#np.zeros((len(mssg),1));
	
	# for imss in range(1): imss==0.
	for imss in range(len(mssg)):
	    ## reset Teq grid for specific Mss, recheck
	    # indices where hydro model stellar masses == 1.0 ([0,1,2,... all indices])
		usg, = np.where(mss == mssg[imss]); # [0,1,...1391]: all hydro grid models have mss=1
		# indices where grid eq. temperatures are within the hydro model eq. temperatures (inside 1100-2000 K)
		Tg1 = Tg0[np.where((Tg0<=max(teq[usg]))*(Tg0>=min(teq[usg])))]; # Tg1 = Tg0 == np.array([1100, 1500, 2000])
		# check if planet of interest eq. temperature is beyond grid
		use = (T_i>max(Tg1) or T_i<min(Tg1)); # T_i = teq from user; use = T or F
		# if so, and if extrapolation is off, shut down
		if(use and wild==0):
			print('Planet parameters are not in the grid range, 1teq');
			ITPL = -1;
			return ITPL; ## moved this inwards to match indent
		
	    ## find closest Teq and set Teq grid
		n_tmp, = np.where(np.abs(Tg1-T_i) == 0);   # index where user teq and grid are equal # empty 
		# if planet eq. temperature and grid temperatures are not the same
		if(len(n_tmp) == 0):
			dtg = np.abs(Tg1-T_i); # difference in grid teq and user teq
			tg = np.array(sorted(Tg1[np.argsort(dtg)][0:2]));# 2 grid teq closest to user teq for interpolation: == array([1100, 1500])
		else: # otherwise set tg to planet eq. temperature
			tg = Tg1[n_tmp]; # array[T_i]

		# mass loss grid in teq: linspace(1,2,2) == array([1., 2.])
		lhy_ii_teq = np.linspace(1, len(tg), len(tg));
	   
	   # for each eq. temperature in tg: iteq == 0,1
		for iteq in range(len(tg)):
			## check euv range
			# indices where hydro model stellar mass == 1., & where hydro model eq. temperature == tg (==[1100,1500])
			usg, = np.where((mss == mssg[imss])*(teq == tg[iteq]));
			# check if EUV within grid bounds
			use = EUV<min(euv[usg]) or EUV>max(euv[usg]);
			if(use and wild==0):
				print('Planet parameters are not in the grid range, 1teq');
				ITPL = -1;
				return ITPL;

			## set euv grid
			# list of unique euv values from hydro model grid, from teq=1100 and teq=1500 models
			# euv[usg] == euv values for models with same teq as iteq
			xuvg1 = sorted(set(euv[usg])); # 1100 K: [11.1,88.1,699.9,27247.5]], 1500: [38.5,304.6,4062.,15501.3,94215.2]

			if((not xuvg1) == False): # true
				xuvg = np.array(xuvg1); # array of euv values from model grid
				
				# mass loss grid for euv values: linspace(1,5,5) == array([1., 2., 3., 4., 5.])
				lhy_ii_xuv = np.linspace(1, len(xuvg), len(xuvg));
				
				# ixuv == 0,1,2,3,4
				# for each xuv in hydro model grid at boundary eq. temperatures
				for ixuv in range(len(xuvg)):
					# indices where hydro model stellar mass == 1., & where hydro model eq. temperature == a boundary teq (1100,1500) & where hydro model xuv the current unique xuv
					usg, = np.where((mss == mssg[imss])*(teq==tg[iteq])*(euv==xuvg[ixuv]));
					#interpol_b0 for bt(usg),rpl(usg)
					# indices where rpl == r_input & mpl = m_input
					n_tmp, = np.where((mss == mssg[imss])*(teq==tg[iteq])*(euv==xuvg[ixuv])*(r_i==rpl)*(m_i==mpl)); # array([]) for interpolation

					if((not n_tmp) and len(xuvg)>1): # if interpolating in R or M
						#print('R or M interpolation \n')
						# hydro model grid unique radii for models with boundary teq and current ixuv
						rplg = np.array(sorted(set(rpl[usg])));
						# mass loss grid for radii data points
						#print('init rplg='+str(rplg))
						rplg = np.asarray(rplg)[np.where(np.asarray(rplg)%1==0)[0]]
						#print('cleaned rplg='+str(rplg))
						lhy_ii_rpl = np.linspace(1, len(rplg), len(rplg)); # can throw errors if radius grid not linear ! 
						#print('initial lhy_rpl='+str(lhy_ii_rpl))
						irpl = np.min(rplg)-1; # == 1.0

						while (irpl < max(rplg) ): # while irpl < 10
							irpl += 1;
							# indices where hydro models have boundary teq & current ixuv & current irpl
							usg1, = np.where((mss == mssg[imss])*(teq==tg[iteq])*(euv==xuvg[ixuv])*(rpl==irpl));
							#print('teq='+str(tg[iteq]))
							#print('xuv='+str(xuvg[ixuv]))
							#print('usg1='+str(usg1))
							if(len(usg1)>1):
								lhy_ar = np.array(lhy[usg1])
								# interpolate mdot from beta of hydro models at usg1
								lhy_ar[np.where(lhy_ar==0.)]=1.e-5
								lhy_ii_rpl[int(irpl-2)] = np.exp(np.interp( np.log(beta_i), np.log(bt[usg1]), np.log(lhy_ar))); ## CHANGED TO -2 NOT -1
								#print(rpl[usg1])
								#print('intp beta;'+str(lhy_ii_rpl[int(irpl-2)]/1e12))
								#print('lhy_rpl now='+str(lhy_ii_rpl[int(irpl-2)]/1e12))
						rplg = rplg[np.where(lhy_ii_rpl>1.e-1)]; # changed 0 to 1e-1
						
						lhy_ii_rpl    = lhy_ii_rpl[np.where(lhy_ii_rpl>1.e-1)];
						# if input planet outside grid radii at current xuv, teq
						if( r_i>max(rplg) or r_i<min(rplg)):
							us, = np.where(np.diff(lhy_ii_rpl) > 0);

							if(r_i<min(rplg)):
								us1 = us[0:2];
							else:
								us1 = (us[-2], us[-1]);
							ind0 = 0*np.linspace(1,len(rplg),len(rplg));
							ind0[us] = 1;

							# indices where mdot increases with radius & where radius increases
							us = np.where((rplg>min(rplg))*(ind0));

							# interp function using only mdot vals that increase with radius
							par = np.polyfit(np.log(rplg[us]), np.log(lhy_ii_rpl[us]), 1);
							par_fit = np.poly1d(par);
							# mass loss grid for euv interpolated as fn of radius
							lhy_ii_xuv[ixuv] = np.exp(par_fit(np.log(r_i)));

						else:
							#print('rplg='+str(rplg)+', lhy_rpl = '+str(lhy_ii_rpl/1e12))
							lhy_ii_xuv[ixuv] = np.exp(np.interp( np.log(r_i), np.log(rplg), np.log(lhy_ii_rpl)));
							#print('intp r:'+str(lhy_ii_xuv[ixuv]/1e12))
					else:
						# if xuvg==0, mass loss grid at xuv = 0
						if(len(n_tmp) == 0):
							lhy_ii_xuv[ixuv]    = 0;
						else:
							# if planet on grid points, record mass loss rate in xuv grid where R,M match at current xuv
							lhy_ii_xuv[ixuv]    = lhy[n_tmp];

			xuvg = xuvg[np.where(lhy_ii_xuv > 0)];
			lhy_ii_xuv = lhy_ii_xuv[np.where(lhy_ii_xuv > 0)];
			# indices where xuv values are closest to input euv
			n_tmp, = np.where(np.abs(xuvg-EUV)<1.);

			## here set par_ii(Tg(iteq))
			# if there are no xuv values close to input euv
			if(len(xuvg)>1 and (not n_tmp)):
				# if EUV in grid at current teq
				if(EUV<=max(xuvg) and EUV>=min(xuvg)):
					# mass loss grid in teq interpolated
					lhy_ii_teq[iteq] = np.exp(np.interp( np.log(EUV), np.log(xuvg), np.log(lhy_ii_xuv)));
					#print('intp xuv:'+str(lhy_ii_teq[iteq]/1e12))
				else:
					#print('xuvs='+str(xuvg)+', ys='+str(lhy_ii_xuv))
					#print(EUV)
					par = np.polyfit(np.log(xuvg), np.log(lhy_ii_xuv), 1);
					par_fit = np.poly1d(par);
					#fn_intp = interpolate.interp1d(xuvg,lhy_ii_xuv,fill_value='extrapolate')
					lhy_ii_teq[iteq] = np.exp(par_fit(np.log(EUV)));
					#lhy_ii_teq[iteq] = fn_intp(EUV);
					#print('intp xuv:'+str(lhy_ii_teq[iteq]/1e12))
					#print('other intp:'+str(np.exp(par_fit(np.log(EUV)))/1e12))
			else:
				# if planet lies between grid xuv values and there are none in xuvg
				if(len(n_tmp) == 0):
					#print('no xuv vals; mdot=0')
					lhy_ii_teq[iteq]    = 0;
				else:
					# if planet is at grid point in xuv
					lhy_ii_teq[iteq]    = lhy_ii_xuv[n_tmp];
		
		## here define par(T_i) : 
		tg = tg[np.where(lhy_ii_teq>0)];
		lhy_ii_teq    = lhy_ii_teq[np.where(lhy_ii_teq>0)];
		# indices where grid teq equal input
		n_tmp, = np.where(tg==T_i);

	    #here set par_ii(Mssg(imss))
	    # if input teq in grid, not on grid point
		if( len(tg)>1 and (not n_tmp)):
			if(len(lhy_ii_teq)>0):
				if( T_i>700 and T_i <= max(tg)):#base on semi-axis
					tg1 = SMAXIS(np.ones((len(tg),1))*mssg[imss], tg); # distances in au
					tg1 = tg1[:,0]; # reshape
					T_i1 = SMAXIS([mssg[imss]], [T_i]);
					T_i1 = T_i1[0,0];
					# interpolate mass loss grid in stellar mass based on teq
					lhy_ii_mss[imss] = np.exp( np.interp( -np.log(T_i1), -np.log(tg1), np.log(lhy_ii_teq)));  
				# if T_i outside grid T range
				else:     
					lhy_ii_mss[imss] = np.exp( np.interp( T_i, tg, np.log(lhy_ii_teq)));
		else:
			# if there are no grid temperatures to interpolate between, and T_i doesn't match any grid points
			if(len(n_tmp) == 0):
				lhy_ii_mss[imss]    = 0;
			else:
				# if T_i matches a grid point just set mdot from grid
				lhy_ii_mss[imss]    = lhy_ii_teq[n_tmp];

	## here define par(Mss)
	# where mass loss grid at teq > 0
	mssg = mssg[np.where(lhy_ii_mss>0)];
	lhy_ii_mss    = lhy_ii_mss[np.where(lhy_ii_mss>0)];
	
	#  M=1; always true
	n_tmp, = np.where(mssg==Mss); 

	if(len(mssg)>1 and (not n_tmp)):
		Lhy_i    = np.interp( Mss, mssg, lhy_ii_mss);
	else:
		if(len(n_tmp) == 0):
			Lhy_i    = 0;
		else:
			Lhy_i    = lhy_ii_mss[n_tmp]; # mass loss grid in teq
	ITPL = Lhy_i;
	
	#print('mdot='+str(ITPL))
	
	#print(beta_i>Lborder(SMAXIS(np.array([1.]),np.array([T_i])),EUV,r_i))
	with open('temp_mdot.txt', 'w') as f:
		f.writelines(str(ITPL[0]))
	return ITPL
#end INTERPOL

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

## astronomical constants	
Lsun=3.9e33;    #solar luminocity	/*erg*s^{-1}*/
Rsun=6.96e10; #solar radius /*cm*/
AU=1.4960e13; #astronomical unit /*cm*/
Re=6.378e8; #earth radius /*cm*/	
Me=5.9722e27; #earth mass /*g*/
kb=1.3807e-16;	#boltzman const /*erg*K^{-1}*/
mh=1.6726e-24;  #hydrogen mass /*g*/
G0=6.6726e-8;   #gravitational constant /* cm^3*g^{-1}*s^{-2}*/

#how to use approximation
# if regime=='analytical_approximation':
	# lb=Lborder(d0,EUV0,Rpl0);
	# if (BETA0(Rpl0,Mpl,T_i)<=lb):
		# Mdot=testCFJ(d0,EUV0,Rpl0,BETA0(Rpl0,Mpl,T_i));#d0
	# else:
		# Mdot=testCF(d0,EUV0,Rpl0,BETA0(Rpl0,Mpl,T_i));#d0
# else:
	# if(regime == 'interpolation'):
		# Mdot = INTERPOL(Mss, EUV0, T_i, Rpl0, Mpl);
	# else:
		# print('Only regimes app and int are supported. Please, set one of those regimes.')
		# return -1

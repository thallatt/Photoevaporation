#
import math
import numpy as np
import os
import shutil
import time
import random
import sys
import inlist_maker as inlist
from lum_evolution import lum_evo
import pandas as pd

#some constants
msun = 1.98855e33
rsun = 6.9598e10
mjup = 1.8986e30
rjup = 7.14e9 #6.9911e9
Lsun = 3.9e33

sigma=5.67e-5
au = 1.496e13
mearth = 5.9722e27
rearth = 6.378e8

###################################################################################################################################

# fiducial parameters
mcore = float(sys.argv[1])*mearth			# total planet mass in gramm
gcr = float(sys.argv[2])			# GCR
mp = mcore + mcore * gcr

rinitial = 3.0 				# inital planet radius in rjup (formation routine)
minitial=1.5316e+29 #                   # initial planet mass in creation phase in gramm (this one fits for planets below   25Mearth; for heavier planets one may use the jupiter mass (0.001Msun.mod) initial model from "mesa-r12115/star/test_suite/irradiated_planet")

fat_0 = (mp-mcore)/mp				# initial atm. mass fraction in planetary masses
rhocore = inlist.ADEN_core(mp/mearth) 	# core density in g/cm^3 (uses rcore=mcore^0.27 approximation for numbers from Rogers et al., 2011)

z = 0.02                                # metallicity of both planet and star
y = 0.24                                # helium fraction of both planet and (initial) star

maxage = 5.e9                          # final age in years
mesh_delta_coeff = 0.5                  # coefficient defining the number of grid points in MESA runs. 0.5 is a good init.guess (smaller the coefficient - more points in grid)

maxEntropy = 9.				# entropy of inflated atmosphere in kb/baryon (init. thermal state; step iv in the paper)
coreLuminosity = 2.0e27			# backup value; ~works for ~4-30 ME; the true value is set in the setS procedure
coolingTime = 1.e7                      # cooling timescale in years

mp_wo_core = mp - mcore			# mass of the atmosphere in gramm; not to change


# parameters having to do with irradiation
ms = 1.# 				# star mass in msun
rs = 1.#				# star radius in rsun
Teff_star = 5780. #			# stellar Teff in K
orb_sep = (float(sys.argv[3])/365.)**(2./3) #				# orbital sepration of the planet in AU
#BA= 0.0		                # planet Bond albedo
period_star = float(sys.argv[4])       # stellar rotation period
thermal_ = bool(str(sys.argv[5])=='True')         # with or without mass loss

initage = 1e7 #                         # initial age when mass loss begins in years (the disc dispersal time)

irrad_col = 300.0			# column depth for depositing stellar radiation as heat
flux_dayside = 0.                       # dayside flux set later in relax_irradiation (Lbol/4./pi/orb_sep[cm]**2)
                                        # (only relevant for relaxing the initial irradiation)

#parameters having to do with atmospheric escape
escape_model = str(sys.argv[6])                # 'HD' stands for hydrodynamic HBA approximation (see Kubyshkina et al., 2018b)
                                   # 'EL' stands for energy limited with eta = 0.15 (to change eta, change x_ctrl(2) in inlist_7_evolve_mass_loss)
                                   # 'ELR' stands for energy limited + radiation-recombination limited escape as described in Chen&Rogers, 2016

##################################################################################################################################

# flags to undertake/skip steps
do_create_planet = False 	#create a ball of gas
do_put_in_core = bool(str(sys.argv[7])=='True') 		#put the core of the right size
do_relaxm = bool(str(sys.argv[8])=='True') 		#remove the atmosphere to the desired value
do_set_entropy = bool(str(sys.argv[9])=='True')		#set an artificial luminocity to inflate the planet
do_cool = bool(str(sys.argv[10])=='True') 			#remove an artificial luminocity and cool the planet for coolingTime
do_relax_irradiation =  bool(str(sys.argv[11])=='True')	#set stellar heating (Teq); 
do_evolve_planet = bool(str(sys.argv[12])=='True') 	#evolve!
##################################################################################################################################


## load MIST track (for irradiation)
Mssg_list = list(range(40,95,5))+list(range(92,132,2));#[40:5:90 92:2:130];    #GRID OF STELLAR MASSES AVAILABLE in MIST ( 0.4..1.3 MSUN)
Mssg = np.array(Mssg_list);

n_tmp, = np.where(abs(ms*100-Mssg)==min(abs(ms*100-Mssg))); #take the closest available mass; these tracks change quite smooth, so it is accurate enough. However, one can also interpolate between two closest tracks (if you need assistance with it, please contact me)
n_tmp  = n_tmp[0];

if Mssg[n_tmp]>=100:
    nameEtrack = 'ev_inp/ei'+str(Mssg[n_tmp])+'.dat';
    nameEtrack1 = 'ev_inp/eho'+str(Mssg[n_tmp])+'.dat';
else:
    nameEtrack = 'ev_inp/ei0'+str(Mssg[n_tmp])+'.dat';
    nameEtrack1 = 'ev_inp/eho0'+str(Mssg[n_tmp])+'.dat';
		
Etrack = np.loadtxt(nameEtrack) #np.loadtxt('ev_inp/ei050.dat') #

#############################################################################################################################
#SETTING NAMES FOR THE MODELS (you may wish to modify them for your convenience)
createmodel = "planet_1_create_" + str(minitial/mearth)[0:5]+ "_ME"  + ".mod"
createmodel_big = "planet_1_create_100.0_ME.mod"
coremodel = "planet_2_core_" + str(mcore/mearth)[0:5] + "_ME_" + str(mp/mearth)[0:6] + "_ME" + ".mod"
relaxedmodel = "planet_3_relxdM_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6]  + ".mod"
entropymodel = "planet_4_setS_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6] + "_" +str(maxEntropy)  + ".mod"
removemodel = "planet_5_removeLc_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6] + "_" +str(maxEntropy)  + ".mod"
relaxirradmodel = "planet_6_relaxirrad_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy)+"_d_"+ str(round(orb_sep,3))+ "_" +str(ms)  + "Msun.mod"
evolvedmodel = "planet_7_evolved_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy) +"_d_"+ str(round(orb_sep,3))+ "_" +str(ms)  + "Msun_" +str(period_star)+ "Pstar_"+str(escape_model)+"Escape.mod"
evolvedmodel_thermal = "planet_7_evolved_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy) +"_d_"+ str(round(orb_sep,3))+ "_" +str(ms)  + "Msun_" +str(period_star)+ "Pstar_Thermal.mod"
#############################################################################################################################
#BEGIN

if (mp - mp * fat_0)/mearth > 25.:
    createmodel=createmodel_big
    print('Mcore > 25. Using alternative createmod.')


shutil.copyfile("src/run_star_extras0.f","src/run_star_extras.f")#otherwise everything goes too slow
os.system('./clean')

# don't need to run this - already have basic model
if do_create_planet:
    inlist1 = "tmp_inlists/inlist_1_create_" + str(mp/mearth)[0:6] + "_ME"
    
    run_time = inlist.create_planet(minitial, y, z, inlist1, createmodel)
    success = True
    if not os.path.exists(createmodel):
        success=False	
    k=open('LOGS/history_1_create','r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "final model number in create=",last[0])
    #print "last[0]==1000",last[0]=="1000"
    if last[0]=="1000":
        success=False
    print( "step 1, Planet created: ", success)

if do_put_in_core:
    inlist2 = "tmp_inlists/inlist_2_core_" + str(mcore/mearth)[0:5] + "_ME_" + str(mp/mearth)[0:6] + "_ME"
    
    run_time = inlist.put_core_in_planet(mesh_delta_coeff, mcore, rhocore, inlist2, createmodel, coremodel)
    success = True
    if not os.path.exists(coremodel):
        success=False	
    k=open('LOGS/history_2_core','r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "final age in core setting =", str(10.**float(last[0]))[0:6], "year")
    #print "last[0]==1000",last[0]=="1000"
    
    print( "step 2, Core inserted: ", success)

if do_relaxm:
    inlist3 = "tmp_inlists/inlist_3_reducem_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6]

    run_time = inlist.relaxm(mesh_delta_coeff, mp,inlist3,coremodel,relaxedmodel)
    success = True
    if not os.path.exists(relaxedmodel):
        success=False	
    k=open('LOGS/history_3_reducemass','r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "final age in fat setting =", str(10.**float(last[0]))[0:6], "year")
    #print "last[0]==1000",last[0]=="1000"
    
    print( "step 3, Atmosphere reduced to fat0 = ", fat_0,": ", success)

if do_set_entropy:	
    with open('LOGS/history_3_reducemass', 'r') as f:
        for line in f:
            pass
        last = line.split()
		
    currentropy= float(last[7])

    if currentropy<float(maxEntropy):
        with open('LOGS/history_3_reducemass', 'r') as f:
            for line in f:
                pass
            last2 = line.split()
		
        if maxEntropy<13.6:
            coreLuminosity= 30*float(last2[1])*3.846e33
        print( coreLuminosity)
        inlist4 = "tmp_inlists/inlist_4_setS_"  + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6] +"_"+ str(maxEntropy)
        run_time = inlist.set_initial_entropy(mesh_delta_coeff, maxEntropy,coreLuminosity,inlist4,relaxedmodel,entropymodel)
        success = True
        if not os.path.exists(entropymodel):
            success=False	
        k=open('LOGS/history_4_setS','r')
        for line in k.readlines():
            pass
        last_temp=line
        last=last_temp.split()
        print( "final entropy set =", last[7])
        #print "last[0]==1000",last[0]=="1000"
    
        print( "step 4, Initial entropy set to = ", maxEntropy,": ", success)
        


if do_cool:
    Lcore_low = mcore*5e-8 #energy generation rate of present Earth; Doesn't matter, will be anyway overtaken by core model in ev.
    inlist5 = "tmp_inlists/inlist_5_removeLc_"  + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6] +"_"+ str(maxEntropy)
    run_time = inlist.remove_Lcore(mesh_delta_coeff, coolingTime, Lcore_low, inlist5,entropymodel,removemodel)
    success = True
    if not os.path.exists(removemodel):
        success=False	
    k=open('LOGS/history_5_removeLc','r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "final age in cooling =", str(10.**float(last[0]))[0:6], "year")
    #print "last[0]==1000",last[0]=="1000"
    
    print( "step 5, Artificial core dissipation removed: ", success)

if do_relax_irradiation:    
    
    Lbol = np.interp(coolingTime, Etrack[:,0], Etrack[:,1]);# [log10(Lsun)]
    Lbol = (10.**(Lbol))*Lsun; #[erg/s]

    Teff = np.interp(coolingTime, Etrack[:,0], Etrack[:,2]);# [log10(K)]
    Teff = 10.**(Teff); #[K]

    Rss  = np.interp(coolingTime, Etrack[:,0], Etrack[:,3]);# [log10(Rsun)]
    Rss  = 10.**(Rss); #[Rsun]

    Teq = Teff*((Rss*rsun/2./orb_sep/au)**0.5); #TEQ(Teff,Rss,d0,0);
    flux_dayside = Lbol/4.0/3.1416/orb_sep/orb_sep/au/au
    relaxage = coolingTime + 1e5
    initage = relaxage
    print('Teq='+str(Teq))
    
    inlist6 = "tmp_inlists/inlist_6_relaxirrad_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy)+"_d_"+ str(orb_sep)+ "_" +str(ms)  + "Msun"
    run_time = inlist.relax_irradiation(mesh_delta_coeff, irrad_col,flux_dayside, relaxage, inlist6,removemodel,relaxirradmodel)

    success = True
    if not os.path.exists(relaxirradmodel):
        success=False	
    k=open('LOGS/history_6_relaxsurfheat','r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "disk dissipation time set to", str(10.**float(last[0]))[0:6], "year")
    #print "last[0]==1000",last[0]=="1000"
    
    print( "step 6, Heating switched on: ", success)
   

if do_evolve_planet:
    
    inlist7 = "tmp_inlists/inlist_7_evolve_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy)+"_d_"+ str(orb_sep)+ "_" +str(ms)  + "Msun"
    
    a_x,b_x,a_euv,b_euv,t_sat = lum_evo(period_star)
    
    run_time = inlist.evolve_planet(mesh_delta_coeff, 1.e7, orb_sep, maxage, nameEtrack1, escape_model, t_sat, a_x, b_x, a_euv, b_euv, inlist7, relaxirradmodel, evolvedmodel)

    success = True
    if not os.path.exists(evolvedmodel):
        success=False	
    k=open('LOGS/history' + evolvedmodel[6:-4],'r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "the planet has evolved for ", str((10.**float(last[0]))/1e9)[0:6], "Gyr")
    print("planet final radius @ "+str(maxage/1.e9)+" Gyr : "+str(last[19][0:6])+" Rearth.")
    print("planet final mass @ "+str(maxage/1.e9)+" Gyr : "+str(round(float(last[18]),3))+" Mearth.")
    print("% of initial atmosphere lost at "+str(maxage/1.e9)+" Gyr : "+str(round(np.abs(100.* ((mcore * gcr / mearth)-(float(last[18])-mcore/mearth))/(mcore * gcr / mearth)),2))+" %.")
    print( "step 7, Evolution complete: ", success)
    print( "atmospheric escape model: ", escape_model)
    print( "Pstar = ", str(period_star), " days")
    print( "evolutionary profiles saved in ", '"LOGS/history' + evolvedmodel[6:-4] + '"')

    if escape_model=='HD':
        dat=pd.read_pickle('mesa_HD.pkl')
        dat.loc[len(dat)]=[mcore/mearth,gcr,float(sys.argv[3]),period_star,last[19][0:6],last[19][0:6],round(float(last[18]),3),round(np.abs(100.* ((mcore * gcr / mearth)-(float(last[18])-mcore/mearth))/(mcore * gcr / mearth)),2)]
        print(dat)
        dat.to_pickle('mesa_HD.pkl')
    elif escape_model=='EL':
        dat=pd.read_pickle('mesa_EL.pkl')
        dat.loc[len(dat)]=[mcore/mearth,gcr,float(sys.argv[3]),period_star,last[19][0:6],last[19][0:6],round(float(last[18]),3),round(np.abs(100.* ((mcore * gcr / mearth)-(float(last[18])-mcore/mearth))/(mcore * gcr / mearth)),2)]
        print(dat)
        dat.to_pickle('mesa_EL.pkl')

if thermal_:
    
    inlist8 = "tmp_inlists/inlist_7_evolve_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy)+"_d_"+ str(orb_sep)+ "_" +str(ms)  + "Msun_Thermal"
    
    a_x,b_x,a_euv,b_euv,t_sat=0.,0.,0.,0.,1.
    print('evolving thermal only.')
    
    run_time = inlist.evolve_planet(mesh_delta_coeff, 1.e7, orb_sep, maxage, nameEtrack1, escape_model, t_sat, a_x, b_x, a_euv, b_euv, inlist8, relaxirradmodel, evolvedmodel_thermal)

    success = True
    if not os.path.exists(evolvedmodel_thermal):
        success=False	
    k=open('LOGS/history' + evolvedmodel_thermal[6:-4],'r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "the planet has evolved for ", str((10.**float(last[0]))/1e9)[0:6], "Gyr")
    print("planet final radius @ "+str(maxage/1.e9)+" Gyr : "+str(last[19][0:6])+" Rearth.")
    print( "step 7, Evolution complete: ", success)
    print( "thermal evolution only.")
    print( "Pstar = ", str(period_star), " days")
    print( "evolutionary profiles saved in ", '"LOGS/history' + evolvedmodel_thermal[6:-4] + '"')
    
    if escape_model=='HD':
        dat=pd.read_pickle('mesa_HD.pkl')
        dat_0 = dat.loc[(dat['coremass']==mcore/mearth)&(dat['gcr']==gcr)&(dat['period']==float(sys.argv[3]))&(dat['Pstar']==period_star)]
        dat['Rf_5Gyr']=dat['Rf_5Gyr'].replace(dat_0['Rf_5Gyr'],(last[19][0:6]))
        print(dat)
        dat.to_pickle('mesa_HD.pkl')
    elif escape_model=='EL':
        dat=pd.read_pickle('mesa_EL.pkl')
        dat_0 = dat.loc[(dat['coremass']==mcore/mearth)&(dat['gcr']==gcr)&(dat['period']==float(sys.argv[3]))&(dat['Pstar']==period_star)]
        dat['Rf_5Gyr']=dat['Rf_5Gyr'].replace(dat_0['Rf_5Gyr'],(last[19][0:6]))
        print(dat)
        dat.to_pickle('mesa_EL.pkl')

#
import math
import numpy as np
import os
import shutil
import time
import random
import sys

msun = 1.98855e33
rsun = 6.9598e10
mjup = 1.8986e30
rjup = 7.14e9 #6.9911e9
Lsun = 3.9e33

sigma=5.67e-5
au = 1.496e13
mearth = 5.9722e27
rearth = 6.378e8

# make the initial planet without the coree

def create_planet(minitial, y, z, inlist1, createmodel):
        start_time = time.time()
        print( "create initial planet")
        f = open('inlist_1_create', 'r')
        g = f.read()
        f.close()
        g = g.replace("<<initial_mass>>", str(minitial))
        g = g.replace("<<z>>",str(z))
        g = g.replace("<<y>>",str(y))
        g = g.replace("<<smwtfname>>", '"' + createmodel + '"')
        h = open(inlist1, 'w')
        h.write(g)
        h.close()
        shutil.copyfile(inlist1,"inlist")
        #os.system('./clean')
        os.system('./mk')
        os.system('./rn')
        run_time = time.time() - start_time
        print( "run time for create_planets in sec = ",run_time)
        return run_time

def put_core_in_planet(mesh_delta_coeff, mcore, rhocore, inlist2, createmodel, coremodel):
        start_time = time.time()
        print( "put core in planet")
        f = open('inlist_2_core', 'r')
        g = f.read()
        f.close()
        g = g.replace("<<loadfile>>",'"' + createmodel + '"')
        g = g.replace("<<smwtfname>>", '"' + coremodel + '"')
        g = g.replace("<<new_core_mass>>", str(mcore/msun))
        g = g.replace("<<core_avg_rho>>", str(rhocore))
        g = g.replace("<<mesh_delta_coeff>>", str(mesh_delta_coeff))
        #g = g.replace("<<core_radius>>",str(coreRadius)) this changes nothing, just the second equal way to intrduce a core
        h = open(inlist2, 'w')
        h.write(g)
        h.close()
        shutil.copyfile(inlist2,"inlist")
        os.system('./mk')
        os.system('./rn')
        run_time = time.time() - start_time
        print( "run time to put in core in sec = ",run_time)
        return run_time

def relaxm(mesh_delta_coeff, mp,inlist3,coremodel,relaxedmodel):
        start_time = time.time()
        print( "reducing atmosphere mass")
        print( "target mass = ",str(mp/mearth), "Mearth")
        f = open('inlist_3_reducemass', 'r')
        g = f.read()
        f.close()
        g = g.replace("<<loadfile>>",'"' + coremodel + '"')
        g = g.replace("<<smwtfname>>", '"' + relaxedmodel + '"')
        g = g.replace("<<new_mass>>",str(mp/msun))
        g = g.replace("<<mesh_delta_coeff>>", str(mesh_delta_coeff))
        h = open(inlist3, 'w')
        h.write(g)
        h.close()
        shutil.copyfile(inlist3,"inlist")
        os.system('./mk')
        os.system('./rn')
        run_time = time.time() - start_time
        print( "run time to remove atmosphere = ",run_time)
        return run_time

def set_initial_entropy(mesh_delta_coeff, maxEntropy,coreLuminosity,inlist4,relaxedmodel,entropymodel):
        start_time = time.time()
        print( "setting initial entropy")
        print( "target entropy=",str(maxEntropy))
        f = open('inlist_4_setS', 'r')
        g = f.read()
        f.close()
        g = g.replace("<<loadfile>>",'"' + relaxedmodel + '"')
        g = g.replace("<<smwtfname>>", '"' + entropymodel + '"')
        g = g.replace("<<Lc>>", str(int(coreLuminosity/1e20)*1e20))
        g = g.replace("<<entropy>>",str(maxEntropy))
        g = g.replace("<<mesh_delta_coeff>>", str(mesh_delta_coeff))
        h = open(inlist4, 'w')
        h.write(g)
        h.close()
        shutil.copyfile(inlist4,"inlist")
        os.system('./mk')
        os.system('./rn')
        run_time = time.time() - start_time
        print( "run time to set entropy = ",run_time)
        return run_time

def remove_Lcore(mesh_delta_coeff, coolingTime, Lcore_low, inlist5,entropymodel,removemodel):
        start_time = time.time()
        print( "remove core dissipation")
        f = open('inlist_5_removeLc', 'r')
        g = f.read()
        f.close()
        g = g.replace("<<loadfile>>",'"' + entropymodel + '"')
        g = g.replace("<<smwtfname>>", '"' + removemodel + '"')
        g = g.replace("<<Lcore_new>>",str(Lcore_low))
        g = g.replace("<<maxage>>",str(coolingTime))
        g = g.replace("<<mesh_delta_coeff>>", str(mesh_delta_coeff))
        #g = g.replace("<<flux_dayside>>", str(flux_dayside) )
        #g = g.replace("<<irrad_col>>", str(irrad_col) )
        #g = g.replace("<<entropy>>",str(lowEntropy))
        h = open(inlist5, 'w')
        h.write(g)
        h.close()
        shutil.copyfile(inlist5,"inlist")
        os.system('./mk')
        os.system('./rn')
        run_time = time.time() - start_time
        print( "run time to evolve in sec = ", run_time)
        return run_time

def relax_irradiation(mesh_delta_coeff, irrad_col, flux_dayside, relaxage, inlist6, removemodel, relaxirradmodel):
        start_time = time.time()
        print( "relax surface heating")
        f = open('inlist_6_relaxsurfheat', 'r')
        g = f.read()
        f.close()
        g = g.replace("<<loadfile>>",'"' + removemodel + '"')
        g = g.replace("<<smwtfname>>", '"' + relaxirradmodel + '"')
        g = g.replace("<<col_depth>>", str(irrad_col))
        g = g.replace("<<irrad_flux>>", str(flux_dayside))
        g = g.replace("<<maxage>>", str(relaxage))
        g = g.replace("<<mesh_delta_coeff>>", str(mesh_delta_coeff))
        #g = g.replace("<<initage>>", '.' + str(initage) + '.')
        h = open(inlist6, 'w')
        h.write(g)
        h.close()
        shutil.copyfile(inlist6,"inlist")
        os.system('./mk')
        os.system('./rn')
        run_time = time.time() - start_time
        print( "run time to evolve in sec = ", run_time)
        return run_time

def relax_irradiation_b(mesh_delta_coeff, irrad_col, flux_dayside, rplfit, inlist6b, relaxedmodel, relaxirradmodel):
        start_time = time.time()
        print( "relax surface heating")
        f = open('inlist_6b_young_planet', 'r')
        g = f.read()
        f.close()
        g = g.replace("<<loadfile>>",'"' + relaxedmodel + '"')
        g = g.replace("<<smwtfname>>", '"' + relaxirradmodel + '"')
        g = g.replace("<<col_depth>>", str(irrad_col))
        g = g.replace("<<irrad_flux>>", str(flux_dayside))
        g = g.replace("<<rplfit>>", str(rplfit))
        g = g.replace("<<mesh_delta_coeff>>", str(mesh_delta_coeff))
        #g = g.replace("<<initage>>", '.' + str(initage) + '.')
        h = open(inlist6b, 'w')
        h.write(g)
        h.close()
        shutil.copyfile(inlist6b,"inlist")
        os.system('./mk')
        os.system('./rn')
        run_time = time.time() - start_time
        print( "run time to evolve in sec = ", run_time)
        return run_time

def evolve_planet(mesh_delta_coeff, initage, orb_sep, maxage, nameEtrack1, escape_model, t_sat, a_x, b_x, a_euv, b_euv, pstar, inlist7, relaxirradmodel, evolvedmodel):
        start_time = time.time()
        print( "evolve escaping atmosphere")
        f = open('inlist_7_evolve_mass_loss', 'r')
        g = f.read()
        f.close()
        g = g.replace("<<loadfile>>",'"' + relaxirradmodel + '"')
        g = g.replace("<<smwtfname>>", '"' + evolvedmodel + '"')
        
        g = g.replace("<<initage>>", str(initage))
        g = g.replace("<<orb_sep>>", str(orb_sep))
        g = g.replace("<<maxage>>", str(maxage))
        g = g.replace("<<history_name>>", '"history' + evolvedmodel[6:-4] + '"')
        g = g.replace("<<mesh_delta_coeff>>", str(mesh_delta_coeff))
        h = open(inlist7, 'w')
        h.write(g)
        h.close()
        shutil.copyfile(inlist7,"inlist")

        f = open('src/run_star_extras_evol.f', 'r')
        k = f.read()
        f.close()
        k = k.replace("<<escape_model>>", escape_model)
        #k = k.replace("<<escape_Cpow>>", str(escape_Cpow))
        #k = k.replace("<<escape_betapow>>", str(escape_betapow))
        #k = k.replace("<<escape_sat>>", str(escape_sat))
        print('a_euv:'+str(a_euv))
        print('tsat:'+str(int(t_sat))+'d6')
        #print('new tsat:'+str(round(t_sat,0)))
        print('new b_euv:'+str(round(b_x*0.58,1)))
        k = k.replace("<<t_sat>>", str(int(t_sat))+'d6')
        k = k.replace("<<a_x>>", str(a_x))
        k = k.replace("<<b_x>>", str(round(b_x,1)))
        k = k.replace("<<a_euv>>", str(round(a_euv,5)))
        k = k.replace("<<b_euv>>", str(round(b_euv,1)))#round(b_x*0.575,1)))
        k = k.replace("<<heating_input>>", '"' +  nameEtrack1 + '"')
        k = k.replace("<<pstar>>", str(pstar).split('.')[0])
        h = open('src/run_star_extras_temp.f', 'w')
        h.write(k)
        h.close()
        shutil.copyfile('src/run_star_extras_temp.f','src/run_star_extras.f')

        os.system('./clean')
        os.system('./mk')
        os.system('./rn')
        run_time = time.time() - start_time
        print( "run time to evolve in sec = ", run_time)
        return run_time

def ADEN_core(mp):
        #Input: Mcore [Mearth]
        #Output: average density in g/cm^3
        rp = mp**0.27
        aveg_den=3.*mp*mearth/(rp*rp*rp*rearth*rearth*rearth*3.1416*4.);
        return aveg_den;


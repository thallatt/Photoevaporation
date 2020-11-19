#!/usr/bin/env python3 

# send submission script to CC for single job

import numpy as np
import os
import os.path
import subprocess

pathset='/home/thallatt/projects/def-evelee/thallatt/kub_mesa_local/'

# parameter grid
coremasses_beluga = [4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.]
coremasses_cedar = [26.,28.,30.,32.,34.,36.,38.,40.,42.]
coremasses_graham = [28.]
gcrs_ = [0.97]
period_star_ = 5.0
periods_ = [1.3, 2.4, 4.2, 7.5, 13.3, 23.6, 42.0, 75.5, 133.0, 239.0]

def write_file(cpu, coremasses, gcrs, period_star, period):
	os.chdir('/Users/Tim/Documents/Work/Msc/code/')
	# switches for creation - no evolution
	do_put_in_core = 'True'
	do_relaxm = 'True'
	do_set_entropy = 'True'
	do_cool = 'True'
	do_relax_irradiation = 'True'
	do_evolve_planet = 'True'
	thermal_ = 'True'

	for itr in range(0,len(coremasses)):		
		for itr_ in range(0,len(gcrs)):

			# file names
			file_name='evolve_CM'+str(coremasses[itr])+'_GCR'+str(gcrs_[itr_])+'_orbper'+str(period)+'_Pstar'+str(period_star)+'.sh'
			completeName = os.path.join('/Users/Tim/Documents/Work/Msc/code/remote/submission_scripts/', file_name)
			
			#file_name_glost='glost_create.txt'
			#completeName_glost = os.path.join('/Users/Tim/Documents/Work/Msc/code/remote/glost_lists/', file_name_glost)
	
			subprocess.run(["rm", str(completeName)]) # delete old version of files
			#subprocess.run(["rm", str(completeName_glost)])
	
			# load files
			fle = open(completeName,'w')	
			#fle_glost = open(completeName_glost,'w')
			
			# set up submission file
			fle.write('#!/bin/bash \n')
			fle.write('#SBATCH --time=02:00:00 \n')
			fle.write('#SBATCH --nodes=1 \n')
			fle.write('#SBATCH --ntasks-per-node=1 \n')
			fle.write('#SBATCH --account=def-evelee \n')
			fle.write('#SBATCH --mem=5120M \n')
			fle.write('module load scipy-stack/2019a \n')
			fle.write('module load StdEnv/2020 \n')
			fle.write('source ~/pyenv/bin/activate \n')

			# current planet
			cm = coremasses[itr]
			GCR = gcrs[itr_]
			#totalmass = coremasses[itr] + gcrs[itr_] * coremasses[itr]
			#fenv = gcrs[itr_] * coremasses[itr] / totalmass

			fle.write('python /home/thallatt/projects/def-evelee/thallatt/kub_mesa_local/bake_your_planet.py '+str(cm)+' '+str(GCR)+' '+str(period)+' '+str(period_star)+' '+str('True')+' '+str('HD')+' '+str(do_put_in_core)+' '+str(do_relaxm)+' '+str(do_set_entropy)+' '+str(do_cool)+' '+str(do_relax_irradiation)+' '+str(do_evolve_planet))

			print(str(cm)+', '+str(GCR))			
			print('submission file: '+str(completeName))
			#print('glost file: '+str(completeName_glost))
			fle.close()
			#fle_glost.close()
			subprocess.run(["scp", str(completeName), str(cpu)+":"+str(pathset)])
			subprocess.run(["scp", str(completeName), str(cpu)+":"+str(pathset)])
			#subprocess.run(["scp", str(completeName_glost), str(cpu)+":"+str(pathset)])
			#subprocess.run(["scp", str(completeName_glost), str(cpu)+":"+str(pathset)]) # do twice since last run doesn't work

# write files
#write_file('cedar', periods_cedar)
write_file('graham', coremasses_graham, gcrs_, period_star_,4.2)
#write_file('beluga', periods_beluga)

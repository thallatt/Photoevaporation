#!/usr/bin/env python
# hard-coded data from MESA calculations

import numpy as np
import consts
import pandas as pd

pss=pd.read_csv('/Users/Tim/Documents/Work/MSc/petigura_saturns.csv')
pss = pss.loc[(pss['x']>=0) & (pss['x']<=200.)]
dps=np.asarray(pss['x'])

kdat=pd.read_csv('/Users/Tim/Documents/Work/MSc/kunimoto_subsaturns.csv')
kdat = kdat.loc[(kdat['x']>=0) & (kdat['x']<=140.)]
kps=np.asarray(kdat['x'])

# critical core masses
# 8 Re
mcrit_1day = [66.,49.,35.,25.5,18.5,13.5,9.5,6.,1.]
mcrit_3day = [60.,42.,30.,22.,16.,11.5,8.0,1.,1.]
mcrit_5day = [55.,38.,27.,20.,14.5,10.5,7.25,1.,1.]
mcrit_10day  = [50.,35.,25.,18.,13.,9.5,6.5,1.,1.]
mcrit_20day = [48.,32.,23.,16.5,12.,8.5,6.,1.,1.]
# 5 Re
mcrit_1day_5  = [71.,54.,42.,31.,22.,17.,13.,1.,1.]
mcrit_3day_5  = [67.,47.,35.,26.,19.,15.,9.75,1.,1.]
mcrit_5day_5  = [60.,42.,31.,23.5,17.,12.,9.,1.,1.]
mcrit_10day_5  = [55.,39.,29.,21.,15.,11.0,7.5,1.,1.]
mcrit_20day_5  = [55.,39.,29.,21.,15.,11.0,7.5,1.,1.] # needs work
# Elim
# 8 Re
mcrit_1day_el = [54.,34.,22.5,16.5,12.05,10.3,8.,1.,1.]
mcrit_1day_5_el = [60.,53.,35.,24.,18.,12.7,8.,1.,1.]
mcrit_5day_el = [33.,22.,16.5,11.5,8.6,6,1.,1.,1.] # needs work
mcrit_5day_5_el = [48.,32.,23.,16.5,18.,12.7,8.,1.,1.]
mcrit_10day_el  = [28.5,18.75,13.48,10.615,8.5,6.,1.,1.,1.]
mcrit_10day_5_el  = [41.,27.,19.,14.,11.5,8.5,1.,1.,1.]
#mcrit_20day = [30.,17.,11.,8.,6.,3.,1.,1.,1.] # needs work

# both
mcrit_1day_ = [54.,34.,22.5,20.,18.5,13.5,9.5,6.,1.]
mcrit_1day_5_ = [60.,53.,35.,24.,18.,17.,13.,1.,1.]
mcrit_5day_ = [33.,23.,20.,18.5,14.5,10.5,7.25,1.,1.]
mcrit_5day_5_ = [48.,32.,31.,23.5,17.,12.,9.,1.,1.]
mcrit_10day_  = [28.5,25.,20.,18.,13.,9.5,6.5,1.,1.]
mcrit_10day_5_  = [41.,35.,29.,21.,15.,11.0,7.5,1.,1.]

# 4 Re @ 5 Gyr
days = np.array([5.,12.,23.])
cms__ = np.array([10.,20.,30.,40.])
gcrs_10_4 = np.array([0.03, 0.08,0.1])
Rpl_10_4 = np.array([4.8,5.15,4.9])
gcrs_20_4 = np.array([0.04,0.075,0.1])
Rpl_20_4 = np.array([4.6,4.7,4.74])
gcrs_30_4 = np.array([0.05,0.06,0.1])
Rpl_30_4 = np.array([4.5,4.4,4.8])
gcrs_40_4 = np.array([0.0375, 0.07,0.075])
Rpl_40_4 = np.array([4.29,4.53,4.46])
# 6 Re @ 5 Gyr
gcrs_10_6 = np.array([0.07, 0.17,0.25])
Rpl_10_6 = np.array([8.2,9.15,9.5])
gcrs_20_6 = np.array([0.15,0.25,0.3])
Rpl_20_6 = np.array([7.24,7.5,7.7])
gcrs_25_6 = np.array([0.15,0.25,0.3])
Rpl_25_6 = np.array([6.8,7.34,7.41])
gcrs_40_6 = np.array([0.15, 0.25,0.3])
Rpl_40_6 = np.array([6.3,6.91,7.11])
# 8 Re @ 5 Gyr
gcrs_10_8 = np.array([0.25, 0.5,0.8])
Rpl_10_8 = np.array([13.0,13.5,12.2])
gcrs_20_8 = np.array([0.5,0.8,1.0])
Rpl_20_8 = np.array([10.4,10.5,10.1])
gcrs_30_8 = np.array([0.4,0.8,1.0])
Rpl_30_8 = np.array([9.1,10.0,9.7])
gcrs_40_8 = np.array([0.45,0.8,1.0])
Rpl_40_8 = np.array([8.9,9.4,9.7])

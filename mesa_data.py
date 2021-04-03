#!/usr/bin/env python
# hard-coded data from MESA calculations

import numpy as np
import consts
import pandas as pd

pss=pd.read_csv('/Users/Tim/Documents/Work/MSc/petigura_saturns.csv')
pss = pss.loc[(pss['x']>=0) & (pss['x']<=200.)]
dps=np.asarray(pss['x'])
print(dps)
kdat=pd.read_csv('/Users/Tim/Documents/Work/MSc/kunimoto_subsaturns.csv')
kdat = kdat.loc[(kdat['x']>=0) & (kdat['x']<=140.)]
kps=np.asarray(kdat['x'])

# critical core masses
# 8 Re
mcrit_1day = [66.,49.,35.,25.5,18.5,13.5,9.5,6.,1.]
mcrit_3day = [60.,42.,30.,22.,16.,11.5,8.0,1.,1.]
mcrit_5day = [55.,38.,27.,20.,14.5,10.5,7.25,1.,1.] # check 27, 10.5
mcrit_10day  = [50.,35.,25.,18.,13.,9.5,6.5,1.,1.]
mcrit_20day = [48.,32.,23.,16.5,12.,8.5,6.,1.,1.]
# 5 Re
mcrit_1day_5  = [71.,54.,42.,31.,22.,17.,13.,8.5,1.] # 75
mcrit_3day_5  = [67.,47.,35.,26.,19.,15.,9.75,6.75,1.]
mcrit_5day_5  = [60.,42.,31.,23.5,17.,12.,9.,6.,1.]
mcrit_10day_5  = [55.,39.,29.,21.,15.,11.0,7.5,1.,1.]
mcrit_20day_5  = [55.,39.,29.,21.,15.,11.0,7.5,1.,1.] # needs work
# Elim
mcrit_1day_el = [54.,34.,22.5,16.5,12.05,10.3,6.5,1.,1.]
mcrit_1day_5_el = [60.,53.,35.,24.,18.,12.7,8.5,1.,1.]
mcrit_5day_el = [33.,22.,16.5,11.5,8.6,7.,5.25,1.,1.] # 23
mcrit_5day_5_el = [48.,32.,23.,16.5,11.75,8.5,6.5,1.,1.]
mcrit_10day_el  = [28.5,18.75,13.48,10.615,8.5,6.,4.,1.,1.]
mcrit_10day_5_el  = [41.,27.,19.,14.,11.5,8.5,6.,1.,1.]

# HD 8 Re; King
mcrit_1day_king = [70.,52.,37.75,27.,20.4,14.25,10.,6.5,1.] # 75 d
mcrit_5day_king = [62.,44.,32.,23.5,16.25,12.,8.7,1.,1.] # 1.3
mcrit_10day_king = [60.,41.,29.,21.,15.,10.9,7.5,1.,1.]
# EL 8 Re; King
mcrit_1day_king_el = [62.,39.,25.5,18.,13.,9.7,7.25,4.,1.]
mcrit_5day_king_el = [45.,28.,19.,13.5,10.,7.5,5.75,1.,1.]
mcrit_10day_king_el = [40.,25.,17.,12.5,9.4,7.25,5.25,1.,1.]

## HD 5 Re; King
mcrit_1day_king_5 = [70.,60.,44.,32.,24.5,18.75,13.5,9.,1.] # 13.3, 23.6, 72 (will need fraction of core mass)
mcrit_5day_king_5 = [60.,50.,36.5,27.,20.,14.5,10.,7.,1.] # 1.3, 4.2
mcrit_10day_king_5 = [59.,47.,35.,25.,18.,13.,9.,6.5,1.]
## EL 5 Re; King
mcrit_1day_king_5_el = [89.,56.,41.,28.,20.,13.,10.,7.,1.] # 56 needs work
mcrit_5day_king_5_el = [67.,42.,29.,21.,14.,10.,8.,1.,1.] # 21 needs work
mcrit_10day_king_5_el = [62.,39.,26.,17.5,12.5,9.,7.,1.,1.] # 1.3 needs work

# HD + EL
mcrit_1day_ = [54.,34.,22.5,20.,18.5,13.5,9.5,6.,1.]
mcrit_1day_5_ = [60.,53.,35.,31.,22.,17.,13.,8.5,1.]
mcrit_5day_ = [33.,23.,20.,18.5,14.5,10.5,7.25,1.,1.]
mcrit_5day_5_ = [48.,32.,31.,23.5,17.,12.,9.,6.,1.]
mcrit_10day_  = [28.5,25.,20.,18.,13.,9.5,6.5,1.,1.]
mcrit_10day_5_  = [41.,35.,29.,21.,15.,11.0,7.5,1.,1.]

# highmass
mcrit_1day_hm = [44.,31.,24.,20.5,17.85,13.5,9.5,6.,1.]
mcrit_1day_5_hm = [46.5,36.5,34.,31.,22.,17.,13.,8.5,1.]# 46.5 model issues
mcrit_5day_hm = [33.5,26.75,21.,18.5,14.5,10.5,7.25,1.,1.] # check 26.75
mcrit_5day_5_hm = [45.,35.25,31.,23.5,17.,12.,9.,6.,1.] # 31 or 32
mcrit_10day_hm = [33.5,26.,20.5,18.,13.,9.5,6.5,1.,1.]
mcrit_10day_5_hm = [41.,35.,28.,21.,15.,11.0,7.5,1.,1.]

mcrit_1day_hm_J = [37.,29.,22.5,20.5,14.75,10.5,7.4,5.,4.] # 2.4 days
mcrit_1day_5_hm_J = [45.5,35.75,33.,25.,18.2,13.,9.,6.,4.]
mcrit_5day_hm_J = [34.5,27.,20.8,17.,12.5,8.75,6.5,4.4,3.5] # 42
mcrit_5day_5_hm_J = [44.,35.25,29.,21.,15.,10.9,7.25,5.0,3.5] # 42
mcrit_10day_hm_J = [34.,26.35,20.35,16.,11.5,8.2,5.75,4.3,3.4] # 42
mcrit_10day_5_hm_J = [41.,35.,27.,20.,14.,9.95,6.9,4.9,3.4] # 42

mcrit_1day_ELP_J = [27.+0.75,17.75+0.75,12.75+0.75,10.+0.75,8.3+0.75,6.75+0.75,5.9,1.,1.] # approximated
mcrit_5day_ELP_J = [27.+0.3,17.75+0.3,13.5,10.3,8.35,7.,5.75,1.,1.] # approximated
mcrit_10day_ELP_J = [27.,17.75,12.75,10.,8.3,6.75,5.65,1.,1.] # 10.7->10 -- local # 8.3 corrected
mcrit_10day_EL_J = [27.,17.75,12.25,9.,6.8,5.5,1.,1.,1.]

mcrit_5day_5_ELP_J = [28.7,19.3,14.5,12.75,9.8,7.8,5.75,1.,1.] # approximated
mcrit_10day_5_ELP_J = [28.,18.5,13.75,11.75,9.,7.5,5.5,1.,1.] # approximated

# for 200d normalization
# 5 Re
# [45.5,35.75,33.,25.,18.2,13.,9.,6.,4.,0.5]
# [44.,35.25,29.,21.,15.,10.9,7.25,5.0,3.5,0.5]
# [41.,35.,27.,20.,14.,9.95,6.9,4.9,3.4,0.5]
# 8 Re
# [37.,29.,22.5,20.5,14.75,10.5,7.4,5.,4.,0.5]
# [34.5,27.,20.8,17.,12.5,8.75,6.5,4.4,3.5,0.5]
# [34.,26.35,20.35,16.,11.5,8.2,5.75,4.3,3.4,0.5]

# 2 Gyr, highmass
# 1,5,10 days -- very similar
mcrit_10day_2gyr = [29.5,22.,14.5,10.,7.5,1.,1.,1.,1.]
mcrit_10day_5_2gyr = [36.,26.,18.5,12.,8.5,1.,1.,1.,1.]

# 2Gyr < 4 days, highmass
mcrit_10day_2gyr_hm = [29.5,22.,14.5,18.,13.,9.5,6.5,1.,1.]
mcrit_10day_5_2gyr_hm = [36.,26.,18.5,21.,15.,11.0,7.5,1.,1.]

# 2 Gyr, EL
mcrit_10day_2gyr_el = [14.,11.,8.,1.,1.,1.,1.,1.,1.] # needs work
mcrit_10day_5_2gyr_el = [22.,15.,11.,8.,1.,1.,1.,1.,1.]

# HD + EL; King
mcrit_1day_king_ = [62.,39.,25.5,19.,17.5,14.25,10.,6.5,1.] # 19, 17.5 give 40 M_Earth
mcrit_1day_5_king_ = [89.,56.,41.,32.,24.5,18.75,13.5,9.,1.]
mcrit_5day_king_ = [45.,28.,20.,19.,16.25,12.,8.7,1.,1.] # 20, 19 give 40 M_Earth
mcrit_5day_5_king_ = [67.,42.,36.5,27.,20.,14.5,10.,7.,1.]
mcrit_10day_king_  = [40.,25.,21.,21.,15.,10.9,7.5,1.,1.] # 21 at 4.2 days gives 40 M_Earth
mcrit_10day_5_king_  = [62.,39.,35.,25.,18.,13.,9.,6.5,1.] # 1 needs work

# highmass king
mcrit_1day_king_hm = [47.5,34.5,26.,21.5,18.,14.25,10.,6.5,1.] # 75 d
mcrit_1day_king_5_hm = [48.,37.,34.25,32.,24.5,18.75,13.5,9.,1.] # check 37
mcrit_5day_king_hm = [43.,31.,23.25,19.5,16.25,12.,8.7,1.,1.] # 31 - 4.3 Re
mcrit_5day_king_5_hm = [48.,36.5,33.5,27.,20.,14.5,10.,7.,1.] # 1.3, 4.2
mcrit_10day_king_hm = [38.,30.,22.5,19.,15.,10.9,7.5,1.,1.] # 30
mcrit_10day_king_5_hm = [47.,36.25,33.,25.,18.,13.,9.,6.5,1.]

# highmass Johnstone; stellar masses; 8 Rearth
mcrit_10day_1msol = [34.,26.35,20.35,16.,11.5,8.2,4.,1.,1.] # 42, 2.4--check 26.5
mcrit_10day_08msol = [25.,19.75,16.5,12.,8.6,6.,1.,1.,1.] # 2.4 & 42
mcrit_10day_12msol = [39.,33.5,27.,20.,13.5,9.5,6.75,1.,1.] # 42

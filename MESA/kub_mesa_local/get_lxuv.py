import numpy as np

def get_lxuv(time_,per_):
	fn=np.load('flxuv_1Msol_Pstar='+str(per_).split('.')[0]+'days.npy',allow_pickle=True)
	fn=fn.item()
	lxuv=float(fn(time_))
	
	with open('temp_lxuv.txt', 'w') as f:
		f.writelines(str(lxuv))

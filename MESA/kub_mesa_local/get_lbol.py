import numpy as np

def get_lbol(time_,per_):
	fn=np.load('flbol_1Msol_Pstar='+str(per_).split('.')[0]+'days.npy',allow_pickle=True)
	fn=fn.item()
	lbol=float(fn(time_))

	with open('temp_lbol.txt', 'w') as f:
		f.writelines(str(lbol))

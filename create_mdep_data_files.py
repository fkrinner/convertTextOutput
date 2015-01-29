from convertTextOutput import readTextFile, getBestFits
import numpy as np
import numpy.linalg as la


def get_raw_data_coma(
			fileName, 
			full_waveset):
	"""Reads a textoutput-file and returns the data and coma correspondign to the waveset"""
	waveSet = []
	for wave in full_waveset:
		waveSet.append(wave.strip())
	rawData = readTextFile(fileName)
	nEvents = rawData[0]['nevents']
	rawRawData = []
	for wave in waveSet:
		for key in rawData[0].iterkeys():
			if wave == key.strip():
				rawRawData.append(rawData[0][key])		

	coma = []
	for i in range(len(waveSet)): # Factor 2 from re/im
		coma.append([0.]*(2*len(waveSet)))
		coma.append([0.]*(2*len(waveSet)))
	data = [0.]*(2*len(waveSet))

	for i in range(len(waveSet)):
		point = rawRawData[i]
		data[2*i  ] = point[0]*(float(nEvents)**.5)
		data[2*i+1] = point[1]*(float(nEvents)**.5)
		iii = point[2]
		for j in range(len(rawRawData)):
			jjj = rawRawData[j][2]
			coma[2*i  ][2*j  ] = rawData[1][2*iii  ][2*jjj  ]*nEvents
			coma[2*i+1][2*j  ] = rawData[1][2*iii+1][2*jjj  ]*nEvents
			coma[2*i  ][2*j+1] = rawData[1][2*iii  ][2*jjj+1]*nEvents
			coma[2*i+1][2*j+1] = rawData[1][2*iii+1][2*jjj+1]*nEvents
	return [data,coma]

def own_pinv(matrix, minev =0.):
	val,vec = la.eig(matrix)
	dim = len(val)
	new = np.zeros((dim,dim))
	for i in range(dim):
		if abs(val[i]) > minev:
			new[i,i] = 1/val[i]
	ret = np.dot(vec,np.dot(new,np.transpose(vec)))
	return ret


def invert_coma(matrix,method = 'pinv'):
	"""Does 'inversion' of covariance matrices"""
	dim = len(matrix)
	if method == 'pinv':
		dim = len(matrix)
		matrix_np = np.zeros((dim,dim))
		for i in range(dim):
			for j in range(dim):
				matrix_np[i,j] = matrix[i][j]
#		matrix_pinv = la.pinv(matrix)
		matrix_pinv = own_pinv(matrix_np,1.E-6)
		print matrix_pinv.shape
		ret = []
		for i in range(dim):
			retLine = []
			for j in range(dim):
				retLine.append(matrix_pinv[i,j])
			ret.append(retLine)
		return ret
	if method == 'diag':
		ret = []
		for i in range(dim):
			ret.append([0.]*dim)
		for i in range(dim):
			if not matrix[i][i] == 0.:
				ret[i][i] = 1./matrix[i][i]
		return ret
	raise ValueError # Unknown method

def get_SDM_data_coma(raw_data_coma):
	"""Creates SpinDensityMatrix data and coma from the amplitudes"""
	raw_data = raw_data_coma[0]
	raw_coma = raw_data_coma[1]
	nWaves = len(raw_data)/2
	SDM = [0.]*nWaves**2
	jacobian = []
	for i in range(2*nWaves):
		jacobian.append([0.]*nWaves**2)

	for i in range(nWaves):
		for j in range(nWaves):
			one_index = i*nWaves+j
			if i==j: # Intensities
				SDM[one_index] = raw_data[2*i  ]**2 + raw_data[2*i+1]**2
				jacobian[2*i  ][one_index] = 2*raw_data[2*i  ]
				jacobian[2*i+1][one_index] = 2*raw_data[2*i+1]
			if i<j: # Real part
				SDM[one_index] = raw_data[2*i  ]*raw_data[2*j  ] + raw_data[2*i+1]*raw_data[2*j+1]
				jacobian[2*i  ][one_index] =   raw_data[2*j  ]
				jacobian[2*j  ][one_index] =   raw_data[2*i  ]
				jacobian[2*i+1][one_index] =   raw_data[2*j+1]
				jacobian[2*j+1][one_index] =   raw_data[2*i+1]
			if i>j: # Imag Part
				SDM[one_index] = - raw_data[2*j  ]*raw_data[2*i+1] + raw_data[2*i  ]*raw_data[2*j+1]
				jacobian[2*i  ][one_index] =   raw_data[2*j+1]
				jacobian[2*j  ][one_index] = - raw_data[2*i+1]
				jacobian[2*i+1][one_index] = - raw_data[2*j  ]
				jacobian[2*j+1][one_index] =   raw_data[2*i  ]
	SDM_coma = []
	for i in range(nWaves**2):
		SDM_coma.append([0.]*nWaves**2)

	for i in range(nWaves**2):
		for j in range(2*nWaves):
			for k in range(2*nWaves):
				for l in range(nWaves**2):
					SDM_coma[i][l]+=jacobian[j][i]*jacobian[k][l]*raw_coma[j][k]
	return [SDM,SDM_coma]



def write_files(
			target_dir, 
			source_dirs, 
			waveset, 
			upper_lims, 
			lower_lims, 
			method = 'pinv',
			name_base_data = 'data_', 
			name_base_coma = 'coma_'):
	"""Writes data and coma files for fit results in the source_dirs (one for each t' bin) for waveset to  target_dir"""
	i=0
	for source_dir in source_dirs:
		write_files_tbin(
					str(i),
					target_dir, 
					source_dir, 
					waveset, 
					upper_lims,
					lower_lims, 
					method = method, 
					name_base_data = name_base_data, 
					name_base_coma = name_base_coma)
		i+=1

def write_files_tbin(
			indentifier,
			target_dir,
			source_dir,
			waveset,
			upper_lims,
			lower_lims,
			method = 'pinv',
			name_base_data = 'data_', 
			name_base_coma = 'coma_'):
	"""Writes data and coma files for one t' bin"""
	file_list = getBestFits(source_dir)
	file_list.sort(key = lambda x:x[1])
	nWaves = len(waveset)	
	nFile = 0
	nWrite= 0
	with open(target_dir+'/'+name_base_data+indentifier,'w') as out_data:
		with open(target_dir+'/'+name_base_coma+indentifier,'w') as out_coma:
			for fn in file_list:
				nFile+=1
				m_low = fn[1]
				m_up  = fn[2]
				m = (m_up+m_low)/2.
				SDM_data_coma = get_SDM_data_coma(get_raw_data_coma(fn[0], waveset))
				for i in range(nWaves): # Apply mass limits
					for j in range(nWaves):
						index = i*nWaves+j
						if upper_lims[i] < m or upper_lims[j] < m or lower_lims[i] > m or lower_lims[j] > m:
							SDM_data_coma[0][index] = 0.
							for index2 in range(nWaves**2):
								SDM_data_coma[1][index ][index2] = 0.
								SDM_data_coma[1][index2][index ] = 0.
				for val in SDM_data_coma[0]:
					out_data.write(str(val)+' ')
				out_data.write('\n')
				inverted = invert_coma(SDM_data_coma[1],method)
				written = 0
				for i in range(nWaves**2):
					for j in range(nWaves**2):
						written +=1
						if not inverted[i][j].imag == 0.:
							print i,j,inverted[i][j],"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
							raise ValueError # Non Real Coma
						out_coma.write(str(inverted[i][j].real)+' ')
						nWrite+=1
				out_coma.write('\n')

if __name__ == "__main__":
	sources = [	"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.100000-0.112853",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.112853-0.127471",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.127471-0.144385",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.144385-0.164401",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.164401-0.188816",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.188816-0.219907",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.219907-0.262177",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.262177-0.326380",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.326380-0.448588",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.448588-0.724294",
			"/nfs/mds/user/fhaas/PWA/fits/2008-all-lH-DT0/PHD/best_11tbins/0.724294-1.000000"]

	waves = [	"1-(1++)0+ rho pi S",
			"1-(1-+)1+ rho pi P",
			"1-(1++)0+ rho pi D",
			"1-(0-+)0+ f0(980) pi S",
			"1-(2++)1+ f2 pi P",
			"1-(2++)1+ rho pi D",
			"1-(2++)2+ rho pi D",
			"1-(2-+)0+ f2 pi D",
			"1-(2-+)0+ f2 pi S",
			"1-(2-+)0+ rho pi F",
			"1-(2-+)1+ f2 pi S",
			"1-(4++)1+ f2 pi F",
			"1-(4++)1+ rho pi G"]

	uppers = [2.3,2.0,2.1,2.3,2.0,2.0,2.0,2.3,2.3,2.1,2.3,2.3,2.3 ]
	lowers = [0.9,0.9,0.9,1.2,1.0,0.9,1.0,1.6,1.4,1.2,1.4,1.4,1.25]


	write_files(	"/nfs/mds/user/fkrinner/data_mdep_fit/testFiles",
			sources,
			waves,
			uppers,
			lowers,
			"diag",
			"testDataOld_",
			"testComaOld_")
	


	

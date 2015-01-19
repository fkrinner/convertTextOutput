from convertTextOutput import getBestLike, readTextFile, getIntegralMatrixAverage
import numpy as np
import numpy.linalg as la
from random import random
import minuit

def getData(path,compares):
	"""calculates the amplitudes with the correct coefficients for the unphysical modes"""
	intFolder = path.replace("/fit/","/integrals/").split("text_fit_")[0]
	bestLikeFile = getBestLike(path)
	rawData = readTextFile(path+'/'+bestLikeFile)
	comparedata = []
	for comparepath in compares:
		comparedata.append(readTextFile(comparepath+'/'+getBestLike(comparepath)))

	mMin = rawData[0]["m3Pi"][0]
	mMax = rawData[0]["m3Pi"][1]
	rawIntegrals = getIntegralMatrixAverage(mMin,mMax,intFolder,False,False)
	isoBases = getDeIsobarredList(rawData)
	for key in isoBases.iterkeys(): # Integral indices are usually the same as fit indices, but keep them anyway to avoid later errors
		for ppt in isoBases[key]:
			ppt.append(rawIntegrals[0][ppt[2]])
	keyList = []
	for key in isoBases.iterkeys():
		keyList.append(key)
	oneRow=[]
	starts=[]
	for key in keyList:
		starts.append(len(oneRow))
		for ppt in isoBases[key]:
			oneRow.append(ppt[:])
			wave = ppt[2]
			for com in comparedata:
				try:
					aa = com[0][wave]
					oneRow[-1].append(aa[0])
					oneRow[-1].append(aa[1])
					break
				except KeyError:
					pass
			if not len(oneRow[-1]) == 9:
				raise IndexError # No compare-vector found.
#			print oneRow[-1]
	goodMatrix=[]

	dim = len(oneRow)
	rowResult=[complex(ppt[3],ppt[4])for ppt in oneRow]	
	comResult=[complex(ppt[7],ppt[8])for ppt in oneRow]
	for i in range(len(oneRow)):
		goodMatrix.append([0.+0.j]*len(oneRow))
	normalizers = []
	global_int_norm = 0.
	GLOBAL_NORM = False
	for i in range(len(oneRow)):
		try:
			normalizers.append(1./rawIntegrals[1][oneRow[i][6]][oneRow[i][6]])
		except ZeroDivisionError:
			normalizers.append(0.+0.j)
		for j in range(len(oneRow)):
			try:
				if not GLOBAL_NORM:
					goodMatrix[i][j] = rawIntegrals[1][oneRow[i][6]][oneRow[j][6]]/(rawIntegrals[1][oneRow[i][6]][oneRow[i][6]]*rawIntegrals[1][oneRow[j][6]][oneRow[j][6]])**.5
				else:
					goodMatrix[i][j] = rawIntegrals[1][oneRow[i][6]][oneRow[j][6]]
			except ZeroDivisionError:
				goodMatrix[i][j] = 0.+0.j

	if GLOBAL_NORM:
		global_int_norm = 0.
		for i in range(len(goodMatrix)):
			global_int_norm+=goodMatrix[i][i]
		for i in range(len(goodMatrix)):
			for j in range(len(goodMatrix)):
				goodMatrix[i][j]/=global_int_norm
	
	eigen = eigenBasis(goodMatrix)
	coeff = decompose(rowResult,eigen)
	normTo = 0 # number of de-isobarred wave to normalize 
	iiii = []
	vecs = []
	DONE = True
	valolo=[]
	for i in range(dim):
#		print coeff[i],
		valolo.append(eigen[i][0].real)
		if eigen[i][0] < 1E-2 and not eigen[i][0] == 0.:
			print eigen[i][0]
#			print "<<<"
			ccc0 = eigen[i][1][starts[normTo]]
			rrr0 = rowResult[starts[normTo]]
#			print rrr0/ccc0,"<<<<<<<<<<<<<<<<<<<<<<<<<<<"
			if not DONE:
				print i,'>>>>(',(rrr0/ccc0).real,',',(rrr0/ccc0).imag,')<<<<',i
				coeff[i] -= rrr0/ccc0
				DONE = True
			vecs.append(eigen[i][1][:])
			iiii.append(i)
		else:
			pass
#			print "   "
#	with open("./de_isobarred_corrected/vec",'w') as outoutout:
#		for i in range(dim):
#			outoutout.write(str(i))
#			for j in range(len(vecs)):
#				outoutout.write(' '+str(vecs[j][i].real))
#			outoutout.write('\n')
	if not len(vecs)==0:
		coeffs = get_physical_coefficients_intens(rowResult, comResult, vecs, starts)
		print coeffs
	#	coaff = [complex(0.,0.)]*len(coeff)
	#	coeff[2] -= complex(0.532412199699,0.33230841749)
	#	coeff[2] -= complex(0.596338098319,0.372208168472)
		for i in range(len(iiii)):
			coeff[iiii[i]]-=coeffs[i]


	reconstruct = [0.+0.j]*dim
	for i in range(dim):
		for j in range(dim):
			reconstruct[i]+=coeff[j]*eigen[j][1][i]
	with open("./de_isobarred_corrected/no_modes",'w') as outoutout:
		for i in range(dim):
			outoutout.write(str(i)+' '+str((abs(rowResult[i])**2*normalizers[i]).real)+' '+str((abs(reconstruct[i])**2*normalizers[i]).real)+' '+str(abs(comResult[i])**2*normalizers[i].real) +'\n')
	valolo.sort()
	with open("de_isobarred_corrected/val",'w') as outoutout:
		for i in range(len(valolo)):
			outoutout.write(str(i)+' '+str(valolo[i])+'\n')
#	print starts


def get_physical_coefficients_intens(result, compare, vectors, start = None): # use start just for interfacing, not needed with intenistites
	"""Fits the coefficients to the compare-points. At the moment a workaraound for 2 or four parameters"""
	if len(vectors) == 1:
		print "One coefficient"
		def f(x,y):
			return chi22(result,compare, vectors,[complex(x,y)],start)
	elif len(vectors) ==2:
		print "Two coefficients"
		def f(x1,y1,x2,y2):
			return chi22(result,compare, vectors,[complex(x1,y1),complex(x2,y2)],start)
	elif len(vectors) > 2:
		raise ValueError # More than two vectors not supported
	m = minuit.Minuit(f)
	m.migrad()
	print "Minuits finished:"
	print "Chi2 =",m.fval
	print m.values
	if len(vectors) ==1:
		return [complex(m.values['x'],m.values['y'])]
	elif len(vectors) ==2:
		return [complex(m.values['x1'],m.values['y1']),complex(m.values['x2'],m.values['y2'])]	



def chi22(result, compare, vectors, coeff, start = None):
	"""Gives a Chi2 for the unphysical modes to minimize"""
	chi2 = 0.
	for i in range(len(result)):
#	for i in range(start[2],len(result)):
		act = result[i]
		for j in range(len(vectors)):
			act -= coeff[j] * vectors[j][i]
		chi2+=(abs(act)**2 - abs(compare[i])**2)**2

	return chi2

	


def get_physical_coefficients(result, compare, vectors, start):
	"""Performs th fit of the parameters"""
	dim  = len(result )
	nComp= len(start)	
	nVec = len(vectors)
	singleComp = []
	for i in range(nComp):
		singleComp.append([complex(0.,0.)]*dim)
	actComp = -1
	for i in range(dim):
		if i in start:
			actComp+=1
		singleComp[actComp][i] = compare[i]
#	for sc in singleComp:
#		print sc


	total_vectors=[]
	for vec in vectors:
		total_vectors.append(vec[:])
	for vec in singleComp:
		total_vectors.append(vec[:])
	coeff = get_coefficients(result,total_vectors)
#	print coeff ,"PPPPPPPPPPPPPPPPPPPPP"
#	coeff = [complex(0.,0.)]*nVec #get_coefficients(result,total_vectors)[:nVec]
	return coeff

def get_coefficients(data, fitvecs, errs = None):
	"""Amplitude fit, not needed anymore"""
	dim = len(data)
	nVec= len(fitvecs)
	if not errs:
		errs=[complex(1.,1.)]*dim
	normchi2 = 0.
	for i in range(dim):
		normchi2+=(data[i].real/errs[i].real)**2+(data[i].imag/errs[i].imag)**2
	coeff = [complex(0.,0.)]*nVec
	coeff[0] = complex(0.532412199699,0.33230841749)
	step_size = 0.0001
	store=coeff[:]
	chi2old = 100000000000.
	nochange =0
	while True:
		for i in range(nVec):
			rrand = 2*random()-1
			irand = 2*random()-1
			coeff[i]+=step_size*complex(rrand,irand)

#			coeff[i]*=complex(2*random()-1,2*random()-1)
#			coeff[i]+=complex(2*random()-1,2*random()-1)
#			coeff[0] = complex(0.532412199699,0.33230841749)
		chi2 = 0.
		for i in range(dim):
			diff = data[i]
			for j in range(nVec):
				diff -= fitvecs[j][i]*coeff[j]
			chi2+=(diff.real/errs[i].real)**2
			chi2+=(diff.imag/errs[i].imag)**2
		if chi2 < chi2old:
			print chi2
			chi2old = chi2
			store = coeff[:]
			nochange =0
		else:
			nochange+=1
			if nochange == 10000:
				break

	print "final chi2:", chi2old, "without",normchi2
	return coeff



def eigenBasis(matrix):
	"""Gets the eigenbasis of the matrix [[val0,[vec0]],[val1,[vec1]],...]"""
	dim = len(matrix)
	aa = np.zeros((dim,dim),dtype = np.complex128)
	eigenBasis=[]
	for i in range(dim):
		for j in range(dim):
			aa[i,j] = matrix[i][j]	

	val,vec = la.eig(aa)
	with open("./de_isobarred_corrected/evs",'w') as outoutout:
		for i in range(dim):
			outoutout.write(str(i)+' '+str(val[i].real)+'\n')

	eigen = []
	for i in range(dim):
		eigen.append([val[i],[vec[j,i] for j in range(dim)]])
	return eigen

def decompose(vals,eigen):
	"""Decomposes vals into the eigenbasis"""
	dim = len(vals)
	coeff = []
	for i in range(dim):
		cc = 0.
		for j in range(dim):
			cc+= vals[j]*eigen[i][1][j]
		coeff.append(cc)
	return coeff


def getDeIsobarredList(rawData):
	"""Gets a list of de-isobarred waves in the data"""
	isoBases = {}
	for key in rawData[0].iterkeys():
		if 'f0_' in key or 'rho_' in key or 'f2_' in key:
			actBase = key.split('_')[0]
			mMinAct = int(key.split('_')[1])
			mMaxAct = int(key.split('_')[2][:4])
			if not isoBases.has_key(actBase):
				isoBases[actBase] = []
			isoBases[actBase].append([mMinAct,mMaxAct,key,rawData[0][key][0],rawData[0][key][1],rawData[0][key][2]])
	for key in isoBases.iterkeys():
		isoBases[key].sort()

	return isoBases
		



if __name__ == "__main__":
	getData("/nfs/mds/user/fkrinner/massIndepententFits/fits/0pp_1mm_2pp_in1pp_MC/fit/0.14077-0.19435/text_fit_0.14077-0.19435.dat_1.50000_1_1500_1540_0141_0194",
		[	"/nfs/mds/user/fkrinner/massIndepententFits/fits/0pp_in1pp_MC/fit/0.14077-0.19435/text_fit_0.14077-0.19435.dat_1.50000_1_1500_1540_0141_0194",
			"/nfs/mds/user/fkrinner/massIndepententFits/fits/1mm_in1pp_MC/fit/0.14077-0.19435/text_fit_0.14077-0.19435.dat_1.50000_1_1500_1540_0141_0194",
			"/nfs/mds/user/fkrinner/massIndepententFits/fits/2pp_in1pp_MC/fit/0.14077-0.19435/text_fit_0.14077-0.19435.dat_1.50000_1_1500_1540_0141_0194"])

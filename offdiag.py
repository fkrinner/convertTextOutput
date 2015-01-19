from convertTextOutput import getIntegralMatrixAverage
import rootpy
import ROOT
import numpy as np
from numpy import linalg as la
import math
from math import pi
rootpy.log.basic_config_colorized()
from rootpy.io import root_open
from rootpy.plotting import Canvas, Hist, Legend
from ROOT import TH1D
from ROOT import TH2D
import numpy
from vector import vector
from vector import remove_vector_from_onb
from convertTextOutput import get2D
from sys import exit

def getAmplitudes(jpc, M, isob, m3pi, path):
	imported = get2D(path,keywave='none',normalizeToIntegrals=False)
	rightpoints=[]	
#	print imported
	for point in imported[0]:
		if point[4]==jpc and point[0] <= m3pi and point[1] > m3pi and point[13]==M and point[14] == isob:
			rightpoints.append(point)
	rightpoints.sort()
#	inte=[]
#	re=[]
#	im=[]
#	errinte=[]
#	errre=[]
#	errim=[]
#	phi=[]
#	errphi=[]
	ampl = []
	for point in rightpoints:
#		print point
#		inte.append(point[5])
#		re.append(point[7])
#		im.append(point[9])
#		phi.append(point[11])
#		errinte.append(point[6])
#		errre.append(point[8])
#		errim.append(point[10])
#		errphi.append(point[12])
		ampl.append(complex(point[7],point[9]))
	return ampl

def getMassDependentIntegrals(wave1, wave2, mMin, mMax, intDir, binWidth = 0.01):
	"""
	Retruns the mass-dependence of the off diagonal integrals for wave1 and wave2
	"""
	nBins = int((mMax-mMin)/binWidth)
	ret=[]
	res=[]
	ims=[]
	for i in range(nBins):
		integrals = getIntegralMatrixAverage(mMin+i*binWidth,mMin+(i+1)*binWidth,intDir, True, False)
		n1 = integrals[0][wave1]
		n2 = integrals[0][wave2]
		integral = integrals[1][n1][n2]
		ret.append([mMin+(i+.5)*binWidth,abs(integral)])
		res.append([mMin+(i+.5)*binWidth,integral.real])
		ims.append([mMin+(i+.5)*binWidth,integral.imag])
	return ret,res,ims

def get2DoffDiagonal(wave1, wave2, mMin, mMax, intDir):
	"""
	Get the two dimensional off-diagonal terms for two de-isobarred waves
	"""
	integrals = getIntegralMatrixAverage(mMin,mMax,intDir,True, False)
	waves1=[]
	waves2=[]
	for key in integrals[0].iterkeys():
		if key.startswith(wave1):
			mev1 = float(key.split('_')[1])/1000
			mev2 = float(key.split('_')[2][:4])/1000
			waves1.append([mev1,mev2,integrals[0][key]])
		if key.startswith(wave2):
			mev1 = float(key.split('_')[1])/1000
			mev2 = float(key.split('_')[2][:4])/1000
			waves2.append([mev1,mev2,integrals[0][key]])
	ints = []
	res=[]
	ims=[]
	for ww1 in waves1:
		for ww2 in waves2:
			ints.append([ww1[0],ww1[1],ww2[0],ww2[1],abs(integrals[1][ww1[2]][ww2[2]])])
			res.append([ww1[0],ww1[1],ww2[0],ww2[1],integrals[1][ww1[2]][ww2[2]].real])
			ims.append([ww1[0],ww1[1],ww2[0],ww2[1],integrals[1][ww1[2]][ww2[2]].imag])
	return ints,res,ims

def compare_diag_offdiag(wave1,wave2,base1,base2,mMin,mMax,intDir):
	"""
	Compares de-isobarred and constant isobar integrals
	"""
	integrals = getIntegralMatrixAverage(mMin,mMax,intDir,False, False)
	nwaves1 = []
	nwaves2 = []
	binning1=[]
	binning2=[]
	for key in integrals[0].iterkeys():
		print key, wave1, wave2
		if key.startswith(base1):
			nwaves1.append(integrals[0][key])
			mup = float(key.split('_')[2][:4])/1000
			mlo = float(key.split('_')[1])/1000
			if not mup in binning1:
				binning1.append(mup)
			if not mlo in binning1:
				binning1.append(mlo)
#			print key
		if key.startswith(wave1):
			nwave1 = integrals[0][key]
			print key,'<---- nwave1'
		if key.startswith(base2):
			nwaves2.append(integrals[0][key])
			mup = float(key.split('_')[2][:4])/1000
			mlo = float(key.split('_')[1])/1000
			if not mup in binning2:
				binning2.append(mup)
			if not mlo in binning2:
				binning2.append(mlo)

#			print key
		if key.startswith(wave2):
			nwave2 = integrals[0][key]
			print key,'<---- nwave2'
	binning1.sort()
	binning2.sort()
	nwaves1.sort()
	nwaves2.sort()
	diag1 = integrals[1][nwave1][nwave1]
	diag2 = integrals[1][nwave2][nwave2]
	offd  = integrals[1][nwave1][nwave2]
	diag1_c = 0.
	for nn in nwaves1:
		for mm in nwaves1:
			diag1_c+=integrals[1][nn][mm]
	diag2_c = 0.
	for nn in nwaves2:
		for mm in nwaves2:
			diag2_c+=integrals[1][nn][mm]
	offd_c = 0.
	for nn in nwaves1:
		for mm in nwaves2:
			offd_c+=integrals[1][nn][mm]
#	print diag1, diag1_c
#	print diag2, diag2_c
#	print offd, offd_c
#	print offd/(diag1*diag2)**.5
#	print binning1,binning2
	binning1 = numpy.asarray(binning1,dtype=numpy.float64)
	binning2 = numpy.asarray(binning2,dtype=numpy.float64)
	hist= TH2D('2d','2d',len(binning1)-1,binning1, len(binning2)-1, binning2)
	mabiwi1 = 0.
	for i in range(len(binning1)-1):
		biwi = binning1[i+1]-binning1[i]
		if biwi > mabiwi1:
			mabiwi1 = biwi
	mabiwi2 = 0.
	for i in range(len(binning2)-1):
		biwi = binning2[i+1]-binning2[i]
		if biwi > mabiwi2:
			mabiwi2=biwi
	print mabiwi1,mabiwi2
	iinntt=0.
	for i in range(len(nwaves1)):
		n = nwaves1[i]
		bw1 = (binning1[i+1]-binning1[i])/mabiwi1
		for j in range(len(nwaves2)):
			bw2 = (binning2[j+1]-binning2[j])/mabiwi2
			m = nwaves2[j]
			integral_value = abs(integrals[1][n][m]/(diag1*diag2)**.5)
			hist.SetBinContent(i+1,j+2,integral_value/bw1/bw2)
			iinntt += integral_value
	print 'iinntt: ',iinntt
	return hist

def get_bad_vector(base1,base2,mMin,mMax,intDir):
	"""
	Gets the 'bad' vector for the two de-siobarred waves
	"""
	integrals = getIntegralMatrixAverage(mMin,mMax,intDir,False, False)
	nwaves1 = []
	nwaves2 = []
	binning1=[]
	binning2=[]
	for key in integrals[0].iterkeys():
		if key.startswith(base1):
			nwaves1.append(integrals[0][key])
			mup = float(key.split('_')[2][:4])/1000
			mlo = float(key.split('_')[1])/1000
			if not mup in binning1:
				binning1.append(mup)
			if not mlo in binning1:
				binning1.append(mlo)
		if key.startswith(base2):
			nwaves2.append(integrals[0][key])
			mup = float(key.split('_')[2][:4])/1000
			mlo = float(key.split('_')[1])/1000
			if not mup in binning2:
				binning2.append(mup)
			if not mlo in binning2:
				binning2.append(mlo)
	binning1.sort()
	binning2.sort()
	nwaves1.sort()
	nwaves2.sort()
	excl_bin=0
	dim1 = len(nwaves1)
	dim2 = len(nwaves2)-excl_bin
	DIM = dim1 + dim2
	aa = np.zeros((DIM,DIM),dtype = np.complex128)
	for i in range(dim1):
		for j in range(dim1):
			aa[i,j] = integrals[1][nwaves1[i]][nwaves1[j]]
	for i in range(dim1):
		for j in range(dim2):
			aa[i,j+dim1] = integrals[1][nwaves1[i]][nwaves2[excl_bin+j]]
			aa[j+dim1,i] = integrals[1][nwaves2[excl_bin+j]][nwaves1[i]]
	for i in range(dim2):
		for j in range(dim2):
			aa[i+dim1,j+dim1] = integrals[1][nwaves2[excl_bin+i]][nwaves2[excl_bin+j]]
	for i in range(DIM):
		for j in range(DIM):
			if i!=j and not (aa[i,i]*aa[j,j]) == 0.+0.j:
				aa[i,j]/=(aa[i,i]*aa[j,j])**.5
			elif i!= j:
				aa[i,j] = 0.+0.j
	for i in range(DIM):
		if aa[i,i] != 0.+0.j:
			aa[i,i]=1.+0.j
	val,vec = la.eig(aa)
	vec = np.transpose(vec)
	mini = 0
	maxi = 0
	minev = abs(val[0]) 	
	maxev = abs(val[0]) 	
	valls = [abs(val[i]) for i in range(DIM)]
	valls.sort()
	outfile = open("GRAM_EV.txt",'w')
	for i in range(DIM):
		outfile.write(str(i)+' '+str(valls[i])+'\n')
	outfile.close()
	for i in range(DIM):
#			print ":::",i,abs(val[i]),minev,maxev
		if abs(val[i]) < minev and not abs(val[i]) == 0.:
			mini = i
			minev = abs(val[i])
		if abs(val[i]) > maxev:
			maxi = i
			maxev = abs(val[i])

	minval = val[mini].real

#	print '------------------------------------------'
#	print minev,maxev,mini,maxi
#	print '------------------------------------------'
#	print vec[mini]
#	print 
#	print vec[maxi]
#	with open("GRAM_VEC",'w') as outfile:
#		for i in range(DIM):
#			outfile.write(str(i)+'  '+str(vec[mini][i].real)+'  '+str(vec[maxi][i].real)+'\n')
	minvec = vector([vec[mini][i].real for i in range(DIM)])
	maxvec = vector([vec[maxi][i].real for i in range(DIM)])
#	return minvec, maxvec
	with open('minvec','w') as minnn:
		for val in minvec:
			minnn.write(str(val)+'  ')
	with open('diag_f0','w') as diagout:
		for nnn in nwaves1:
			diagout.write(str(integrals[1][nnn][nnn].real)+'  ')

	with open('diag_rho','w') as diagout:
		for nnn in nwaves2:
			diagout.write(str(integrals[1][nnn][nnn].real)+'  ')

	return minval

if __name__ == "__main__":
	canv = Canvas(width=1000,height=1000)
#	with open('min_egenvalues_over_mass_bin.txt','w') as oouutt:
#		for i in range(200):
#			minval = get_bad_vector('1-(1++)0+ f0_','1-(1++)0+ rho_',.50+i*.01,.51+i*.01,'/nfs/mds/user/fkrinner/massIndepententFits/fits/tst/integrals/0.10000-0.14077')	
#			oouutt.write(str(i)+'  '+str(minval)+'\n')

	get_bad_vector('1-(1++)0+ f0_','1-(1++)0+ rho_',1.50,1.54,'/nfs/mds/user/fkrinner/massIndepententFits/fits/tst/integrals/0.10000-0.14077')	

#	minvec.normalize()
	
	hist = compare_diag_offdiag('1-(1++)0+ (pipi)_S pi P','1-(1++)0+ rho pi S','1-(1++)0+ f0_','1-(1++)0+ rho_',1.50,1.54,'/nfs/mds/user/fkrinner/massIndepententFits/fits/tst/integrals/0.10000-0.14077')
	hist.Draw('col')
	exit(1)##################################


#	print minvec
	dat_f0 = getAmplitudes('1++', '0', 'f0_',  1.52,'/nfs/mds/user/fkrinner/massIndepententFits/fits/0pp_1mm_2pp_in1pp_MC/fit/0.14077-0.19435')
	dat_rho= getAmplitudes('1++', '0', 'rho_',1.52,'/nfs/mds/user/fkrinner/massIndepententFits/fits/0pp_1mm_2pp_in1pp_MC/fit/0.14077-0.19435')
	dat = vector(dat_f0 + dat_rho)
	coeff = minvec*dat
	dat_new = dat - minvec*coeff
	print coeff
	with open('outfile.txt','w') as out:
		for i in range(len(dat)):
			out.write(str(i)+"  "+str(abs(dat[i])**2)+"  "+str(abs(dat_new[i])**2)+"  "+str(minvec[i])+'\n')
	n_f0 =len(dat_f0)
	n_rho = len(dat_rho)

	bad_f0 = minvec[:n_f0]
	bad_f0.normalize()
	basis_f0=[]
	for i in range(n_f0):
		basis_f0.append(vector([0.]*n_f0))
		basis_f0[i][i] = 1.
	
	bad_rho = minvec[n_f0:]
	bad_f0.normalize()
	basis_rho =[]
	for i in range(n_rho):
		basis_rho.append(vector([0.]*n_rho))
		basis_rho[i][i] = 1.


	rembas_f0=remove_vector_from_onb(basis_f0,bad_f0)
	rembas_rho=remove_vector_from_onb(basis_rho,bad_rho)
	with open("outfile.txt",'w') as out:
		dim = len(rembas_f0[0])
		for i in range(dim):
			out.write(str(i))
			for bv in rembas_f0:
				out.write("  "+str(bv[i]))
			out.write("\n")

	dim=len(rembas_f0)
	for i in range(dim):
		sp = rembas_f0[i]*bad_f0
		for j in range(dim):
			sp = rembas_f0[i]*rembas_f0[j]

	for i in range(len(rembas_f0)):
		with open('./ONB_f0/'+str(i),'w') as fill:
			rembas_f0[i] = vector([0.]*len(rembas_f0[i]))
			rembas_f0[i][i] = 1.
			for val in rembas_f0[i]:
				fill.write(str(val)+"  ")

	for i in range(len(rembas_rho)):
		with open('./ONB_rho/'+str(i),'w') as fill:
			rembas_rho[i] = vector([0.]*len(rembas_rho[i]))
			rembas_rho[i][i] = 1.
			for val in rembas_rho[i]:
				fill.write(str(val)+"  ")


	print len(rembas_rho)

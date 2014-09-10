import rootpy
import ROOT
import numpy
from numpy import linalg as la
import math
from math import pi
rootpy.log.basic_config_colorized()
from rootpy.io import root_open
from rootpy.plotting import Canvas, Hist, Legend
from ROOT import TH1D
from ROOT import TH2D
import os  
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re

import resource

"""Collection of more or less useful functions to handle the PWA text-output"""

#####################################################################################################################################
##								GLOBAL DEFINITIONS						   ##
#####################################################################################################################################

keyWaves=[
		'1-(1++)0+ rho pi S                                          ',
		'1-(1++)1+ rho pi S                                          ',# Wavename has to match name in the card (Length 60)
		'1-(4++)1+ rho pi G                                          ',
		'1-(1++)0+ f2 pi P                                           ',
		'1-(2++)1+ rho pi D                                          ']

keyWave=keyWaves[0]

integralsDefault='/nfs/mds/user/fkrinner/massIndepententFits/integrals/pipiS/0.10000-0.14077'

binning_f0=[.278,.320,.360, .400, .440, .480, .520, .560, .600, .640, .680, .720, .760, .800, .840, .880, .920, .930, .940, .950, .960, .970, .980, .990, 1.000, 1.010, 1.020, 1.030, 1.040, 1.050, 1.060, 1.070, 1.080, 1.120, 1.160, 1.200, 1.240, 1.280, 1.320, 1.360, 1.400, 1.440, 1.480, 1.520, 1.560, 1.600, 1.640, 1.680, 1.720, 1.760, 1.800, 1.840, 1.880, 1.920, 1.960, 2.000, 2.040, 2.080, 2.120, 2.160, 2.200, 2.240, 2.280]

f0_isobars={'(pipi)_S':0,'f0(980)':1,'f0(1500)':2} # To generate pseudoplots, the order of the isobars has to match the order in the tabulated functions file.

mPi=0.13957018

deisobared_f0_waves={	'0-+':[	'1-(0-+)0+ (pipi)_S pi S                                     ' ,
				'1-(0-+)0+ f0(980) pi S                                      ' ,
				'1-(0-+)0+ f0(1500) pi S                                     '],
			'1++':[	'1-(1++)0+ (pipi)_S pi P                                     ' ,
				'1-(1++)0+ f0(980) pi P                                      '],
			'2-+':[	'1-(2-+)0+ (pipi)_S pi D                                     ' ,
				'1-(2-+)0+ f0(980) pi D                                      ']}
#####################################################################################################################################
##								READ FILES							   ##
#####################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------
def getLike(inFile):
	"""Returns the log-likelihood of 'inFile'"""
	count_calls('getLike')
	nextIsLike=False
	nextIsNevents=False
	data=open(inFile,'r')
	for line in data.readlines():
		if nextIsLike:
			like=float(line)
			if nevents >0:
				data.close()
				return like
			else:
				data.close()
				return -1.E10
		if nextIsNevents:
			nevents = int(line)
			nextIsNevents=False
		if 'likelihood' in line:
			nextIsLike=True
		if 'Number of events' in line:
			nextIsNevents=True
	data.close()
	return -1.E10
#------------------------------------------------------------------------------------------------------------------------------------
def getFitStatus(inFile): 
	"""Retuns [mMin,mMax,nevents,likelihood,t'Min,t'Max]"""
	count_calls('getFitStatus')
	nextIsLike=False
	nextIsNevents=False
	nextIsMass=False
	nextIsT=False
	data=open(inFile,'r')
	nevents=0
	like=0.
	mMin=0.
	mMax=0.
	tMin=0.
	tMax=0.
	for line in data.readlines():
		if nextIsLike:
			like=float(line)
			nextIsLike=False
		if nextIsNevents:
			nevents = int(line)
			nextIsNevents=False
		if nextIsMass:
			mMin=float(line.split(';')[0])
			mMax=float(line.split(';')[1])
			nextIsMass=False
		if nextIsT:
			tMin=float(line.split(';')[0])
			tMax=float(line.split(';')[1])
			nextIsT=False
		if 'likelihood' in line:
			nextIsLike=True
		if 'Number of events' in line:
			nextIsNevents=True
		if 'Mass bin [GeV]:' in line:
			nextIsMass=True
		if "// t' bin" in line:
			nextIsT =True
		if not nevents==0 and not like ==0. and not mMin==0. and not mMax==0. and not tMin ==0. and not tMax==0.:
			break
	return [mMin,mMax,nevents,like,tMin,tMax]			
#------------------------------------------------------------------------------------------------------------------------------------
def getEvents(inFile): 
	"""Returns the number of events of 'inFile'"""
	count_calls('getEvents')
	nextIsN=False
	data=open(inFile,'r')
	for line in data.readlines():
		if nextIsN:
			n=int(line)
			return n
		if 'Number of events' in line:
			nextIsN=True
	return 0
#------------------------------------------------------------------------------------------------------------------------------------
def getBestLike(direct='.'): 
	"""Returns the file with the bes log-likelihood in the directory 'direct'"""
	count_calls('getBestLike')
	maxLike=-1.E10
	bestFile=''
	for fn in os.listdir(direct):
    		if os.path.isfile(direct+'/'+fn):
			likeNew=getLike(direct+'/'+fn)
			print "  Found file with likelihood: "+str(likeNew)
			if maxLike<likeNew:
				maxLike=likeNew
				bestFile=fn
	print "Best file is: '"+bestFile+"'"
	return bestFile
#------------------------------------------------------------------------------------------------------------------------------------
def getTotalFitStatus(direct='./'):
	"""Return the total fit status of 'direct'"""
	count_calls('getTotalFitStatus')
	data=[]
	for fn in os.listdir(direct):
		if 'text_fit_' in fn:
			for fn2 in os.listdir(direct+'/'+fn):
				if os.path.isfile(direct+'/'+fn+'/'+fn2):
					print " -  Read file .../"+fn2
					data.append(getFitStatus(direct+'/'+fn+'/'+fn2))
	return data
#------------------------------------------------------------------------------------------------------------------------------------
def getBestFits(inDirect='./'): 
	"""Returns a list of files with the best log-likelihoods"""
	count_calls('getBestFits')
	direct=inDirect
	while True: # Remove // from path so it can be found in the 'bestFits.txt' even if //, /// or ... is given in the path
		directOld=direct
		direct=direct.replace('//','/')
		if direct == directOld:
			break
	fits=[]
	foundDirect=False
	foundDirectTot=False
	if not os.path.isfile('bestFits.txt'): #Store results, so if the best likelihoods are determined once, they can be found more easyly
		open('bestFits.txt', 'a').close()
	read=open('bestFits.txt','r')
	for line in read.readlines():
#		print line
		chunks=line.split()
		if chunks[0]=='DIRECTORY:' and foundDirect:
			foundDirect=False
		if foundDirect:
			fits.append([chunks[0],float(chunks[1]),float(chunks[2])])
		if chunks[0]=='DIRECTORY:':
			if chunks[1]==direct:
				print "Found directory in 'bestFits.txt': Do not scan the folders, but take the results from the file."
				foundDirect=True	
				foundDirectTot=True	
	if not foundDirectTot:
		for fn in os.listdir(direct):
			if 'text_fit_' in fn:
				chunks=fn.split('_')
				m3PiMin=float(chunks[5])/1000
				m3PiMax=float(chunks[6])/1000	
				bestFile=getBestLike(direct+'/'+fn)
				fits.append([direct+'/'+fn+'/'+bestFile,m3PiMin,m3PiMax])
		store=open('bestFits.txt','a')
		store.write('DIRECTORY: '+direct+'\n')
		for i in range(0,len(fits)):
			store.write(fits[i][0]+'  '+str(fits[i][1])+'  '+str(fits[i][2])+'\n')
		store.close()
	return fits		
#------------------------------------------------------------------------------------------------------------------------------------
def deleteDirectory(inDirect='./'):
	"""Deletes a directory from the 'bestFits.txt' file (e.g., if more attempts are made)"""
	count_calls('deleteDirectory')
	direct=inDirect
	while True:
		directOld=direct
		direct=direct.replace('//','/')
		if direct == directOld:
			break
	if os.path.isfile('bestFits.txt'):
		inFile = open('bestFits.txt','r')
		lines = inFile.readlines()
		inFile.close()
		outFile=open('bestFits.txt','w')
		writeNext=True
		for line in lines:
			if line.startswith('DIRECTORY:'):
				chunks=line.split()
				if chunks[1]==direct:
					writeNext=False
				else:
					writeNext=True
			if writeNext:
				outFile.write(line)
		outFile.close()
	else:
		print "Directory not found in 'bestFits.txt'."
#------------------------------------------------------------------------------------------------------------------------------------
def updateDirectory(direct='./'):
	"""Updates a directory in the 'bestFits.txt' file"""
	count_calls('updateDirectory')
	deleteDirectory(direct)
	getBestFits(direct)
#------------------------------------------------------------------------------------------------------------------------------------
def get_values(string):
	"""Extracts the values from '(...,...)(...,...)(...,...) shaped string, where every ... is a float"""
	count_calls('get_values')
	valst = string.strip()[1:-1] # remove (...) at the end
	chunks = valst.split(')(')
	values = []
	for chunk in chunks:
		values.append([float(bit) for bit in chunk.split(',')])
	return values

#------------------------------------------------------------------------------------------------------------------------------------
def readTextFile(inFile): 
	"""Reads the text-output 'inFile' and returns [{wave: ...},[[COMA]]]. Works for arbitrary rank."""
	count_calls('readTextFile')
	waves={}
	print "  OPEN: "+inFile
	data=open(inFile,'r')
	print 'Reading data file: \n   - '+inFile
	enterTheMatrix=False
	nextisNevents=False
	nextisMasses=False
	nextisTprime=False
	covarianceMatrix=[]
	prevWaves=[]
	prevRank=1
	for line in data.readlines():
		if line.startswith("'"):
			wave=line[1:61]
			values = get_values(line.split("'")[-1])
			re=values[0][0]
			im=values[0][1]
			if prevRank == len(values)-1:# Addd omitted entries for ranks != 1 (0.000,0.000)
				rankstr ='R0'+str(len(values))
				for i in range(1,len(values)):
					prevWave = prevWaves[-i]
					waves[prevWave[:60-len(rankstr)]+rankstr] = [0.,0.]	
				if prevRank == 1: # Remove the rankless entry if necessary
					waves[prevWave[:57]+'R01'] = waves[prevWave]	
					waves.pop(prevWave)
			if len(values)== 1:	
				waves[wave]=[re,im]		
			else:									
				for i in range(len(values)):
					rankstr='R0'+str(i+1)
					waves[wave[:60-len(rankstr)]+rankstr] = [values[i][0],values[i][1]]
			prevRank = len(values)
			prevWaves.append(wave)
		if enterTheMatrix:
			modLine=line.replace('-',' -') # Separate enries, where the '-' takes the whitespace
			modLine=modLine.replace('E -','E-') # undo the separation from the line before, if an exponent was affected
			modLine=modLine.strip()
			covarianceMatrixLine = [float(chunk) for chunk in modLine.split()]
			covarianceMatrix.append(covarianceMatrixLine)
		if nextisTprime:
			nextisTprime=False
			waves['tprime']=[float(ttt) for ttt in line.split(';')]
		if nextisMasses:
			nextisMasses=False
			waves['m3Pi']=[float(mmm) for mmm in line.split(';')]
		if nextisNevents:
			nextisNevents=False
			waves['nevents']=int(line)
		if 'missing rank entries' in line:
			enterTheMatrix=True
		if "t' bin" in line:
			nextisTprime=True
		if 'Mass bin ' in line:
			nextisMasses=True
		if 'Number of events' in line:
			nextisNevents=True
	nWave = 0
	for wave in prevWaves: #Count the waves..
		try:
			waves[wave].append(nWave)
			nWave+=1
		except KeyError:
			rank =1
			while True:
				try:
					rankstr='R0'+str(rank)
					waves[wave[:60-len(rankstr)]+rankstr].append(nWave)
					nWave+=1
					rank+=1
				except KeyError:
					break
	return [waves,covarianceMatrix]	
#------------------------------------------------------------------------------------------------------------------------------------
def getIntegrals(
			inFile	): 
	"""Returns diagonal integrals in 'inFile' (real) as {wave: integral}"""
	count_calls('getIntegrals')
	enterTheMatrix=False
	data=open(inFile,'r')
	integrals=[]
	waves=[]
	intMap={}
	nWave=0
	for line in data.readlines():
		if line.startswith("'"):
			waves.append(line[1:61])
		if enterTheMatrix:
			nWave+=1
			chunks=line.split('(')
			vals=chunks[nWave]
			reIm=vals.replace(')',' ').split(',')
			integral=[float(reIm[0]),float(reIm[1])]
			integrals.append(integral)
		if 'Matrix' in line:
			enterTheMatrix=True
	for i in range(0,len(waves)):
		intMap[waves[i]]=integrals[i]
	return intMap
#------------------------------------------------------------------------------------------------------------------------------------
def getIntegralMatrix(
			inFile	): 
	"""Returns the integrals in 'inFile' as complex numbers in the form [{wave: ID},[[INTEGRALS]]], where ID gives the position in the matrix"""
	count_calls('getIntegralMatrix')
	enterTheMatrix=False
	data=open(inFile,'r')
	integrals=[]
	waves=[]
	intMap={}
	nWave=0
	for line in data.readlines():
		if line.startswith("'"):
			waves.append(line[1:61])
		if enterTheMatrix:
			integralLine=[]
			chunks=line.split('(')
			for i in range(1,len(chunks)):
				vals=chunks[i]
				reIm=vals.replace(')',' ').split(',')
				integral=float(reIm[0])+1j*float(reIm[1])
				integralLine.append(integral)
			integrals.append(integralLine)
		if 'Matrix' in line:
			enterTheMatrix=True
	for i in range(0,len(waves)):
		intMap[waves[i]]=i
	if not isHermitian(integrals): # There was an error in the textoutput. The matrix was not hermitian. Fix this by hand.
		hermitianize(integrals)
	return [intMap,integrals]
#------------------------------------------------------------------------------------------------------------------------------------
def getIntegralAverage(		
			m3PiMin,m3PiMax,
			intDir=integralsDefault,
			acceptanceCorrected=False,
			normalizeToDiag=False		): # == COMPENSATE_AMP 0
	"""Returns the diagonal integrals (real) averaged from 'm3PiMin' to 'm3PiMax'. Formatted as in 'getIntegrals(...)'"""
	count_calls('getIntegralAverage')
	fileString='PWANormIntegralsNAcc'
	if acceptanceCorrected:
		fileString='PWANormIntegralsAcc'
	filesInRange=[]
	print 'Reading integral files: '
	for fn in os.listdir(intDir):
		if fileString in fn:
			chunks=fn.split('_')
			mFileMin=float(chunks[1])/1000
			mFileMax=float(chunks[2])/1000
			if mFileMax-0.005>=m3PiMin and mFileMin<=m3PiMax-0.005: #the 0.005 are there, to avoid rounding errors.
				filesInRange.append(intDir+'/'+fn)
				print '   - '+intDir+'/'+fn
	ints=[]
	for intFile in filesInRange:
		ints.append(getIntegrals(intFile))
	finInts=ints[0]
	for wave in finInts.iterkeys():
		sumRe=0
		sumIm=0
		for i in range(0,len(ints)):
			sumRe+=ints[i][wave][0]
			sumIm+=ints[i][wave][1]
		sumRe/=len(ints)
		sumIm/=len(ints)
		finInts[wave]=[sumRe,sumIm]
		if normalizeToDiag:
			finInts[wave]=[1.,0.]
	return finInts
#------------------------------------------------------------------------------------------------------------------------------------
def getIntegralMatrixAverage(m3PiMin,m3PiMax,intDir=integralsDefault,normalizeToDiag=False,acceptanceCorrected=True): 
	"""Returns the integral-matrix averaged from 'm3PiMin' to 'm3PiMax'. Formatted as in 'getIntegralMatrix(...)'"""
	count_calls('getIntegralMatrixAverage')
	filesInRange=[]
	fileString='PWANormIntegralsNAcc'
	if acceptanceCorrected:
		fileString='PWANormIntegralsAcc'
	print 'Reading integral files: '
	for fn in os.listdir(intDir):
		if fileString in fn:
			chunks=fn.split('_')
			mFileMin=float(chunks[1])/1000
			mFileMax=float(chunks[2])/1000
			if mFileMax-0.005>=m3PiMin and mFileMin<=m3PiMax-0.005: #the 0.005 are there, to avoid rounding errors.
				filesInRange.append(intDir+'/'+fn)
				print '   - '+intDir+'/'+fn
	ints=[]
	for intFile in filesInRange:
#		print "::reading::\n"+intFile
		ints.append(getIntegralMatrix(intFile))
#		print "::done::"
	finInts=ints[0]
	dim = len(finInts[1])
	for i in range(0,dim):
		for j in range(0,dim):
			summ=0.+0.j
			for actInt in ints:
				summ+=actInt[1][i][j]
			summ/=len(ints)
			finInts[1][i][j]=summ
	if normalizeToDiag:
		for i in range(dim):
			for j in range(dim):
				if finInts[1][i][i]!= 0.+0.j and finInts[1][j][j]!= 0.+0.j:
					if not i==j:
						finInts[1][i][j]/= (finInts[1][i][i]*finInts[1][j][j])**(.5)
				else:
					finInts[1][i][j] =  0.+0.j
		for i in range(dim):
			try:
				finInts[1][i][i]/=finInts[1][i][i]
			except ZeroDivisionError:
				finInts[1][i][i]=0.+0.j
	return finInts
#------------------------------------------------------------------------------------------------------------------------------------
def getWholeFit( direct ):
	"""Returns the whole data of a fit"""
	count_calls('getWholeFit')
	fileList=getBestFits(direct)
	fitData=[]
	print fileList
	for i in fileList:
		actData=readTextFile(i[0])
		fitData.append(actData)
	return fitData
#####################################################################################################################################
##								CALCULATE STUFF							   ##
#####################################################################################################################################
def binIntensities(binData):
	"""Returns sorted list of intensities with the respective wave names"""
	count_calls('binIntensities')
	nevents = binData[0]['nevents']
	ints = []
	for key in binData[0].iterkeys():
		if len(key) > 30: # Other keys (t', nevents ... ) are not so long, wave have len(key) == 60.
			intens = (binData[0][key][0]**2. + binData[0][key][1]**2.)*nevents
			ints.append([intens,key.strip()])
	ints.sort()
	total=0.
	for intens in ints:
		total+=intens[0]
	for intens in ints:
		intens.append(intens[0]/total)
	return ints

#------------------------------------------------------------------------------------------------------------------------------------
def getRelevantData(waves,direct):
	"""Gives the T and covariance matrix for all 'waves[i]' in 'direct'"""
	count_calls('getRelevantData')
	points =[]
	wavesStrip=[]
	for wave in waves:
		wavesStrip.append(wave.strip())
	fileList=getBestFits(direct)
	waveNumbers={}
	dataZero = readTextFile(fileList[0][0])
	for wave in dataZero[0].iterkeys():
		if wave.strip() in wavesStrip:
			waveNumbers[wave.strip()] = [wave,dataZero[0][wave][2]]
	for fil in fileList:
		fitData = readTextFile(fil[0])
		mMin = fil[1]
		mMax = fil[2]
		Ts=[]
		coma =[]
		for wave in wavesStrip:
			nevents = fitData[0]['nevents']
			re = fitData[0][waveNumbers[wave][0]][0]*nevents**.5
			im = fitData[0][waveNumbers[wave][0]][1]*nevents**.5
			nn = waveNumbers[wave][1]
			comaLine1=[]
			comaLine2=[]
			for wave2 in wavesStrip:
				mm = waveNumbers[wave2][1]
				comaLine1.append(fitData[1][2*nn][2*mm]*nevents)
				comaLine1.append(fitData[1][2*nn][2*mm+1]*nevents)
				comaLine2.append(fitData[1][2*nn+1][2*mm]*nevents)
				comaLine2.append(fitData[1][2*nn+1][2*mm+1]*nevents)
			Ts.append(re)
			Ts.append(im)
			coma.append(comaLine1)
			coma.append(comaLine2)
		points.append([mMin,mMax,Ts,coma])
	points.sort()
	return points
#------------------------------------------------------------------------------------------------------------------------------------
def getComaAmp(	waves,				# List of waves
		up,				# Upper mass limits for waves
		low, 				# Lower mass limits for waves
		direct,				# Directory
		CONJUGATE=True		):	# Conjugate Amplitudes (Required for some fits, issue in PWA)
	"""Prepares the data to be used by chi2amp.LoadFit(...) in 'chi2ampextended.py'"""
	count_calls('getComaAmp')
	COMPARE=False
	SHOW = False # Shows jacobian and its numerical computed counterpart
	raw_data = getRelevantData(waves,direct)
	if CONJUGATE:
		for p in range(len(raw_data)):
			for i in range(len(raw_data[p][2])):
				if i%2==1:
					raw_data[p][2][i]*=-1
				for j in range(len(raw_data[p][2])):
					if (i+j)%2==1:
						raw_data[p][3][i][j]*=-1
	nWaves = len(waves)
	nBins  = len(raw_data)
	phases = []
	dcdR =[] # dc / dR
	dcdI =[] # dc / dI
	dsdR =[] # ds / dR
	dsdI =[] # ds / dI
	comps=[]
	binning=[raw_data[0][0]]
	for point in raw_data:
		comp=[]
		binning.append(point[1])
		for i in range(len(point[2])/2):
			comp.append(point[2][2*i]+1.j*point[2][2*i+1])
		comps.append(comp)	
	for comp in comps:
		if not abs(point[2][0]==0.):
			phases.append(comp[0]/abs(comp[0])) # Get the phases from the first wave = anchor wave
			re = comp[0].real
			im = comp[0].imag
			ab = (re**2+im**2)**.5
			#			# cos = re/sqrt(re^2 + im^2)
			#			# sin = im/sqrt(re^2 + im^2)
			dcdR.append(1/ab - re**2/ab**3)
			dcdI.append(-re*im/ab**3)
			dsdR.append(-re*im/ab**3)
			dsdI.append(1/ab - im**2/ab**3)
		else:
			phases.append(1.+0.j)
			dcdR.append(0.)
			dcdI.append(0.)
			dsdR.append(0.)
			dsdI.append(0.)
#		print "::::"+str(phases)+"::::"
	for i in range(len(comps)):
		for j in range(len(comps[i])):
			comps[i][j]/=phases[i]
#	for comp in comps:
#		print comp
	comas=[]
	for i in range(len(raw_data)):
		jac = []
		c = phases[i].real
		s = phases[i].imag
		for j in range(int(len(raw_data[i][2])/2)):
			jac_line1 = []		
			jac_line2 = []
			for k in range(int(len(raw_data[i][2])/2)):
				if j == k:
					# c/exp(i phi) = c/(cos(phi) + i sin(phi)) = c*(cos(phi)-i sin(phi))
					jac_line1.append( c) # dR'/dR
					jac_line1.append( s) # dR'/dI
					jac_line2.append(-s) # dI'/dR
					jac_line2.append( c) # dI'/dI
				else:
					jac_line1.append(0.)
					jac_line1.append(0.)
					jac_line2.append(0.)
					jac_line2.append(0.)
			jac.append(jac_line1)
			jac.append(jac_line2)
		for j in range(int(len(raw_data[i][2])/2)):
			re = raw_data[i][2][2*j  ]
			im = raw_data[i][2][2*j+1]
			jac[2*j  ][0]+= (re*dcdR[i] + im*dsdR[i]) # dR'/dRa
			jac[2*j  ][1]+= (re*dcdI[i] + im*dsdI[i]) # dR'/dIa
			jac[2*j+1][0]+= (im*dcdR[i] - re*dsdR[i]) # dI'/dRa
			jac[2*j+1][1]+= (im*dcdI[i] - re*dsdI[i]) # dI'/dIa
		if SHOW:
			numericalJac(raw_data[i][2])
			print
			prettyPrint(jac)
			raw_input()
		coma_rot = (numpy.matrix(jac)*numpy.matrix(raw_data[i][3])*numpy.transpose(numpy.matrix(jac))).tolist()	
#		coma_rot = (numpy.transpose(numpy.matrix(jac))*numpy.matrix(raw_data[i][3])*numpy.matrix(jac)).tolist()	
#		raw_data[i][3]= coma_rot
		comas.append(coma_rot)
	points=[]
	for comp in comps:
		point=[]
		for dat in comp:
			point.append(dat.real)
			point.append(dat.imag)
		points.append(point)
	if COMPARE:
		for i in range(len(raw_data)):
#			raw_T = raw_data[i][2]
#			rot_T = points[i]
#			coma_raw = raw_data[i][3]
#			coma_rot = comas[i]
#			print
#			print rot_T
#			print
#			print coma_rot
			print "--------------"
			print "rotated: "+str(numpy.matrix(rot_T)*numpy.matrix(coma_rot)*numpy.transpose(numpy.matrix(rot_T)))
			print "Raw    : "+str(numpy.matrix(raw_T)*numpy.matrix(coma_raw)*numpy.transpose(numpy.matrix(raw_T)))
			print "T_rot: "+str(rot_T)
			print "T_raw: "+str(raw_T)
	comas_inv=[]
	boi = -1 #BinOfInterest
	for bin in range(len(comas)):
		coma = comas[bin]
		bc = binning[bin]+binning[bin+1]
		bc/=2.
		for i in range(len(coma)):
			upi = up[i/2]
			lowi = low[i/2]
			for j in range(len(coma)):
				upj=up[j/2]
				lowj=low[j/2]
				if bin == boi and i==7 and j==7:
					print '-------'
					print upi
					print lowi
					print upj
					print lowj
					print bc
				if upi < bc or lowi > bc or upj < bc or lowj > bc:
					coma[i][j]=0.
		if bin == boi:
			print (binning[boi]+binning[boi+1])/2.
			print binning
			prettyPrint(coma)
		coma_inv= la.pinv(numpy.matrix(coma)).tolist() # Since Im_anc == 0, the coma is singular. Use pinv instead of inv. Should give the same
		if bin == boi:
			prettyPrint(coma_inv)
			raw_input('<enter>')
		comas_inv.append(coma_inv)
	return [points,comas_inv]
#------------------------------------------------------------------------------------------------------------------------------------
def getComaData(	waves,		# List of waves
			up,		# Upper limits
			low,		# Lower Limirs for the 'waves'
			direct, 	# Directory
			flagg ='PINV', 	# Mode of matrix-inversion (default: PSEUDO-INVERSE)
			eps=1.E-3, 
			CONJUGATE = True		):
	"""Prepares the data to be used by chi2coma.LoadFit() defined in 'chi2comaextended.py'"""
	count_calls('getComaData')
	raw_data = getRelevantData(waves,direct)
	nWaves = len(waves)
	nBins  = len(raw_data)
	final_comas_inv=[]
	data_points=[]
	if CONJUGATE:
		conj=-1
	else:
		conj=1
	for point in raw_data:
		masss = (point[0]+point[1])/2
		wave_on=[]						# Check, which waves are active in the current bin
		for i in range(len(up)):				# Set all other elements to zero
			if masss > low[i] and masss < up[i]:		# So that the wrong correlations won't be taken into account
				wave_on.append(1.)
			else:
				wave_on.append(0.)
		for i in range(len(wave_on)):
			point[2][2*i  ]*=wave_on[i]
			point[2][2*i+1]*=wave_on[i]	
			for j in range(len(wave_on)):
				point[3][2*i  ][2*j  ]*=wave_on[i]*wave_on[j]			
				point[3][2*i+1][2*j  ]*=wave_on[i]*wave_on[j]			
				point[3][2*i  ][2*j+1]*=wave_on[i]*wave_on[j]			
				point[3][2*i+1][2*j+1]*=wave_on[i]*wave_on[j]		
#		prettyPrint(point[3])
		jacobian = []
		data_point=[]
		raw_coma=numpy.matrix(point[3])
		for i in range(nWaves*2):
			jacLine =[0.]*2*nWaves
			ii = int(i/nWaves)
			jj = i - nWaves*ii
			if ii == jj: # Intensities
				addpoint = point[2][2*ii]**2 + point[2][2*ii+1]**2
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  = 2*point[2][2*ii]   # dInt_i / dRe_i = 2*Re_i
					jacLine[2*ii+1]= 2*point[2][2*ii+1] # dInt_i / dIm_i = 2*Im_i
			elif ii > jj: # Real Part
				addpoint = point[2][2*ii]*point[2][2*jj] + point[2][2*ii+1]*point[2][2*jj+1]
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  = point[2][2*jj]	    # dRe_ij / dRe_i = Re_j			
					jacLine[2*jj]  = point[2][2*ii]	    # dRe_ij / dRe_j = Re_i
					jacLine[2*ii+1]= point[2][2*jj+1]   # dRe_ij / dIm_i = Im_j
					jacLine[2*jj+1]= point[2][2*ii+1]   # dRe_ij / dIm_j = Im_i
			else: # Imaginary Part		# POTENTIAL MINUS SIGN HERE ... INVESTIGATE
			#			R_j      *   I_i           -     R_i       *     I_j
				addpoint = point[2][2*jj]*point[2][2*ii+1] - point[2][2*ii]*point[2][2*jj+1]
				addpoint*=conj
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  =-point[2][2*jj+1] *conj  # dIm_ij / dRe_i =-Im_j				
					jacLine[2*jj]  = point[2][2*ii+1] *conj  # dIm_ij / dRe_j = Im_i
					jacLine[2*ii+1]= point[2][2*jj]	  *conj  # dIm_ij / dIm_i = Re_j
					jacLine[2*jj+1]=-point[2][2*ii]   *conj  # dIm_ij / dIm_j =-Re_i
			jacobian.append(jacLine)
		jacobian = numpy.matrix(jacobian)
		final_coma = jacobian*raw_coma*jacobian.T
#		print data_point
		final_coma_inv = invert_right_submatrix(final_coma,flagg, eps)
		final_comas_inv.append(final_coma_inv)
		data_points.append(data_point)
	return [data_points,final_comas_inv]
#------------------------------------------------------------------------------------------------------------------------------------
def getComaDataNon(waves,direct):
	"""Returns the data like 'getComaData(...)', but with noninverted coma"""
	count_calls('getComaDataNon')
	raw_data = getRelevantData(waves,direct)
	nWaves = len(waves)
	nBins  = len(raw_data)
	final_comas_inv=[]
	data_points=[]
	for point in raw_data:
		jacobian = []
		data_point=[]
		raw_coma=numpy.matrix(point[3])
		for i in range(nWaves**2):
			jacLine = [0.]*2*nwaves
			ii = int(i/nWaves)
			jj = i - nWaves*ii
			if ii == jj: # Intensities
				data_point.append(point[2][2*ii]**2 + point[2][2*ii+1]**2)
				jacLine[2*ii]  = 2*point[2][2*ii]   # dInt_i / dRe_i = 2*Re_i
				jacLine[2*ii+1]= 2*point[2][2*ii+1] # dInt_i / dIm_i = 2*Im_i
			elif ii > jj: # Real Part
				data_point.append(point[2][2*ii]*point[2][2*jj] + point[2][2*ii+1]*point[2][2*jj+1])
				jacLine[2*ii]  = point[2][2*jj]	    # dRe_ij / dRe_i = Re_j			
				jacLine[2*jj]  = point[2][2*ii]	    # dRe_ij / dRe_j = Re_i
				jacLine[2*ii+1]= point[2][2*jj+1]   # dRe_ij / dIm_i = Im_j
				jacLine[2*jj+1]= point[2][2*ii+1]   # dRe_ij / dIm_j = Im_i
			else: # Imaginary Part
				data_point.append(-point[2][2*jj]*point[2][2*ii+1] + point[2][2*ii]*point[2][2*jj+1])
				jacLine[2*ii]  = point[2][2*jj+1]   # dIm_ij / dRe_i =-Im_j				
				jacLine[2*jj]  =-point[2][2*ii+1]   # dIm_ij / dRe_j = Im_i
				jacLine[2*ii+1]=-point[2][2*jj]	    # dIm_ij / dIm_i = Re_j
				jacLine[2*jj+1]= point[2][2*ii]     # dIm_ij / dIm_j =-Re_i
			jacobian.append(jacLine)
		jacobian = numpy.matrix(jacobian)
		final_coma = jacobian*raw_coma*jacobian.T
		final_comas_inv.append(final_coma)
		data_points.append(data_point)
	return [data_points,final_comas_inv]
#------------------------------------------------------------------------------------------------------------------------------------
def deisobarredRatio(jpc,direct, intdir='',acceptanceCorrected=True,normalizeToDiag=True):
	"""Gives the ratio of de-isobarred to total wave"""
	count_calls('deisobarredRatio')
	if intdir =='':
		tprime=filter(lambda a: a != '', direct.split('/'))[-1]
		intdir=direct+'/../../integrals/'+tprime+'/'		
	data=getWholeFit( direct )
	if jpc == '0++':
		deiso = 'f0_'
		iso = ['(pipi)', 'f0(']
	elif jpc == '1--':
		deiso = 'rho_'
		iso = ['rho ']
	elif jpc=='2++':
		deiso = 'f2_'
		iso = ['f2 ']
	deisoWaves=[]
	gesWaves=[]
	for wave in data[0][0].iterkeys():
		if len(wave) > 50 and not wave[-3] == 'R':
			if deiso in wave:
				deisoWaves.append(wave)
				gesWaves.append(wave)
			nIn=0
			for isobar in iso:
				if isobar in wave:
					nIn+=1
			if nIn:
				gesWaves.append(wave)
	# one could write out all single intensities here, the data is available, but it is not done by now.
	deisoPlot=getTotal(direct,deisoWaves,intDir=intdir,normalizeToDiag=normalizeToDiag,acceptanceCorrected=acceptanceCorrected)
	gesPlot=getTotal(direct,gesWaves,intDir=intdir,normalizeToDiag=normalizeToDiag,acceptanceCorrected=acceptanceCorrected)
	return {'total vs deisobarred':[deisoPlot,gesPlot]}
#------------------------------------------------------------------------------------------------------------------------------------
def compareFits(direct1, direct2, intdir1='', intdir2='',acceptanceCorrected=True,normalizeToDiag=True):
	"""Compares two fits by waves, if a wave does not appear in both, if compares the corresponding spin-totals"""
	count_calls('compareFits')
	if intdir1 =='':
		tprime=filter(lambda a: a != '', direct1.split('/'))[-1]
		intdir1=direct1+'/../../integrals/'+tprime+'/'
	if intdir2 =='':
		tprime=filter(lambda a: a != '', direct2.split('/'))[-1]
		intdir2=direct2+'/../../integrals/'+tprime+'/'
	data1=getWholeFit( direct1 )
	data2=getWholeFit( direct2 )
	matched={}
	unmatched1=[]
	unmatched2=[]
	for wave in data1[0][0].iterkeys():
		if len(wave) > 50 and not wave[-3] == 'R': # Do not take additional information, such as nevents, tbin m3Pi..., dot not take negative refelctivity
			if data2[0][0].has_key(wave):
				matched[wave.strip()]=[getSingleIntensity(wave,data1),getSingleIntensity(wave,data2)]
			else:
				unmatched1.append(wave)
	for wave in data2[0][0].iterkeys():
		if len(wave) > 50 and not wave[-3]=='R':
			if not matched.has_key(wave.strip()):
				unmatched2.append(wave)
	unmatchedJPC1={}
	for wave in unmatched1:
		jpc = getJPCfromWave(wave)
		if not unmatchedJPC1.has_key(jpc):
			unmatchedJPC1[jpc]=[]
		unmatchedJPC1[jpc].append(wave)
	unmatchedJPC2={}
	for wave in unmatched2:
		jpc = getJPCfromWave(wave)
		if not unmatchedJPC2.has_key(jpc):
			unmatchedJPC2[jpc]=[]
		unmatchedJPC2[jpc].append(wave)
	for jpc in unmatchedJPC1.iterkeys():
		if unmatchedJPC2.has_key(jpc):
			matched[jpc]=[getTotal(direct1,unmatchedJPC1[jpc],intDir=intdir1,normalizeToDiag=normalizeToDiag,acceptanceCorrected=acceptanceCorrected),getTotal(direct2,unmatchedJPC2[jpc],intDir=intdir2,normalizeToDiag=True,acceptanceCorrected=acceptanceCorrected)]
		else:
			print jpc +' is still unmatched.'
	for jpc in unmatchedJPC2.iterkeys():
		if not unmatchedJPC1.has_key(jpc):
			print jpc +' is still unmatched.'
	unnat1=[]
	unnat2=[]
	for wave in data1[0][0].iterkeys():
		if wave[-3]=='R':
			unnat1.append(wave)
	for wave in data2[0][0].iterkeys():
		if wave[-3]=='R':
			unnat2.append(wave)
	matched['unnatural']=[getTotal(direct1,unnat1,intDir=intdir1,normalizeToDiag=normalizeToDiag,acceptanceCorrected=acceptanceCorrected),getTotal(direct2,unnat2,intDir=intdir2,normalizeToDiag=normalizeToDiag,acceptanceCorrected=acceptanceCorrected)]
	return matched
#------------------------------------------------------------------------------------------------------------------------------------
def getSingleIntensity( wave, fitData ):
	"""Returns the intensity of a single wave"""
	count_calls('getSingleIntensity')
	points=[]
	for dataPoint in fitData:
		mMin=dataPoint[0]['m3Pi'][0]
		mMax=dataPoint[0]['m3Pi'][1]
		data=dataPoint[0][wave]
		nEvents=dataPoint[0]['nevents']
		re=data[0]*nEvents**.5
		im=data[1]*nEvents**.5
		intens=re**2+im**2
		nR = 2*data[2]
		nI = nR+1
		coma=[ 	[	dataPoint[1][nR][nR],	dataPoint[1][nR][nI]	],
			[	dataPoint[1][nI][nR],	dataPoint[1][nI][nI]	]	]
		rescale(coma,[nEvents**.5,nEvents**.5])
		jacobian = [2*re,2*im]
		err=vectorMatrixVector(jacobian,coma)**.5
		points.append([mMin,mMax,intens,err])
	points.sort()
	return points
#------------------------------------------------------------------------------------------------------------------------------------
def getFitCoefficient(
			inFile,
			wave1,
			wave2,
			intDir=integralsDefault,
			acceptanceCorrected=False,
			normalizeToDiag=False			): 
	"""Returns the fit coefficients from 'inFile' between the two waves."""
	count_calls('getFitCoefficient')
	fitData=readTextFile(inFile)
	mMin=fitData[0]['m3Pi'][0]
	mMax=fitData[0]['m3Pi'][1]
	integrals=getIntegralAverage(mMin,mMax,intDir,acceptanceCorrected=acceptanceCorrected,normalizeToDiag=normalizeToDiag)
	nEvents=fitData[0]['nevents']
	data1=fitData[0][wave1]
	data2=fitData[0][wave2]
	Int1=integrals[wave1][0]**(.5)
	Int2=integrals[wave2][0]**(.5)
	if Int1 !=0:
		re1=data1[0]*nEvents**(.5)#/Int1 ## Do not normalize the usual plots to the integrals...
		im1=data1[1]*nEvents**(.5)#/Int1 ## Do not normalize the usual plots to the integrals...
	else:
		re1=0
		im1=0
	if Int2!=0:			
		re2=data2[0]*nEvents**(.5)#/Int2 ## Do not normalize the usual plots to the integrals...
		im2=data2[1]*nEvents**(.5)#/Int2 ## Do not normalize the usual plots to the integrals...
	else:
		re2=0
		im2=0
	nR1=2*data1[2]
	nI1=nR1+1
	nR2=2*data2[2]
	nI2=nR2+1
	coma=[	[fitData[1][nR1][nR1],	fitData[1][nR1][nI1],	fitData[1][nR1][nR2],	fitData[1][nR1][nI2]	],
		[fitData[1][nI1][nR1],	fitData[1][nI1][nI1],	fitData[1][nI1][nR2],	fitData[1][nI1][nI2]	],
		[fitData[1][nR2][nR1],	fitData[1][nR2][nI1],	fitData[1][nR2][nR2],	fitData[1][nR2][nI2]	],
		[fitData[1][nI2][nR1],	fitData[1][nI2][nI1],	fitData[1][nI2][nR2],	fitData[1][nI2][nI2]	]]
#	rescale(coma,[nEvents**(.5)/Int1,nEvents**(.5)/Int1,nEvents**(.5)/Int2,nEvents**(.5)/Int2]) ## Do not normalize the usual plots to the integrals...
	rescale(coma,[nEvents**(.5),nEvents**(.5),nEvents**(.5),nEvents**(.5)])
	int1=re1**2+im1**2
	int2=re2**2+im2**2
	re=re1*re2+im1*im2
	im=im1*re2-im2*re1
	phase=math.atan2(im,re)
	int1Jac=[2*re1,2*im1,0,0]
	int2Jac=[0,0,2*re2,2*im2]
	reJac=[re2,im2,re1,im1]
	imJac=[-im2,re2,im1,-re1]
	phaseJac=[0,0,0,0]
	for i in range(0,4):
		try:
#			phaseJac[i]=1./(1.+re**2./im**2.)*(1./re*imJac[i]-im/re**2.*reJac[i]) # Looks like there was something wrong here
			phaseJac[i]=1./(1.+im**2./re**2.)*(1./re*imJac[i]-im/re**2.*reJac[i])
		except:
			phaseJac[i]=0.
	int1err=vectorMatrixVector(int1Jac,coma)**(.5)
	int2err=vectorMatrixVector(int2Jac,coma)**(.5)
	reerr=vectorMatrixVector(reJac,coma)**(.5)
	imerr=vectorMatrixVector(imJac,coma)**(.5)
	phaseerr=vectorMatrixVector(phaseJac,coma)**(.5)
	return [[mMin,mMax,int1,int1err,int2,int2err,re,reerr,im,imerr,phase,phaseerr],{'nevents':nEvents,'tprime':fitData[0]['tprime'],'wave1':wave1,'wave2':wave2}] ## Gives [[FitData (int,re,im,...)],{Fit Info}]
#------------------------------------------------------------------------------------------------------------------------------------
def get2PiSDM(
		inFile,
		intDir,
		waves,
		normalizeToIntegrals= True,
		acceptanceCorrected = False,
		normalizeToDiag=False			 ):
	"""
	Returns the spin-density-matrix for the 2-pion subsystem
	The de-isbarred waves have to be given as follows:
	[{'id':'anything','jpc':'0-+','iso','f0_','M':'0'},{...},...]
	"""
	count_calls('PiSDM')
	if normalizeToIntegrals:
		integralExponent = 1
	else:
		integralExponent = 0
	fitData=readTextFile(inFile)		
	mMin=fitData[0]['m3Pi'][0]
	mMax=fitData[0]['m3Pi'][1]
	integrals=getIntegralAverage(mMin,mMax,intDir,acceptanceCorrected=acceptanceCorrected,normalizeToDiag=normalizeToDiag)
	nEvents=fitData[0]['nevents']
	sdm =[]
	for wave1 in waves:
		binning = []
		for wave in fitData[0].iterkeys():
			if wave[10:14].startswith(wave1['iso']) and wave[3:6]==wave1['jpc'] and wave[7]==wave1['M']:
				if wave1['iso']=='f0_' or wave1['iso']=='f2_':
					m2Min=float(wave[13:17])/1000
					m2Max=float(wave[18:22])/1000
				elif wave1['iso']=='rho_':
					m2Min=float(wave[14:18])/1000
					m2Max=float(wave[19:23])/1000
				if not m2Min in binning:
					binning.append(m2Min)
				if not m2Max in binning:
					binning.append(m2Max)
		binning.sort()
		for wave2 in waves:
			dat=[]
			for bin in range(1,len(binning)):

				m = (binning[bin-1]+binning[bin])/2
				for wave in fitData[0].iterkeys():
					if wave[10:14].startswith(wave1['iso']) and wave[3:6]==wave1['jpc'] and wave[7]==wave1['M']:
						if wave1['iso']=='f0_' or wave1['iso']=='f2_':
							m2Min=float(wave[13:17])/1000
							m2Max=float(wave[18:22])/1000
						elif wave1['iso']=='rho_':
							m2Min=float(wave[14:18])/1000
							m2Max=float(wave[19:23])/1000
						if m2Min <= m and m2Max > m:
							integral1 = integrals[wave][0]**.5
							if integral1 != 0.:
								re1 = fitData[0][wave][0]*nEvents**.5/integral1**integralExponent
								im1 = fitData[0][wave][1]*nEvents**.5/integral1**integralExponent
							else:
									re1=0.
									im1=0.
							nRe1=2*fitData[0][wave][2]
							nIm1=nRe1+1
					if wave[10:14].startswith(wave2['iso']) and wave[3:6]==wave2['jpc'] and wave[7]==wave2['M']:
						if wave2['iso']=='f0_' or wave2['iso']=='f2_':
							m2Min=float(wave[13:17])/1000
							m2Max=float(wave[18:22])/1000
						elif wave2['iso']=='rho_':
							m2Min=float(wave[14:18])/1000
							m2Max=float(wave[19:23])/1000
						if m2Min <= m and m2Max > m:
							integral2 = integrals[wave][0]**.5
							if integral2 != 0.:
								re2 = fitData[0][wave][0]*nEvents**.5/integral2**integralExponent
								im2 = fitData[0][wave][1]*nEvents**.5/integral2**integralExponent
							else:
								re2=0.
								im2=0.
							nRe2=2*fitData[0][wave][2]
							nIm2=nRe2+1
				coma=[	[fitData[1][nRe1][nRe1],fitData[1][nRe1][nIm1],	fitData[1][nRe1][nRe2],	fitData[1][nRe1][nIm2]	],
					[fitData[1][nIm1][nRe1],fitData[1][nIm1][nIm1],	fitData[1][nIm1][nRe2],	fitData[1][nIm1][nIm2]	],
					[fitData[1][nRe2][nRe1],fitData[1][nRe2][nIm1],	fitData[1][nRe2][nRe2],	fitData[1][nRe2][nIm2]	],
					[fitData[1][nIm2][nRe1],fitData[1][nIm2][nIm1],	fitData[1][nIm2][nRe2],	fitData[1][nIm2][nIm2]	]]
				if integral1 !=0. and integral2 !=0.:
					rescale(coma,[nEvents**(.5)/integral1**integralExponent,nEvents**(.5)/integral1**integralExponent,nEvents**(.5)/integral2**integralExponent,nEvents**(.5)/integral2**integralExponent])
				elif integral1 != 0.:
					rescale(coma,[nEvents**(.5)/integral1**integralExponent,nEvents**(.5)/integral1**integralExponent,0.,0.])
				elif integral2 != 0.:
					rescale(coma,[0.,0.,nEvents**(.5)/integral2**integralExponent,nEvents**(.5)/integral2**integralExponent])
				else:
					rescale(coma,[0.,0.,0.,0.])
				Re=re1*re2+im1*im2
				Im=re1*im2-im1*re2
				phase=math.atan2(Im,Re)
				reJac=[re2,im2,re1,im1]
				imJac=[im2,-re2,-im1,re1]
				phaseJac=[0,0,0,0]
				for i in range(0,4):
					try:
						phaseJac[i]=1./(1.+Im**2./Re**2.)*(1./Re*imJac[i]-Im/Re**2.*reJac[i])
					except:
						phaseJac[i]=0.
				errRe=vectorMatrixVector(reJac,coma)**(.5)
				if not wave1 == wave2:
					errIm=vectorMatrixVector(imJac,coma)**(.5)
					errPhase=vectorMatrixVector(phaseJac,coma)**(.5)
				else:
					errIm =0.
					errPhase=0.
				dat.append([binning[bin-1],binning[bin],Re,errRe,Im,errIm,phase,errPhase,wave1['id'],wave2['id']])
			sdm.append(dat)
	return sdm
#------------------------------------------------------------------------------------------------------------------------------------
def getFit2Pi(
		inFile,
		intDir=integralsDefault,
		normalizeToIntegrals=True,
		acceptanceCorrected=True,
		deisobarredWaves=['f0_','rho_','f2_'],
		normalizeToDiag=False,				
		normalizeToBinWidth = True,
		keywave = keyWave		): # If not nomrlized to integrals, the points are too high. NOrmalize the by BindWidth to obtain smooth curves 
	"""Returns the 2-Pi coefficients from 'inFile' as: [[coefficients],{info}]"""
	count_calls('getFit')
	if normalizeToIntegrals: ## Points are divided by int**integralExponent --> nothing happens, if integralExponent == 0.
		integralExponent=1
	else:
		integralExponent=0
	fitData=readTextFile(inFile)
	mMin=fitData[0]['m3Pi'][0]
	mMax=fitData[0]['m3Pi'][1]
	integrals=getIntegralAverage(mMin,mMax,intDir,acceptanceCorrected=acceptanceCorrected,normalizeToDiag=normalizeToDiag)
	nEvents=fitData[0]['nevents']
	useKey = True
	if keywave != 'none':
		keyData=fitData[0][keywave]
		intKey=integrals[keywave][0]**(.5)
		isKey=True
		if not intKey ==0.:
			keyRe=keyData[0]*nEvents**(.5)/intKey**integralExponent
			keyIm=keyData[1]*nEvents**(.5)/intKey**integralExponent
		else:
			keyRe=0.
			keyIm=0.
			isKey=False
		kRe=2*keyData[2]
		kIm=kRe+1
	else:
		useKey = False
		keyData = 1.+0.j
		intKey = 1.
		isKey = True
		keyRe = 1.
		keyIm = 0.
		kRe = -1
		kIm = -1
	Data2Pi=[]
	for deisobarredWave in deisobarredWaves:
		for wave in fitData[0].iterkeys():
			if wave[10:14].startswith(deisobarredWave):
				integral=integrals[wave][0]**(.5)
				if integral != 0:
					re=fitData[0][wave][0]*nEvents**(.5)/integral**integralExponent
					im=fitData[0][wave][1]*nEvents**(.5)/integral**integralExponent
				else:
					re=0.
					im=0.		
				nRe=2*fitData[0][wave][2]
				jpc=wave[3:6]
				M=wave[7]
				nIm=nRe+1
				if useKey:
					coma=[	[fitData[1][nRe][nRe],	fitData[1][nRe][nIm],	fitData[1][nRe][kRe],	fitData[1][nRe][kIm]	],
						[fitData[1][nIm][nRe],	fitData[1][nIm][nIm],	fitData[1][nIm][kRe],	fitData[1][nIm][kIm]	],
						[fitData[1][kRe][nRe],	fitData[1][kRe][nIm],	fitData[1][kRe][kRe],	fitData[1][kRe][kIm]	],
						[fitData[1][kIm][nRe],	fitData[1][kIm][nIm],	fitData[1][kIm][kRe],	fitData[1][kIm][kIm]	]]
				else:
					coma=[	[fitData[1][nRe][nRe],	fitData[1][nRe][nIm],	0. 		    ,	0.			],
						[fitData[1][nIm][nRe],	fitData[1][nIm][nIm],	0. 		    ,	0.			],
						[0. 		     ,	0. 		    ,	0. 		    ,	0.			],
						[0. 		     ,	0. 		    ,	0. 		    ,	0.			]]
				if integral!=0:
					rescale(coma,[nEvents**(.5)/integral**integralExponent,nEvents**(.5)/integral**integralExponent,nEvents**(.5)/intKey**integralExponent,nEvents**(.5)/intKey**integralExponent])
				elif isKey:
					rescale(coma,[0,0,nEvents**(.5)/intKey**integralExponent,nEvents**(.5)/intKey**integralExponent])
				else:
					rescale(coma,[0,0,0,0])
				if deisobarredWave=='f0_' or deisobarredWave=='f2_':
					m2Min=float(wave[13:17])/1000
					m2Max=float(wave[18:22])/1000
				elif deisobarredWave=='rho_':
					m2Min=float(wave[14:18])/1000
					m2Max=float(wave[19:23])/1000
				intens=re**2+im**2
				Re=re*keyRe+im*keyIm
				Im=re*keyIm-im*keyRe
				phase=math.atan2(re*keyIm-im*keyRe,re*keyRe+im*keyIm)
				if re*keyRe+im*keyIm != 0:
					phaseJacFac=(1+((re*keyIm-im*keyRe)/(re*keyRe+im*keyIm))**2)**(-1)
				else:
					phaseJacFac=0
				intensJac=[2*re,2*im,0.,0.]
				reJac=[keyRe,keyIm,re,im]
				imJac=[keyIm,-keyRe,-im,re]
				if re*keyRe+im*keyIm != 0: ## Formula is hardcoded here... in 'getFitCoefficient(...)', a better definition is used. Nevertheless, both should be mathematically equivalent.
					phaseJac=[	phaseJacFac*( keyIm/(re*keyRe+im*keyIm)-keyRe*(re*keyIm-im*keyRe)/(re*keyRe+im*keyIm)**2),
							phaseJacFac*(-keyRe/(re*keyRe+im*keyIm)-keyIm*(re*keyIm-im*keyRe)/(re*keyRe+im*keyIm)**2),
							phaseJacFac*(   -im/(re*keyRe+im*keyIm)-   re*(re*keyIm-im*keyRe)/(re*keyRe+im*keyIm)**2),
							phaseJacFac*(    re/(re*keyRe+im*keyIm)-   im*(re*keyIm-im*keyRe)/(re*keyRe+im*keyIm)**2)]
				else:
					phaseJac=[0,0,0,0]
				errIntens=vectorMatrixVector(intensJac,coma)**(.5)
				try:
					errRe=vectorMatrixVector(reJac,coma)**(.5)
				except:
					if abs(vectorMatrixVector(reJac,coma))<1.E-6:
						print vectorMatrixVector(reJac,coma)
						errRe = (-vectorMatrixVector(reJac,coma))**(.5)
					else:
						raise Exception('Negative error of the real part')
				try:
					errIm=vectorMatrixVector(imJac,coma)**(.5)
				except:
					if abs(vectorMatrixVector(imJac,coma))<1.E-6:
						print vectorMatrixVector(imJac,coma)
						errIm=(-vectorMatrixVector(imJac,coma))**(.5)
					else:
						raise Exception('Negative error of the imag part')
#					raw_input()
				try:
					errPhase=vectorMatrixVector(phaseJac,coma)**(.5)
				except:
					if abs(vectorMatrixVector(phaseJac,coma)) < 1.E-6:
						print vectorMatrixVector(phaseJac,coma)
						errPhase=(-vectorMatrixVector(phaseJac,coma))**(.5)
					else:
						raise Exception('Negative error of the phase')
#					raw_input()
				if not normalizeToIntegrals and normalizeToBinWidth:
					width = m2Max - m2Min
					intens/=width
					errIntens/=width
					Re/=width**.5
					errRe/=width**.5
					Im/=width**.5
					errIm/=width**.5
				Data2Pi.append([m2Min,m2Max,jpc,intens,errIntens,Re,errRe,Im,errIm,phase,errPhase,M,deisobarredWave])
	Data2Pi.sort()
	return [Data2Pi,{'nevents':nEvents,'m3Pi':fitData[0]['m3Pi'],'tprime':fitData[0]['tprime']}]
#------------------------------------------------------------------------------------------------------------------------------------
def getFits(
		wave1,
		wave2,
		direct,
		intDir='',
		acceptanceCorrected=False,
		normalizeToDiag=False		): 
	"""Retuns the 3-Pi dependence of the fit in 'direct'"""
	count_calls('getFits')
	if intDir =='': 					# If no intDir is given, it look for the integrals at the usual place
		tprime=filter(lambda a: a != '', direct.split('/'))[-1]
		intDir=direct+'/../../integrals/'+tprime+'/'
	fileList=getBestFits(direct)
	fitData=[]
	for i in range(0,len(fileList)):
		actData=getFitCoefficient(fileList[i][0],wave1,wave2,intDir,acceptanceCorrected=acceptanceCorrected,normalizeToDiag=normalizeToDiag)
		fitData.append(actData)
	fitData.sort()
	return fitData
#------------------------------------------------------------------------------------------------------------------------------------
def get2D(	direct,
		intDir='',
		normalizeToIntegrals=True,
		acceptanceCorrected=False,
		normalizeToDiag=False,
		keywave = keyWave		): 
	"""Returns a 2-dimensional spectrum (m3Pi-m2Pi)"""
	count_calls('get2D')
	if intDir =='':
		tprime=filter(lambda a: a != '', direct.split('/'))[-1]
		intDir=direct+'/../../integrals/'+tprime+'/'
	fileList=getBestFits(direct)
	fitData=[]
	for i in range(len(fileList)):
		actData=getFit2Pi(fileList[i][0],intDir=intDir,normalizeToIntegrals=normalizeToIntegrals,acceptanceCorrected=acceptanceCorrected,normalizeToDiag=normalizeToDiag,keywave=keywave)
		for i in range(len(actData[0])):
			actData[0][i].reverse()
			actData[0][i].append(actData[1]['m3Pi'][1])
			actData[0][i].append(actData[1]['m3Pi'][0])
#			actData[0][i].append(actData[1]['jpc'])
			actData[0][i].reverse()
			fitData.append(actData[0][i])
		plotInfo=actData[1]
	del plotInfo['m3Pi']
#	del plotInfo['jpc']
	fitData.sort()
	return [fitData,plotInfo]
#------------------------------------------------------------------------------------------------------------------------------------
def getSDM2D(		
			direct,
			intDir='',
			waves=[{'id':'f0','jpc':'0-+','M':'0','iso':'f0_'},{'id':'rho','jpc':'0-+','M':'0','iso':'rho_'},{'id':'f2','jpc':'0-+','M':'0','iso':'f2_'}],
			normalizeToIntegrals= True,
			acceptanceCorrected = False,
			normalizeToDiag=False		 ):
	"""Returns the 2 dimensional 2pi spin-density-matrices"""
	count_calls('getSDM2D')
	if intDir =='':
		tprime=filter(lambda a: a != '', direct.split('/'))[-1]
		intDir=direct+'/../../integrals/'+tprime+'/'
	fileList=getBestFits(direct)
	fitData=[]	
	for i in range(len(fileList)):
		actData = get2PiSDM( fileList[i][0], intDir, waves, normalizeToIntegrals= normalizeToIntegrals,	acceptanceCorrected = acceptanceCorrected,normalizeToDiag=normalizeToDiag)
		fitData.append([fileList[i][1],fileList[i][2],actData])
	return fitData
#------------------------------------------------------------------------------------------------------------------------------------
def makeArgands2Pi(data, slices, jpc):
	"""
	Creates 2pi argand diagrams from 'data' at 'slices'
	'data' is the output of 'get2D(...)'
	"""
	count_calls('makeArgands2Pi')
	argands_raw = {}
	for slic in slices:
		argands_raw[slic] = []
		for dat in data[0]:	
			if dat[0] <= slic and dat[1] > slic and dat[4]== jpc:
				argands_raw[slic].append(dat)
	argands_fine={}
	for mass in argands_raw.iterkeys():
		argands_raw[mass].sort()
		argands_fine[mass]=[[],[],[],[]]
		for point in argands_raw[mass]:
			argands_fine[mass][0].append(point[ 7])
			argands_fine[mass][1].append(point[ 8])
			argands_fine[mass][2].append(-point[ 9])
			argands_fine[mass][3].append(point[10])
	return argands_fine		
#------------------------------------------------------------------------------------------------------------------------------------
def getRelevantMatrices(inFile, 
			isobar, 
			intDir=integralsDefault,
			normalizeToDiag=True,
			acceptanceCorrected=False	):
	"""Returns the relevant matrices to calculate the spin totals for all waves in 'isobar'"""
	count_calls('getRelevantMatrices')
	fitData=readTextFile(inFile)
	m3PiMin=fitData[0]['m3Pi'][0]
	m3PiMax=fitData[0]['m3Pi'][1]
	integrals=getIntegralMatrixAverage(m3PiMin,m3PiMax,intDir,normalizeToDiag,acceptanceCorrected=acceptanceCorrected)
	takeData=[]
	for wave in fitData[0].iterkeys():
		if wave in isobar and not wave[-3]=='R': # Do not use rank 2 negative refelctivity waves (hardcoded)
#			print "Use wave: "+wave.strip()
			data=fitData[0][wave]
			takeData.append([wave.strip(),data[0]+1j*data[1],data[2],integrals[0][wave],wave[-1]])
		elif wave in isobar and wave[-3]=='R':
			waveRep=wave.replace('R01','   ').replace('R02','   ')
			data=fitData[0][wave]
			takeData.append([waveRep.strip()+'R0'+wave[-1],data[0]+1j*data[1],data[2],integrals[0][waveRep],wave[-1]])
#		else:
#			print "Reject wave: "+wave.strip()
	T=[]
	Iij=[]
	coma=[]
	for i in range(0,len(takeData)):
		T.append(takeData[i][1])
		iInt=takeData[i][3]
		iComa=takeData[i][2]
		intLine=[]
		comaLine1=[]
		comaLine2=[]
		for j in range(0,len(takeData)):
			jInt=takeData[j][3]
			if takeData[i][4] == takeData[j][4]:
				intLine.append(integrals[1][iInt][jInt])
#				intLine.append('(1 '+takeData[i][4]+takeData[j][4]+')')
			else:
				intLine.append(0.)
#				intLine.append('(0 '+takeData[i][4]+takeData[j][4]+')')
			jComa=takeData[j][2]
			comaLine1.append(fitData[1][ 2*iComa ][ 2*jComa ])
			comaLine1.append(fitData[1][ 2*iComa ][2*jComa+1])
			comaLine2.append(fitData[1][2*iComa+1][ 2*jComa ])
			comaLine2.append(fitData[1][2*iComa+1][2*jComa+1])
		Iij.append(intLine)
		coma.append(comaLine1)
		coma.append(comaLine2)
	return [T,Iij,coma,{'m3Pi':fitData[0]['m3Pi'],'tprime':fitData[0]['tprime'],'nevents':fitData[0]['nevents']}]
#------------------------------------------------------------------------------------------------------------------------------------
def getTotalPoint(
			inFile,	
			isobar, 
			intDir=integralsDefault,
			normalizeToDiag=True,
			acceptanceCorrected=False,
			lookAtInput=False,
			interference_only = False	):
	"""Calculates one point of the spin totals"""
	count_calls('getTotalPoint')
	data=getRelevantMatrices(inFile,isobar, intDir,normalizeToDiag=normalizeToDiag,acceptanceCorrected=acceptanceCorrected)
	if lookAtInput:
		print "[T,I,coma,info]:"
		print
		print data #### PRINT JUST FOR TEST CASE
		print
		raw_input("Press <enter> to continue")
	value=0.+0.j
	T=data[0]
	Iij=data[1]
	nevents=data[3]['nevents']
	if interference_only:
		for i in range(len(Iij)):
			Iij[i][i] = 0.+0.j
	for i in range(0,len(T)):
		for j in range(0,len(T)):
#			value+=T[j]*Iij[i][j]*T[i].conjugate() ### EXCHANGE i <--> j, to check, if sth. changes.
			value+=T[i]*Iij[i][j]*T[j].conjugate() ### OLD VERSION (SEEMS TO BE RIGHT)
	if abs(value.imag)>10.E-10:
		print "Error: Number of events is not real."
	jacobian=[]
#	if value.real == 0.:
#		for i in range(0,len(T)):
#			for j in range(0,len(T)):
#				print "value+= "+str(T[i])+' * '+str(Iij[i][j])+' * '+str(T[j].conjugate())
	for i in range(0,len(T)):
		currJac=0.+0.j
		for j in range(0,len(T)):
			currJac+=Iij[i][j]*T[j].conjugate()
		jacobian.append(2*currJac.real)
		jacobian.append(-2*currJac.imag)
	err=0.
	coma=data[2]
	for i in range(0,2*len(T)):
		for j in range(0,2*len(T)):
			err+= jacobian[i]*coma[i][j]*jacobian[j]
#			if jacobian[i]*coma[i][j]*jacobian[j] > value.real**2/2:
#				print jacobian[i]
#				print jacobian[j]
#				print coma[i][j]
#				raw_input() 
	err**=.5
	err*=nevents
	nRes=nevents*value.real
#	print '\n'+str(nRes)+' = '+str(nevents)+' * '+str(value.real)+'\n'
	m3PiMin=data[3]['m3Pi'][0]
	m3PiMax=data[3]['m3Pi'][1]
	return [m3PiMin,m3PiMax,nRes,err,data[3]['tprime'][0],data[3]['tprime'][1]]	
#------------------------------------------------------------------------------------------------------------------------------------
def getTotal(
		direct,
		isobar,
		intDir='',
		normalizeToDiag=True,
		acceptanceCorrected=False,
		interference_only = False	):
	"""Gets spin totals over the whole mass range"""
	if intDir =='':
		tprime=filter(lambda a: a != '', direct.split('/'))[-1]
		intDir=direct+'/../../integrals/'+tprime+'/'
	count_calls('getTotal')
	files=getBestFits(direct)
	points=[]
	for fileIn in files:
		points.append(getTotalPoint(fileIn[0],isobar, intDir,normalizeToDiag=normalizeToDiag,acceptanceCorrected=acceptanceCorrected,interference_only=interference_only))
	points.sort()
	return points
#------------------------------------------------------------------------------------------------------------------------------------
def getPercentages(
		direct, 
		intDir='',	# Result will be normalized to the integrals, if son directory is set.
		weightNevents=True		):
	"""Gets the percentages of the waves intensities"""
	count_calls('getPercentages')
	if intDir =='':
		tprime=filter(lambda a: a != '', direct.split('/'))[-1]
		intDir=direct+'/../../integrals/'+tprime+'/'
	fits=getBestFits(direct)
	sizeMap={}
	nEvents=1
	for fit in fits:
		actFit=readTextFile(fit[0])
		if not intDir=='':
			actInt=getIntegralAverage(actFit[0]['m3Pi'][0],actFit[0]['m3Pi'][1],intDir)
		intt=1
		if weightNevents:
			nEvents=actFit[0]['nevents']	
		for wave in actFit[0].iterkeys():
			if len(wave)>50: #Do not use the other entries, that are no waves
				try:
					sizeMap[wave.strip()]
				except:
					sizeMap[wave.strip()]=0.
				if not intDir=='':
					intt=actInt[wave.replace('R02','   ').replace('R01','   ')][0]
				if not intt ==0.:	
					sizeMap[wave.strip()]+=	(actFit[0][wave][0]**2+actFit[0][wave][1]**2)*nEvents#/intt
	#				print '('+str(actFit[0][wave][0])+'^2 + '+str(actFit[0][wave][1])+'^2) * '+str(nEvents)+' / '+str(intt)+' :: '+wave.replace('R02','   ').replace('R01','   ').strip()
	total=0.
	for wave in sizeMap.iterkeys():
		total+=sizeMap[wave]
	perces=[]
	for wave in sizeMap.iterkeys():
		perces.append([sizeMap[wave]/total,wave])
	perces.sort()
	return perces
#------------------------------------------------------------------------------------------------------------------------------------
def getNevents(direct):
	"""Return the number of events over m3Pi"""
	count_calls('getNevents')
	files=getBestFits(direct)
	events=[]
	for fileIn in files:
		events.append([fileIn[1],fileIn[2],getEvents(fileIn[0]),0.])
	events.sort()
	return events
#####################################################################################################################################
##								PRINT ROUTINES							   ##
#####################################################################################################################################
def writeSDM2DtoRoot(inData, fileName):
	"""Writes the two-dimensional 2pi spin-density-matrix to a ROOT file """
	count_calls('writeSDM2DtoRoot')
	outROOT = root_open('./ROOT/'+fileName, mode = "RECREATE")
	binning3=[]
	for dataset in inData:
		if not dataset[0] in binning3:
			binning3.append(dataset[0])
		if not dataset[1] in binning3:
			binning3.append(dataset[1])
	binning3.sort()
	binning3 = numpy.asarray(binning3,dtype=numpy.float64)		
	meta = inData[0][2]
	for i in range(len(meta)):
		plot = meta[i]
		binning2=[]
		id2 = plot[0][-1]
		id1 = plot[0][-2]
		for point in plot:
			if not point[0] in binning2:
				binning2.append(point[0])
			if not point[1] in binning2:
				binning2.append(point[1])		
		binning2.sort()
		binning2 = numpy.asarray(binning2,dtype=numpy.float64)
		hi = TH2D(id1+" vs. "+id2,id1+" vs. "+id2,len(binning3)-1,binning3, len(binning2)-1,binning2)
		for dataset in inData:
			m3 = (dataset[0]+dataset[1])/2
			bin3 = hi.GetXaxis().FindBin(m3)
			plot = dataset[2][i]
			for point in plot:
				m2 = (point[0]+point[1])/2
				bin2 = hi.GetYaxis().FindBin(m2)
				if id1 == id2:
					v = point[2]
					e = point[3]
				else:
					v = point[6]
					e = point[7]
				hi.SetBinContent(bin3,bin2,v)
				hi.SetBinError(bin3,bin2,e)
		if not id1 == id2:
			hi = removePhaseAmbiguities(hi)
		hi.Write()
	outROOT.Close()
#------------------------------------------------------------------------------------------------------------------------------------
def plotArgands(data):
	"""Plots the argand-diagrams"""
	count_calls('plotArgands')
	for mass in data.iterkeys():
		plt.errorbar(data[mass][0],data[mass][2],xerr=data[mass][1],yerr=data[mass][3])
		plt.show()
#------------------------------------------------------------------------------------------------------------------------------------
def plotComparison(dataSet,key='',label1='',label2='',stat='show'):
	"""Plots the comparison of two fits"""
	count_calls('plotComparison')
	if not len(dataSet[0]) == len(dataSet[1]):
		print 'Number of bins does not match'
	for i in range(1,len(dataSet[0])):
		if not dataSet[0][i][0] == dataSet[0][i-1][1]:
			print 'Binning does not match1'
		if not dataSet[1][i][0] == dataSet[1][i-1][1]:
			print 'Binning does not match2'
	for i in range(len(dataSet[0])):
		if not dataSet[0][i][0] == dataSet[1][i][0]:
			print 'Binning does not match3'
	x0=[]
	x1=[]
	ex0=[]
	ex1=[]
	y0=[]
	y1=[]
	ey0=[]
	ey1=[]
	for data in dataSet[0]:
		x0.append((data[1]+data[0])/2)
		ex0.append((data[1]-data[0])/2)
		y0.append(data[2])
		ey0.append(data[3])
	plt.errorbar(x0,y0,xerr=ex0,yerr=ey0,linestyle='None',label=label1)
	for data in dataSet[1]:
		x1.append((data[1]+data[0])/2)
		ex1.append((data[1]-data[0])/2)
		y1.append(data[2])
		ey1.append(data[3])
	plt.title(key)
	plt.errorbar(x1,y1,xerr=ex1,yerr=ey1,linestyle='None',label=label2)
	plt.legend()
	if stat=='show':
		plt.show()
#------------------------------------------------------------------------------------------------------------------------------------
def doComparison(data,name='comparison.pdf',label1='',label2=''):
	"""Writes the comparison of two fits to a .pdf file"""
	count_calls('doComparison')
	pdf_pages = PdfPages(name)
	keys=[]
	for key in data.iterkeys():
		keys.append(key)
	keys.sort()
	for key in keys:
		plotComparison(data[key],key,label1,label2,'')
#		if not folder=='':
#			if not os.path.isdir(folder):
#				os.makedirs(folder)
		pdf_pages.savefig()
		plt.clf()
	pdf_pages.close()
#------------------------------------------------------------------------------------------------------------------------------------
def print2DtoFiles(
			dataSet,
			outFolder,
			jpc='0-+',
			isobar='rho_',
			M='0'			):
	"""
	Creates a folder with files, that contain the whole two dimensional fit
	The outout format can be used by 'Chi2.LoadDataFile(...)'
	"""
	count_calls('print2DtoFiles')
	if not os.path.isdir("./"+outFolder):
		os.makedirs("./"+outFolder)
	dat = dataSet[0]
	ms  = {}
	for point in dat:
		mmed = (point[0] + point[1])/2
		if not ms.has_key(mmed):
			str1 = str(point[0])
			str2 = str(point[1])
			while len(str1) < 4:
				str1 = str1 + '0'
			while len(str2) < 4:
				str2 = str2 + '0'
			ms[mmed]=(str1+"_"+str2).replace('.','')
	for m in ms.iterkeys():
		print2PiToFile(dataSet,m,'./'+outFolder+'/'+ms[m]+'.dat',jpc,isobar,M)

#------------------------------------------------------------------------------------------------------------------------------------
def print2PiToFile(
			dataSet,
			m3Pi,
			outfile='',
			jpc='0-+',			
			isobar='rho_',
			M='0'			):
	"""Prints data with 'jpc' and isobar to a file to be used by mdep-fit"""
	count_calls('print2PiToFile')
	if outfile =='':
		outfile = '2pi_'+isobar+jpc.replace('-','m').replace('+','p')+'.txt'
	fileOut=open(outfile,'w')
	dat = dataSet[0]	
	accepted = []
	for point in dat:
		if point[0] <= m3Pi and point[1]>m3Pi  and point[4] == jpc and point[14]==isobar and point[13] == M:
			accepted.append(point)	
	accepted.sort()
	fileOut.write('0/0/b')
	fileOut.write('/'+str(accepted[0][2]))
	for point in accepted:
		fileOut.write('/'+str(point[3]))
	fileOut.write('\n')
	fileOut.write('0/0/v')
	for point in accepted:
		reim  = (1.+0.j)*point[7] + (0.+1.j)*point[9]	
		inte  = point[5]
		if inte==0.:
			inte = 1.
		int0 = abs(reim)**2/inte
		fileOut.write('/'+str(int0))
	fileOut.write('\n')
	fileOut.write('0/0/e')
	for point in accepted:
		fileOut.write('/1.') # Constant error, but does not matter here, because this will be matched exactly
	fileOut.write('\n')
	fileOut.write('1/1/v')
	for point in accepted:
		fileOut.write('/'+str(point[5]))
	fileOut.write('\n')
	fileOut.write('1/1/e')
	for point in accepted:
		fileOut.write('/'+str(point[6]))
	fileOut.write('\n')
	fileOut.write('0/1/v')
	for point in accepted:
		fileOut.write('/'+str(point[7]))
	fileOut.write('\n')
	fileOut.write('0/1/e')
	for point in accepted:
		fileOut.write('/'+str(point[8]))
	fileOut.write('\n')
	fileOut.write('1/0/v')
	for point in accepted:
		fileOut.write('/'+str(point[9]))
	fileOut.write('\n')
	fileOut.write('1/0/e')
	for point in accepted:
		fileOut.write('/'+str(point[10]))
	fileOut.write('\n')
	fileOut.close()
#------------------------------------------------------------------------------------------------------------------------------------
def print2PiToRoot(dataSet,outFile='fit2PiOut.root'):
	"""Prints 2Pi data to a .root file"""
	count_calls('print2PiToRoot')
	listOfJpc=['0-+','1++','2-+']
	listOfM=['0','1','2']
	mass = str((dataSet[1]['m3Pi'][0]+dataSet[1]['m3Pi'][1])/2)
	outROOT = root_open(outFile, mode = "RECREATE")
	for jpc in listOfJpc:	
		for M in listOfM:
			data=[]
			for i in range(0,len(dataSet[0])):	
				if dataSet[0][i][2]==jpc and dataSet[0][i][11]==M:
					data.append([dataSet[0][i][0],dataSet[0][i][1],dataSet[0][i][3],dataSet[0][i][4],dataSet[0][i][5],dataSet[0][i][6],dataSet[0][i][7],dataSet[0][i][8],dataSet[0][i][9],dataSet[0][i][10]])
			data.sort()
			binning=[data[0][0]]
			for i in range(0,len(data)):
				binning.append(data[i][1])
				if binning[-2] != data[i][0]:
					print "Error in binning, lower and upper edges of bins do not match."
			binning = numpy.asarray(binning,dtype=numpy.float64)	
			hiInt=TH1D('Intensity_of_'+jpc+'_'+M+'_at_m(3#pi)='+mass,'Intensity_of_'+jpc+'_'+M+'_at_m(3#pi)='+mass,len(data),binning)
			hiInt.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}} System (GeV/#it{c}^{2})")
			hiInt.GetYaxis().SetTitle("BW-Intensity")
			hiRe=TH1D('Real_part_of_'+jpc+'_'+M+'_at_m(3#pi)='+mass,'Real_part_of_'+jpc+'_'+M+'_at_m(3#pi)='+mass,len(data),binning)
			hiRe.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}} System (GeV/#it{c}^{2})")
			hiRe.GetYaxis().SetTitle("Real part")
			hiIm=TH1D('Imaginary_part_of_'+jpc+'_'+M+'_at_m(3#pi)='+mass,'Imaginary_part_of_'+jpc+'_'+M+'_at_m(3#pi)='+mass,len(data),binning)
			hiIm.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}} System (GeV/#it{c}^{2})")
			hiIm.GetYaxis().SetTitle("Imaginary part")
			hiPh=TH1D('Phase_of_'+jpc+'_'+M+'_at_m(3#pi)='+mass,'Phase_of_'+jpc+'_'+M+'_at_m(3#pi)='+mass,len(data),binning)
			hiPh.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}} System (GeV/#it{c}^{2})")
			hiPh.GetYaxis().SetTitle("Phase")
			for i in range(0,len(data)):
				hiInt.SetBinContent(i+1,data[i][2])
				hiInt.SetBinError(i+1,data[i][3])
				hiRe.SetBinContent(i+1,data[i][4])
				hiRe.SetBinError(i+1,data[i][5])
				hiIm.SetBinContent(i+1,data[i][6])
				hiIm.SetBinError(i+1,data[i][7])
				hiPh.SetBinContent(i+1,data[i][8])
				hiPh.SetBinError(i+1,data[i][9])
			hiInt.Write()
			hiRe.Write()
			hiIm.Write()
			hiPh.Write()
	outROOT.Close()
#------------------------------------------------------------------------------------------------------------------------------------
def print3PiToRoot(dataSet,outFile='fit3PiOut.root'):
	"""Prints three Pi data to a .root file"""
	count_calls('print3PiToRoot')
	outROOT = root_open(outFile, mode = "RECREATE")	
	binning=[dataSet[0][0][0]]
	for i in range(0,len(dataSet)):
		binning.append(dataSet[i][0][1])
		if binning[-2] != dataSet[i][0][0]:
			print "Error in binning, lower and upper edges of bins do not match."
	binning = numpy.asarray(binning,dtype=numpy.float64)	
	hiInt1=TH1D("Intensity_of_"+dataSet[0][1]['wave1'].strip(),"Intensity_of_"+dataSet[0][1]['wave1'].strip(),len(dataSet),binning)
	hiInt1.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	hiInt1.GetYaxis().SetTitle("Intensity")
	hiInt2=TH1D("Intensity_of_"+dataSet[0][1]['wave2'].strip(),"Intensity_of_"+dataSet[0][1]['wave2'].strip(),len(dataSet),binning)
	hiInt2.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	hiInt2.GetYaxis().SetTitle("Intensity")
	hiRe=TH1D("Real_part_of_"+dataSet[0][1]['wave1'].strip()+'-'+dataSet[0][1]['wave2'].strip(),"Real_part_of_"+dataSet[0][1]['wave1'].strip()+'-'+dataSet[0][1]['wave2'].strip(),len(dataSet),binning)
	hiRe.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	hiRe.GetYaxis().SetTitle("Real part")
	hiIm=TH1D("Imaginary_part_of_"+dataSet[0][1]['wave1'].strip()+'-'+dataSet[0][1]['wave2'].strip(),"Imaginary_part_of_"+dataSet[0][1]['wave1'].strip()+'-'+dataSet[0][1]['wave2'].strip(),len(dataSet),binning)
	hiIm.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	hiIm.GetYaxis().SetTitle("Imaginary part")
	hiPh=TH1D("Phase_of_"+dataSet[0][1]['wave1'].strip()+'-'+dataSet[0][1]['wave2'].strip(),"Phase_of_"+dataSet[0][1]['wave1'].strip()+'-'+dataSet[0][1]['wave2'].strip(),len(dataSet),binning)
	hiPh.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	hiPh.GetYaxis().SetTitle("Phase")
	for i in range(0,len(dataSet)):
		hiInt1.SetBinContent(i+1,dataSet[i][0][2])
		hiInt2.SetBinContent(i+1,dataSet[i][0][4])
		hiRe.SetBinContent(i+1,dataSet[i][0][6])
		hiIm.SetBinContent(i+1,dataSet[i][0][8])
		hiPh.SetBinContent(i+1,dataSet[i][0][10])
		hiInt1.SetBinError(i+1,dataSet[i][0][3])
		hiInt2.SetBinError(i+1,dataSet[i][0][5])
		hiRe.SetBinError(i+1,dataSet[i][0][7])
		hiIm.SetBinError(i+1,dataSet[i][0][9])
		hiPh.SetBinError(i+1,dataSet[i][0][11])
	hiInt1.Write()
	hiInt2.Write()
	hiRe.Write()
	hiIm.Write()
	hiPh.Write()
	outROOT.Close()
#------------------------------------------------------------------------------------------------------------------------------------
def print2DtoRoot(
			dataSet,
			outFile='fit2D.root',
			cuts={'0-+':[{'cut':'m3pi','location':1.80}],'1++':[],'2-+':[]},		
			jpcs=['0-+','1-+','1++','2-+','2++','4++'],
			Ms=['0','1','2'],
			isobarrrs=['f0_','rho_','f2_'],		
			suppressSigma=0		):
 	"""Prints two dimensional data to a .root file"""
	count_calls('print2DtoRoot')
	outROOT=root_open('./ROOT/'+outFile,mode="RECREATE")
	for isobarrr in isobarrrs:
		for jpc in jpcs:
			for M in Ms:
				addString = '_'+isobarrr
				bins2Pi=[]
				bins3Pi=[]
				for i in range(0,len(dataSet[0])):
					if dataSet[0][i][4]==jpc and dataSet[0][i][13] == M and isobarrr in dataSet[0][i][14]:
						bins3Pi.append(dataSet[0][i][0])
						bins3Pi.append(dataSet[0][i][1])
						bins2Pi.append(dataSet[0][i][2])
						bins2Pi.append(dataSet[0][i][3])
				bins3Pi.sort()
				bins2Pi.sort()
				if not len(bins3Pi) + len(bins2Pi) == 0:
					binning2Pi=[bins2Pi[0]]
					binning3Pi=[bins3Pi[0]]
					for i in range(1,len(bins3Pi)-1):
						if binning3Pi[-1] != bins3Pi[i]:
							binning3Pi.append(bins3Pi[i])
							if bins3Pi[i] != bins3Pi[i+1]:
								print "Warning: Binning in m(3Pi) is wrong."
					binning3Pi.append(bins3Pi[-1])
					for i in range(1,len(bins2Pi)-1):
						if binning2Pi[-1] != bins2Pi[i]:
							binning2Pi.append(bins2Pi[i])
							if bins2Pi[i] != bins2Pi[i+1]:
								print "Warning: Binning in m(2Pi) is wrong."
					binning2Pi.append(bins2Pi[-1])
					binning2Pi= numpy.asarray(binning2Pi,dtype=numpy.float64)
					binning3Pi= numpy.asarray(binning3Pi,dtype=numpy.float64)
					histIn  = TH2D("Intensity of "+jpc+'_'+M+addString,"Intensity of "+jpc+addString,len(binning3Pi)-1,binning3Pi,len(binning2Pi)-1,binning2Pi)
					histIn.SetDrawOption('col')
					histIn.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histIn.GetYaxis().SetTitle("Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histIn.GetZaxis().SetTitle("Intensity of "+jpc+addString)
					histRe  = TH2D("Real part of "+jpc+'_'+M+addString,"Real part of "+jpc+addString,len(binning3Pi)-1,binning3Pi,len(binning2Pi)-1,binning2Pi)
					histRe.SetDrawOption('col')
					histRe.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histRe.GetYaxis().SetTitle("Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histRe.GetZaxis().SetTitle("Real part of "+jpc+addString)
					histIm  = TH2D("Imag part of "+jpc+'_'+M+addString,"Imag part of "+jpc+addString,len(binning3Pi)-1,binning3Pi,len(binning2Pi)-1,binning2Pi)
					histIm.SetDrawOption('col')
					histIm.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histIm.GetYaxis().SetTitle("Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histIm.GetZaxis().SetTitle("Imag part of "+jpc+addString)
					histPh  = TH2D("Phase of "+jpc+'_'+M+addString,"Phase of "+jpc+addString,len(binning3Pi)-1,binning3Pi,len(binning2Pi)-1,binning2Pi)
					histPh.SetDrawOption('col')
					histPh.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histPh.GetYaxis().SetTitle("Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histPh.GetZaxis().SetTitle("Phase of "+'_'+M+jpc+addString)
					for i in range(0,len(dataSet[0])):
						if dataSet[0][i][4] == jpc and dataSet[0][i][13] == M and isobarrr in dataSet[0][i][14]:
							m2Center = (dataSet[0][i][3] + dataSet[0][i][2])/2
							m3Center = (dataSet[0][i][0] + dataSet[0][i][1])/2
							n2 = histIn.GetYaxis().FindBin(m2Center)
							n3 = histIn.GetXaxis().FindBin(m3Center)	
							histIn.SetBinContent(n3,n2,dataSet[0][i][5])
							histRe.SetBinContent(n3,n2,dataSet[0][i][7])
							histIm.SetBinContent(n3,n2,dataSet[0][i][9])
							histPh.SetBinContent(n3,n2,dataSet[0][i][11])
							histIn.SetBinError(n3,n2,dataSet[0][i][6])
							histRe.SetBinError(n3,n2,dataSet[0][i][8])
							histIm.SetBinError(n3,n2,dataSet[0][i][10])
							histPh.SetBinError(n3,n2,dataSet[0][i][12])
							if histIn.GetBinContent(n3,n2) < suppressSigma * histIn.GetBinError(n3,n2):
								histIn.SetBinContent(n3,n2,0.)
								histIn.SetBinError(n3,n2,0.)
					histPh=removePhaseAmbiguities(histPh)
					histIn.Write()
					histRe.Write()
					histIm.Write()
					histPh.Write()
					try:
						for cut in cuts[jpc]:
							direct=cut['cut']
							mass=cut['location']
							if direct == 'm3pi':
								nBin=histIn.GetXaxis().FindBin(mass)
								cutHist=histIn.ProjectionY("m(3#pi)="+str(mass),nBin,nBin,"e")
								cutHist.Write()
							elif direct == 'm2pi':
								nBin=histIn.GetYaxis().FindBin(mass)
								cutHist=histIn.ProjectionX("m(2#pi)="+str(mass),nBin,nBin,"e")
								cutHist.Write()
					except:
						print "No cuts for JPC = "+jpc+'_'+M
				else:
					print "Nothing for "+jpc+"_"+M
	outROOT.close()
#------------------------------------------------------------------------------------------------------------------------------------
def printStatusToRoot(dataSet,ROOTname='fit_status.root'):
	"""Prints a fit status to .root"""
	count_calls('printStatusToRoot')
	binning=[]
	tString="t' = "+str(dataSet[0][4])+'-'+str(dataSet[0][5])+'(GeV/#it{c})^{2}'
	for data in dataSet:
		if not data[0] in binning:
			binning.append(data[0])
		if not data[1] in binning:
			binning.append(data[1])
	binning.sort()
	binning = numpy.asarray(binning,dtype=numpy.float64)
	eventsHist=TH1D("Number of events",tString,len(binning)-1,binning)
	eventsHist.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	eventsHist.GetYaxis().SetTitle("nEvents")
	maxLikeHist=TH1D("Maximum Likelihood",tString,len(binning)-1,binning)
	maxLikeHist.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	maxLikeHist.GetYaxis().SetTitle("max log like")
	maxLikePerEventHist=TH1D("Maximum Likelihood per event",tString,len(binning)-1,binning)
	maxLikePerEventHist.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	maxLikePerEventHist.GetYaxis().SetTitle("max log like per event")
	attemptsHist=TH1D("Number of attempts",tString,len(binning)-1,binning)
	attemptsHist.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	attemptsHist.GetYaxis().SetTitle("fit attempts")
	for data in dataSet:
		bin=eventsHist.GetXaxis().FindBin((data[0]+data[1])/2)
		if not data[2] ==0:	
			eventsHist.SetBinContent(bin,data[2])
			if data[3] > maxLikeHist.GetBinContent(bin):
				maxLikeHist.SetBinContent(bin,data[3])
				maxLikePerEventHist.SetBinContent(bin,float(data[3])/data[2])
			attemptsHist.SetBinContent(bin,attemptsHist.GetBinContent(bin)+1)
	likeDistHist=TH2D('Likelihood distribution',tString,len(binning)-1,binning,100,0,100)
	likeDistHist.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
	likeDistHist.GetYaxis().SetTitle("Likelihood difference")
	for data in dataSet:
		bin=eventsHist.GetXaxis().FindBin((data[0]+data[1])/2)
		if not data[2] == 0:	
			dist=-int(math.ceil(data[3] - maxLikeHist.GetBinContent(bin)))+1
			if dist < 100:
				likeDistHist.SetBinContent(bin,dist,likeDistHist.GetBinContent(bin,dist)+1)
	outROOT = root_open(ROOTname, mode = "RECREATE")	
	eventsHist.Write()	
	maxLikeHist.Write()
	maxLikePerEventHist.Write()
	attemptsHist.Write()
	likeDistHist.Write()
	outROOT.Close()
#------------------------------------------------------------------------------------------------------------------------------------
def prettyPrint(matrix):
	"""Prints 'matrix' in a nice way"""
	count_calls('prettyPrint')
	s = [[str(e) for e in row] for row in matrix]
	lens = [max(map(len, col)) for col in zip(*s)]
	fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
	table = [fmt.format(*row) for row in s]
	print '\n'.join(table)
#####################################################################################################################################
##								SOME FUNCTIONS							   ##
#####################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------
def init_array(n,zero = 0.):
	"""Initializes an n x n array with 'zero' in every entry"""
	count_calls('init_array')
	arr = []
	for i in range(n):
		line =[]
		for j in range(n):
			line.append(zero)
		arr.append(line)
	return arr
#------------------------------------------------------------------------------------------------------------------------------------
def isHermitian(matrix):
	"""Checks if the matrix is hermitian"""
	count_calls('isHermitian')
	for i in range(0,len(matrix)):
		for j in range(0,len(matrix)):
			if matrix[i][j] != matrix[j][i].conjugate():
				return False
	return True
#------------------------------------------------------------------------------------------------------------------------------------
def hermitianize(matrix):
	"""Makes the matrix hermitian"""
	count_calls('hermitianize')
	for i in range(0,len(matrix)):
		for j in range(0,i):
			matrix[i][j]=matrix[j][i].conjugate()
	return matrix
#------------------------------------------------------------------------------------------------------------------------------------
def isSymmetric(matrix):
	"""Checks if the matrix is symmetric"""
	count_calls('isSymmetric')
	for i in range(0,len(matrix)):
		for j in range(0,len(matrix)):
			if abs(matrix[i][j] - matrix[j][i])>1.E-10:
				print str(i)+'   '+str(j)+'   '+str(matrix[i][j])+' - '+str(matrix[j][i])+' = '+str(matrix[j][i]-matrix[i][j]) 
				return False
	return True
#------------------------------------------------------------------------------------------------------------------------------------
def vectorMatrixVector(vector,matrix):
	"""Returns vector^T.matrix.vector"""
	count_calls('vectorMatrixVector')
	if len(matrix) != len(matrix[0]) or len(vector) != len(matrix):
		return 'ERROR'
	else:
		val=0
		for i in range(0,len(vector)):
			for j in range(0,len(vector)):
				val+=vector[i]*matrix[i][j]*vector[j]
		return val
#------------------------------------------------------------------------------------------------------------------------------------
def rescale(matrix,vector):
	"""Returns matrix[i][j] -> matrix[i][j]*vector[i]*vector[j]"""
	count_calls('rescale')
	for i in range(0,len(matrix)):
		for j in range(0,len(matrix[0])):
			matrix[i][j]*=vector[i]*vector[j]
#------------------------------------------------------------------------------------------------------------------------------------
def getFourDigits(intIn):
	"""Returns a string of the intIn, that has at least four digits (1 -> '0001')"""
	count_calls('getFourDigits')
	string = str(intIn)
	while len(string)<4:
		string='0'+string
	return string
#------------------------------------------------------------------------------------------------------------------------------------
def getListOfSteplikes(wave='0-+',deisobarredWave='f0_'):
	"""Returns a list of steplike functions to be used in getTotal(...)"""
	count_calls('getListOfSteplikes')
	if deisobarredWave=='f0_':
		prefixes={'0-+':'1-(0-+)0+ f0_','1++':'1-(1++)0+ f0_','2-+':'1-(2-+)0+ f0_'}
		suffixes={'0-+':' pi S                                 ','1++':' pi P                                 ','2-+':' pi D                                 '}
		stepList=[]
		massStep=10
		massStep2=40
		massStart=320
		actMass=massStart
		massStop=2500
		while(actMass<massStop): 
			stepList.append(prefixes[wave]+getFourDigits(actMass)+'_'+getFourDigits(actMass+massStep)+suffixes[wave])
			actMass+=massStep
		actMass=massStart
		while(actMass<massStop): 
			stepList.append(prefixes[wave]+getFourDigits(actMass)+'_'+getFourDigits(actMass+massStep2)+suffixes[wave])
			actMass+=massStep2
		stepList.append(prefixes[wave]+getFourDigits(278)+'_'+getFourDigits(320)+suffixes[wave])
		return stepList
#------------------------------------------------------------------------------------------------------------------------------------
def getListOfQuantumNumbers(j='',p='',c='',m=''):
	"""Returns a list of waves with certain quantum numbers"""
	count_calls('getListOfQuantumNumbers')
	wavesData=open('waveData.txt','r')
	waves=[]
	for wave in wavesData.readlines():
		waveStrip=wave.rstrip('\n')
		if (j=='' or j==waveStrip[3]) and (p=='' or p==waveStrip[4]) and (c=='' or c==waveStrip[5]) and (m=='' or m==waveStrip[7]) and wave[8]=='+':
			waves.append(waveStrip)
	wavesData.close()
	return waves

#------------------------------------------------------------------------------------------------------------------------------------
def writeWaveList(inFile):
	"""Writes the list of all waves to a file"""
	count_calls('writeWaveList')
	inFile=open(inFile,'r')
	outFile=open('waveData.txt','w')
	for line in inFile.readlines():
		if line.startswith("'"):
			outFile.write(line[1:61]+'\n')
	inFile.close()
	outFile.close()
#------------------------------------------------------------------------------------------------------------------------------------
def removePhaseAmbiguities(histIn):
	"""Removes phase ambiguities from TH1D and TH2D"""
	count_calls('removePhaseAmbiguities')
	hist=histIn
	className=hist.ClassName()
	if className.startswith('TH1'):
		for i in range(2,hist.GetNbinsX()+1):
			while True:
				act=hist.GetBinContent(i)
				err=hist.GetBinError(i)
				old=hist.GetBinContent(i-1)
				if act > old + pi:
					hist.SetBinContent(i,act -2*pi)
				elif act < old -pi:
					hist.SetBinContent(i,act+2*pi)
				else:
					break
		if hist.GetBinError(i) ==0.:
			hist.SetBinContent(i,0.)
	if className.startswith('TH2'):
		for j in range(1,hist.GetNbinsX()+1):
			for i in range(2,hist.GetNbinsY()+1):
				while True:
					act=hist.GetBinContent(j,i)
					err=hist.GetBinError(j,i)
					old=hist.GetBinContent(j,i-1)
					if act > old + pi and act !=0. and err !=0.:
						hist.SetBinContent(j,i,act -2*pi)
					elif act < old -pi and act !=0. and err !=0.:
						hist.SetBinContent(j,i,act+2*pi)
					else:
						break
				if hist.GetBinError(j,i)==0.:
					hist.SetBinContent(j,i,0.)
		for i in range(2,hist.GetNbinsX()+1):
			while True:
				act=hist.GetBinContent(i,1)
				err=hist.GetBinError(i,1)
				old=hist.GetBinContent(i-1,1)
				if act > old + pi and act !=0. and err !=0.:
					for j in range(1,hist.GetNbinsY()+1):
						actSub=hist.GetBinContent(i,j)
						errSub=hist.GetBinContent(i,j)
						if actSub != 0. and errSub != 0.:
							hist.SetBinContent(i,j,actSub - 2*pi)
				elif act < old -pi and act !=0. and err !=0.:
					for j in range(1,hist.GetNbinsY()+1):
						actSub=hist.GetBinContent(i,j)
						errSub=hist.GetBinContent(i,j)
						if actSub != 0. and errSub != 0.:
							hist.SetBinContent(i,j,actSub + 2*pi)
				else:
					break
	return hist
#------------------------------------------------------------------------------------------------------------------------------------
def getJPCfromWave(wave): 
	""" Get the corresponding JPC from a wave (Works only for regular waves (not FLAT or reflectivity = -1))"""
	count_calls('getJPCfromWave')
	jpc=wave[:9]
	if 'f0' in wave or '(pipi)' in wave:
		jpc = jpc +'(0++)'
	if 'f2' in wave:
		jpc = jpc+'(2++)'
	if 'rho' in wave and not 'rho3' in wave:
		jpc = jpc+'(1--)'
	if 'rho3' in wave:
		jpc=jpc+'(3--)'
	for i in range(1,len(wave)):
		if not wave[-i] == ' ':
			jpc = jpc+'pi'+wave[-i]
			break
	return jpc
#####################################################################################################################################
##							DIFFERENT MATRIX INVERSIONS						   ##
#####################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------
#def invert_right_submatrix(full_coma,data_point):
#	"""Returns the inverse of the submatrix of 'full_coma', where corresponding entries in 'data_point' are nonzero"""
#	count_calls('invert_right_submatrix')
#	subComa = []
#	index_map={}
#	nzin = 0 #NonZeroIndexNumber
#	for i in range(len(data_point)):
#		if not data_point[i] == 0.:
#			index_map[i] = nzin
#			nzin+=1
#		else:
#			index_map[i] =  -1
#	for i in range(len(data_point)):
#		if not data_point[i] == 0.:
#			subComaLine=[]
#			for j in range(len(data_point)):
#				if not data_point[j] == 0.:
#					subComaLine.append(full_coma[i,j])		
#			subComa.append(subComaLine)
#	subComa = numpy.matrix(subComa)
#	print subComa
#	subComa_inv = la.inv(subComa)
#	final_coma_inv=[]
#	for i in range(len(data_point)):
#		final_coma_inv_line=[]
#		for j in range(len(data_point)):
#			if index_map[i] < 0 or index_map[j] <0: # data_point is zero -> leave point out
#				final_coma_inv_line.append(0.)
#			else:
#				final_coma_inv_line.append(subComa_inv[index_map[i],index_map[j]])
#		final_coma_inv.append(final_coma_inv_line)
#	return final_coma_inv
#------------------------------------------------------------------------------------------------------------------------------------
def invert_right_submatrix(coma, flag = 'PINV', epsilon = 1.E-3):
	"""
	Selects different inversion methods for covariance matrices:
	'EPSILON': Daigonalizes the matrix, sets every EV smaller than epsilon to epsilon and inverts the matrix
	'EPSILON_NONZERO': Does the same as 'EPSILON' but leaves out zero lines and columns first and the puts them back in
	'OMIT_ZERO_EV': Diagonalizes the matrix, inverts every nonzero EV and rotates back
	'INTENS_REAL': Inverts only the submatrix for intensities and real parts, sets the rest to zero
	'DIAG_ERRORS': Just sets the diagonal elements to 1/coma[i][i]. If this is set, chi2coma.GetComaFromErrors() has to be used in addition
	'ANCHOR_FIRST': Inverts only the submatrix, where the first wave appears, sets all other elements to zero.
	"""
	count_calls('invert_right_submatrix')
	if flag == 'EPSILON':
		return invert_epsilon(coma,epsilon)
	if flag == 'EPSILON_NONZERO':
		return invert_epsilon_nonzero(coma,epsilon)
	if flag == 'OMIT_ZERO_EV':
		return invert_epsilon(coma,0.)
	if flag == 'INTENS_REAL':
		return invert_intens_real(coma)
	if flag == 'DIAG_ERRORS':
		return invert_epsilon(coma,0.,False)
	if flag == 'ANCHOR_FIRST':
		return invert_anchor(coma)
	if flag == 'PINV':
		return la.pinv(coma).tolist()
	print "No valid inversion-method specified"
	return coma
#------------------------------------------------------------------------------------------------------------------------------------
def invert_epsilon(full_coma, epsilon =0., TRUE_INVERT=True):
	"""Returns the invertes coma, without modes with eigenvalue < 'minEV'"""
	count_calls('invert_epsilon')
	# if True: Inverts COMA with regularized Eigenvalues, else: Gives back COMA with inverted diagonal elements (DO NOT CHANGE UNLESS YOU KNOW WHAT YOUR'RE DOING)
	limm = epsilon # Upper limit for an EV to count as zero
	if limm==0.:
		limm=1.E-8
	if TRUE_INVERT:
		exponent = -1.
	else:
		exponent = 1.
	val, vec = la.eig(full_coma)
	inverted_vals=[]
	for i in range(len(vec)):
		if abs(val[i]) > limm:
			inverted_vals.append(val[i]**exponent)
		else:
			if exponent == 0.:
				inverted_vals.append(0.)
			else:
				inverted_vals.append(epsilon**exponent)
	diag = numpy.zeros((len(val),len(val)))
	for i in range(len(val)):
		diag[i][i] = inverted_vals[i]
	diag = numpy.matrix(diag)
	if TRUE_INVERT:
		return (vec*diag*la.inv(vec)).tolist()
	else:
		licoma = (vec*diag*la.inv(vec)).tolist()
		print "WARNING: For testing, not the real inverted COMA is returned, but just the normal COMA with an inverted diagonal"
		for i in range(len(licoma)):
			if not licoma[i][i]==0.:
				licoma[i][i]**=-1
		return licoma
#------------------------------------------------------------------------------------------------------------------------------------
def invert_epsilon_nonzero(full_coma_with_zeros, epsilon = 1.E-3):
	"""Returns the invertes coma, without modes with eigenvalue < 'minEV' omitting complete zero lines"""
	count_calls('invert_epsilon_nonzero')
	fcwz = full_coma_with_zeros.tolist()
	nWaves = len(fcwz)
	zerolines =[]
	for i in range(nWaves):			
		nnzero=0			
		for j in range(nWaves):		
			if not fcwz[i][j] == 0:	
				nnzero+=1	
				break		
		if nnzero == 0:			
			zerolines.append(i)	
	fcnz = []				
	for i in range(nWaves):			
		if not i in zerolines:
			fcnz_line = []		
			for j in range(nWaves):					
				if not j in zerolines:
					fcnz_line.append(fcwz[i][j])
			fcnz.append(fcnz_line)
	fcnz = numpy.matrix(fcnz)
	fcnz_inv = invert_epsilon(fcnz,epsilon)
	fcwz_inv=[]
	line_count=0
	zero_line=[0.]*nWaves
	for i in range(nWaves):	
		if i in zerolines:
			fcwz_inv.append(zero_line)
		else:
			col_count =0
			fcwz_inv_line=[]
			for j in range(nWaves):
				if j in zerolines:
					fcwz_inv_line.append(0.)
				else:
					fcwz_inv_line.append(fcnz_inv[line_count][col_count])
					col_count+=1
			fcwz_inv.append(fcwz_inv_line)
			line_count+=1
	if len(fcnz)>1:
		val, vec = la.eig(fcnz)
		print val
		print "DET: "+str(la.det(fcnz))
		raw_input()
	return fcwz_inv

#------------------------------------------------------------------------------------------------------------------------------------
def invert_intens_real(full_coma):
	"""Inverts the submatrix corresponding only to intensities and Re(T_i T_{i+1}^*)"""
	count_calls('invert_intens_real')
	nWaves = int(len(full_coma)**.5)
	intReIndices=[]
	for i in range(nWaves):
		intReIndices.append(i*(nWaves+1))
		if i*(nWaves+1)+1 < nWaves**2: # Does not look for Re(T_{n-1} T_n^*) T_n does not exist
			intReIndices.append(i*(nWaves+1)+1)
	subcoma=[[0.]*len(intReIndices) for iiii in range(len(intReIndices))]
	for i in range(len(intReIndices)):
		ii = intReIndices[i]
		for j in range(len(intReIndices)):
			jj = intReIndices[j]
			subcoma[i][j] = full_coma[ii,jj]
	subcoma = numpy.matrix(subcoma)
	subcoma_inv = invert_nonzero(subcoma)
	full_coma_inv=numpy.zeros((len(full_coma),len(full_coma)))
	full_coma_inv=numpy.matrix(full_coma_inv)
	for i in range(len(intReIndices)):
		ii = intReIndices[i]
		for j in range(len(intReIndices)):
			jj = intReIndices[j]
			full_coma_inv[ii,jj] = subcoma_inv[i,j]
	return full_coma_inv.tolist()
#------------------------------------------------------------------------------------------------------------------------------------
def invert_anchor(full_coma):
	"""Inverts the submatrix corresponding only to interferences with the first wave"""	
	count_calls('invert_anchor')
	nWaves = int(len(full_coma)**.5)
	intReIndices=[]
	for i in range(len(full_coma)):
		ii = int(i/nWaves)
		jj = i - ii*nWaves
		if ii ==0 or jj == 0:
			intReIndices.append(i)
	subcoma=[[0.]*len(intReIndices) for iiii in range(len(intReIndices))]
	for i in range(len(intReIndices)):
		ii = intReIndices[i]
		for j in range(len(intReIndices)):
			jj = intReIndices[j]
			subcoma[i][j] = full_coma[ii,jj]
	subcoma = numpy.matrix(subcoma)
	subcoma_inv = invert_nonzero(subcoma)
	full_coma_inv=numpy.zeros((len(full_coma),len(full_coma)))
	full_coma_inv=numpy.matrix(full_coma_inv)
	for i in range(len(intReIndices)):
		ii = intReIndices[i]
		for j in range(len(intReIndices)):
			jj = intReIndices[j]
			full_coma_inv[ii,jj] = subcoma_inv[i,j]
	return full_coma_inv.tolist()
	
#------------------------------------------------------------------------------------------------------------------------------------
def invert_nonzero(matrix):
	"""Inverts matrix, leaving out zero-lines, that would render the matrix singular"""
	count_calls('invert_nonzero')
	n = len(matrix)
	arr = matrix.tolist()
	zeros = []
	for i in range(n):
		line = arr[i]
		nnz = 0
		for val in line:
			if not val ==0.:
				nnz+=1
				break
		if nnz == 0:
			zeros.append(i)
	nonzero_matrix = []
	for i in range(n):
		if not i in zeros:
			nonzero_line=[]
			for j in range(n):
				if not j in zeros:
					nonzero_line.append(matrix[i,j])
			nonzero_matrix.append(nonzero_line)
	if not nonzero_matrix==[]:
		nonzero_matrix=numpy.matrix(nonzero_matrix)
		nonzero_inv = la.inv(nonzero_matrix)
	else:
		nonzero_inv = numpy.matrix(nonzero_matrix)
	full_inv=numpy.zeros((n,n))
	full_inv=numpy.matrix(full_inv)
	inonz=0
	for i in range(n):
		if not i in zeros:
			jnonz=0
			for j in range(n):
				if not j in zeros:
					full_inv[i,j] = nonzero_inv[inonz,jnonz]
					jnonz+=1
			inonz+=1
	return full_inv

#####################################################################################################################################
##							Scripts									   ##
#####################################################################################################################################

def do2D(name):
	count_calls('do2D')
	path = "/nfs/mds/user/fkrinner/massIndepententFits/fits/"+name+"/fit/0.14077-0.19435/"
	for kw in keyWaves:
		try:
			qw = get2D(path, keywave = kw)
			break
		except:
			print kw.strip()+" not in the waveset. Try next"
			pass
	print2DtoRoot(qw,name+'.root')


#####################################################################################################################################
##							COMA COMA TUEDILUE							   ##
#####################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------
def getSDM_complex(SDM_v):
	"""Builds the complex spin-density-matrix out of the real-vector from used by the program"""
	count_calls('getSDM_complex')
	nWaves = int(len(SDM_v)**.5)
	SDM_c = init_array(nWaves,0.+0.j)
	for i in range(nWaves**2):
		ii = int(i/nWaves)
		jj = i - ii*nWaves
		SDM_c[ii][jj] += SDM_v[i]	# This gets only the j < i entries right
		SDM_c[jj][ii] += SDM_v[i]*1.j	#
	for i in range(nWaves):
		for j in range(nWaves):
			if i==j: # The diagonal elements are wrong, fis this here
				SDM_c[i][j] = SDM_c[i][j].real+0.j
			if j > i: # the j > i entries are set wrong before, fix this here. 
				SDM_c[i][j] = SDM_c[j][i].conjugate()
	return SDM_c
#------------------------------------------------------------------------------------------------------------------------------------
def getSubSDM_c(SDM_c):
	"""Returns SDM with only [i][i] and [i][i+1] entries nonzero"""
	count_calls('getSubSDM_c')
	nWaves=len(SDM_c)
	SDM_s = init_array(nWaves,0.+0.j)
	for i in range(nWaves):
		SDM_s[i][i] = SDM_c[i][i]
		if not i == nWaves-1:
			SDM_s[i][i+1] = SDM_c[i][i+1]
	return SDM_s
#------------------------------------------------------------------------------------------------------------------------------------
def getSubSDM_v(SDM_v):
	"""Returns vector SDM with only [i][i] and [i][i+1] entries nonzero"""
	count_calls('getSubSDM_v')
	nWaves = int(len(SDM_v)**.5)
	SDM_s=SDM_v[:]
	for i in range(len(SDM_v)):
		ii = int(i/nWaves)
		jj = i-ii*nWaves
		if abs(ii-jj)>1:
			SDM_s[i]=0.
	return SDM_s	
#------------------------------------------------------------------------------------------------------------------------------------
def reconstrucSDM_c(SDM_s):
	"""Reconscructs the full SDM from the getSubSDM(...) output"""
	count_calls('reconstrucSDM_c')
	nWaves = len(SDM_s)	
	SDM_c = init_array(nWaves,0.+0.j)
	for i in range(nWaves):
		SDM_c[i][i] = SDM_s[i][i]
		if not i == nWaves-1:
			SDM_c[i][i+1] = SDM_s[i][i+1]
	for n in range(2,nWaves):
		for i in range(nWaves):
			if not i+n > nWaves-1:
				SDM_c[i][i+n] = inferSDM(i,i+n,SDM_c)
	for i in range(nWaves):
		for j in range(i):
			SDM_c[i][j] = SDM_c[j][i].conjugate()	
	return SDM_c
#------------------------------------------------------------------------------------------------------------------------------------
def getSubComa(coma,SDM_vr):
	count_calls('getSubComa')
	if not len(coma) == len(SDM_vr):
		print "Dimensions differ. Abort"
		return
	coma_s = init_array(len(coma),0.)
	for i in range(len(coma)):
		if not SDM_vr[i]==0.:
			for j in range(len(coma)):
				if not SDM_vr[j] ==0.:
					coma_s[i][j] = coma[i][j]
	return coma_s
#------------------------------------------------------------------------------------------------------------------------------------
def inferSDM(i,j,SDM):
	"""Calculates: SDM[i][j] = SDM[i+1][j] SDM[i][j-1] / SDM[i+1][j-1]^* """
	count_calls('inferSDM')
	return SDM[i+1][j]*SDM[i][j-1] / (SDM[i+1][j-1])
#------------------------------------------------------------------------------------------------------------------------------------
def rotate_phi(phi, vec):
	"""Returns vec = {re, im, re, im,...}, rotated in the complex plane"""
	count_calls('rotate_phi')
	vec_r = vec[:]
	for i in range(len(vec)/2):
		vec_r[2*i] = math.cos(phi)*vec[2*i] - math.sin(phi)*vec[2*i+1]
		vec_r[2*i+1] = math.cos(phi)*vec[2*i+1] + math.sin(phi)*vec[2*i]
	return vec_r
#------------------------------------------------------------------------------------------------------------------------------------
def rotateT(T):
	count_calls('rotateT')
	ab = (T[0]**2+T[1]**2)**.5
	c = T[0]/ab
	s = T[1]/ab
	Tr = []
	for i in range(len(T)/2):
		Tr.append(c*T[2*i]+s*T[2*i+1])
		Tr.append(-s*T[2*i]+c*T[2*i+1])
	return Tr
#------------------------------------------------------------------------------------------------------------------------------------
def numericalJac(T):
	count_calls('numericalJac')
	epsilon=1.E-10
	jac=[]
	Tr = rotateT(T)
	for i in range(len(T)):
		Teps = T[:]
		Teps[i]+=epsilon
		Treps = rotateT(Teps)
		for j in range(len(T)):
			Treps[j]-=Tr[j]
			Treps[j]/=epsilon
		jac.append(Treps)
	prettyPrint(jac)
#------------------------------------------------------------------------------------------------------------------------------------
def count_calls(name):
	"""Since there are many methods in this file, and most of them are unused, count their calls to see, which ones are needed"""
	try:
		inf = open('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput/count_function_calls/'+name,'r')
		try:
			count = int(inf.read())
		except ValueError:
			count =0
		inf.close()
	except IOError:
		count = 0
	outf = open('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput/count_function_calls/'+name,'w')
	outf.write(str(count+1))
	outf.close()
#------------------------------------------------------------------------------------------------------------------------------------

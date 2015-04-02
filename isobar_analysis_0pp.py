from convertTextOutput import get2D
from convertTextOutput import removePhaseAmbiguities
from convertTextOutput import getTotal
import numpy
import rootpy
import ROOT
rootpy.log.basic_config_colorized()
from rootpy.io import root_open
from rootpy.plotting import Canvas, Hist, Legend
from ROOT import TH1D
from ROOT import TH2D
from ROOT import TGraphErrors
from ROOT import TMultiGraph

canv = Canvas(name = "canv0", title = "PWA")

def isobarred_analysis_0pp():
	jpcs = ['0-+','1++','2-+']
	tbins=['0.10000-0.14077','0.14077-0.19435','0.19435-0.32617','0.32617-1.00000']
	root_name='isobarred_fit.root'
	outROOT=root_open('./ROOT/'+root_name,mode="RECREATE")
	totallists={
		'0-+':['1-(0-+)0+ (pipi)_S pi S                                     ','1-(0-+)0+ f0(980) pi S                                      ','1-(0-+)0+ f0(1500) pi S                                     '],
		'1++':['1-(1++)0+ (pipi)_S pi P                                     ','1-(1++)0+ f0(980) pi P                                      '],
		'2-+':['1-(2-+)0+ (pipi)_S pi D                                     ','1-(2-+)0+ f0(980) pi D                                      ']}
	sums = {}
	for tbin in tbins:
		for jpc in jpcs:
			data = getTotal('/nfs/mds/user/fkrinner/massIndepententFits/fit/isobared/'+tbin,totallists[jpc],'/nfs/mds/user/fkrinner/massIndepententFits/integrals/pipiS/'+tbin,normalizeToDiag=True)
			binning = [data[0][0]]
			for point in data:
				binning.append(point[1])
			binning= numpy.asarray(binning,dtype=numpy.float64)
			mmin = binning[0]
			mmax = binning[-1]
			hist = TH1D('isobarred_0++_total_of_'+jpc+"_"+tbin,'isobarred_0++_total_of_'+jpc+"_"+tbin,len(binning)-1,binning)
			for i in range(len(data)):
				hist.SetBinContent(i+1,data[i][2])
				hist.SetBinError(i+1,data[i][3])
			hist.Write()
			hist.Draw()
			canv.Print('./pdfs/isobarred_'+jpc+"_"+tbin+".pdf")
			name = 'sum_over_tbins_'+jpc
			print 'name:',name
			if not sums.has_key(name):
				sums[name] = hist
				sums[name].SetName('isobarred_sum_over_tbins_'+jpc)
			else:
				sums[name].Add(hist)
	for key in sums.iterkeys():
		sums[key].Write()
		sums[key].Draw()
		canv.Print('./pdfs/'+sums[key].GetName()+".pdf")
	outROOT.close()	

def to_60(string):
	while len(string) < 60:
		string += ' '
	return string


def make_f0_wavelist(mmin,mmax,prefix,suffix):
	binns = ['0278', '0320', '0360', '0400', '0440', '0480', '0520', '0560', '0600', '0640', '0680', '0720', '0760', '0800', '0840', '0880', '0920', '0930', '0940', '0950', '0960', '0970', '0980', '0990', '1000', '1010', '1020', '1030', '1040', '1050', '1060', '1070', '1080', '1120', '1160', '1200', '1240', '1280', '1320', '1360', '1400', '1440', '1480', '1520', '1560', '1600', '1640', '1680', '1720', '1760', '1800', '1840', '1880', '1920', '1960', '2000', '2040', '2080', '2120', '2160', '2200', '2240', '2280']
	retwaves = []
	for i in range(len(binns)-1):
		up = float(binns[i+1])/1000
		low= float(binns[i])/1000
		if up >= mmin and low <= mmax:
			retwaves.append(to_60(prefix+binns[i]+'_'+binns[i+1]+suffix))
	return retwaves

def isobar_analysis_0pp():
	"""Does the 2D isobar analysis, based on 'print2DtoRoot', with all cuts, argands and what not"""
	root_name='isobar_analysis.root'
	outROOT=root_open('./ROOT/'+root_name,mode="RECREATE")
	isobar = 'f0_'
	jpcs = ['0-+','1++','2-+']
#	jpcs = ['0-+']
	M='0'
	iso_slices = { # do not use bin borders in definitions, a bin will be used, if any part of the defined interval overlaps with the bin #hence the *.**1 and *.**9 at the end of each definition
		'0-+':[[1.661,1.699,'below_resonance'],[1.781,1.819,'on_resonance'],[1.901,1.939,'above_resonance']],
		'1++':[[1.261,1.299,'below_resonance'],[1.381,1.419,'on_resonance'],[1.501,1.539,'above_resonance']],
		'2-+':[[1.781,1.819,'below_resonance'],[1.901,1.939,'on_resonance'],[2.021,2.059,'above_resonance']]
	}
	prefixes = {'0-+':'1-(0-+)0+ f0_','1++':'1-(1++)0+ f0_','2-+':'1-(2-+)0+ f0_'}
        suffixes = {'0-+':' pi S','1++':' pi P','2-+':' pi D'}
	X_slices = [[0.961,0.999,'f_0(980)'],[1.401,1.559,'f_0(1500)'],[0.2781,2.279,'Incoherent_sum']]
	suppressSigma=0
	tbins=['0.10000-0.14077','0.14077-0.19435','0.19435-0.32617','0.32617-1.00000']
#	tbins=['0.10000-0.14077']
	sumintens={}
	for tbin in tbins:
		dataSet = get2D('/nfs/mds/user/fkrinner/massIndepententFits/fit/pipiS/'+tbin, '/nfs/mds/user/fkrinner/massIndepententFits/integrals/pipiS/'+tbin ,normalizeToIntegrals = False)
		for jpc in jpcs:
			addString = '_'+isobar
			bins2Pi=[]
			bins3Pi=[]
			for i in range(0,len(dataSet[0])):
				if dataSet[0][i][4]==jpc and dataSet[0][i][13] == M and isobar in dataSet[0][i][14]:
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
				histIn  = TH2D("Intensity of "+jpc+'_'+M+addString+"_"+tbin,"Intensity of "+jpc+addString,len(binning3Pi)-1,binning3Pi,len(binning2Pi)-1,binning2Pi)
				histIn.SetDrawOption('col')
				histIn.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
				histIn.GetYaxis().SetTitle("Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
				histIn.GetZaxis().SetTitle("Intensity of "+jpc+addString)
				histRe  = TH2D("Real part of "+jpc+'_'+M+addString+"_"+tbin,"Real part of "+jpc+addString,len(binning3Pi)-1,binning3Pi,len(binning2Pi)-1,binning2Pi)
				histRe.SetDrawOption('col')
				histRe.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
				histRe.GetYaxis().SetTitle("Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
				histRe.GetZaxis().SetTitle("Real part of "+jpc+addString)
				histIm  = TH2D("Imag part of "+jpc+'_'+M+addString+"_"+tbin,"Imag part of "+jpc+addString,len(binning3Pi)-1,binning3Pi,len(binning2Pi)-1,binning2Pi)
				histIm.SetDrawOption('col')
				histIm.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
				histIm.GetYaxis().SetTitle("Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
				histIm.GetZaxis().SetTitle("Imag part of "+jpc+addString)
				histPh  = TH2D("Phase of "+jpc+'_'+M+addString+"_"+tbin,"Phase of "+jpc+addString,len(binning3Pi)-1,binning3Pi,len(binning2Pi)-1,binning2Pi)
				histPh.SetDrawOption('col')
				histPh.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
				histPh.GetYaxis().SetTitle("Mass of the #pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
				histPh.GetZaxis().SetTitle("Phase of "+'_'+M+jpc+addString)
				for i in range(0,len(dataSet[0])):
					if dataSet[0][i][4] == jpc and dataSet[0][i][13] == M and isobar in dataSet[0][i][14]:
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
				histIn.Draw("col")
				canv.Print('./pdfs/'+histIn.GetName()+".pdf")
				for slic in iso_slices[jpc]:
					minm = slic[0]
					maxm = slic[1]
					name = slic[2]
					for i in range(len(binning3Pi)-1):
						bul = binning3Pi[i+1]
						bll = binning3Pi[i]
						if bul > maxm and bll < maxm:
							maxb = i+1
						if bul > minm and bll < minm:
							minb = i+1
					massstring = 'm3Pi='+str(binning3Pi[minb-1])+'-'+str(binning3Pi[maxb])
					slcInt=histIn.ProjectionY('slice_of_'+jpc+'_intens_'+name+"_"+tbin   ,minb,maxb)
					slcRe =histRe.ProjectionY('slice_of_'+jpc+'_real_part_'+name+"_"+tbin,minb,maxb)
					slcIm =histIm.ProjectionY('slice_of_'+jpc+'_imag_part_'+name+"_"+tbin,minb,maxb)
					slcPha=histPh.ProjectionY('slice_of_'+jpc+'_phase_'+name+"_"+tbin    ,minb,maxb)
					slcInt.Draw()
					canv.Print('./pdfs/'+slcInt.GetName()+".pdf")
					re = []
					im = []
					ree= []
					ime= []
					for i in range(1,slcRe.GetNbinsX()+1):
						e = slcRe.GetBinError(i)+ slcIm.GetBinError(i)
						if not e==0.:
							re.append(slcRe.GetBinContent(i))
							im.append(slcIm.GetBinContent(i))
							ree.append(slcRe.GetBinError(i))
							ime.append(slcIm.GetBinError(i))			
					if len(re)>0 and len(re) == len(im) and len(re) == len(ree) and len(re) == len(ime):
						while re[-1] ==0. and im[-1]==0.: # Kill the last (zero) point
							re = re[:-1]
							im = im[:-1]
							ree=ree[:-1]
							ime=ime[:-1]
						re= numpy.asarray(re,dtype=numpy.float64)
						im= numpy.asarray(im,dtype=numpy.float64)
						ree= numpy.asarray(ree,dtype=numpy.float64)
						ime= numpy.asarray(ime,dtype=numpy.float64)
						argand = TGraphErrors(len(re),re,im,ree,ime)
						argand.SetName('slice_of_'+jpc+'_argand_'+name+"_"+tbin)
#						argand.Write()
						argand.Draw('apl')
						canv.Print('./pdfs/'+argand.GetName()+".pdf")
						argand_wrapper = TMultiGraph()
						argand_wrapper.SetName(argand.GetName())
						argand_wrapper.Add(argand)
						argand_wrapper.Write()
					slcInt.SetTitle(massstring)
					slcRe.SetTitle(massstring)
					slcIm.SetTitle(massstring)
					slcPha.SetTitle(massstring)
					slcInt.Write()
					slcRe.Write()
					slcIm.Write()
					slcPha.Write()
				for slic in X_slices:
					minm = slic[0]
					maxm = slic[1]
					name = slic[2]
					total_list = make_f0_wavelist(minm,maxm,prefixes[jpc],suffixes[jpc])
#					if name == 'f_0(1500)':
#						print total_list
#						raise Exception
					total_data = getTotal('/nfs/mds/user/fkrinner/massIndepententFits/fit/pipiS/'+tbin,total_list, '/nfs/mds/user/fkrinner/massIndepententFits/integrals/pipiS/'+tbin,normalizeToDiag=True)
					total_binning = [total_data[0][0]]
					for total_point in total_data:
						total_binning.append(total_point[1])
					total_binning= numpy.asarray(total_binning,dtype=numpy.float64)
					total_mmin = total_binning[0]
					total_mmax = total_binning[-1]
					total_hist = TH1D('coherent_sum_of_'+jpc+'_'+name+'_'+tbin,'coherent_sum_of_'+jpc+'_'+name+'_'+tbin,len(total_binning)-1,total_binning)
					total_hist.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}}#pi^{#font[122]{-}} System (GeV/#it{c}^{2})")
					histIn.GetZaxis().SetTitle("Intensity of "+jpc+addString)
					for i in range(len(total_data)):
						total_hist.SetBinContent(i+1,total_data[i][2])
						total_hist.SetBinError(i+1,total_data[i][3])
					for i in range(len(binning2Pi)-1):
						bul = binning2Pi[i+1]
						bll = binning2Pi[i]
						if bul > maxm and bll < maxm:
							maxb = i+1
						if bul > minm and bll < minm:
							minb = i+1
					massstring = 'm2Pi='+str(binning2Pi[minb-1])+'-'+str(binning2Pi[maxb])
					slcInt=histIn.ProjectionX('slice_of_'+jpc+'_intens_'+name+"_"+tbin   ,minb,maxb)
					slcRe =histRe.ProjectionX('slice_of_'+jpc+'_real_part_'+name+"_"+tbin,minb,maxb)
					slcIm =histIm.ProjectionX('slice_of_'+jpc+'_imag_part_'+name+"_"+tbin,minb,maxb)
					slcPha=histPh.ProjectionX('slice_of_'+jpc+'_phase_'+name+"_"+tbin    ,minb,maxb)
					slcInt.SetTitle(massstring)
					slcRe.SetTitle(massstring)
					slcIm.SetTitle(massstring)
					slcPha.SetTitle(massstring)
					total_hist.SetTitle(massstring)
					slcInt.Write()
					slcRe.Write()
					slcIm.Write()
					slcPha.Write()
					total_hist.Write()
					slcInt.Draw()
					canv.Print('./pdfs/'+slcInt.GetName()+".pdf")
					if not sumintens.has_key(name+"_"+jpc):
						sumintens[name+"_"+jpc] = slcInt
						sumintens[name+"_"+jpc].SetName('incoherent_sum_'+jpc+'_'+name)
					else:
						sumintens[name+"_"+jpc].Add(slcInt)
	for key in sumintens.iterkeys():
		sumintens[key].Write()
		sumintens[key].Draw()
		canv.Print('./pdfs/'+sumintens[key].GetName()+".pdf")
	outROOT.close()
	print "ran with no exceptions"

if __name__ == "__main__":
	
	isobar_analysis_0pp()

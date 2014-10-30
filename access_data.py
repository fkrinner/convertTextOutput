import rootpy
rootpy.log.basic_config_colorized()
from rootpy.io import root_open
from rootpy.tree import Tree
import ROOT
from ROOT import TH2D
from rootpy.plotting import Canvas, Hist, Legend
import sys
import os
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything')
from random import random

from four_vector import four_vector
from offdiag import compare_diag_offdiag

m_Pi= 0.1395
m_p = 0.938272046


def get_events(mMin, mMax, directory = '/nfs/nas/data/compass/hadron/2008/comSkim/2008-binned/all/skim_2012/', n = 10000, treeName = 'USR52mb',preSelect = False):
	useFiles=[]
	for fn in os.listdir(directory):
		try:
			mmin = float(fn[:4])/1000.
			mmax = float(fn[5:9])/1000.
		except ValueError:
			continue
		if mmin >= mMin and mmax <= mMax:
			useFiles.append(fn)
	evtPerFile = n/len(useFiles)
	print useFiles
	print len(useFiles),"in the mass range -> take",evtPerFile,"events per file"
	returnEvents=[]
	for fn in useFiles:
		print "open file:",directory+"/"+fn
		ROOT_file = root_open(directory+'/'+fn,'READ')
		tree = ROOT_file.Get(treeName)
		evtCount = 0
		for event in tree:
			if preSelect:
				if not event.IsTriggered:
					continue
				if not event.IsInTarget:
					continue
				if not event.IsExclusive:
					continue
				if event.CentralProdVeto:
					continue
				if not event.IsPlanar_extended:
					continue
				if not event.CorrectNbrRPDTracks:
					continue
			evtData = [	event.run,
					event.E_p0,
					event.px0,
					event.py0,
					event.pz0,
					event.px1,
					event.py1,
					event.pz1,
					event.px2,
					event.py2,
					event.pz2,
					event.px3,
					event.py3,
					event.pz3,
					event.xVtx,
					event.yVtx,
					event.zVtx,
					event.massX,
					event.IsTriggered,
					event.IsInTarget,
					event.IsExclusive,
					event.IsInT,
					event.CentralProdVeto,
					event.IsInBeamTime,
					event.RICH_Veto,
					event.CEDAR_Veto,
					event.IsInDeltaPhi,
					event.CorrectNbrRPDTracks,
					event.IsPlanar,
					event.IsPlanar_extended]
			returnEvents.append(evtData)
			if evtCount >= evtPerFile:
				break
			evtCount+=1
	return returnEvents

def normToBinWidth(hist):
	for i in range(hist.GetNbinsX()):
		for j in range(hist.GetNbinsY()):
			widthX = hist.GetXaxis().GetBinWidth(i+1)
			widthY = hist.GetYaxis().GetBinWidth(j+1)
			hist.SetBinContent(i+1,j+1,hist.GetBinContent(i+1,j+1)/widthX/widthY)

if __name__ == "__main__":
	canv = Canvas(name = "canv0", title = "PWA",height=1000,width=1000)
	odhist = compare_diag_offdiag('1-(1++)0+ (pipi)_S','1-(1++)0+ rho pi S','1-(1++)0+ f0_','1-(1++)0+ rho_',1.50,1.54,'/nfs/mds/user/fkrinner/massIndepententFits/fits/tst/integrals/0.10000-0.14077')
	prhist = TH2D(odhist)
	for i in range(prhist.GetNbinsX()):
		for j in range(prhist.GetNbinsY()):
			prhist.SetBinContent(i+1,j+1,0.)
			prhist.SetBinError(i+1,j+1,0.)	
	events = get_events(1.50,1.54,n=100000,preSelect = True)
	mpis = []
	for event in events:
		p1=[event[5],event[6],event[7]]
		p2=[event[8],event[9],event[10]]
		p3=[event[11],event[12],event[13]]

		p1 = four_vector([(p1[0]**2+p1[1]**2+p1[2]**2+m_Pi**2)**.5,p1[0],p1[1],p1[2]])
		p2 = four_vector([(p2[0]**2+p2[1]**2+p2[2]**2+m_Pi**2)**.5,p2[0],p2[1],p2[2]])
		p3 = four_vector([(p3[0]**2+p3[1]**2+p3[2]**2+m_Pi**2)**.5,p3[0],p3[1],p3[2]])

		m13 = (p1+p3).mass()
		m23 = (p2+p3).mass()

		prhist.Fill(m13,m23)
	normToBinWidth(prhist)
	
	prhist.Draw()
	canv.Update()
	raw_input()
	odhist.Draw()
	canv.Update()
	raw_input()


	for i in range(prhist.GetNbinsX()):
		for j in range(prhist.GetNbinsY()):
			if prhist.GetBinContent(i+1,j+1) != 0.:
				odhist.SetBinContent(i+1,j+1,odhist.GetBinContent(i+1,j+1)/prhist.GetBinContent(i+1,j+1))
			else:
				odhist.SetBinContent(i+1,j+1,0.)

	odhist.Draw()		








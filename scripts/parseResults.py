import csv
import pprint
import numpy as np
from ROOT import TGraph
from ROOT import TGraphErrors
from ROOT import TCanvas
from array import array
from ROOT import TH1D
from ROOT import gROOT
from ROOT import TFile
import ROOT

c0 = TCanvas("c0","c0",0,0,600,400)
c0.Print("results.pdf(","pdf")
hRes = TH1D("hRes","hRes",500,0,50)
gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")

elossGraph = (TFile.Open("../info/eloss.root")).Get("elossVsBeamE");
elossGraph.Draw("ALP")
c0.Print("elossCurve.pdf","pdf")


scan2_3a = { 
	"start" : 10397,
	"stop"  : 10406}


scan2_3b = { 
	"start" : 10437,
	"stop"  : 10444}
	

efficiency = {
	0 : 2.04e-4,
	1 : 2.25e-4,
	2 : 1.85e-4,
	3 : 1.16e-4,
	4 : 1.63e-4,		
	5 : 1.75e-4,
	6 : 1.65e-4,
	7 : 1.61e-4,
	8 : 1.48e-4,	
	9 : 1.84e-4,
	10: 9.64e-5,
	11: 2.13e-4
}	
name = {
	0 : "0_0",
 	1 : "0_1",
	2 : "0_2",
	3 : "0_3",
	4 : "0_4",
	5 : "0_5",
	6 : "0_6",
	7 : "0_7",
	8 : "1_0",
	9 : "1_1",
	10 : "1_4",
	11 : "1_5"
}

angle = {
	0 : 150,
	1 : 135,
	2 : 105,
	3 : 90,
	4 : 75,
	5 : 60,
	6 : 45,
	7 : 30,
	8 : 15,
	9 : 0,
	10 : 90,#270
	11 : 120
}


with open('../info/17Oang_FixedRunSheet.dat','r') as csvfile:
	runInfo = csv.DictReader(csvfile,delimiter='\t')
	infoList = list(runInfo)
	#print(readerList.index("10154"))


for det in [9,8,7,6,5,4,3,10,2,11,1,0]:
	with open('../output/n1Results{}.txt'.format(det),'r') as result:
		rNum = array( 'd') 
		y    = array( 'd')
		yE   = array( 'd')
		gain = array( 'd')
		offset = array( 'd')
		resolution = array( 'd')
		yPerQ = array( 'd')
		yEPerQ = array( 'd')
		beamE = array( 'd')
		labAngle = array( 'd')
		
		scan2_3aY = array( 'd')
		scan2_3aYE = array( 'd')
		scan2_3aEb = array( 'd')
		scan2_3bY = array( 'd')
		scan2_3bYE = array( 'd')
		scan2_3bEb = array( 'd')
		
		
		for row in result:
			d = row.split()
			#grab info from runInfo

			for entry in infoList:
#				print(entry['Run'],d[0])
				if entry['Run']==d[0]:
					break				
						
			if entry['Ea (keV)']=='':
				#if det==0:
				#	print(entry['Run'],"Empty Beam Energy")
				continue
			if entry['Q(target)']=='':
				#if det==0:
				#	print(entry['Run'],"No Q")
				continue	
			if float(entry['Ea (keV)'])>4500:
				continue	
			if entry['good?']=='#':
				#if det==0:
				#	print(entry['Run'],entry['Ea (keV)'],"Run Commented Out")
				continue
			if float(d[5])>25:
				print(d[0],name[det])
			if float(d[2])>=float(d[1]):
				continue		
			rNum.append(int(d[0]))
			y.append(float(d[1]))
			yE.append(float(d[2]))
			gain.append(float(d[3]))
			offset.append(float(d[4]))
			resolution.append(float(d[5]))
			if int(entry['Run'])<10681:
				yPerQ.append(float(d[1])/float(entry['Q(target)'])/efficiency[det])
				yEPerQ.append(float(d[2])/float(entry['Q(target)'])/efficiency[det])
			else:	
				yPerQ.append(2.0*float(d[1])/float(entry['Q(target)'])/efficiency[det])
				yEPerQ.append(2.0*float(d[2])/float(entry['Q(target)'])/efficiency[det])
			
			inBeamE = float(entry['Ea (keV)'])
			BeamEloss = elossGraph.Eval(inBeamE/1000)
			#print(BeamEloss)
			beamE.append(inBeamE/1000-BeamEloss)
			labAngle.append(angle[det])		
			hRes.Fill(float(d[5]))
		
			print(entry['Run'])
			if int(entry['Run'])>=int(scan2_3a['start']) and int(entry['Run'])<=int(scan2_3a['stop']):
				scan2_3aY.append(yPerQ[-1])
				scan2_3aYE.append(yEPerQ[-1])
				scan2_3aEb.append(beamE[-1])

				# these runs were omitted due to suprsession voltage being off. Go back and check this laters
			if int(entry['Run'])>=int(scan2_3b['start']) and int(entry['Run'])<=int(scan2_3b['stop']):
				scan2_3bY.append(yPerQ[-1])
				scan2_3bYE.append(yEPerQ[-1])
				scan2_3bEb.append(beamE[-1])


			

	c0.cd()
	c0.SetLogy(0)
	
	grRes = TGraph(len(rNum),rNum,resolution)
	grRes.SetTitle("Detector {} Resolution".format(name[det]))
	hRes.SetTitle("Detector {} Resolution".format(name[det]))
	grRes.Draw("AL*")
	c0.Update()
	c0.Print("results.pdf","pdf")	
	hRes.Draw()
	c0.Update()
	c0.Print("results.pdf","pdf")
	hRes.Reset()	
		
	grGain = TGraph(len(rNum),rNum,gain)
	grGain.SetTitle("Detector {} Gain".format(name[det]))
	grGain.Draw("AL*")
	c0.Update()
	c0.Print("results.pdf","pdf")		
		
	grOffset = TGraph(len(rNum),rNum,offset)
	grOffset.SetTitle("Detector {} Offset".format(name[det]))
	grOffset.Draw("AL*")
	c0.Update()	
	c0.Print("results.pdf","pdf")
	
	npBeamE = np.asarray(beamE)
	npYperQ = np.asarray(yPerQ)
	npYEPerQ    = np.asarray(yEPerQ)
	iE = np.argsort(npBeamE)
	npLabAngle = np.asarray(labAngle)
	
	grYield = TGraphErrors(len(npBeamE),npBeamE[iE],npYperQ[iE],np.zeros(len(npYEPerQ)),npYEPerQ[iE])
	if det==10:
		with open('2017o17n1_90b.dat','w') as az:
			writer = csv.writer(az,delimiter='\t')
			writer.writerows(np.column_stack((npBeamE[iE],npLabAngle[iE],npYperQ[iE],npYEPerQ[iE])))
	else:
		with open('2017o17n1.dat', 'a') as az:
			writer = csv.writer(az,delimiter='\t')
			writer.writerows(np.column_stack((npBeamE[iE],npLabAngle[iE],npYperQ[iE],npYEPerQ[iE])))
	grYield.SetTitle("Detector {} Yield {}".format(name[det],labAngle[det]))
	grYield.Draw("APL")
	grYield.GetXaxis().SetRangeUser(0,5);

	c0.SetLogy(1)
	c0.Update()
	c0.Print("results.pdf","pdf")
	
	
	grYield.GetXaxis().SetRangeUser(3.5,4.0);
	c0.SetLogy(0)
	
	
	#npScan2_3bEb = np.asarray
	#print(type(scan2_3aEb[0]))
	#print(type(scan2_3aY[0]))
	#print(type(scan2_3aYE[0]))
	#print(len(scan2_3aEb),len(scan2_3aY),len(scan2_3aYE))
	#print(type(np.zeros(len(scan2_3aEb))[0]))
	#print(type(rNum[0]),type(offset[0]))

#	grScan2_3a = TGraph(len(scan2_3aEb),scan2_3aEb,scan2_3aY)#,np.zeros(len(scan2_3aEb)),scan2_3aYE)
	#grScan2_3a = TGraph(len(scan2_3aEb),scan2_3aEb,scan2_3aY)
	
#	grScan2_3b = TGraph(len(scan2_3bEb),scan2_3bEb,scan2_3bY)#,np.zeros(len(scan2_3bEb)),scan2_3bYE)
	
#	grScan2_3a.SetLineColor(ROOT.kAzure)
#	grScan2_3a.Draw();
#	grScan2_3b.Draw("SAME")
#	c0.Update()
#	c0.Print("results.pdf","pdf")

		
c0.Print("results.pdf)","pdf")	




		
			


	



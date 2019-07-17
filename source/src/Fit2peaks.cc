//
// Fit2peaks.cc
// --------------
//
// Class containing the method to fit a single Gaussian peak on 
// top of the detector background. The background is used as an
// internal energy calibration. The energy of the additional 
// peak to be fit is fixed relative to the pre-determined energy
// calibration of the background spectra.
//
// K. T. Macon (macon.kevin@gmail.com)
//

#include "Fit2peaks.hh"
#include "TAxis.h"

#include <iomanip>

Fit2peaks::Fit2peaks(TFile* bgf) : Egamma{1778.97,2838.3}, Eedge{1555.56,2603.9}, backgroundFile(bgf), backgroundHist(NULL), dataFile(NULL), dataTree(NULL), spectrum(NULL), subtracted(NULL), calibrated(NULL), fitFunction(NULL), bgFunction(NULL), reject(false), rejectMore(false), gain(0.0), offset(0.0), peakArea{0.0,0.0}, peakUnc{0.0,0.0}
{
	fitFunction = new TF1("fitFunction",this,0,16384,11);
	bgFunction  = new TF1("bgFunction",this,0,16384,11);
	bgFunction->SetLineColor(1);
	bgFunction->SetNpx(1e5);
	fitFunction->SetNpx(1e5);
	fitFunction->SetParNames("N_{1}","#sigma_{1}","offset","gain","quad","scale","H_{1}","breadth",
		"N_{2}","#sigma_{2}","H_{2}");
	bgFunction->SetParNames("N_{1}","#sigma_{1}","offset","gain","quad","scale","H_{1}","breadth",
		"N_{2}","#sigma_{2}","H_{2}");
	spectrum = new TH1D("htemp","^{25}Mg(#alpha,n_{1,2})",16000,0,16000);
	subtracted = new TH1D("htemp2","Background Subtracted ^{25}Mg(#alpha,n_{1})",16384,0,16384);
	calibrated = new TH1D("htemp3","Background Subtracted ^{25}Mg(#alpha,n_{1})",16384,0,16384);
}

Fit2peaks::~Fit2peaks()
{
	delete fitFunction;
	delete bgFunction;
	delete spectrum;
	delete subtracted;
}

void Fit2peaks::SetEgamma(int id,double eg)
{
	Egamma[id] = eg;
	Eedge[id] = eg - eg/(1.0+2.0*eg/511);
	fitFunction->SetParName((id==0)?0:8,Form("N_{%.0f}",Egamma[id]));
	fitFunction->SetParName((id==0)?1:9,Form("#sigma_{%.0f}",Egamma[id]));
}

void Fit2peaks::Fit()
{
 	spectrum->Fit(fitFunction,"R");
 	peakArea[0] = fitFunction->GetParameter(0);
 	sigma = fitFunction->GetParameter(1);
 	offset = fitFunction->GetParameter(2);
 	gain = fitFunction->GetParameter(3);
}


void Fit2peaks::Fit(TString opt1)
{
	spectrum->Fit(fitFunction,opt1);
	peakArea[0] = fitFunction->GetParameter(0);
 	sigma = fitFunction->GetParameter(1);
 	offset = fitFunction->GetParameter(2);
 	gain = fitFunction->GetParameter(3);
}

void Fit2peaks::Fit(TString opt1,double x1,double x2)
{
	spectrum->Fit(fitFunction,opt1,"",x1,x2);
	peakArea[0] = fitFunction->GetParameter(0);
 	sigma = fitFunction->GetParameter(1);
 	offset = fitFunction->GetParameter(2);
 	gain = fitFunction->GetParameter(3);
}

void Fit2peaks::Rebin(int factor)
{
	spectrum->Rebin(factor);
	subtracted->Rebin(factor);
	calibrated->Rebin(factor);
}

void Fit2peaks::LoadBackground(int bnum,int cnum)
{
	LoadBackground(Form("bgcal%i_%i",bnum,cnum));
}

void Fit2peaks::LoadBackground(TString hname)
{
	backgroundHist = static_cast<TH1D*>(backgroundFile->Get(hname));
	if( backgroundHist ==NULL)
	{
		cerr << "Background Histogram "+hname+" Not Found!" << endl;
	}
	return;
}

void Fit2peaks::CreateSpectrum(TFile *dFile, int bnum,int cnum,TString spectrumName)
{
	dataFile = dFile;
	dataTree = static_cast<TTree*>(dataFile->Get("Data"));
	spectrum->SetName(spectrumName);
	dataTree->Project(spectrumName,"lgate","lgate>0");	
	return;
}


void Fit2peaks::CreateSpectrum(TTree *dTree, int bnum,int cnum,TString spectrumName)
{
	dataTree = dTree;
	spectrum->SetName(spectrumName);
	dataTree->Project(spectrumName,"Energy",Form("Energy>0 && Board==%d && Channel==%d",bnum,cnum));	
	return;
}


void Fit2peaks::Subtract()
{
	bgFunction->SetParameters(fitFunction->GetParameters());
	bgFunction->SetParameter(0,0.0);
	bgFunction->SetParameter(8,0.0);
	bool tempReject = reject;
	reject = false;
	if( spectrum == NULL || backgroundHist == NULL)
	{
		cerr << "Background or Data spectrum missing" << endl;
	}
	subtracted->SetName(static_cast<TString>(spectrum->GetName())+"_subtracted");
	calibrated->SetName(static_cast<TString>(spectrum->GetName())+"_calibrated");
	calibrated->GetXaxis()->SetLimits(0.0*gain+offset,16384*gain+offset);
	
	
	for(int i = 0; i<spectrum->GetNbinsX(); i++)
	{
		int counts = spectrum->GetBinContent(i);
		double background = bgFunction->Eval(spectrum->GetBinCenter(i));

		
		// DOUBLE CHECK ERROR CALCULATION
		subtracted->SetBinContent(i,counts-background);
		subtracted->SetBinError(i,TMath::Sqrt(counts+background));
		calibrated->SetBinContent(i,counts-background);
		calibrated->SetBinError(i,TMath::Sqrt(counts+background));
	}
	reject = tempReject;
	int il = calibrated->GetXaxis()->FindBin(Egamma[0]-4.0*sigma);
	double delta = 4.0*sigma/calibrated->GetBinWidth(7);
	int ir = il+delta+delta+delta+delta-1;
	//int ill = il-delta;
	//int irr = ir+delta+delta+delta;
	peakArea[0] = calibrated->IntegralAndError(il,ir,peakUnc[0]);  //[ill, il-1 = ill+delta-1] [ il , ir=il+4delta-1 ] [ir+1,irr=ir+delta]
	
	return;
}

void Fit2peaks::SetParameters(double p0,double p1, double p2, double p3, double p4,
	double p5,double p6,double p7, double p8, double p9, double p10)
{
	fitFunction->SetParameters(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);
	SetLimits();
	return;
}

void Fit2peaks::SetParameters(double* p)
{
	fitFunction->SetParameters(p);
	SetLimits();
	return;
}

void Fit2peaks::SetLimits()
{
	fitFunction->SetParLimits(0,0,1000000);
	fitFunction->SetParLimits(1,1,25);
	fitFunction->SetParLimits(3,.5*fitFunction->GetParameter(3)
		,1.5*fitFunction->GetParameter(3));
	fitFunction->SetParLimits(4,0,0.01);
	//fitFunction->SetParLimits(7,2,500);
	fitFunction->SetParLimits(8,0,1000000);
	fitFunction->SetParLimits(9,1,25);
	fitFunction->SetParLimits(10,0.0001,1000);
	return;
}

double Fit2peaks::operator() (double *x, double *par)
{
	double norm1   = par[0]; // Gaussian normalization (peak area)
	double sigma1  = par[1]; // Energy resolution parameter
	double offset  = par[2]; // spectra energy offset in keV 
	double gain    = par[3]; // spectra gain in keV per channel
	double quad		= par[4]; // quadratic term in the energy calibration
	double scale   = par[5]; // relative amplitude of background (ratio of livetimes)
	double ratio1  = par[6]; // height of multiple scattering step relative to peak area
	//double broaden = par[7]; // ``resolution'' of multiple scattering events should be worse, right?
	double norm2   = par[8]; // second peak
	double constBG = par[9]; // second energy resolution parameter
	double ratio2  = par[10];
	
	double sigma2   = TMath::Sqrt(Egamma[1]/Egamma[0])*sigma1;
	double sigmaCE1 = TMath::Sqrt(Eedge[0]/Egamma[0])*sigma1;
	double sigmaCE2 = TMath::Sqrt(Eedge[1]/Egamma[0])*sigma1;
	
	double xch1 = (-gain + TMath::Sqrt(gain*gain-4*quad*(offset-Egamma[0])))/2/quad;
	double xch2 = (-gain + TMath::Sqrt(gain*gain-4*quad*(offset-Egamma[1])))/2/quad;
	double xchCE1 = (-gain + TMath::Sqrt(gain*gain-4*quad*(offset-Eedge[0])))/2/quad;
	double xchCE2 = (-gain + TMath::Sqrt(gain*gain-4*quad*(offset-Eedge[1])))/2/quad;
	double result = norm1*TMath::Gaus(x[0],xch1,sigma1/gain,true);
	result += 0.5*norm1*ratio1*TMath::Erfc((x[0]-xchCE1)/TMath::Sqrt(2)/sigmaCE1/gain);
	result += norm2*TMath::Gaus(x[0],xch2,sigma2/gain,true);
	result += 0.5*norm2*ratio2*TMath::Erfc((x[0]-xchCE2)/TMath::Sqrt(2)/sigmaCE2/gain);


	//scaled background histogram
	// variable bin widths here for a non-linear energy calibration
	result += scale*(backgroundHist->Interpolate(x[0]*x[0]*quad+x[0]*gain+offset))*(gain + 2.0*quad*x[0]);

	
	bool roi[2];
	roi[0] = (Eedge[0]<(x[0]+4*sigma1)*(x[0]+4*sigma1)*quad+(x[0]+4*sigma1)*gain+offset) 
		&& (Egamma[0]>(x[0]-3.0*sigma1)*(x[0]-3.0*sigma1)*quad+(x[0]-3.0*sigma1)*gain+offset);
	roi[1] = (Eedge[1]<(x[0]+1.0*sigma2)*(x[0]+1.0*sigma2)*quad+(x[0]+1.0*sigma2)*gain+offset) 
		&& (Egamma[1]>(x[0]-3.0*sigma2)*(x[0]-3.0*sigma2)*quad+(x[0]-3.0*sigma2)*gain+offset);
	
	if(reject && !(roi[0] || roi[1]))
	{
		TF1::RejectPoint(true);
		return 0;
	}
	else if(rejectMore && (x[0]*gain+offset<1575 || x[0]*gain+offset>2700))
	{//Custom fit ROI with Min/Max given an input calibration to keV.
		TF1::RejectPoint(true);
		return 0;
	}
	
	return result;
	
}







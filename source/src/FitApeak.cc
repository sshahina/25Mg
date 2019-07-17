//
// FitApeak.cc
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

#include "FitApeak.hh"
#include "TAxis.h"

#include <iomanip>

FitApeak::FitApeak(TFile* bgf) : Egamma(1778.97), Eedge(1555.56), backgroundFile(bgf), backgroundHist(NULL), dataFile(NULL), dataTree(NULL), spectrum(NULL), subtracted(NULL), calibrated(NULL), fitFunction(NULL), bgFunction(NULL), reject(false), rejectMore(false), gain(0.0), offset(0.0), peakArea(0.0), peakUnc(0.0)
{
	fitFunction = new TF1("fitFunction",this,0,16384,8);
	bgFunction  = new TF1("bgFunction",this,0,16384,8);
	bgFunction->SetLineColor(1);
	bgFunction->SetNpx(1e5);
	fitFunction->SetNpx(1e5);
	fitFunction->SetParNames("Norm","#sigma","offset","gain","scale","adjustment","height","breadth");
	bgFunction->SetParNames("Norm","#sigma","offset","gain","scale","adjustment","height","breadth");
	spectrum = new TH1D("htemp","^{17}O(#alpha,n_{1})",16000,0,16000);
	subtracted = new TH1D("htemp2","Background Subtracted ^{17}O(#alpha,n_{1})",16384,0,16384);
	calibrated = new TH1D("htemp3","Background Subtracted ^{17}O(#alpha,n_{1})",16384,0,16384);
}

FitApeak::~FitApeak()
{
	delete fitFunction;
	delete bgFunction;
	delete spectrum;
	delete subtracted;
}

void FitApeak::SetEgamma(double eg)
{
	Egamma = eg;
	Eedge = Egamma - Egamma/(1.0+2.0*Egamma/511);
}

void FitApeak::Fit()
{
 	spectrum->Fit(fitFunction,"R");
 	peakArea = fitFunction->GetParameter(0);
 	sigma = fitFunction->GetParameter(1);
 	offset = fitFunction->GetParameter(2);
 	gain = fitFunction->GetParameter(3);
}


void FitApeak::Fit(TString opt1)
{
	spectrum->Fit(fitFunction,opt1);
	peakArea = fitFunction->GetParameter(0);
 	sigma = fitFunction->GetParameter(1);
 	offset = fitFunction->GetParameter(2);
 	gain = fitFunction->GetParameter(3);
}

void FitApeak::Fit(TString opt1,double x1,double x2)
{
	spectrum->Fit(fitFunction,opt1,"",x1,x2);
	peakArea = fitFunction->GetParameter(0);
 	sigma = fitFunction->GetParameter(1);
 	offset = fitFunction->GetParameter(2);
 	gain = fitFunction->GetParameter(3);
}

void FitApeak::Rebin(int factor)
{
	spectrum->Rebin(factor);
	subtracted->Rebin(factor);
	calibrated->Rebin(factor);
}

void FitApeak::LoadBackground(int bnum,int cnum)
{
	LoadBackground(Form("bgcal%i_%i",bnum,cnum));
}

void FitApeak::LoadBackground(TString hname)
{
	backgroundHist = static_cast<TH1D*>(backgroundFile->Get(hname));
	if( backgroundHist ==NULL)
	{
		cerr << "Background Histogram "+hname+" Not Found!" << endl;
	}
	return;
}

void FitApeak::CreateSpectrum(TFile *dFile, int bnum,int cnum,TString spectrumName)
{
	dataFile = dFile;
	dataTree = static_cast<TTree*>(dataFile->Get("Data"));
	spectrum->SetName(spectrumName);
	dataTree->Project(spectrumName,"lgate","lgate>0");	
	return;
}


void FitApeak::CreateSpectrum(TTree *dTree, int bnum,int cnum,TString spectrumName)
{
	dataTree = dTree;
	spectrum->SetName(spectrumName);
	dataTree->Project(spectrumName,"Energy",Form("Energy>0 && Board==%d && Channel==%d",bnum,cnum));	
	return;
}


void FitApeak::Subtract()
{
	bgFunction->SetParameters(fitFunction->GetParameters());
	bgFunction->SetParameter(0,0);
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
	int il = calibrated->GetXaxis()->FindBin(Egamma-4.0*sigma);
	double delta = 4.0*sigma/calibrated->GetBinWidth(7);
	int ir = il+delta+delta+delta+delta-1;
	//int ill = il-delta;
	//int irr = ir+delta+delta+delta;
	peakArea = calibrated->IntegralAndError(il,ir,peakUnc);  //[ill, il-1 = ill+delta-1] [ il , ir=il+4delta-1 ] [ir+1,irr=ir+delta]
	
	return;
}

void FitApeak::SetParameters(double p0,double p1, double p2, double p3, double p4,double p5,double p6,double p7)
{
	fitFunction->SetParameters(p0,p1,p2,p3,p4,p5,p6,p7);
	SetLimits();
	return;
}

void FitApeak::SetParameters(double* p)
{
	fitFunction->SetParameters(p);
	SetLimits();
	return;
}

void FitApeak::SetLimits()
{
	fitFunction->SetParLimits(0,0,1000000);
	fitFunction->SetParLimits(1,1,25);
	fitFunction->SetParLimits(3,.5*fitFunction->GetParameter(3)
		,1.5*fitFunction->GetParameter(3));
	fitFunction->SetParLimits(4,0.0001,1000);
	fitFunction->SetParLimits(5,.998,1.002);
	fitFunction->SetParLimits(7,2,500);
	return;
}

double FitApeak::operator() (double *x, double *par)
{
	double norm   = par[0]; // Gaussian peak normalizatio (peak area)
	double sigma  = par[1]; // Energy resolution parameter
	double offset = par[2]; // spectra energy offset in keV 
	double gain   = par[3]; // spectra gain in keV per channel
	double scale  = par[4]; // relative amplitude of background (ratio of livetimes)
	double adj    = par[5];
	double ratio  = par[6]; // height of multiple scattering step relative to peak area
	double broaden = par[7]; // ``resolution'' of multiple scattering events should be worse, right?
	
	double result = norm*TMath::Gaus(x[0],(Egamma-offset)/gain*adj,sigma/gain,true);
	result += scale*(backgroundHist->Interpolate(x[0]*gain+offset));
	result += 0.5*norm*ratio*TMath::Erfc((x[0]-(Egamma-offset)/gain)/TMath::Sqrt(2)/sigma/gain/broaden);
	if(reject && ((x[0]-sigma)*gain+offset<Eedge || x[0]*gain+offset>2700))
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







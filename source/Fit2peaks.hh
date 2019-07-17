//
// Fit2peaks.hh
// --------------
//
// Class containing the method to fit a two Gaussian peaks on 
// top of the detector background + flat compton continuum background. 
// The background is used as an internal energy calibration. The energy of 
// the additional peak to be fit is fixed relative to the pre-determined 
// energy calibration of the background spectra.
//
// K. T. Macon (macon.kevin@gmail.com)
//

#ifndef Fit2peaks_hh
#define Fit2peaks_hh

//ROOT Stuff
#include <TROOT.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>


#include <iostream>
using std::cout;
using std::ends;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <string>
using std::string;

class Fit2peaks
{
	double Egamma[2];
	double Eedge[2];
//	double mean;// Peak location determined by calibration of background histogram
// which is then fit to the spectrum under analysis. Hopefully the quenching factor is robust for the alpha background peaks!!
	TFile *backgroundFile;  // Background histogram calibrated in units of keV.
	TH1D  *backgroundHist;
	TFile *dataFile;
	TTree *dataTree;
	TH1D  *spectrum;
	TH1D	*subtracted;
	TH1D	*calibrated;
	TF1   *fitFunction; // 
	TF1	*bgFunction;
	bool reject = false;
	bool rejectMore = false;
	double gain;
	double offset;
	double sigma;
	double peakArea[2];
	double peakUnc[2];
	public:
		Fit2peaks(TFile*);
		~Fit2peaks();
		void LoadBackground(TString);
		void LoadBackground(int,int);
		void LoadCalibration(TString);
		void SaveCalibration(TString);
		void CreateSpectrum(TFile*,int,int,TString);
		void CreateSpectrum(TTree*,int,int,TString);
		void SetParameters(double*);
		void SetParameters(double,double,double,double,double,double,double,double,double,double,double);
		void SetParameter(int id, double val) {fitFunction->SetParameter(id,val);}
		void FixParameter(int id, double val) {fitFunction->FixParameter(id,val);}
		void ReleaseParameter(int id) 		  {fitFunction->ReleaseParameter(id);}
		void FixLiveTimeRatio(double val)     {fitFunction->FixParameter(5,val);}
		void SetLimits();
		void SetEgamma(int,double);
		void Reject(bool flag) {reject = flag;}
		void RejectMore(bool flag) {rejectMore = flag;}
		void Fit();
		void Fit(TString opt1);
		void Fit(TString opt1,double x1,double x2);
		void Subtract();
		void Rebin(int);
		double operator() (double*,double*);
		double EnergyToChannel(double,double*);
		double GetPeakArea(int id)   {return peakArea[id];}
		double GetPeakUnc(int id)    {return peakUnc[id];}
		double Gain()          {return gain;}
		double Offset()        {return offset;}
		TH1D *GetBackground()  {return backgroundHist;}
		TH1D *GetSpectrum()    {return spectrum;}
		TH1D *GetSubtracted()  {return subtracted;}
		TH1D *GetCalibrated()  {return calibrated;}
		TF1  *GetFitFunction() {return fitFunction;}
		TF1  *GetBgFunction()  {return bgFunction;}
		double GetParameter(int pnum)  {return fitFunction->GetParameter(pnum);}
};

#endif // Fit2peaks_hh


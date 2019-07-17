#include <iostream>
using std::cout;
using std::endl;

#include<unistd.h>
#include <string>
using std::string;
//#include<bits/stdc++.h>

#include <vector>
using std::vector;

#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TError.h"

// Fitting routines
#include "calibration/fitFunctions.h"

// Gain match routine
#include "calibration/gainMatch.h"

// Silicon line from first run, these get updated run by run with the new ranges
double sLow[] = {430,558,413,524,551,381,480,565,495,331,356,430,482,500};
double sHigh[] = {452,579,440,550,580,400,510,591,515,356,376,458,510,530};

// global array for gain match constants
double _a_gain[13];
double _b_gain[13];

// global array for calibration constants
double _a_calibrator[13];
double _b_calibrator[13];
void peakFitter(const char *fileName, const char *fileBack, const char *detector,
	const vector < double > &peak1BackLow, const vector < double > &peak1BackHigh,
	const vector < double > &peak2BackLow, const vector < double > &peak2BackHigh,
	vector < double > &peak1SpecLow, vector < double > &peak1SpecHigh,
	vector < double > &peak2SpecLow, vector < double > &peak2SpecHigh,
	double low, double high, int detLoop){


	TFile *fyield = new TFile(fileName);
	TFile *fbackground = new TFile(fileBack);	// bg spectra

	// Get runtime info from root file
	TVectorD *run_t0 = static_cast<TVectorD*>(fyield->Get("lastTimeStampSeconds-0"));
	TVectorD *back_t0 = static_cast<TVectorD*>(fbackground->Get("lastTimeStampSeconds-0"));
	//run_t0->Print();
	//back_t0->Print();
	// Run time and background time for each board
	double runTime0 = (*run_t0)[0];
	double backTime0 = (*back_t0)[0];

	// Scale background spectra to ratio of run time
	double scale0 = runTime0/backTime0;

	// Get histograms from root file
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	TH1D *HBACK = static_cast<TH1D*>(fbackground->Get(detector));

	// Scale the background
	HBACK->Scale(scale0);
	HBACK->SetDirectory(0);

	// Get integrated charge information
	TH1D *qcharge = static_cast<TH1D*>(fyield->Get("h_0_7"));
	double charge = qcharge->GetEntries();
	gStyle->SetOptFit(1111);

	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",600,400);
	c0->Divide(1,2);
	c0->Update();

	c0->cd(1);

	// Prepare bg subtracted histogram and other histograms
	TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));
	TH1D *h2 = new TH1D(Form("h2 - det-%i",detLoop),Form("h2 - det-%i",detLoop),10000,0,2000);
	TH1D *h3 = new TH1D(Form("h3 - det-%i",detLoop),Form("h3 - ECAL - det-%i",detLoop),10000,0,2000);

	double area;
	double area_err;
	double chi2NDF;
	double sig1;
	double sig2;

	double a = 0;
	double b = 0;
	double a_gm = 0;
	double b_gm = 0;
	double a_cal = 0;
	double b_cal = 0;

	double linear = 0;
	double offset = 0;

	vector < double > p1Peak;


	// Only gain match longer runs where background peaks appear
	/*if(runTime0 > 850){

		/////////////////////////
		///    GAIN MATCH     ///
		/////////////////////////

		// Find BG peak positions
		vector < double > peak1Back;
		vector < double > peak2Back;

		peak1Back = iterative_single_gauss_peak(peak1BackLow[detLoop],peak1BackHigh[detLoop],HBACK);
		peak2Back = iterative_double_gauss_peak_same_width(peak2BackLow[detLoop],peak2BackHigh[detLoop],HBACK);

		// Find the BG peak positions in runs*/
		vector < double > peak1Position;
		vector < double > peak2Position;

		peak1Position = iterative_single_gauss_peak(peak1SpecLow[detLoop],peak1SpecHigh[detLoop],hyield);
		peak2Position = iterative_double_gauss_peak_same_width(peak2SpecLow[detLoop],peak2SpecHigh[detLoop],hyield);

		/*// Add new ranges for each detector
		peak1SpecLow[detLoop] = peak1Position[1];
		peak1SpecHigh[detLoop] = peak1Position[2];
		peak2SpecLow[detLoop] = peak2Position[1];
		peak2SpecHigh[detLoop] = peak2Position[2];

		// New gain matched BG histogram -> h2
		vector < double > gain;
		gain = gain_match(peak1Position[0],peak2Position[0],peak1Back[0],peak2Back[0],h2,HBACK);
		a = gain[0];
		b = gain[1];
		a_gm = gain[0];
		b_gm = gain[1];

		_a_gain[detLoop] = a_gm;
		_b_gain[detLoop] = b_gm;*/

		// Draw results
		hyield->Draw();
		hyield->GetXaxis()->SetRangeUser(0,2000);
		hyield->SetStats(kFALSE);


		// Recalculate errors manually
		for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
			double yval = ysubtracted->GetBinContent(i);
			double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
			double yerr = ysubtracted->GetBinError(i);
			double yerr2 = h2->GetBinError(i);	// Takes the error from the gain matched bg spectrum

			ysubtracted->SetBinContent(i,yval-yval2);
			ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
		}

		ysubtracted->SetLineColor(kRed);
		ysubtracted->Draw("SAME");

		c0->cd(2);

		// Now I've gainmatched and subtracted the BG -> ysubtracted, find flourine line position
		vector < double > sPosition;
		sPosition = single_gauss_peak(sLow[detLoop],sHigh[detLoop],ysubtracted);

		// Found the lines update their ranges run by run
		sLow[detLoop] = sPosition[1];
		sLow[detLoop] = sPosition[2];

		vector < double > calibrators;

		calibrators = calibrate(1468,1778.969,peak1Position[0],sPosition[0],h3,ysubtracted);
		a_cal = calibrators[0];
		b_cal = calibrators[1];

		// Use these calibrators for when the runtime is too short and no gainmatch fitting will be performed
		_a_calibrator[detLoop] = a_cal;
		_b_calibrator[detLoop] = b_cal;
	//}
	/*else{
		 cout << Form("Runtime too short: %f seconds. Not performing gain match",runTime0);
		a_cal = _a_calibrator[detLoop];
		b_cal = _b_calibrator[detLoop];
		a_gm = _a_gain[detLoop];
		b_gm = _b_gain[detLoop];

		h2 = gain_match2(a_gm,b_gm,h2,HBACK);

		// Draw results
		hyield->Draw();
		hyield->GetXaxis()->SetRangeUser(0,2000);
		hyield->SetStats(kFALSE);

		// Recalculate errors manually*/
		for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
			double yval = ysubtracted->GetBinContent(i);
			double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
			double yerr = ysubtracted->GetBinError(i);
			double yerr2 = h2->GetBinError(i);	// Takes the error from the gain matched bg spectrum

			ysubtracted->SetBinContent(i,yval-yval2);
			ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
		}

		ysubtracted->SetLineColor(kRed);
		ysubtracted->Draw("SAME");

		c0->cd(2);

		h3 = calibrate2(a_cal,b_cal,h3,ysubtracted);
	//}


	h3->GetXaxis()->SetRangeUser(0,2000);
	h3->Draw();
	h3->SetStats(kFALSE);

	p1Peak = single_gauss_area(1778.969,h3); // 843.76

	area = p1Peak[0];
	area_err = p1Peak[1];
	chi2NDF = p1Peak[2];
	sig1 = p1Peak[3];
	linear = p1Peak[4];
	offset = p1Peak[5];

	// The charge is integrated charge of alpha which is 2+
	double yield = area/(charge);
	double yield_err = area_err/(charge);

	double goodFit;
	if (chi2NDF <= 1.4 && chi2NDF >=.6 ) goodFit = 0;
	else goodFit = 1;

	string runNum = fileName;
	//if (loc==1) runNum = runNum.substr(9,3);
	 runNum = runNum.substr(78,3);

	string detNum = detector;

	c0->SaveAs(Form("mkdir -p Yields/P1/run324/det_%s_Fit.png",detNum.c_str()));
	c0->SaveAs(Form("mkdir -p Yields/P1/det-%i/run324_Fit.png",detLoop));


	ofstream myfile;
	myfile.open (" mkdir -p Yields/P1/_P1.csv",std::ios::app);
	myfile<<Form("run%s",runNum.c_str())<<","<< Form("det_%s",detNum.c_str())<<","<<
				yield<<","<<yield_err<<","<<area<<","<<area_err<<","<<runTime0<<","<<
				goodFit<<","<<a<<","<<b<<","<<sig1<<","<<chi2NDF<<","<<linear<<","<<
				offset<<","<<charge<<"\n";
	myfile.close();


	c0->Clear();
	fyield->Close();
	fbackground->Close();
	delete c0;

	gROOT->Reset();
}

/*==============================MAIN=========================================*/
void mg25Yields_new(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	const char *path = "/afs/crc.nd.edu/group/nsl/ast/projects/25Mgang/n1analysis/scripts_new";
	chdir(path);

	// Background Spectrum file
	const char *fileBackground = "/afs/crc.nd.edu/group/nsl/ast/projects/25Mgang/n1analysis/scripts_new/runs/run244.root";

	const char *detect;
	const char *files;
	const char *outFile;

	// Prepare structure of data output in CSV file
	ofstream myfile;
	myfile.open ("mkdir -p Yields/P1/_P1.csv",std::ios::out);
	myfile<<"Run"<<","<<"Detector"<<","<<"Yield"<<","<<"Yield err"<<","<<"Area"<<","
			<<"Area err"<<","<<"Time"<<","<<"Fit Status"<<","<<"a"<<","<<"b"<<","
			<<"sig1"<<","<<"X2NDF"<<","<<"Linear"<<","<<"Offset"<<","<<"Q_int"<<"\n";
	myfile.close();

	// BG spectra peak ranges
	const vector < double >  peak1BackLow {193,248,185,239,247,171,219,255,218,150,160,194,217,226};
	const vector < double >  peak1BackHigh {207,268,199,254,263,185,229,275,235,161,173,205,231,237};
	const vector < double >  peak2BackLow {350,450,331,410,450,306,395,460,400,270,290,350,390,405};
	const vector < double >  peak2BackHigh {380,500,372,460,495,340,436,500,440,295,316,385,426,440};

	// BG peak ranges starting from run0159 \
		These get updated as runs progress
	vector < double >  peak1SpecLow {191,251,189,231,246,173,220,253,220,148,159,194,218,224};
	vector < double >  peak1SpecHigh {207,268,199,252,266,186,237,279,238,163,171,206,233,237};
	vector < double >  peak2SpecLow {350,454,336,420,451,310,390,460,400,272,288,350,391,406};
	vector < double >  peak2SpecHigh {377,485,369,460,484,332,420,490,430,297,315,380,426,437};


	// Make directory to visually inspect the fits
	for(int ii = 0; ii<13; ii++){
		try {
			gSystem->Exec(Form("mkdir -p Yields/P1/det-%i",ii));
		}catch(...){}
	}


	// Loop through runs: 159-410
	cout << "\nBEGINNING PEAK FITTING:" << endl;
	//int fileNum = 1;

	int upToRun;
	//if (loc==1) upToRun = 175;
	/*else*/ upToRun = 434;
	for(int i=324;i<upToRun;i++){

		// Skip bad runs
		if(i==331) continue;
		else if(i==334) continue;
		else if(i==337) continue;
		else if(i==339) continue;
		else if(i==341) continue;
		else if(i==343) continue;
		else if(i==345) continue;
		else if(i==347) continue;
		else if(i==349) continue;
		else if(i==351) continue;
		else if(i==353) continue;
		else if(i==356) continue;
		else if(i==358) continue;
		else if(i==360) continue;
		else if(i==362) continue;
		else if(i==364) continue;
		else if(i==367) continue;
		else if(i==369) continue;
		else if(i==374) continue;
		else if(i==376) continue;
		else if(i==378) continue;
		else if(i==380) continue;
		else if(i==382) continue;
		else if(i==384) continue;
		else if(i==387) continue;
		else if(i==403) continue;
		else if(i==404) continue;
		else if(i==412) continue;
		else if(i==417) continue;
		else if(i==430) continue;
		else if(i==434) continue;
		try {
			gSystem->Exec(Form("mkdir -p Yields/P1/run%d",i));
		}catch(...){}

		double p1;

		// Loop through detectors on board 1 (0-7) and board 2 (8-12)
		for(int j=0;j<13;j++){ // 13

			files = Form("/afs/crc.nd.edu/group/nsl/ast/projects/25Mgang/n1analysis/scripts_new/runs/run324.root");

			if (j<8){
				detect = Form("h_0_%d",j);
			}
			else{
				detect = Form("h_1_%d",j-8);
			}

			// Estimated peak positions from calibration
			double peakPos[] = {442.775,569.12,430.36,537.476,567.59,390.537,494,577.402,506.42,344.102,367.343,444.597,496.303,514.243};
			p1 = peakPos[j];

			// Perform peak fitting
			peakFitter(files,fileBackground,detect,peak1BackLow,peak1BackHigh,
				peak2BackLow,peak2BackHigh,peak1SpecLow,peak1SpecHigh,
				peak2SpecLow,peak2SpecHigh,p1-40,p1+90,j);
			
		}

		//cout << Form("Fitting  complete",fileNum) << endl;
		//fileNum+=1;
   }
	
}



{
	const int NFILES = 2;
	gStyle->SetOptFit(1111);
	TFile *f[NFILES];
	TH1D  *h[NFILES];
	TTree *t[NFILES];
	TFile *f[0] = new TFile("../rootfiles/run_244/UNFILTERED/compass_run_244.root");
	TFile *f[1] = new TFile("../rootfiles/run_324/UNFILTERED/compass_run_324.root");
	f[0]->ls();
	
	double energies[8] = {1460.820,1473.0,1778.969, 511,583.19,609.31, 788.74,2838.29};//Energies in keV
	double channelsErrors[8] = {4,2,2, 3,3,2, 4,2};
	double channels[8] = {464.673,473.301,568.302,163.265,190.099,196.528,258,897.061};
	double energyErrors[8] = {0,0,0,0,0,0,0,0};
	
	TF1 *fitFunc = new TF1("fitFunc","[0]*TMath::Gaus(x[0],[1],[2])+[3]*TMath::Gaus(x[0],[4],[2])",0,5000);
	fitFunc->SetNpx(1e5);
	fitFunc->SetParNames("N1","mean","sigma","N2","mean2");
	fitFunc->SetParameters(2000,channels[0],3,3600,channels[1],20);
	fitFunc->SetParLimits(0,0,10000);	
	fitFunc->SetParLimits(1,channels[0]-5,channels[0]+5);
	fitFunc->SetParLimits(2,2,5);
	fitFunc->SetParLimits(3,0,10000);
	fitFunc->SetParLimits(4,channels[1]-5,channels[1]+5);
	
	for(int i = 0; i<NFILES; i++)
	{
		TTree *t[i] = static_cast<TTree*>(f[i]->Get("Data"));
		TString name = Form("h_%i",i);
		TH1D *h[i] = new TH1D(name,name,10000,0,5000);
		t[i]->Project(name,"Energy+rndm()-0.5","Energy>0 && Board==0 && Channel==1");
	//	h[i]->Rebin(2);
	}
 	
	TCanvas *c0 = new TCanvas("c0","c0",0,0,600,400);
	c0->cd();
	h[0]->Draw("histo");
	h[1]->SetLineColor(1);
	//double xa = h[0]->FindBin(450);
	//double xb = h[0]->FindBin(650);
	//double i0 = h[0]->Integral(xa,xb);//0,16000);
	//double i1 = h[1]->Integral(xa,xb);//0,16000);
	double liveTime0 = 3987;//235;//seconds
	double liveTime1 = 235;//3987;//seconds
	cout << "liveTime0 "<<liveTime0 << "livetime1 " << liveTime1 << " " << liveTime0/liveTime1 << endl;
	h[0]->Sumw2();
	h[1]->Sumw2();
	h[1]->Scale(liveTime0/liveTime1);
	h[1]->Draw("SAME");
	h[0]->GetXaxis()->SetRangeUser(450,485);//Fit range for 1.4 MeV peak
	h[0]->Fit("fitFunc","Br");
	fitFunc->Draw("SAME");
	h[0]->GetXaxis()->SetRangeUser(0,5000);
	c0->Update();
	///  Calibration of background histogram
	
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,400);
	c1->cd();
	
	channels[0] = fitFunc->GetParameter(1);
	channels[1] = fitFunc->GetParameter(4);
      
	TGraphErrors *gr = new TGraphErrors(8,channels,energies,channelsErrors,energyErrors);
	gr->Draw("ap*");
	
	TF1 *f2 = new TF1("mypol2","pol2",0,6000);
	f2->SetParameter(1,0.5);
	f2->SetParameter(0,0);
	gr->Fit("mypol2","B");
	cout << fitFunc->GetChisquare() << " " << fitFunc->GetNDF() << endl;
	
	
	
	
}











{
	const int NFILES = 2;
	gStyle->SetOptFit(1111);
	TFile *f[NFILES];
	TH1D  *h[NFILES];
	TTree *t[NFILES];
	TFile *f[0] = new TFile("../rootfiles/run10427.root");
	TFile *f[1] = new TFile("../rootfiles/run10421.root");
	
	
	TF1 *fitFunc = new TF1("fitFunc","[0]*TMath::Gaus(x[0],[1],[2])+[3]*TMath::Gaus(x[0],[4],[2])",0,5000);
	fitFunc->SetNpx(1e5);
	fitFunc->SetParNames("N1","mean","sigma","N2","mean2");
	fitFunc->SetParameters(1000,3150,20,1000,3225,20);
	fitFunc->SetParLimits(0,0,5000);	
	fitFunc->SetParLimits(1,3100,3190);
	//fitFunc->SetParLimits(2,30,32);
	fitFunc->SetParLimits(3,0,5000);
	fitFunc->SetParLimits(4,3190,3300);
	
	for(int i = 0; i<NFILES; i++)
	{
		TTree *t[i] = static_cast<TTree*>(f[i]->Get("pt0_0"));
		TString name = Form("h_%i",i);
		TH1D *h[i] = new TH1D(name,name,16000,0,16000);
		t[i]->Project(name,"lgate","lgate>0");
		h[i]->Rebin(2);
	}
	TCanvas *c0 = new TCanvas("c0","c0",0,0,600,400);
	c0->cd();
	h[0]->Draw("histo");
	h[1]->SetLineColor(kBlack);
	double xa = h[0]->FindBin(4000);
	double xb = h[0]->FindBin(10000);
	double i0 = h[0]->Integral(xa,xb);//0,16000);
	double i1 = h[1]->Integral(xa,xb);//0,16000);
	cout << i0 << " " << i1 << " " << i0/i1 << endl;
	h[1]->Scale(i0/i1);
	h[1]->Draw("SAME");
	h[1]->GetXaxis()->SetRangeUser(3100,3290);
	h[1]->Fit("fitFunc","Br");
	fitFunc->Draw("SAME");
	h[0]->GetXaxis()->SetRangeUser(0,6000);
	c0->Update();
	///  Calibration of background histogram
	
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,400);
	c1->cd();
	double energies[5] = {511,605,1460.820,1473.0,1633.602};
	double channelsErrors[5] = {20,20,50,15,10};
	double channels[5] = {1133,1351,3151,3219,3551};
	double energyErrors[5] = {0,0,0,0,0};
	channels[2] = fitFunc->GetParameter(1);
	channels[3] = fitFunc->GetParameter(4);
	TGraphErrors *gr = new TGraphErrors(5,channels,energies,channelsErrors,energyErrors);
	gr->Draw("alp");
	
	TF1 *f2 = new TF1("mypol1","pol1",0,6000);
	f2->SetParameter(1,0.5);
	f2->FixParameter(0,0);
	gr->Fit("mypol1","B");
	cout << fitFunc->GetChisquare() << " " << fitFunc->GetNDF() << endl;
	
	
	
	
}











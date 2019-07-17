{
	int RUNNUMBER = 244;
	const int Ndetectors = 13;
	TFile *ifile = new TFile(Form("../rootfiles/run_%d/UNFILTERED/compass_run_%d.root",RUNNUMBER,RUNNUMBER));
	//TFile *ifile = new TFile(Form("../rootfiles/run_%d/FILTERED/compass_run_%d.root",RUNNUMBER,RUNNUMBER));
	TFile *ofile = new TFile("calibratedBackgroundHists.root","RECREATE");
	
	
	const int bnum[Ndetectors] = {0,0,0,0,0,0,0,  1,1, 1,1,1,1 };//
	const int cnum[Ndetectors] = {0,1,2,3,4,5,6,  0,1, 3,4,5,6 };//
	double gain[Ndetectors];
	double offset[Ndetectors];
	double quad[Ndetectors];
	double blah;
	TH1F *h[Ndetectors];
	TCanvas *c0 = new TCanvas("c0","c0",0,0,600,400);
	ifstream ecal("../cal/calDefault.txt");
	if(ecal.is_open())
	{
		for(int i = 0; i<Ndetectors;i++)
		{
			TString cname;
			
			ecal >> cname;
			ecal >> offset[i];
			ecal >> gain[i];
			ecal >> blah;
			TString check = Form("b%i_c%i",bnum[i],cnum[i]);
			TString hname = Form("bgcal%i_%i",bnum[i],cnum[i]);
			
			if(cname != check)
			{
				cout << "Channel mismatch in calibration file (" << cname 
					<<  "!=" << cname << ")" << endl;
			}
			ofile->cd();
			h[i] = new TH1F(hname,hname,10000,0,10000);
			TTree *t = static_cast<TTree*>(ifile->Get("Data"));
			t->Project(hname,Form("(Energy+rndm()-0.5)*%g+%g",gain[i],offset[i]),
				Form("Energy>0 && Board==%d && Channel==%d",bnum[i],cnum[i]));
			//h[i]->SetLineColor(i+1);
			if(i==0)
			{
				h[i]->SetLineColor(kRed);
				h[i]->Draw();	
			}
			else
			{	
				h[i]->Draw("SAME");				

			}
			h[i]->Write();
			c0->cd();

		}
		c0->SaveAs("c0.root");
		ofile->Close();
	}
	else
	{
		cout << " could not open background calibration parameter file " << endl;
	}
	
	
}

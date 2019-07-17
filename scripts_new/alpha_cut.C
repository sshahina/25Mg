


int alpha_cut()
{
	TFile *f = new TFile("../rootfiles/run_244/UNFILTERED/compass_run_244.root");
	TFile *cutFile =new TFile("../scripts/alphaCut_b0c6.root");
	TCutG *alphas = dynamic_cast<TCutG*>(cutFile->Get("alphas"));
	TTree *t = static_cast<TTree*>(f->Get("Data"));
	TCutG *alphas = dynamic_cast<TCutG*>(cutFile->Get("alphas"));
        
	for(int i = 0; i<7; i++)
	{
		
		TString name = Form("h_%i",i);
		TH1D *h[i] = new TH1D(name,name,5000,0,5000);
		t->Project(name,"Energy+rndm()-0.5","Energy>0 && Board==0 && Channel==%d &&!alphas",i);
		TCanvas *c[i] = new TCanvas("c0","c0",0,0,600,400);
		h[i]->Draw();
	}
 	




return 0;
}

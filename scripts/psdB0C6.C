 int psdB0C6(){
	
	gStyle->SetPalette(53);//
	TFile *f = new TFile("../rootfiles/run_244/UNFILTERED/compass_run_244.root");
        TFile *cutFile = new TFile("alphaCut_b0c6.root");
	TCutG *alphas = dynamic_cast<TCutG*>(cutFile->Get("alphas"));
	TCanvas *c0 = new TCanvas("c0","c0");
	c0->cd();
        //alphas->SetName("alphas");
	dataTree->Draw("(EnergyShort+rndm(1)-0.5)/(Energy+rndm(7)-0.5):Energy+rndm(3)-0.5>>hist1(1500,0,1500,50,0.9,1.0)","Energy>0 &&     Board==0 && Channel==0 && Energy>50 && !alphas","colz");
        //alphas->Draw("same");
	



     


	TCanvas *c1 = new TCanvas("c1","c1");
	c1->cd();
	
	dataTree->Draw("Energy>>hist2(5000,0,5000)","!alphas && Board==0 && Channel==6");
	
}

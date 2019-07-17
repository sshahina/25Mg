{
	TFile *f = new TFile("../rootfiles/run_324/UNFILTERED/compass_run_324.root");
	TTree *datatree = static_cast<TTree*>(f->Get("Data"));
	datatree->Scan("Channel:Timestamp:Board:Energy","Board==0 && Channel==1 && Energy>0");
	


}

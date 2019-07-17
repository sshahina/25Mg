{
//========= Macro generated from object: alphas/Graph
//========= by ROOT version5.34/19
   
   TCutG *cutg = new TCutG("alphas",9);
   cutg->SetVarX("Energy+rndm(3)-0.5");
   cutg->SetVarY("(EnergyShort+rndm(1)-0.5)/(Energy+rndm(7)-0.5)");
   cutg->SetTitle("Graph");
   cutg->SetFillColor(1);
   cutg->SetPoint(0,463.06,0.955807);
   cutg->SetPoint(1,482.134,0.946252);
   cutg->SetPoint(2,525.048,0.942057);
   cutg->SetPoint(3,602.931,0.940775);
   cutg->SetPoint(4,779.359,0.940775);
   cutg->SetPoint(5,833.4,0.946835);
   cutg->SetPoint(6,823.863,0.958487);
   cutg->SetPoint(7,634.72,0.959186);
   cutg->SetPoint(8,463.06,0.955807);
   cutg->Draw("");
}

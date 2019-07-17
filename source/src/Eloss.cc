#include "Eloss.hh"

void Eloss::Init(const char* LossFileName, double EMax)
{
 ifstream lossfile;
 lossfile.open(LossFileName);
 string garbage;
 cout << "Using Energy Loss File: " << LossFileName << endl;
 if (lossfile.is_open()) {
    while (garbage!="Straggling") {
    	lossfile >> garbage;
	//	 cout << garbage << endl;
    }
    getline(lossfile,garbage);
    getline(lossfile,garbage);
// cout << dashline <<endl;
 
    int i;
    string units;
    units = "keV";
 
    i=0;

    Energy[0]  = 0;
    dEdXE[0]   = 0;
    dEdXN[0]   = 0;
 
    while (Energy[i++]<EMax) {
       lossfile >> Energy[i] >> units >> dEdXE[i] >> dEdXN[i];
       getline(lossfile,garbage);
       if (units=="keV") {Energy[i] = 0.001*Energy[i];}
    }
    Npoints = i-1;
 } else {
 	 cout << "Cannot open energy loss file " << LossFileName << endl; 
 }
}

double Eloss::GetLoss(double Path, double FinalEnergy)
{
	int k;
	double IntEnergy,DE;
	IntEnergy = FinalEnergy;
	
	for (int i=0;i<NPSteps;i++) {
	   k = 0;	
	   while (IntEnergy>Energy[k++] && k<Npoints);
	   k--;
	   DE = ((dEdXN[k]+dEdXE[k])-(dEdXN[k-1]+dEdXE[k-1]))*(IntEnergy-Energy[k-1])/
	         (Energy[k]-Energy[k-1]) + dEdXN[k-1]+dEdXE[k-1];
//	   cout << "DE=" << DE << endl;
           DE = DE*Path/float(NPSteps);
 //          cout << "DE 2 =" << DE << "Path " << Path <<  endl;
	   IntEnergy += DE; 
	}
//	cout << " " << FinalEnergy << " " << IntEnergy << endl;
	return IntEnergy - FinalEnergy;
}


double Eloss::GetLossFromInit(double Path, double InitialEnergy)
{
	int k;
	double IntEnergy,DE;
	IntEnergy = InitialEnergy;
	
	for (int i=0;i<NPSteps;i++) {
	   k = 0;	
	   while (IntEnergy>Energy[k++] && k<Npoints);
	   k--;
	   DE = ((dEdXN[k]+dEdXE[k])-(dEdXN[k-1]+dEdXE[k-1]))*(IntEnergy-Energy[k-1])/
	         (Energy[k]-Energy[k-1]) + dEdXN[k-1]+dEdXE[k-1];
//	   cout << "DE=" << DE << endl;
           DE = DE*Path/float(NPSteps);
 //          cout << "DE 2 =" << DE << "Path " << Path <<  endl;
	   IntEnergy -= DE; 
	}
//	cout << " " << FinalEnergy << " " << IntEnergy << endl;
	return InitialEnergy - IntEnergy;
}





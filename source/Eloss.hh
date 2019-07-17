//
// Eloss.hh
// --------------
//
// Class containing the method to calculate energy loss using 
// SRIM stopping power table file
//
// K. T. Macon (macon.kevin@gmail.com)
//

#ifndef Eloss_hh
#define Eloss_hh


#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#define NumLoss 1000
#define NPSteps 1000

class Eloss 
{	
	double Energy[NumLoss],dEdXE[NumLoss],dEdXN[NumLoss];
	int Npoints;
public:	
	void Init(const char* LossFileName, double EMax);
	double GetLoss(double Path, double FinalEnergy);
	double GetLossFromInit(double Path, double InitialEnergy);
	virtual ~Eloss() {};
};
#endif // Eloss_hh



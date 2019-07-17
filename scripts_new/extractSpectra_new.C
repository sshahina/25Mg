//
// `root-config --cxx --cflags` -O2 -Wall -Wextra -o extractSpectra extractSpectra.C `root-config --libs`
// ./extractSpectra
//

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TVectorD.h"

#include <iostream>
using std::cout;
using std::endl;

#include <map>
using std::map;

#include <utility>
using std::pair;
using std::make_pair;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <string>
using std::string;

void extractSpectra_new(TString inDirName, TString outFileName)
{
    // const UShort_t chargeBoard   = 1;
    // const UShort_t chargeChannel = 7;
    
    TChain inChain("Data");
    inChain.Add(inDirName + "*.root");

    UShort_t board, channel, energy;
    ULong64_t timestamp;

    inChain.SetBranchStatus("*", 0);

    inChain.SetBranchStatus("Board", 1);
    inChain.SetBranchAddress("Board", &board);
    
    inChain.SetBranchStatus("Channel", 1);
    inChain.SetBranchAddress("Channel", &channel);
    
    inChain.SetBranchStatus("Energy", 1);
    inChain.SetBranchAddress("Energy", &energy);

    inChain.SetBranchStatus("Timestamp", 1);
    inChain.SetBranchAddress("Timestamp", &timestamp);
    
    map <UShort_t, map <UShort_t, TH1D*> > histograms;
    
    map <UShort_t, map <UShort_t, ULong_t > > firstEventTimeStamp;
    map <UShort_t, map <UShort_t, ULong_t > > lastEventTimeStamp;

    TFile outFile(outFileName, "RECREATE");
    
    for (int iBoard = 0; iBoard < 2; iBoard++) {
        for (int iChannel = 0; iChannel < 8; iChannel++) {
            TString histoName = TString::Format("h_%i_%i", iBoard, iChannel);
            histograms[iBoard][iChannel] = new TH1D(histoName, histoName, 8192, 0, 8192);

            firstEventTimeStamp[iBoard][iChannel] = 0;
            lastEventTimeStamp[iBoard][iChannel] = 0;
        }
    }
    
    const size_t nEntries = inChain.GetEntries();
    
    cout << "Filling spectra." << endl;
    for (size_t iEntry = 0; iEntry < nEntries; iEntry++) {
        inChain.GetEntry(iEntry);
	if ((board > 1) || (channel > 7)) continue; // just a precaution
        histograms[board][channel]->Fill(energy);
        
        if (firstEventTimeStamp[board][channel] == 0) {
            firstEventTimeStamp[board][channel] = timestamp;
        }
        lastEventTimeStamp[board][channel] = timestamp;
    }
    cout << endl;
 
    cout << "Writing spectra to file:" << endl; 
    outFile.cd();
    for (int iBoard = 0; iBoard < 2; iBoard++) {
        ULong_t lastTimeStampForThisBoard = 0;
        for (int iChannel = 0; iChannel < 8; iChannel++) {
            cout << "  * " << histograms[iBoard][iChannel]->GetName() << ": " << histograms[iBoard][iChannel]->GetEntries() << " entries" << endl;
            histograms[iBoard][iChannel]->Write("", TObject::kOverwrite);

            if (lastTimeStampForThisBoard < lastEventTimeStamp[iBoard][iChannel]) {
                lastTimeStampForThisBoard = lastEventTimeStamp[iBoard][iChannel];
            }
        }
        TVectorD lastTimeStamp(1);
        lastTimeStamp[0] = 1e-12*lastTimeStampForThisBoard;
        lastTimeStamp.Write(TString::Format("lastTimeStampSeconds-%i", iBoard));
    }
    cout << endl;
}


int extractSpectra_new()
{
	 int upToRun = 434;
	for(int i=324;i<upToRun;i++)
	{
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
		TString inDirName = Form("../rootfiles/run_%d/UNFILTERED/compass_run_%d",i,i);
		// inDirName = ".";
		TString OutFileName = Form("run%d.root",i);
		extractSpectra_new(inDirName,OutFileName);

	}
return 0;
}


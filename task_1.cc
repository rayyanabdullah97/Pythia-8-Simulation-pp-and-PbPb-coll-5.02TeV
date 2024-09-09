//pT distribution of pions, kaons and protons at 5.02 TeV

#include <iostream>
#include "Pythia8/Pythia.h" // Include Pythia headers.
#include "TH1.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
using namespace Pythia8; // Let Pythia8:: be implicit.

int main(int argc, char* argv[]) { // Begin main program.

TApplication theApp("hist", &argc, argv);

// Set up generation.
Pythia pythia; // Declare a Pythia object
//pythia.readString("Top:gg2ttbar = on"); // Switch on process.
pythia.readString("Beams:idA = 2212");
pythia.readString("Beams:idB = 2212");
pythia.readString("Beams:eCM = 5020."); 
pythia.readString("HardQCD:all = on");
pythia.readString("SoftQCD:all = on");
//pythia.readString("Top:qqbar2ttbar = on");
//pythia.readString("Next:numberShowEvent = 5");
pythia.init(); // Initialize; incoming pp beams is default.

TFile* outFile = new TFile("pThist.root", "RECREATE");

  // Book histogram.
TH1F *pT_pion = new TH1F("pT","pions pT distribution", 300, 0, 3.);
TH1F *pT_kaon = new TH1F("pT","kaons pT distribution", 300, 0, 3.);
TH1F *pT_proton = new TH1F("pT","proton pT distribution", 300, 0, 3.);

TCanvas *c1 = new TCanvas("c1", "pT distribution", 800, 600);
    c1->Divide(3, 1); // Divide the canvas into two pads.
// Generate event(s).

//pythia.next(); // Generate an(other) event. Fill event record.
for (int iEvent = 0; iEvent < 10000; ++iEvent) {
pythia.next();
    int iPion =0;
    int iKaon = 0;
    int iProton = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {

        if (pythia.event[i].id() == 111 || pythia.event[i].id() == 211 || pythia.event[i].id() == -211) {
        iPion = i;
        pT_pion->Fill(pythia.event[iPion].pT());
    }
        if (pythia.event[i].id() == 130 || pythia.event[i].id() == 310 || pythia.event[i].id() == 311 || pythia.event[i].id() == 321) {
        iKaon = i;
        pT_kaon->Fill(pythia.event[iKaon].pT());
    }
        if (pythia.event[i].id() == 2212 || pythia.event[i].id() == -2212) {
        iProton = i;
        pT_proton->Fill(pythia.event[iProton].pT());
    }
    }
}
//pythia.stat();
// After creating the histograms
pT_pion->SetXTitle("p_{T} [GeV/c]");  
pT_pion->SetYTitle("Number of Events");  

pT_kaon->SetXTitle("p_{T} [GeV/c]");  
pT_kaon->SetYTitle("Number of Events");  

pT_proton->SetXTitle("p_{T} [GeV/c]");  
pT_proton->SetYTitle("Number of Events"); 

c1->cd(1); // Switch to the first pad.
pT_pion->Draw(); // Draw the pT histogram.

c1->cd(2); // Switch to the second pad.
pT_kaon->Draw(); // Draw the eta histogram.

c1->cd(3);
pT_proton->Draw();
std::cout << "\nDouble click on the histogram window to quit.\n";
gPad->WaitPrimitive();

  // Save histogram on file and close file.
pT_pion->Write();
pT_kaon->Write();
pT_proton->Write();
delete outFile;

c1->SaveAs("pT_distribution.png");
return 0;
} 
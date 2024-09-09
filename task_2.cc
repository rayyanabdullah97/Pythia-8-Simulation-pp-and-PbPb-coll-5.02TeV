#include<iostream>
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TFile.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {
    TApplication theApp("hist", &argc, argv);
    Pythia pythia;

    pythia.readString("Beams:idA = 1000822080");
    pythia.readString("Beams:idB = 1000822080");
    pythia.readString("Beams:eCM = 5020.");
    //pythia.readString("HeavyIon:mode = 1");
    //pythia.readString("HeavyIon:Angantyr = on");
    pythia.readString("Beams:frameType = 1");

      // Initialize the Angantyr model to fit the total and semi-includive
    // cross sections in Pythia within some tolerance.
    pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
    // These parameters are typicall suitable for sqrt(S_NN)=5TeV
    pythia.readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
    // A simple genetic algorithm is run for 20 generations to fit the
    // parameters.
    double genlim[] = {2979.4, 2400.1, 1587.5, 1028.8, 669.9,
                     397.4, 220.3, 116.3, 54.5};
    // If you change any parameters these should also be changed.

    // The upper edge of the correponding percentiles:
    double pclim[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

    pythia.readString("HeavyIon:SigFitNGen = 20");

    pythia.init();

    TH1F *pT_pion = new TH1F("pT_pion", "pT distribution of Pions", 100, 0, 10.);
    TH1F *pT_kaon = new TH1F("pT_kaon", "pT distribution of Kaons", 100, 0, 10.);
    TH1F *pT_proton = new TH1F("pT_proton", "pT distribution of Protons", 100, 0, 10.);

    for (int iEvent = 0; iEvent < 1000; ++iEvent) {
        if (!pythia.next()) continue;

        for (int i = 0; i < pythia.event.size(); ++i) {
            int id = pythia.event[i].id();
            double pT = pythia.event[i].pT();
            
            if (id == 111 || id == 211 || id == -211) { // Pions
                pT_pion->Fill(pT);
            }
            if (id == 130 || id == 310 || id == 311 || id == 321) { // Kaons
                pT_kaon->Fill(pT);
            }
            if (id == 2212 || id == -2212) { // Protons
                pT_proton->Fill(pT);

            }
        }
    }

    TCanvas *c1 = new TCanvas("c1", "pT Distributions in Pb-Pb Collisions at 5.02 TeV", 800, 600);
    pT_pion->SetLineColor(kRed);
    pT_kaon->SetLineColor(kBlue);
    pT_proton->SetLineColor(kGreen);
    
    pT_pion->SetTitle("pT Distribution in Pb-Pb Collisions at 5.02 TeV; pT (GeV); dN/dp_T;");
    pT_pion->Draw();
    pT_kaon->Draw("SAME");
    pT_proton->Draw("SAME");

    // Add a legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(pT_pion, "Pions", "l");
    leg->AddEntry(pT_kaon, "Kaons", "l");
    leg->AddEntry(pT_proton, "Protons", "l");
    leg->Draw();

    TFile *outFile3 = new TFile("pT_distributions_Pb_Pb.root", "RECREATE");
    c1->Write("mycanvas3");
    outFile3->Close();

    pythia.stat();

    std::cout << "\nDouble-click on the canvas to quit." << std::endl;
    c1->WaitPrimitive();
    return 0;
}
//Invariant diff cross section vs pT for Pb-Pb Coll at 5.02TeV CoM Energy
#include<iostream>
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
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

    //Book histogram for invariant differential cross section
    TH1F *cross_pion = new TH1F("cross_pion", "Invariant differential cross section of Pions; pT (GeV); d/sigma/dp_T", 100, 0, 10.);
    TH1F *cross_kaon = new TH1F("cross_kaon", "Invariant differential cross section of Kaons; pT (GeV); d/sigma/dp_T", 100, 0, 10.);
    TH1F *cross_proton = new TH1F("cross_proton", "Invariant differential cross section of Protons; pT (GeV); d/sigma/dp_T", 100, 0, 10.);

    // Book 2D histogram for (pT, y)
    TH2F *h2D_pion = new TH2F("h2D_pion", "pT vs y for Pions; pT (GeV); y", 100, 0, 10., 100, -5, 5);
    TH2F *h2D_kaon = new TH2F("h2D_kaon", "pT vs y for Kaons; pT (GeV); y", 100, 0, 10., 100, -5, 5);
    TH2F *h2D_proton = new TH2F("h2D_proton", "pT vs y for Protons; pT (GeV); y", 100, 0, 10., 100, -5, 5);

    TH1F *nCharged = new TH1F("nCharged","Charged Multiplicity; Number of Charged Particles; Number of Events", 1000, 0, 2000.);

    for (int iEvent = 0; iEvent < 10000; ++iEvent) {
        if (!pythia.next()) continue;
        int multiplicity = 0;
        for (int i = 0; i < pythia.event.size(); ++i) {
            int id = pythia.event[i].id();
            double pT = pythia.event[i].pT();
            double y = pythia.event[i].y();
            
            if (pythia.event[i].isCharged()){
            ++multiplicity;}

            if ( id == 211 || id == -211 /*id == 111 ||*/) { // Pions
                pT_pion->Fill(pT);
                h2D_pion->Fill(pT, y);
            }
            if ( id == 321 || id == -321 /*id == 130 || id == 310 || id == 311 ||*/) { // Kaons
                pT_kaon->Fill(pT);
                h2D_kaon->Fill(pT, y);
                
            }
            if (id == 2212 || id == -2212) { // Protons
                pT_proton->Fill(pT);
                h2D_proton->Fill(pT, y);

            }
        }
        nCharged->Fill(multiplicity);
    }

    // Calculate and plot invariant differential cross section
    TGraph *graph_pion = new TGraph();
    TGraph *graph_kaon = new TGraph();
    TGraph *graph_proton = new TGraph();

    int pointIndex = 0;
    //invariant differential cross section for Pions
    for (int i = 1; i <= h2D_pion->GetNbinsX(); ++i) {
        double pT_pi = h2D_pion->GetXaxis()->GetBinCenter(i);
        double Pionsum_d2N_dpT_dy = 0.0;

        for (int j = 1; j <= h2D_pion->GetNbinsY(); ++j) {
            Pionsum_d2N_dpT_dy += h2D_pion->GetBinContent(i, j);
        }

        double inv_diff_cross_section_pi = Pionsum_d2N_dpT_dy / (2 * M_PI * pT_pi);
        graph_pion->SetPoint(pointIndex++, pT_pi, inv_diff_cross_section_pi);
        cross_pion->Fill(pT_pi, inv_diff_cross_section_pi);
    }
    //invariant differential cross section for Kaons
    for (int i = 1; i <= h2D_kaon->GetNbinsX(); ++i) {
        double pT_ka = h2D_kaon->GetXaxis()->GetBinCenter(i);
        double Kaonsum_d2N_dpT_dy = 0.0;

        for (int j = 1; j <= h2D_kaon->GetNbinsY(); ++j) {
            Kaonsum_d2N_dpT_dy += h2D_kaon->GetBinContent(i, j);
        }

        double inv_diff_cross_section_ka = Kaonsum_d2N_dpT_dy / (2 * M_PI * pT_ka);
        graph_kaon->SetPoint(pointIndex++, pT_ka, inv_diff_cross_section_ka);
        cross_kaon->Fill(pT_ka, inv_diff_cross_section_ka);
    }
    //invariant differential cross section for Protons
    for (int i = 1; i <= h2D_proton->GetNbinsX(); ++i) {
        double pT_pr = h2D_proton->GetXaxis()->GetBinCenter(i);
        double Protonsum_d2N_dpT_dy = 0.0;

        for (int j = 1; j <= h2D_proton->GetNbinsY(); ++j) {
            Protonsum_d2N_dpT_dy += h2D_proton->GetBinContent(i, j);
        }

        double inv_diff_cross_section_pr = Protonsum_d2N_dpT_dy / (2 * M_PI * pT_pr);
        graph_proton->SetPoint(pointIndex++, pT_pr, inv_diff_cross_section_pr);
        cross_proton->Fill(pT_pr, inv_diff_cross_section_pr);
    }
    
    TCanvas *c3 = new TCanvas("c3", "Charged Multiplicity in Pb-Pb Coll at 5.02 TeV", 800, 600);
    c3->SetLogy();
    nCharged->SetTitle("Charged Mutiplicity in Pb-Pb Collisions at 5.02 TeV");

    TCanvas *c1 = new TCanvas("c1", "pT Distributions in Pb-Pb Collisions at 5.02 TeV", 800, 600);
    c1->SetLogy();
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
    pT_pion->Write();
    pT_kaon->Write();
    pT_proton->Write(); 
    nCharged->Write();   
    //h2D_pion->Write();
    //h2D_kaon->Write();
    //h2D_proton->Write();
    outFile3->Close();

    // Plot the graph
    TCanvas *c2 = new TCanvas("c2", "Invariant Differential Cross Section vs pT", 800, 600);
    c2->SetLogy();

    // Set different colors for each graph
    graph_pion->SetMarkerColor(kRed);    // Set pion points to red
    graph_kaon->SetMarkerColor(kBlue);   // Set kaon points to blue
    graph_proton->SetMarkerColor(kGreen);// Set proton points to green

    graph_pion->SetTitle("Invariant Differential Cross Section vs pT for Pb-Pb Coll at 5.02TeV; pT (GeV); d^2#sigma/dp_{T}dy (GeV^{-2})");
    graph_pion->SetMarkerStyle(20);
    graph_kaon->SetMarkerStyle(21);
    graph_proton->SetMarkerStyle(22);

    graph_pion->Draw("AP");
    graph_kaon->Draw("P SAME");
    graph_proton->Draw("P SAME");

    TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(graph_pion, "Pions", "p");
    legend->AddEntry(graph_kaon, "Kaons", "p");
    legend->AddEntry(graph_proton, "Protons", "p");
    legend->Draw();

    // Save the graph to a file
    TFile *outFile2 = new TFile("InvariantDiffCrossSections_Pb_Pb.root", "RECREATE");
    c2->Write("myCanvas4");
    //graph_pion->Write("InvariantDiffCrossSection");
    //graph_kaon->Write();
    //graph_proton->Write();
    h2D_pion->Write();
    h2D_kaon->Write();
    h2D_proton->Write();
    cross_pion->Write();
    cross_kaon->Write();
    cross_proton->Write();
    outFile2->Close();

    pythia.stat();

    std::cout << "\nDouble-click on the canvas to quit." << std::endl;
    c1->WaitPrimitive();
    return 0;
}

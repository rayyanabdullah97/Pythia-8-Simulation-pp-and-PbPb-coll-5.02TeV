//Invariant diff cross section vs pT for pp Coll at 5.02TeV CoM Energy
#include <iostream>
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TFile.h"
#include <TStopwatch.h>
#include "TLatex.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {
    TApplication theApp("hist", &argc, argv);
    TStopwatch timer;
    Pythia pythia;
    
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 5020.");
    pythia.readString("HardQCD:all = on");
    pythia.init();

    TH1F *pT_pion = new TH1F("pT_pion", "pT distribution of Pions", 1000, 0, 10);
    TH1F *pT_kaon = new TH1F("pT_kaon", "pT distribution of Kaons", 1000, 0, 10);
    TH1F *pT_proton = new TH1F("pT_proton", "pT distribution of Protons", 1000, 0, 10);
    
    //Book histogram for invariant differential cross section
    TH1F *cross_pion = new TH1F("cross_pion", "Invariant differential cross section of Pions; pT (GeV); d/sigma/dp_T", 100, 0, 10.);
    TH1F *cross_kaon = new TH1F("cross_kaon", "Invariant differential cross section of Kaons; pT (GeV); d/sigma/dp_T", 100, 0, 10.);
    TH1F *cross_proton = new TH1F("cross_proton", "Invariant differential cross section of Protons; pT (GeV); d/sigma/dp_T", 100, 0, 10.);
    
    // Book 2D histogram for (pT, y)
    TH2F *h2D_pion = new TH2F("h2D_pion", "pT vs y for Pions; pT (GeV); y", 1000, 0, 10., 1000, -5, 5);
    TH2F *h2D_kaon = new TH2F("h2D_kaon", "pT vs y for Kaons; pT (GeV); y", 1000, 0, 10., 1000, -5, 5);
    TH2F *h2D_proton = new TH2F("h2D_proton", "pT vs y for Protons; pT (GeV); y", 1000, 0, 10., 1000, -5, 5);
    
    TH1F *nCharged = new TH1F("nCharged","Charged Multiplicity; Number of Charged Particles; Number of Events", 1000, 0, 1000.);

        // Calculate and plot invariant differential cross section
    TGraph *graph_pion = new TGraph();
    TGraph *graph_kaon = new TGraph();
    TGraph *graph_proton = new TGraph();

    timer.Start();
    // Generate events and fill histograms
    int nEvents = 10000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        int multiplicity = 0;
        for (int i = 0; i < pythia.event.size(); ++i) {
            int id = pythia.event[i].id();
            double pT = pythia.event[i].pT();
            double y = pythia.event[i].y();
            if (pythia.event[i].isCharged()){
                ++multiplicity;}
            if (/*id == 111 ||*/ id == 211 || id == -211) { // Pions
                pT_pion->Fill(pT);
                h2D_pion->Fill(pT, y);
            }
            if (/*id == 130 || id == 310 || id == 311 ||*/ id == 321 || id == -321) { // Kaons
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
    timer.Stop();

    TCanvas *c3 = new TCanvas("c3", "Charged Multiplicity in pp Coll at 5.02 TeV", 800, 600);
    c3->SetLogy();
    nCharged->SetTitle("Charged Mutiplicity in pp Collisions at 5.02 TeV");
    nCharged->Draw();

    TCanvas *c1 = new TCanvas("c1", "pT Distributions", 800, 600);
    c1->SetLogy();
    pT_pion->SetLineColor(kRed);
    pT_kaon->SetLineColor(kBlue);
    pT_proton->SetLineColor(kGreen);
    
    pT_pion->SetTitle("pT Distribution; pT (GeV); dN/dp_T;");
    pT_pion->Draw();
    pT_kaon->Draw("SAME");
    pT_proton->Draw("SAME");

    // Add a legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(pT_pion, "Pions", "l");
    leg->AddEntry(pT_kaon, "Kaons", "l");
    leg->AddEntry(pT_proton, "Protons", "l");
    leg->Draw();

    // Add text annotations using TLatex
    TLatex *latex = new TLatex();
    latex->SetNDC(); // Use normalized device coordinates (NDC) for positioning

    // Set text properties (optional)
    latex->SetTextSize(0.03);
    latex->SetTextAlign(23); // Align at top-left corner of the text box

    // Add the text annotations
    latex->DrawLatex(0.6, 0.85, "#sqrt{s} = 5.02 TeV");  // Center-of-mass energy
    latex->DrawLatex(0.6, 0.80, "pp Collisions");         // Type of collision
    latex->DrawLatex(0.6, 0.75, Form("Events: %d", nEvents)); 

    c1->Update();
    // Save the histograms to a file
    TFile *outFile1 = new TFile("pT_distributions.root", "RECREATE");
    c1->Write("mycanvas1");
    pT_pion->Write();
    pT_kaon->Write();
    pT_proton->Write();
    nCharged->Write(); 
    //h2D_pion->Write();
    //h2D_kaon->Write();
    //h2D_proton->Write();
    outFile1->Close();

    // Plot the graph
    TCanvas *c2 = new TCanvas("c2", "Invariant Differential Cross Section vs pT", 800, 600);
    c2->SetLogy();

// Set names for the graphs explicitly
    graph_pion->SetName("graph_pion");
    graph_kaon->SetName("graph_kaon");
    graph_proton->SetName("graph_proton");
    // Set different colors for each graph
    graph_pion->SetMarkerColor(kRed);    // Set pion points to red
    graph_kaon->SetMarkerColor(kBlue);   // Set kaon points to blue
    graph_proton->SetMarkerColor(kGreen);// Set proton points to green

    graph_pion->SetMarkerSize(0.5);  // Increase pion point size to 1.5
    graph_kaon->SetMarkerSize(0.5);  // Increase kaon point size to 1.5
    graph_proton->SetMarkerSize(0.5);// Increase proton point size to 1.5
    graph_pion->SetTitle("Invariant Differential Cross Section vs pT; pT (GeV); d^2#sigma/dp_{T}dy (GeV^{-2})");
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

    // Add text annotations using TLatex
    TLatex *latex2 = new TLatex();
    latex2->SetNDC(); // Use normalized device coordinates (NDC) for positioning

    // Set text properties (optional)
    latex2->SetTextSize(0.03);
    latex2->SetTextAlign(23); // Align at top-left corner of the text box

    // Add the text annotations
    latex2->DrawLatex(0.6, 0.85, "#sqrt{s} = 5.02 TeV");  // Center-of-mass energy
    latex2->DrawLatex(0.6, 0.80, "pp Collisions");         // Type of collision
    latex2->DrawLatex(0.6, 0.75, Form("Events: %d", nEvents)); 
    // Save the graph to a file
    TFile *outFile2 = new TFile("InvariantDiffCrossSections.root", "RECREATE");
    c2->Write("myCanvas2");
    h2D_pion->Write();
    h2D_kaon->Write();
    h2D_proton->Write();
    graph_pion->Write();
    graph_kaon->Write();
    graph_proton->Write();
    cross_pion->Write();
    cross_kaon->Write();
    cross_proton->Write();
    outFile2->Close();

    //c1->SaveAs("InvariantDiffCrossSection_pions.png");
    //std::cout << "Real time: " << timer.RealTime() << " seconds" << std::endl;
    //std::cout << "CPU time: " << timer.CpuTime() << " seconds" << std::endl;
    timer.Print();
    std::cout << "\nDouble-click on the canvas to quit." << std::endl;
    c1->WaitPrimitive();
    return 0;
}

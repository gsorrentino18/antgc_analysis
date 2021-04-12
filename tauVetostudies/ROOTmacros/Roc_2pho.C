#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <THStack.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TMarker.h>

using namespace std;

void Roc_2pho(){

   TCanvas * c1 = new TCanvas("c1","c1");
   c1->cd(); 
   gPad->SetGrid();

   const Int_t npoints = 6;
   Double_t bkg_rej_20[npoints] = {1-0.6676, 1-0.9438, 1-0.9639, 1-0.9705, 1-0.9632, 1-0.9516};
   Double_t sig_eff_20[npoints] = {0.869496, 0.9923, 0.9939, 0.9950, 0.9938, 0.9927};

   Double_t bkg_rej_15[npoints] = {1-0.5958, 1-0.9150, 1-0.9411, 1-0.9519, 1-0.9397, 1-0.9201};
   Double_t sig_eff_15[npoints] = {0.8365, 0.9898, 0.9922, 0.9936, 0.9919, 0.9899};

   TMultiGraph *roc = new TMultiGraph();
    
   TGraph* roc_20 = new TGraph(npoints, sig_eff_20, bkg_rej_20);
   TGraph* roc_15 = new TGraph(npoints, sig_eff_15, bkg_rej_15);

   roc_15->SetMarkerStyle(20);
   roc_20->SetMarkerStyle(22);

   roc->SetTitle("");
   roc->GetXaxis()->SetTitle("Signal Efficiency (%)");
   roc->GetYaxis()->SetTitle("Background Rejection (%)");
   roc->Add(roc_15, "AP");
   roc->Add(roc_20, "AP");
   roc->Draw("A");

   roc->GetHistogram()->GetXaxis()->SetLimits(.82, .996);
   roc->GetHistogram()->SetMinimum(.82);
   roc->GetHistogram()->SetMaximum(.996);
   roc->GetHistogram()->GetYaxis()->SetRangeUser(0.02,0.43);
   gPad->Modified();
   gPad->Update();
   
   Int_t colors[npoints] = {kOrange, kBlack, kRed, kMagenta, kGreen, kBlue};
   TMarker *m_20[npoints];
   TMarker *m_15[npoints];
   for (Int_t i=0; i<npoints; i++) {
      m_20[i] = new TMarker(sig_eff_20[i], bkg_rej_20[i], 22);
      m_15[i] = new TMarker(sig_eff_15[i], bkg_rej_15[i], 20);
      m_20[i]->SetMarkerColor(colors[i]);
      m_20[i]->SetMarkerSize(1.3);
      m_20[i]->Draw();
      m_15[i]->SetMarkerColor(colors[i]);
      m_15[i]->SetMarkerSize(1.3);
      m_15[i]->Draw();
   }

   //TLegend* leg = new TLegend(0.1, 0.7, 0.4, 0.89); //x1, y1, x2, y2
   TLegend* leg = new TLegend(0.60, 0.6, 0.9, 0.891);
   leg->SetBorderSize(0);
   leg->SetEntrySeparation(0.01);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(m_20[0],"pT > 20 GeV","p");
   leg->AddEntry(m_20[1],"pT > 20 GeV + LooseID","p");
   leg->AddEntry(m_20[2],"pT > 20 GeV + MediumID","p");
   leg->AddEntry(m_20[3],"pT > 20 GeV + TightID","p");
   leg->AddEntry(m_20[4],"pT > 20 GeV + WP80ID","p");
   leg->AddEntry(m_20[5],"pT > 20 GeV + WP90ID","p");
   leg->AddEntry(m_15[0],"pT > 15 GeV","p");
   leg->AddEntry(m_15[1],"pT > 15 GeV + LooseID","p");
   leg->AddEntry(m_15[2],"pT > 15 GeV + MediumID","p");
   leg->AddEntry(m_15[3],"pT > 15 GeV + TightID","p");
   leg->AddEntry(m_15[4],"pT > 15 GeV + WP80ID","p");
   leg->AddEntry(m_15[5],"pT > 15 GeV + WP90ID","p");
   leg->Draw("SAME");
   c1->SaveAs("roc.png");

   Double_t punzi_20[npoints+1] = {0};
   Double_t punzi_15[npoints+1] = {0};
   Double_t a = 2.;
   Double_t b = 5.;
   Double_t bkg_20[npoints+1] = {428.578, 605.915, 618.802, 623.04, 618.393, 610.909, 642.005};
   Double_t bkg_15[npoints+1] = {382.481, 587.414, 604.168, 611.088, 603.275, 590.71, 642.005};
   for (Int_t j=0; j<npoints+1; j++) {
      punzi_20[j] = pow(b,2)/2. + a*sqrt(bkg_20[j]) + (b/2.)*sqrt(b*b +4*a*sqrt(bkg_20[j]) + 4*bkg_20[j]);
      punzi_15[j] = pow(b,2)/2. + a*sqrt(bkg_15[j]) + (b/2.)*sqrt(b*b +4*a*sqrt(bkg_15[j]) + 4*bkg_15[j]);
      //punzi[j] = a*a/8. +9*b*b/13 + a*sqrt(bkg[j]) + (b/2.)*sqrt(b*b +4*a*sqrt(bkg[j]) + 4*bkg[j]);
   }

   cout << "pT > 20 Gev -->   Punzi Smin/SignalEff " << endl;
   cout << "noID Punzi: " << punzi_20[0] << " " << sig_eff_20[0] << " " << punzi_20[0]/sig_eff_20[0] << endl;
   cout << "Loose Punzi: " << punzi_20[1] << " " << sig_eff_20[1] << " " << punzi_20[1]/sig_eff_20[1] << endl;
   cout << "Medium Punzi: " << punzi_20[2] << " "  << sig_eff_20[2] << " " << punzi_20[2]/sig_eff_20[2] << endl;
   cout << "Tight Punzi: " << punzi_20[3] << " "  << sig_eff_20[3] << " " << punzi_20[3]/sig_eff_20[3] << endl;
   cout << "WP80 Punzi: " << punzi_20[4] << " " << sig_eff_20[4] << " " << punzi_20[4]/sig_eff_20[4] << endl;
   cout << "WP90 Punzi: " << punzi_20[5] << " " << sig_eff_20[5] << " " << punzi_20[5]/sig_eff_20[5] << endl;
   cout << "no 2nd photon veto Punzi: " << punzi_20[6] << " " << "1." << " " << punzi_20[6]/1. << endl;

   cout << "pT > 15 Gev -->   Punzi Smin/SignalEff " << endl;
   cout << "noID Punzi: " << punzi_15[0] << " " << sig_eff_15[0] << " " << punzi_15[0]/sig_eff_15[0] << endl;
   cout << "Loose Punzi: " << punzi_15[1] << " " << sig_eff_15[1] << " " << punzi_15[1]/sig_eff_15[1] << endl;
   cout << "Medium Punzi: " << punzi_15[2] << " "  << sig_eff_15[2] << " " << punzi_15[2]/sig_eff_15[2] << endl;
   cout << "Tight Punzi: " << punzi_15[3] << " "  << sig_eff_15[3] << " " << punzi_15[3]/sig_eff_15[3] << endl;
   cout << "WP80 Punzi: " << punzi_15[4] << " " << sig_eff_15[4] << " " << punzi_15[4]/sig_eff_15[4] << endl;
   cout << "WP90 Punzi: " << punzi_15[5] << " " << sig_eff_15[5] << " " << punzi_15[5]/sig_eff_15[5] << endl;
   cout << "no 2nd photon veto Punzi: " << punzi_15[6] << " " << "1." << " " << punzi_15[6]/1. << endl;
  
}


#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <THStack.h>
#include "CMS_lumi.C"
#include <TStyle.h>

using namespace std;

void SecondPhoVeto(){

   TFile *zgamma_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_15GeV_noID/zgamma.root");
   TFile *wgamma_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_15GeV_noID/wgamma.root");

   TTree *zgamma_tree = (TTree*)zgamma_file->Get("tnpPhoIDs/fitter_tree");
   TTree *wgamma_tree = (TTree*)wgamma_file->Get("tnpPhoIDs/fitter_tree");
   
   Float_t ph_etZG;
   Float_t ph_etaZG;
   Float_t ph_phiZG;
   Float_t puWeightZG;
   //Float_t genWeightZG;
   Int_t noIDZG;
   Int_t looseZG;
   Int_t mediumZG;
   Int_t tightZG;
   Int_t WP80ZG;
   Int_t WP90ZG;
   
   Float_t ph_etWG;
   Float_t ph_etaWG;
   Float_t ph_phiWG;
   Float_t puWeightWG;
   //Float_t genWeightWG;
   Int_t noIDWG;
   Int_t looseWG;
   Int_t mediumWG;
   Int_t tightWG;
   Int_t WP80WG;
   Int_t WP90WG;

   zgamma_tree->SetBranchAddress("ph_et",&ph_etZG);
   zgamma_tree->SetBranchAddress("ph_eta",&ph_etaZG);
   zgamma_tree->SetBranchAddress("ph_phi", &ph_phiZG);
   zgamma_tree->SetBranchAddress("puWeight", &puWeightZG);
   zgamma_tree->SetBranchAddress("noID", &noIDZG);
   zgamma_tree->SetBranchAddress("loose", &looseZG);
   zgamma_tree->SetBranchAddress("medium", &mediumZG);
   zgamma_tree->SetBranchAddress("tight", &tightZG);
   zgamma_tree->SetBranchAddress("WP80", &WP80ZG);
   zgamma_tree->SetBranchAddress("WP90", &WP90ZG);
   //zgamma_tree->SetBranchAddress("genWeight", &genWeightZG);
   
   wgamma_tree->SetBranchAddress("ph_et",&ph_etWG);
   wgamma_tree->SetBranchAddress("ph_eta",&ph_etaWG);
   wgamma_tree->SetBranchAddress("ph_phi", &ph_phiWG);
   wgamma_tree->SetBranchAddress("puWeight", &puWeightWG);
   wgamma_tree->SetBranchAddress("noID", &noIDWG);
   wgamma_tree->SetBranchAddress("loose", &looseWG);
   wgamma_tree->SetBranchAddress("medium", &mediumWG);
   wgamma_tree->SetBranchAddress("tight", &tightWG);
   wgamma_tree->SetBranchAddress("WP80", &WP80WG);
   wgamma_tree->SetBranchAddress("WP90", &WP90WG);
   //wgamma_tree->SetBranchAddress("genWeight", &genWeightWG);

   const Int_t nhist=3;
   TH1F* histzgamma[nhist];
   TH1F* histwgamma[nhist];
 
   THStack *hs[nhist];
 
   gStyle->SetOptStat(false);

   histzgamma[0]   = new TH1F("phoEtZG","phoEtZG",40,200,1000);
   histzgamma[1]   = new TH1F ("phoEtaZG","phoEtaZG",20,-1.4442,1.4442);
   histzgamma[2]   = new TH1F ("phoPhiZG", "phoPhiZG", 20, -3.14, 3.14);
 
   histwgamma[0]   = new TH1F("phoEtWG","phoEtWG",40,200,1000);
   histwgamma[1]   = new TH1F("phoEtaWG","phoEtaWG",20,-1.4442,1.4442);
   histwgamma[2]   = new TH1F ("phoPhiWG", "phoPhiWG", 20, -3.14, 3.14);

   Int_t zgamma_nentries = (Int_t)zgamma_tree->GetEntries();
   Int_t wgamma_nentries = (Int_t)wgamma_tree->GetEntries();

   //Float_t sumWeightsZ = 0;
   //Float_t sumWeightsW = 0;
   Float_t zgamma_noID = 0;
   Float_t zgamma_loose = 0;
   Float_t zgamma_medium = 0;
   Float_t zgamma_tight = 0;
   Float_t zgamma_WP80 = 0;
   Float_t zgamma_WP90 = 0;

   Float_t wgamma_noID = 0;
   Float_t wgamma_loose = 0;
   Float_t wgamma_medium = 0;
   Float_t wgamma_tight = 0;
   Float_t wgamma_WP80 = 0;
   Float_t wgamma_WP90 = 0;
 
   for (Int_t i=0; i<zgamma_nentries; i++) {
      zgamma_tree->GetEntry(i);
      //std::cout << looseZG << std::endl;
      //std::cout << "zgamma_loose: " << zgamma_loose << std::endl;
      if (noIDZG==0) zgamma_noID++;
      if (looseZG==0) zgamma_loose++;
      if (mediumZG==0) zgamma_medium++; 
      if (tightZG==0) zgamma_tight++;
      if (WP80ZG==0) zgamma_WP80++;
      if (WP90ZG==0) zgamma_WP90++;
      //sumWeightsZ += genWeightZG;
      histzgamma[0]->Fill(ph_etZG, puWeightZG);
      histzgamma[1]->Fill(ph_etaZG, puWeightZG);
      histzgamma[2]->Fill(ph_phiZG, puWeightZG);
   }

   for (Int_t i=0; i<wgamma_nentries; i++) {
      wgamma_tree->GetEntry(i);
      if (noIDWG==0) wgamma_noID++;
      if (looseWG==0) wgamma_loose++;
      if (mediumWG==0) wgamma_medium++;
      if (tightWG==0) wgamma_tight++;
      if (WP80WG==0) wgamma_WP80++;
      if (WP90WG==0) wgamma_WP90++;
      //sumWeightsW += genWeightWG;
      histwgamma[0]->Fill(ph_etWG, puWeightWG);
      histwgamma[1]->Fill(ph_etaWG, puWeightWG);
      histwgamma[2]->Fill(ph_phiWG, puWeightWG);
   }

   Float_t normZG = 0.2038 * 1000 * 41.56/5110900.;
   Float_t normWG = 0.7158 * 1000 * 41.56/4801400;

   Float_t zw = (zgamma_nentries*normZG*1.0)/(wgamma_nentries*normWG);
   std::cout << "zgamma/wgamma baseline:  " << zgamma_nentries*normZG << "/" << wgamma_nentries*normWG << " = " << zw << std::endl;

   Float_t zw_noID = (zgamma_noID*normZG*1.0)/(wgamma_noID*normWG);
   Float_t zw_loose = (zgamma_loose*normZG*1.0)/(wgamma_loose*normWG);
   Float_t zw_medium = (zgamma_medium*normZG*1.0)/(wgamma_medium*normWG);
   Float_t zw_tight = (zgamma_tight*normZG*1.0)/(wgamma_tight*normWG);
   Float_t zw_WP80 = (zgamma_WP80*normZG*1.0)/(wgamma_WP80*normWG);
   Float_t zw_WP90 = (zgamma_WP90*normZG*1.0)/(wgamma_WP90*normWG);

   std::cout << "zgamma/wgamma noID:  " << zgamma_noID*normZG << "/" << wgamma_noID*normWG << " = " << zw_noID << std::endl;
   std::cout << "zgamma/wgamma LOOSE:  " << zgamma_loose*normZG << "/" << wgamma_loose*normWG << " = " << zw_loose << std::endl;
   std::cout << "zgamma/wgamma MEDIUM:  " << zgamma_medium*normZG << "/" << wgamma_medium*normWG << " = " << zw_medium << std::endl;
   std::cout << "zgamma/wgamma TIGHT:  " << zgamma_tight*normZG << "/" << wgamma_tight*normWG << " = " << zw_tight << std::endl;
   std::cout << "zgamma/wgamma WP80:  " << zgamma_WP80*normZG << "/" << wgamma_WP80*normWG << " = " << zw_WP80 << std::endl;
   std::cout << "zgamma/wgamma WP90:  " << zgamma_WP90*normZG << "/" << wgamma_WP90*normWG << " = " << zw_WP90 << std::endl;

   std::cout << "signal eff noID:  " << zgamma_noID/zgamma_nentries << std::endl;
   std::cout << "signal eff LOOSE:  " << zgamma_loose/zgamma_nentries << std::endl;
   std::cout << "signal eff MEDIUM:  " << zgamma_medium/zgamma_nentries << std::endl;
   std::cout << "signal eff TIGHT:  " << zgamma_tight/zgamma_nentries << std::endl;
   std::cout << "signal eff WP80:  " << zgamma_WP80/zgamma_nentries << std::endl; 
   std::cout << "signal eff WP90:  " << zgamma_WP90/zgamma_nentries << std::endl; 

   std::cout << "background eff noID:  " << wgamma_noID/wgamma_nentries << std::endl;
   std::cout << "background eff LOOSE:  " << wgamma_loose/wgamma_nentries << std::endl;
   std::cout << "background eff MEDIUM:  " << wgamma_medium/wgamma_nentries << std::endl;
   std::cout << "background eff TIGHT:  " << wgamma_tight/wgamma_nentries << std::endl;
   std::cout << "background eff WP80:  " << wgamma_WP80/wgamma_nentries << std::endl;
   std::cout << "background eff WP90:  " << wgamma_WP90/wgamma_nentries << std::endl;

   TString title[nhist] = {"ph_et", "ph_eta", "ph_phi"};
   TString title_axis[nhist] = {"p_{T}^{#gamma}", "#eta", "#phi"}; 
   
   for (int i=0; i<nhist; i++) {

      TCanvas *c2 = new TCanvas("c2", "c2");
      c2->cd();

      TPad* pad11 = new TPad("pad11", "pad11", 0.0, 0.3, 1.0, 1.0);
      pad11->SetBottomMargin(0);
      pad11->Draw();
      pad11->cd();

      hs[i] = new THStack(title[i]+"_stack", "");

      histzgamma[i]->SetFillColor(kSpring+9);
      histzgamma[i]->SetLineColor(kSpring+9);

      histwgamma[i]->SetFillColor(kAzure+1);
      histwgamma[i]->SetLineColor(kAzure+1);
      std::cout << "test" << std::endl;

      hs[i]->Add(histzgamma[i]);
      hs[i]->Add(histwgamma[i]);

      hs[i]->Draw("hist");
      hs[i]->GetYaxis()->SetTitle("Events");

      TLegend* leg_stack = new TLegend();
      leg_stack = new TLegend(0.79, 0.5, 0.9, 0.891);
      leg_stack->SetBorderSize(0);
      leg_stack->SetEntrySeparation(0.01);
      leg_stack->SetFillColor(0);
      leg_stack->SetFillStyle(0);

      leg_stack->AddEntry(histwgamma[i], "WG", "f");
      leg_stack->AddEntry(histzgamma[i], "ZG", "f");
      leg_stack->Draw("SAME");
      
      pad11->SetLogy();

      CMS_lumi(pad11, 17, 0);
      pad11->Update();
      c2->Update();
      c2->cd();

      c2->SaveAs(title[i]+"stack.png");
   }

}

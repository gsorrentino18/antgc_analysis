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

void tauveto(){

   TFile *zgamma_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/tauVetostudies/withpixelveto/Tight_final/zgamma.root");
   TFile *wgamma_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/tauVetostudies/withpixelveto/Tight_tau/wgamma.root");

   TTree *zgamma_tree = (TTree*)zgamma_file->Get("tnpPhoIDs/fitter_tree");
   TTree *wgamma_tree = (TTree*)wgamma_file->Get("tnpPhoIDs/fitter_tree");
   
   Float_t ph_etZG;
   Float_t ph_etaZG;
   Float_t ph_phiZG;
   Float_t puWeightZG;
   Char_t lepVetoZG;
   Int_t tauIsGenMatchedZG;    

   Float_t ph_etWG;
   Float_t ph_etaWG;
   Float_t ph_phiWG;
   Float_t puWeightWG;
   Char_t lepVetoWG;
   Int_t tauIsGenMatchedWG;
   Int_t hasPromptTauWG;

   zgamma_tree->SetBranchAddress("ph_et",&ph_etZG);
   zgamma_tree->SetBranchAddress("ph_eta",&ph_etaZG);
   zgamma_tree->SetBranchAddress("ph_phi", &ph_phiZG);
   zgamma_tree->SetBranchAddress("puWeight", &puWeightZG);
   zgamma_tree->SetBranchAddress("lepVeto", &lepVetoZG);  
   zgamma_tree->SetBranchAddress("tauIsGenMatched", &tauIsGenMatchedZG);  

   wgamma_tree->SetBranchAddress("ph_et",&ph_etWG);
   wgamma_tree->SetBranchAddress("ph_eta",&ph_etaWG);
   wgamma_tree->SetBranchAddress("ph_phi", &ph_phiWG);
   wgamma_tree->SetBranchAddress("puWeight", &puWeightWG);
   wgamma_tree->SetBranchAddress("lepVeto", &lepVetoWG);  
   wgamma_tree->SetBranchAddress("tauIsGenMatched", &tauIsGenMatchedWG);
   wgamma_tree->SetBranchAddress("hasPromptTau", &hasPromptTauWG);

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

   Float_t tau_zgamma_loose = 0;
   Float_t tau_zgamma_loose_gm = 0;
   Float_t tau_wgamma_loose = 0;
   Float_t tau_wgamma_loose_gm = 0;
   Float_t wgamma = 0;

   for (Int_t i=0; i<zgamma_nentries; i++) {
      zgamma_tree->GetEntry(i);
      if (lepVetoZG==(char)4) {
        tau_zgamma_loose++;
        //if (tauIsGenMatchedZG) {
          //tau_zgamma_loose_gm++;
        //}
      }
      histzgamma[0]->Fill(ph_etZG, puWeightZG);
      histzgamma[1]->Fill(ph_etaZG, puWeightZG);
      histzgamma[2]->Fill(ph_phiZG, puWeightZG);
   }

   for (Int_t i=0; i<wgamma_nentries; i++) {
      wgamma_tree->GetEntry(i);
      //if (!hasPromptTauWG) continue;
      wgamma++;
      if (lepVetoWG==(char)4) {
        tau_wgamma_loose++;
        if (tauIsGenMatchedWG) {
          tau_wgamma_loose_gm++;
        }
      }
      histwgamma[0]->Fill(ph_etWG, puWeightWG);
      histwgamma[1]->Fill(ph_etaWG, puWeightWG);
      histwgamma[2]->Fill(ph_phiWG, puWeightWG);
   }

   Float_t normZG = 0.2038 * 1000 * 41.56/5110900.;
   Float_t normWG = 0.7158 * 1000 * 41.56/4801400;

   Float_t zw = (zgamma_nentries*normZG)/(wgamma*normWG);

   Float_t zgamma_loose = zgamma_nentries - tau_zgamma_loose;
   Float_t zgamma_loose_gm = zgamma_nentries - tau_zgamma_loose_gm;
   Float_t wgamma_loose = wgamma - tau_wgamma_loose;
   Float_t wgamma_loose_gm = wgamma - tau_wgamma_loose_gm;

   std::cout << "#wgamma evts: " << wgamma*normWG << std::endl;
   std::cout << "#tau in wgamma: " << tau_wgamma_loose*normWG << "   #tau gen matched in wgamma: " << tau_wgamma_loose_gm*normWG << std::endl;

   std::cout << "zgamma/wgamma baseline:  " << zgamma_nentries*normZG << "/" << wgamma*normWG << " = " << zw << std::endl;
   std::cout << "zgamma/wgamma:  " << zgamma_loose*normZG << "/" << wgamma_loose*normWG << " = " << (zgamma_loose*normZG)/(wgamma_loose*normWG) << std::endl;
   std::cout << "zgamma/wgamma only genmatched tau:  " << zgamma_loose_gm*normZG << "/" << wgamma_loose_gm*normWG << " = " << (zgamma_loose_gm*normZG)/(wgamma_loose_gm*normWG) << std::endl; 

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

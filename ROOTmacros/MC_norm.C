#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>
#include <iostream>

#include "CMS_lumi.C"

using namespace std;

void MC_norm(Int_t index) {
   
   //double norm = xsec * 1000. * lumi2017 / ngen; 
   Double_t lumi = 41.55;
   Double_t norm100to200 = 93.66 * 1000. * lumi / 2230077;
   Double_t norm200to400 = 4.135 * 1000. * lumi / 688400;
   Double_t norm400to800 = 0.2457 * 1000. * lumi / 955396;
   Double_t norm800to2000 = 0.009936 * 1000. * lumi / 769820;

   TFile *MC100to200_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/test2/dy100to200.root"); 
   TFile *MC200to400_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/test2/dy200to400.root");
   TFile *MC400to800_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/test2/dy400to800.root");
   TFile *MC800to2000_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/test2/dy800to2000.root");
   TFile *data_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/testdata2/data.root");

   TCanvas *c1 = new TCanvas("c1", "c1");
   c1->cd();

   TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
   pad1->SetBottomMargin(0.001);
   pad1->Draw();
   pad1->cd();

   gStyle->SetOptStat(false);

   const Int_t nhist=8;
   TH1F* histMC100to200[nhist];
   TH1F* histMC200to400[nhist];
   TH1F* histMC400to800[nhist];
   TH1F* histMC800to2000[nhist];
   TH1F* histdata[nhist];

   TString title[nhist] = {"phoBDT", "invmass", "deltaPhi", "deltaEta", "deltaPt", "deltaRs", "elePt", "phoPt"};

   for (Int_t i=0; i<nhist; i++) { 
      histMC100to200[i]=(TH1F*)MC100to200_file->Get(title[i]);
      histMC200to400[i]=(TH1F*)MC200to400_file->Get(title[i]);
      histMC400to800[i]=(TH1F*)MC400to800_file->Get(title[i]);
      histMC800to2000[i]=(TH1F*)MC800to2000_file->Get(title[i]);
      histdata[i]=(TH1F*)data_file->Get(title[i]);
   }

   for (Int_t i=0; i<nhist; i++) { 
      histMC100to200[i]->Scale(norm100to200);
      histMC200to400[i]->Scale(norm200to400);
      histMC400to800[i]->Scale(norm400to800);
      histMC800to2000[i]->Scale(norm800to2000);
   }

   TH1F* histMC = (TH1F*)histMC100to200[index]->Clone("histMC");
   histMC->Add(histMC200to400[index]);
   histMC->Add(histMC400to800[index]);
   histMC->Add(histMC800to2000[index]);

   TH1F* h_ratio = (TH1F*)histdata[index]->Clone("h_ratio");
   h_ratio->Divide(histMC);
   
   histMC->SetFillColor(kRed-9);
   //h_MC_phoBDTpred->SetLineColor(kRed-7);
   histdata[index]->SetMarkerStyle(20);

   histMC->GetYaxis()->SetTitle("Events");
   histMC->GetYaxis()->SetTitleSize(0.05);
   histMC->GetYaxis()->SetTitleOffset(0.95);
   histMC->GetYaxis()->SetLabelSize(0.045);

   TLegend* leg = new TLegend();
   //leg = new TLegend(0.1289,0.7436,0.5100,0.8823);
   leg = new TLegend(0.669, 0.754, 0.986, 0.891);
   leg->SetBorderSize(0);
   leg->SetEntrySeparation(0.01);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);

   leg->AddEntry(histMC, "DYJetsToEE", "f");
   leg->AddEntry(histdata[index], "Data", "ep");
   c1->Update();
  
   histMC->SetTitle("");
   histdata[index]->SetTitle(""); 

   //histMC->GetYaxis()->SetRangeUser(histMC->GetMinimum()+0.00001,histdata[index]->GetMaximum()+10000);
   histMC->Draw("hist");
   histdata[index]->Draw("ep, SAME");
   //pad1->Update();

   leg->Draw("SAME");
   if ((title[index] == "phoBDT") || (title[index] == "elePt") || (title[index] == "phoPt") || (title[index] == "deltaPhi") || (title[index] == "deltaEta") || (title[index] == "deltaPt")) {
      pad1->SetLogy();
   }
   CMS_lumi(pad1, 17, 0);
   pad1->Update();
   c1->Update();
   c1->cd();

   TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
   pad2->SetTopMargin(0.01);
   pad2->SetBottomMargin(0.2);
   pad2->Draw();
   pad2->cd();


   if (title[index] == "phoBDT") {
      h_ratio->GetXaxis()->SetTitle("Photon BDT score");
   }

   if (title[index] == "invmass") {
      h_ratio->GetXaxis()->SetTitle("e#gamma invariant mass");
   }
   if (title[index] == "deltaPhi") {
      h_ratio->GetXaxis()->SetTitle("(#Delta#ph)i_{e#gamma}");
   }
   if (title[index] == "deltaEta") {
      h_ratio->GetXaxis()->SetTitle("(#Delta#eta)_{e#gamma}");
   }
   if (title[index] == "deltaPt") {
      h_ratio->GetXaxis()->SetTitle("(#DeltaPt)_{e#gamma}");
   }
   if (title[index] == "deltaRs") {
      h_ratio->GetXaxis()->SetTitle("(#DeltaR)_{e#gamma}");
   }
   if (title[index] == "elePt") {
      h_ratio->GetXaxis()->SetTitle("p_{T}^{e}");
   }
   if (title[index] == "phoPt") {
      h_ratio->GetXaxis()->SetTitle("p_{T}^{#gamma}");
   }
   
   h_ratio->SetTitle("");
   h_ratio->GetXaxis()->SetTitleFont(42);
   h_ratio->GetXaxis()->SetTitleSize(0.11);
   h_ratio->GetXaxis()->SetTitleOffset(0.8);
   h_ratio->GetXaxis()->SetLabelFont(42);
   h_ratio->GetXaxis()->SetLabelSize(0.1);
   h_ratio->GetYaxis()->SetTitle("Data/MC");
   h_ratio->GetYaxis()->SetTitleSize(0.11);
   h_ratio->GetYaxis()->SetTitleOffset(0.43);
   h_ratio->GetYaxis()->SetLabelSize(0.1);
   h_ratio->GetYaxis()->SetLabelOffset(0.01);
   h_ratio->GetYaxis()->SetNdivisions(505);
   h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
   h_ratio->SetMarkerStyle(20);
  
   TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
   line->SetLineColor(kRed);
   line->SetLineWidth(1);
   h_ratio->Draw();
   line->Draw("SAME");
   c1->Update();
   
   c1->SaveAs(title[index]+"_norm.png");


}

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <THStack.h>

using namespace std;

void BDTValidation(){

   TFile *signalMC_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/testPLOTMC/signalMC.root"); 
   TFile *qcdMC_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/testPLOTMC/qcdMC.root");
   TFile *data_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/testPLOTDATA/data.root");

   TTree *signalMC_tree = (TTree*)signalMC_file->Get("fullEB/fullEBTree");
   TTree *qcdMC_tree = (TTree*)qcdMC_file->Get("fullEB/fullEBTree");
   TTree *data_tree = (TTree*)data_file->Get("fullEB/fullEBTree"); 

   Float_t signalMC_phoBDTpred;
   Float_t qcdMC_phoBDTpred;
   Float_t data_phoBDTpred;

   signalMC_tree->SetBranchAddress("phoBDTpred",&signalMC_phoBDTpred);
   qcdMC_tree->SetBranchAddress("phoBDTpred",&qcdMC_phoBDTpred);
   data_tree->SetBranchAddress("phoBDTpred",&data_phoBDTpred);

   TCanvas *c1 = new TCanvas("c1", "c1");
   c1->cd();

   THStack *hs = new THStack("hs","");
   //hs->GetXaxis()->SetTitle("BDT score");
   //hs->GetYaxis()->SetTitle("Events");

   TH1F *h_signalMC_phoBDTpred   = new TH1F("h_signalMC_phoBDTpred","signalMC BDT distribution",100,0,1);
   TH1F *h_qcdMC_phoBDTpred   = new TH1F("h_qcdMC_phoBDTpred","qcdMC BDT distribution",100,0,1);
   TH1F *h_data_phoBDTpred   = new TH1F("h_data_phoBDTpred","data BDT distribution",100,0,1);
   
   Int_t signalMC_nentries = (Int_t)signalMC_tree->GetEntries();
   Int_t qcdMC_nentries = (Int_t)qcdMC_tree->GetEntries();
   Int_t data_nentries = (Int_t)data_tree->GetEntries();

   for (Int_t i=0; i<signalMC_nentries; i++) {
      signalMC_tree->GetEntry(i);
      h_signalMC_phoBDTpred->Fill(signalMC_phoBDTpred);
   }
   for (Int_t i=0; i<qcdMC_nentries; i++) {
      qcdMC_tree->GetEntry(i);
      h_qcdMC_phoBDTpred->Fill(qcdMC_phoBDTpred);
   }
   for (Int_t i=0; i<data_nentries; i++) {
      data_tree->GetEntry(i);
      h_data_phoBDTpred->Fill(data_phoBDTpred);
   }   
   
   //Double_t factor = 1.0;
   //h_signalMC_phoBDTpred->Scale(factor/h_signalMC_phoBDTpred->Integral());
   //h_qcdMC_phoBDTpred->Scale(factor/h_qcdMC_phoBDTpred->Integral());
   //h_data_phoBDTpred->Scale(factor/h_data_phoBDTpred->Integral());

   h_signalMC_phoBDTpred->SetFillColor(kBlue-4);
   h_qcdMC_phoBDTpred->SetFillColor(kRed-7);
   //h_data_phoBDTpred->SetFillColor(kGreen);
  
   hs->Add(h_signalMC_phoBDTpred);
   hs->Add(h_qcdMC_phoBDTpred);
   //hs->Add(h_data_phoBDTpred);
   //hs->Draw();
   h_data_phoBDTpred->Draw("ep");
   
   hs->GetXaxis()->SetTitle("BDT score");
   hs->GetYaxis()->SetTitle("Events");
   c1->Update();
  
   TLegend* leg = new TLegend();
   leg = new TLegend(0.1289,0.7436,0.5100,0.8823);
   leg->SetBorderSize(0);
   leg->SetEntrySeparation(0.01);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);

   leg->AddEntry(h_signalMC_phoBDTpred, "DY", "f");
   leg->AddEntry(h_qcdMC_phoBDTpred, "QCD", "f");

   c1->Update();
   hs->Draw();
   leg->Draw("SAME");
   c1->SaveAs("BDTPlot.png");
}

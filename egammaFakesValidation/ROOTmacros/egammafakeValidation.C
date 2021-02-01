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

void egammafakeValidation(){

   TFile *data_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/egammaFakesValidation/data_01FebV2.root");
   TFile *signalMC_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/egammaFakesValidation/MCall_01FebV2.root");
   /*TFile *wjets_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/wgjets_18Jan.root");
   TFile *wlnu_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/wglnu_18Jan.root");
   TFile *gjets_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/gjet_18Jan.root");
   TFile *qcd_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/qcdht_18Jan.root");
   TFile *ttgjets_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/TTgjets_18Jan.root");
   TFile *ww_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/ww_18Jan.root");*/
   TFile *wele_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/egammaFakesValidation/wele_01FebV2.root");  //wele_28Jan_bis.root
   TFile *wz_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/egammaFakesValidation/wz_01FebV2.root");
   TFile *zz_file = new TFile("/local/cms/user/gsorrent/antgc_analysis/egammaFakesValidation/zz_01FebV2.root");

   TTree *data_tree = (TTree*)data_file->Get("tnpPhoIDs/fitter_tree"); 
   TTree *signalMC_tree = (TTree*)signalMC_file->Get("tnpPhoIDs/fitter_tree");
   /*TTree *gjets_tree = (TTree*)gjets_file->Get("tnpPhoIDs/fitter_tree");
   TTree *wjets_tree = (TTree*)wjets_file->Get("tnpPhoIDs/fitter_tree");
   TTree *wlnu_tree = (TTree*)wlnu_file->Get("tnpPhoIDs/fitter_tree");
   TTree *qcd_tree = (TTree*)qcd_file->Get("tnpPhoIDs/fitter_tree");
   TTree *ttgjets_tree = (TTree*)ttgjets_file->Get("tnpPhoIDs/fitter_tree");
   TTree *ww_tree = (TTree*)ww_file->Get("tnpPhoIDs/fitter_tree");*/
   TTree *wele_tree = (TTree*)wele_file->Get("tnpPhoIDs/fitter_tree");
   TTree *wz_tree = (TTree*)wz_file->Get("tnpPhoIDs/fitter_tree");
   TTree *zz_tree = (TTree*)zz_file->Get("tnpPhoIDs/fitter_tree");

   Float_t puWeight;
   Float_t ph_pFECALClusIsoCorr;
   Float_t ph_pFHCALClusIsoCorr;
   Float_t ph_TrkIsoCorr;
   Float_t ph_HoverE;
   Float_t ph_HoverECorr;
   Float_t ph_BDTpred;
   Float_t ph_et;
   Float_t ph_eta;
   Float_t ph_R9Full5x5;
   Float_t ph_S4Full5x5;
   Float_t ph_EmaxOESCrFull5x5;
   Float_t ph_E2ndOESCrFull5x5;
   Float_t ph_E1x3OESCrFull5x5;
   Float_t ph_E2x5OESCrFull5x5;
   Float_t ph_E5x5OESCrFull5x5;
   Float_t ph_2x2OE3x3Full5x5;
   Float_t ph_SigmaIEtaIEta;
   Float_t ph_SigmaIEtaIPhi;
   Float_t ph_SigmaIPhiIPhi;
   Float_t ph_SieieOSipipFull5x5;
   Float_t ph_EtaWidth;
   Float_t ph_PhiWidth;
   Float_t ph_EtaWOPhiWFull5x5;
   Bool_t pass95;
   Float_t FRweight;
   Float_t met;

   Float_t puWeightMC;
   Float_t ph_pFECALClusIsoCorrMC;
   Float_t ph_pFHCALClusIsoCorrMC;
   Float_t ph_TrkIsoCorrMC;
   Float_t ph_HoverEMC;
   Float_t ph_HoverECorrMC;
   Float_t ph_BDTpredMC;
   Float_t ph_etMC;
   Float_t ph_etaMC;
   Float_t ph_R9Full5x5MC;
   Float_t ph_S4Full5x5MC;
   Float_t ph_EmaxOESCrFull5x5MC;
   Float_t ph_E2ndOESCrFull5x5MC;
   Float_t ph_E1x3OESCrFull5x5MC;
   Float_t ph_E2x5OESCrFull5x5MC;
   Float_t ph_E5x5OESCrFull5x5MC;
   Float_t ph_2x2OE3x3Full5x5MC;
   Float_t ph_SigmaIEtaIEtaMC;
   Float_t ph_SigmaIEtaIPhiMC;
   Float_t ph_SigmaIPhiIPhiMC;
   Float_t ph_SieieOSipipFull5x5MC;
   Float_t ph_EtaWidthMC;
   Float_t ph_PhiWidthMC;
   Float_t ph_EtaWOPhiWFull5x5MC;
   Bool_t pass95MC;
   Float_t FRweightMC;
   Float_t metMC;
           
   /*Float_t puWeightGj;
   Float_t ph_pFECALClusIsoCorrGj;
   Float_t ph_pFHCALClusIsoCorrGj;
   Float_t ph_TrkIsoCorrGj;
   Float_t ph_HoverEGj;
   Float_t ph_HoverECorrGj;
   Float_t ph_BDTpredGj;
   Float_t ph_etGj;
   Float_t ph_R9Full5x5Gj;
   Float_t ph_S4Full5x5Gj;
   Float_t ph_EmaxOESCrFull5x5Gj;
   Float_t ph_E2ndOESCrFull5x5Gj;
   Float_t ph_E1x3OESCrFull5x5Gj;
   Float_t ph_E2x5OESCrFull5x5Gj;
   Float_t ph_E5x5OESCrFull5x5Gj;
   Float_t ph_2x2OE3x3Full5x5Gj;
   Float_t ph_SigmaIEtaIEtaGj;
   Float_t ph_SigmaIEtaIPhiGj;
   Float_t ph_SigmaIPhiIPhiGj;
   Float_t ph_SieieOSipipFull5x5Gj;
   Float_t ph_EtaWidthGj;
   Float_t ph_PhiWidthGj;
   Float_t ph_EtaWOPhiWFull5x5Gj;
   Bool_t pass95Gj;

   Float_t puWeightWj;
   Float_t ph_pFECALClusIsoCorrWj;
   Float_t ph_pFHCALClusIsoCorrWj;
   Float_t ph_TrkIsoCorrWj;
   Float_t ph_HoverEWj;
   Float_t ph_HoverECorrWj;
   Float_t ph_BDTpredWj;
   Float_t ph_etWj;
   Float_t ph_R9Full5x5Wj;
   Float_t ph_S4Full5x5Wj;
   Float_t ph_EmaxOESCrFull5x5Wj;
   Float_t ph_E2ndOESCrFull5x5Wj;
   Float_t ph_E1x3OESCrFull5x5Wj;
   Float_t ph_E2x5OESCrFull5x5Wj;
   Float_t ph_E5x5OESCrFull5x5Wj;
   Float_t ph_2x2OE3x3Full5x5Wj;
   Float_t ph_SigmaIEtaIEtaWj;
   Float_t ph_SigmaIEtaIPhiWj;
   Float_t ph_SigmaIPhiIPhiWj;
   Float_t ph_SieieOSipipFull5x5Wj;
   Float_t ph_EtaWidthWj;
   Float_t ph_PhiWidthWj;
   Float_t ph_EtaWOPhiWFull5x5Wj;
   Bool_t pass95Wj;

   Float_t puWeightWlnu;
   Float_t ph_pFECALClusIsoCorrWlnu;
   Float_t ph_pFHCALClusIsoCorrWlnu;
   Float_t ph_TrkIsoCorrWlnu;
   Float_t ph_HoverEWlnu;
   Float_t ph_HoverECorrWlnu;
   Float_t ph_BDTpredWlnu;
   Float_t ph_etWlnu;
   Float_t ph_R9Full5x5Wlnu;
   Float_t ph_S4Full5x5Wlnu;
   Float_t ph_EmaxOESCrFull5x5Wlnu;
   Float_t ph_E2ndOESCrFull5x5Wlnu;
   Float_t ph_E1x3OESCrFull5x5Wlnu;
   Float_t ph_E2x5OESCrFull5x5Wlnu;
   Float_t ph_E5x5OESCrFull5x5Wlnu;
   Float_t ph_2x2OE3x3Full5x5Wlnu;
   Float_t ph_SigmaIEtaIEtaWlnu;
   Float_t ph_SigmaIEtaIPhiWlnu;
   Float_t ph_SigmaIPhiIPhiWlnu;
   Float_t ph_SieieOSipipFull5x5Wlnu;
   Float_t ph_EtaWidthWlnu;
   Float_t ph_PhiWidthWlnu;
   Float_t ph_EtaWOPhiWFull5x5Wlnu;
   Bool_t pass95Wlnu;

   Float_t puWeightQcd;
   Float_t ph_pFECALClusIsoCorrQcd;
   Float_t ph_pFHCALClusIsoCorrQcd;
   Float_t ph_TrkIsoCorrQcd;
   Float_t ph_HoverEQcd;
   Float_t ph_HoverECorrQcd;
   Float_t ph_BDTpredQcd;
   Float_t ph_etQcd;
   Float_t ph_R9Full5x5Qcd;
   Float_t ph_S4Full5x5Qcd;
   Float_t ph_EmaxOESCrFull5x5Qcd;
   Float_t ph_E2ndOESCrFull5x5Qcd;
   Float_t ph_E1x3OESCrFull5x5Qcd;
   Float_t ph_E2x5OESCrFull5x5Qcd;
   Float_t ph_E5x5OESCrFull5x5Qcd;
   Float_t ph_2x2OE3x3Full5x5Qcd;
   Float_t ph_SigmaIEtaIEtaQcd;
   Float_t ph_SigmaIEtaIPhiQcd;
   Float_t ph_SigmaIPhiIPhiQcd;
   Float_t ph_SieieOSipipFull5x5Qcd;
   Float_t ph_EtaWidthQcd;
   Float_t ph_PhiWidthQcd;
   Float_t ph_EtaWOPhiWFull5x5Qcd;
   Bool_t pass95Qcd;

   Float_t puWeightTTGj;
   Float_t ph_pFECALClusIsoCorrTTGj;
   Float_t ph_pFHCALClusIsoCorrTTGj;
   Float_t ph_TrkIsoCorrTTGj;
   Float_t ph_HoverETTGj;
   Float_t ph_HoverECorrTTGj;
   Float_t ph_BDTpredTTGj;
   Float_t ph_etTTGj;
   Float_t ph_R9Full5x5TTGj;
   Float_t ph_S4Full5x5TTGj;
   Float_t ph_EmaxOESCrFull5x5TTGj;
   Float_t ph_E2ndOESCrFull5x5TTGj;
   Float_t ph_E1x3OESCrFull5x5TTGj;
   Float_t ph_E2x5OESCrFull5x5TTGj;
   Float_t ph_E5x5OESCrFull5x5TTGj;
   Float_t ph_2x2OE3x3Full5x5TTGj;
   Float_t ph_SigmaIEtaIEtaTTGj;
   Float_t ph_SigmaIEtaIPhiTTGj;
   Float_t ph_SigmaIPhiIPhiTTGj;
   Float_t ph_SieieOSipipFull5x5TTGj;
   Float_t ph_EtaWidthTTGj;
   Float_t ph_PhiWidthTTGj;
   Float_t ph_EtaWOPhiWFull5x5TTGj;
   Bool_t pass95TTGj;

   Float_t puWeightWW;
   Float_t ph_pFECALClusIsoCorrWW;
   Float_t ph_pFHCALClusIsoCorrWW;
   Float_t ph_TrkIsoCorrWW;
   Float_t ph_HoverEWW;
   Float_t ph_HoverECorrWW;
   Float_t ph_BDTpredWW;
   Float_t ph_etWW;
   Float_t ph_R9Full5x5WW;
   Float_t ph_S4Full5x5WW;
   Float_t ph_EmaxOESCrFull5x5WW;
   Float_t ph_E2ndOESCrFull5x5WW;
   Float_t ph_E1x3OESCrFull5x5WW;
   Float_t ph_E2x5OESCrFull5x5WW;
   Float_t ph_E5x5OESCrFull5x5WW;
   Float_t ph_2x2OE3x3Full5x5WW;
   Float_t ph_SigmaIEtaIEtaWW;
   Float_t ph_SigmaIEtaIPhiWW;
   Float_t ph_SigmaIPhiIPhiWW;
   Float_t ph_SieieOSipipFull5x5WW;
   Float_t ph_EtaWidthWW;
   Float_t ph_PhiWidthWW;
   Float_t ph_EtaWOPhiWFull5x5WW;
   Bool_t pass95WW;*/

   Float_t puWeightWele;
   Float_t ph_pFECALClusIsoCorrWele;
   Float_t ph_pFHCALClusIsoCorrWele;
   Float_t ph_TrkIsoCorrWele;
   Float_t ph_HoverEWele;
   Float_t ph_HoverECorrWele;
   Float_t ph_BDTpredWele;
   Float_t ph_etWele;
   Float_t ph_etaWele;
   Float_t ph_R9Full5x5Wele;
   Float_t ph_S4Full5x5Wele;
   Float_t ph_EmaxOESCrFull5x5Wele;
   Float_t ph_E2ndOESCrFull5x5Wele;
   Float_t ph_E1x3OESCrFull5x5Wele;
   Float_t ph_E2x5OESCrFull5x5Wele;
   Float_t ph_E5x5OESCrFull5x5Wele;
   Float_t ph_2x2OE3x3Full5x5Wele;
   Float_t ph_SigmaIEtaIEtaWele;
   Float_t ph_SigmaIEtaIPhiWele;
   Float_t ph_SigmaIPhiIPhiWele;
   Float_t ph_SieieOSipipFull5x5Wele;
   Float_t ph_EtaWidthWele;
   Float_t ph_PhiWidthWele;
   Float_t ph_EtaWOPhiWFull5x5Wele;
   Bool_t pass95Wele;
   Float_t metWele;
   Float_t FRweightWele;
  
   Float_t puWeightWZ;
   Float_t ph_pFECALClusIsoCorrWZ;
   Float_t ph_pFHCALClusIsoCorrWZ;
   Float_t ph_TrkIsoCorrWZ;
   Float_t ph_HoverEWZ;
   Float_t ph_HoverECorrWZ;
   Float_t ph_BDTpredWZ;
   Float_t ph_etWZ;
   Float_t ph_etaWZ;
   Float_t ph_R9Full5x5WZ;
   Float_t ph_S4Full5x5WZ;
   Float_t ph_EmaxOESCrFull5x5WZ;
   Float_t ph_E2ndOESCrFull5x5WZ;
   Float_t ph_E1x3OESCrFull5x5WZ;
   Float_t ph_E2x5OESCrFull5x5WZ;
   Float_t ph_E5x5OESCrFull5x5WZ;
   Float_t ph_2x2OE3x3Full5x5WZ;
   Float_t ph_SigmaIEtaIEtaWZ;
   Float_t ph_SigmaIEtaIPhiWZ;
   Float_t ph_SigmaIPhiIPhiWZ;
   Float_t ph_SieieOSipipFull5x5WZ;
   Float_t ph_EtaWidthWZ;
   Float_t ph_PhiWidthWZ;
   Float_t ph_EtaWOPhiWFull5x5WZ;
   Bool_t pass95WZ;
   Float_t metWZ;
   Float_t FRweightWZ;

   Float_t puWeightZZ;
   Float_t ph_pFECALClusIsoCorrZZ;
   Float_t ph_pFHCALClusIsoCorrZZ;
   Float_t ph_TrkIsoCorrZZ;
   Float_t ph_HoverEZZ;
   Float_t ph_HoverECorrZZ;
   Float_t ph_BDTpredZZ;
   Float_t ph_etZZ;
   Float_t ph_etaZZ;
   Float_t ph_R9Full5x5ZZ;
   Float_t ph_S4Full5x5ZZ;
   Float_t ph_EmaxOESCrFull5x5ZZ;
   Float_t ph_E2ndOESCrFull5x5ZZ;
   Float_t ph_E1x3OESCrFull5x5ZZ;
   Float_t ph_E2x5OESCrFull5x5ZZ;
   Float_t ph_E5x5OESCrFull5x5ZZ;
   Float_t ph_2x2OE3x3Full5x5ZZ;
   Float_t ph_SigmaIEtaIEtaZZ;
   Float_t ph_SigmaIEtaIPhiZZ;
   Float_t ph_SigmaIPhiIPhiZZ;
   Float_t ph_SieieOSipipFull5x5ZZ;
   Float_t ph_EtaWidthZZ;
   Float_t ph_PhiWidthZZ;
   Float_t ph_EtaWOPhiWFull5x5ZZ;
   Bool_t pass95ZZ;
   Float_t metZZ;
   Float_t FRweightZZ;

   data_tree->SetBranchAddress("puWeight",&puWeight);
   data_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorr);
   data_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorr);
   data_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorr);
   data_tree->SetBranchAddress("ph_hoe",&ph_HoverE);
   data_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorr);
   data_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpred);
   data_tree->SetBranchAddress("ph_et",&ph_et);
   data_tree->SetBranchAddress("ph_eta",&ph_eta);
   data_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5);
   data_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5);
   data_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5);
   data_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5);
   data_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5);
   data_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5);
   data_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5);
   data_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5);
   data_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEta);
   data_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhi);
   data_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhi);
   data_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5);
   data_tree->SetBranchAddress("ph_se", &ph_EtaWidth);
   data_tree->SetBranchAddress("ph_sp", &ph_PhiWidth);
   data_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5);
   data_tree->SetBranchAddress("pass95", &pass95);
   data_tree->SetBranchAddress("FRweight", &FRweight);
   data_tree->SetBranchAddress("met", &met);

   signalMC_tree->SetBranchAddress("puWeight",&puWeightMC);
   signalMC_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrMC);
   signalMC_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrMC);
   signalMC_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrMC);
   signalMC_tree->SetBranchAddress("ph_hoe",&ph_HoverEMC);
   signalMC_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrMC);
   signalMC_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredMC); 
   signalMC_tree->SetBranchAddress("ph_et",&ph_etMC);
   signalMC_tree->SetBranchAddress("ph_eta",&ph_etaMC);
   signalMC_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5MC);
   signalMC_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5MC);
   signalMC_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5MC);
   signalMC_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5MC);
   signalMC_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5MC);
   signalMC_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5MC);
   signalMC_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5MC);
   signalMC_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5MC);
   signalMC_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaMC);
   signalMC_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiMC);
   signalMC_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiMC);
   signalMC_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5MC);
   signalMC_tree->SetBranchAddress("ph_se", &ph_EtaWidthMC);
   signalMC_tree->SetBranchAddress("ph_sp", &ph_PhiWidthMC);
   signalMC_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5MC);
   signalMC_tree->SetBranchAddress("pass95", &pass95MC);
   signalMC_tree->SetBranchAddress("FRweight", &FRweightMC);
   signalMC_tree->SetBranchAddress("met", &metMC);

   /*gjets_tree->SetBranchAddress("puWeight",&puWeightGj);
   gjets_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrGj);
   gjets_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrGj);
   gjets_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrGj);
   gjets_tree->SetBranchAddress("ph_hoe",&ph_HoverEGj);
   gjets_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrGj);
   gjets_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredGj);
   gjets_tree->SetBranchAddress("ph_et",&ph_etGj);
   gjets_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5Gj);
   gjets_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5Gj);
   gjets_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5Gj);
   gjets_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5Gj);
   gjets_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5Gj);
   gjets_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5Gj);
   gjets_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5Gj);
   gjets_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5Gj);
   gjets_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaGj);
   gjets_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiGj);
   gjets_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiGj);
   gjets_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5Gj);
   gjets_tree->SetBranchAddress("ph_se", &ph_EtaWidthGj);
   gjets_tree->SetBranchAddress("ph_sp", &ph_PhiWidthGj);
   gjets_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5Gj);
   gjets_tree->SetBranchAddress("pass95", &pass95Gj);

   wjets_tree->SetBranchAddress("puWeight",&puWeightWj);
   wjets_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrWj);
   wjets_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrWj);
   wjets_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrWj);
   wjets_tree->SetBranchAddress("ph_hoe",&ph_HoverEWj);
   wjets_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrWj);
   wjets_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredWj);
   wjets_tree->SetBranchAddress("ph_et",&ph_etWj);
   wjets_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5Wj);
   wjets_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5Wj);
   wjets_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5Wj);
   wjets_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5Wj);
   wjets_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5Wj);
   wjets_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5Wj);
   wjets_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5Wj);
   wjets_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5Wj);
   wjets_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaWj);
   wjets_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiWj);
   wjets_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiWj);
   wjets_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5Wj);
   wjets_tree->SetBranchAddress("ph_se", &ph_EtaWidthWj);
   wjets_tree->SetBranchAddress("ph_sp", &ph_PhiWidthWj);
   wjets_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5Wj);
   wjets_tree->SetBranchAddress("pass95", &pass95Wj);

   wlnu_tree->SetBranchAddress("puWeight",&puWeightWlnu);
   wlnu_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrWlnu);
   wlnu_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrWlnu);
   wlnu_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrWlnu);
   wlnu_tree->SetBranchAddress("ph_hoe",&ph_HoverEWlnu);
   wlnu_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrWlnu);
   wlnu_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredWlnu);
   wlnu_tree->SetBranchAddress("ph_et",&ph_etWlnu);
   wlnu_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaWlnu);
   wlnu_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiWlnu);
   wlnu_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiWlnu);
   wlnu_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5Wlnu);
   wlnu_tree->SetBranchAddress("ph_se", &ph_EtaWidthWlnu);
   wlnu_tree->SetBranchAddress("ph_sp", &ph_PhiWidthWlnu);
   wlnu_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5Wlnu);
   wlnu_tree->SetBranchAddress("pass95", &pass95Wlnu);

   qcd_tree->SetBranchAddress("puWeight",&puWeightQcd);
   qcd_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrQcd);
   qcd_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrQcd);
   qcd_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrQcd);
   qcd_tree->SetBranchAddress("ph_hoe",&ph_HoverEQcd);
   qcd_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrQcd);
   qcd_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredQcd);
   qcd_tree->SetBranchAddress("ph_et",&ph_etQcd);
   qcd_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5Qcd);
   qcd_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5Qcd);
   qcd_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5Qcd);
   qcd_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5Qcd);
   qcd_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5Qcd);
   qcd_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5Qcd);
   qcd_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5Qcd);
   qcd_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5Qcd);
   qcd_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaQcd);
   qcd_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiQcd);
   qcd_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiQcd);
   qcd_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5Qcd);
   qcd_tree->SetBranchAddress("ph_se", &ph_EtaWidthQcd);
   qcd_tree->SetBranchAddress("ph_sp", &ph_PhiWidthQcd);
   qcd_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5Qcd);
   qcd_tree->SetBranchAddress("pass95", &pass95Qcd);

   ttgjets_tree->SetBranchAddress("puWeight",&puWeightTTGj);
   ttgjets_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrTTGj);
   ttgjets_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrTTGj);
   ttgjets_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrTTGj);
   ttgjets_tree->SetBranchAddress("ph_hoe",&ph_HoverETTGj);
   ttgjets_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrTTGj);
   ttgjets_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredTTGj);
   ttgjets_tree->SetBranchAddress("ph_et",&ph_etTTGj);
   ttgjets_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaTTGj);
   ttgjets_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiTTGj);
   ttgjets_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiTTGj);
   ttgjets_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5TTGj);
   ttgjets_tree->SetBranchAddress("ph_se", &ph_EtaWidthTTGj);
   ttgjets_tree->SetBranchAddress("ph_sp", &ph_PhiWidthTTGj);
   ttgjets_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5TTGj);
   ttgjets_tree->SetBranchAddress("pass95", &pass95TTGj);

   ww_tree->SetBranchAddress("puWeight",&puWeightWW);
   ww_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrWW);
   ww_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrWW);
   ww_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrWW);
   ww_tree->SetBranchAddress("ph_hoe",&ph_HoverEWW);
   ww_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrWW);
   ww_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredWW);
   ww_tree->SetBranchAddress("ph_et",&ph_etWW);
   ww_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5WW);
   ww_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5WW);
   ww_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5WW);
   ww_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5WW);
   ww_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5WW);
   ww_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5WW);
   ww_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5WW);
   ww_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5WW);
   ww_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaWW);
   ww_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiWW);
   ww_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiWW);
   ww_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5WW);
   ww_tree->SetBranchAddress("ph_se", &ph_EtaWidthWW);
   ww_tree->SetBranchAddress("ph_sp", &ph_PhiWidthWW);
   ww_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5WW);
   ww_tree->SetBranchAddress("pass95", &pass95WW);*/

   wele_tree->SetBranchAddress("puWeight",&puWeightWele);
   wele_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrWele);
   wele_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrWele);
   wele_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrWele);
   wele_tree->SetBranchAddress("ph_hoe",&ph_HoverEWele);
   wele_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrWele);
   wele_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredWele);
   wele_tree->SetBranchAddress("ph_et",&ph_etWele);
   wele_tree->SetBranchAddress("ph_eta",&ph_etaWele);
   wele_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5Wele);
   wele_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5Wele);
   wele_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5Wele);
   wele_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5Wele);
   wele_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5Wele);
   wele_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5Wele);
   wele_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5Wele);
   wele_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5Wele);
   wele_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaWele);
   wele_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiWele);
   wele_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiWele);
   wele_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5Wele);
   wele_tree->SetBranchAddress("ph_se", &ph_EtaWidthWele);
   wele_tree->SetBranchAddress("ph_sp", &ph_PhiWidthWele);
   wele_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5Wele);
   wele_tree->SetBranchAddress("pass95", &pass95Wele);
   wele_tree->SetBranchAddress("met", &metWele);
   wele_tree->SetBranchAddress("FRweight", &FRweightWele);
 
   wz_tree->SetBranchAddress("puWeight",&puWeightWZ);
   wz_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrWZ);
   wz_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrWZ);
   wz_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrWZ);
   wz_tree->SetBranchAddress("ph_hoe",&ph_HoverEWZ);
   wz_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrWZ);
   wz_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredWZ);
   wz_tree->SetBranchAddress("ph_et",&ph_etWZ);
   wz_tree->SetBranchAddress("ph_eta",&ph_etaWZ);
   wz_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5WZ);
   wz_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5WZ);
   wz_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5WZ);
   wz_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5WZ);
   wz_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5WZ);
   wz_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5WZ);
   wz_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5WZ);
   wz_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5WZ);
   wz_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaWZ);
   wz_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiWZ);
   wz_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiWZ);
   wz_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5WZ);
   wz_tree->SetBranchAddress("ph_se", &ph_EtaWidthWZ);
   wz_tree->SetBranchAddress("ph_sp", &ph_PhiWidthWZ);
   wz_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5WZ);
   wz_tree->SetBranchAddress("pass95", &pass95WZ);
   wz_tree->SetBranchAddress("met", &metWZ);
   wz_tree->SetBranchAddress("FRweight", &FRweightWZ);
 
   zz_tree->SetBranchAddress("puWeight",&puWeightZZ);
   zz_tree->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_pFECALClusIsoCorrZZ);
   zz_tree->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_pFHCALClusIsoCorrZZ);
   zz_tree->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorrZZ);
   zz_tree->SetBranchAddress("ph_hoe",&ph_HoverEZZ);
   zz_tree->SetBranchAddress("ph_hoeCorr",&ph_HoverECorrZZ);
   zz_tree->SetBranchAddress("ph_BDTpred",&ph_BDTpredZZ);
   zz_tree->SetBranchAddress("ph_et",&ph_etZZ);
   zz_tree->SetBranchAddress("ph_eta",&ph_etaZZ);
   zz_tree->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5ZZ);
   zz_tree->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5ZZ);
   zz_tree->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5ZZ);
   zz_tree->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5ZZ);
   zz_tree->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5ZZ);
   zz_tree->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5ZZ);
   zz_tree->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5ZZ);
   zz_tree->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5ZZ);
   zz_tree->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEtaZZ);
   zz_tree->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhiZZ);
   zz_tree->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhiZZ);
   zz_tree->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5ZZ);
   zz_tree->SetBranchAddress("ph_se", &ph_EtaWidthZZ);
   zz_tree->SetBranchAddress("ph_sp", &ph_PhiWidthZZ);
   zz_tree->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5ZZ);
   zz_tree->SetBranchAddress("pass95", &pass95ZZ);
   zz_tree->SetBranchAddress("met", &metZZ);
   zz_tree->SetBranchAddress("FRweight", &FRweightZZ);

   /*TCanvas *c1 = new TCanvas("c1", "c1");
   c1->cd();

   TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
   pad1->SetBottomMargin(0);
   pad1->Draw();
   pad1->cd();*/

   const Int_t nhist=24;
   TH1F* histdata[nhist];
   TH1F* histMC[nhist];
   /*TH1F* histgjets[nhist];
   TH1F* histwjets[nhist];
   TH1F* histwlnu[nhist];
   TH1F* histqcd[nhist];
   TH1F* histttgjets[nhist];
   TH1F* histww[nhist];*/
   TH1F* histwele[nhist];
   TH1F* histwz[nhist]; 
   TH1F* histzz[nhist];
 
   THStack *hs[nhist]; 
   
   histdata[0]   = new TH1F("PFECALClusIsoCorr","PFECALClusIsoCorr",15,0,30);
   histdata[1]   = new TH1F("PFHCALClusIsoCorr","PFHCALClusIsoCorr",20,0,40);
   histdata[2]   = new TH1F("TrkIsoCorr","TrkIsoCorr",10,0,20);
   histdata[3]   = new TH1F("phoBDT","phoBDT",25,0,1); //original plot 50, 0, 1
   histdata[4]   = new TH1F("phoEt","phoEt",10,200,1000);
   histdata[5]   = new TH1F("phoR9Full5x5", "phoR9Full5x5",100, 0.1, 1.5 );
   histdata[6]   = new TH1F ("phoS4Full5x5", "phoS4Full5x5", 100, 0.4, 1.1);
   histdata[7]   = new TH1F ("phoEmaxOESCrFull5x5", "phoEmaxOESCrFull5x5", 150, 0, 1);
   histdata[8]   = new TH1F ("phoE2ndOESCrFull5x5", "phoE2ndOESCrFull5x5", 150, 0, 0.6);
   histdata[9]   = new TH1F ("phoE1x3OESCrFull5x5", "phoE1x3OESCrFull5x5", 150, 0.2, 1.1);
   histdata[10]   = new TH1F ("phoE2x5OESCrFull5x5", "phoE2x5OESCrFull5x5", 150, 0.6, 1.2);
   histdata[11]   = new TH1F ("phoE5x5OESCrFull5x5", "phoE5x5OESCrFull5x5", 150, 0.6, 1.2);
   histdata[12]   = new TH1F ("pho2x2OE3x3Full5x5", "pho2x2OE3x3Full5x5", 100, 0.5, 1);
   histdata[13]   = new TH1F ("phoSigmaIEtaIEta", "phoSigmaIEtaIEta", 100, 0, 0.04);
   histdata[14]   = new TH1F ("phoSigmaIEtaIPhi", "phoSigmaIEtaIPhi", 100, -0.001, 0.001);
   histdata[15]   = new TH1F ("phoSigmaIPhiIPhi", "phoSigmaIPhiIPhi", 100, 0, 0.07);
   histdata[16]   = new TH1F ("phoSieieOSipipFull5x5", "phoSieieOSipipFull5x5",150, 0, 1.5);
   histdata[17]   = new TH1F ("phoPhiWidth", "phoPhiWidth", 100, 0., 0.06);
   histdata[18]   = new TH1F ("phoEtaWidth", "phoEtaWidth", 100, 0, 0.02);
   histdata[19]   = new TH1F ("phoEtaWOPhiWFull5x5", "phoEtaWOPhiWFull5x5", 150, 0, 1.5);
   histdata[20]   = new TH1F ("phoHovE", "phoHovE", 20, 0, 0.05);
   histdata[21]   = new TH1F ("phoHovECorr", "phoHovECorr", 20, 0, 0.05);
   histdata[22]   = new TH1F ("met", "met", 10, 150, 1000);
   histdata[23]   = new TH1F ("phoEta","phoEta",10,-1.4442,1.4442);

   gStyle->SetOptStat(false);
   histMC[0]   = new TH1F("PFECALClusIsoCorrMC","PFECALClusIsoCorrMC",15,0,30);
   histMC[1]   = new TH1F("PFHCALClusIsoCorrMC","PFHCALClusIsoCorrMC",20,0,40);
   histMC[2]   = new TH1F("TrkIsoCorrMC","TrkIsoCorrMC",10,0,20);
   histMC[3]   = new TH1F("phoBDTMC","phoBDTMC",25,0,1); //original plot 50, 0, 1
   histMC[4]   = new TH1F("phoEtMC","phoEtMC",10,200,1000);
   histMC[5]   = new TH1F("phoR9Full5x5MC", "phoR9Full5x5MC",100, 0.1, 1.5 );
   histMC[6]   = new TH1F ("phoS4Full5x5MC", "phoS4Full5x5MC", 100, 0.4, 1.1);
   histMC[7]   = new TH1F ("phoEmaxOESCrFull5x5MC", "phoEmaxOESCrFull5x5MC", 150, 0, 1);
   histMC[8]   = new TH1F ("phoE2ndOESCrFull5x5MC", "phoE2ndOESCrFull5x5MC", 150, 0, 0.6);
   histMC[9]   = new TH1F ("phoE1x3OESCrFull5x5MC", "phoE1x3OESCrFull5x5MC", 150, 0.2, 1.1);
   histMC[10]   = new TH1F ("phoE2x5OESCrFull5x5MC", "phoE2x5OESCrFull5x5MC", 150, 0.6, 1.2);
   histMC[11]   = new TH1F ("phoE5x5OESCrFull5x5MC", "phoE5x5OESCrFull5x5MC", 150, 0.6, 1.2);
   histMC[12]   = new TH1F ("pho2x2OE3x3Full5x5MC", "pho2x2OE3x3Full5x5MC", 100, 0.5, 1);
   histMC[13]   = new TH1F ("phoSigmaIEtaIEtaMC", "phoSigmaIEtaIEtaMC", 100, 0, 0.04);
   histMC[14]   = new TH1F ("phoSigmaIEtaIPhiMC", "phoSigmaIEtaIPhiMC", 100, -0.001, 0.001);
   histMC[15]   = new TH1F ("phoSigmaIPhiIPhiMC", "phoSigmaIPhiIPhiMC", 100, 0, 0.07);
   histMC[16]   = new TH1F ("phoSieieOSipipFull5x5MC", "phoSieieOSipipFull5x5MC",150, 0, 1.5);
   histMC[17]   = new TH1F ("phoPhiWidthMC", "phoPhiWidthMC", 100, 0., 0.06);
   histMC[18]   = new TH1F ("phoEtaWidthMC", "phoEtaWidthMC", 100, 0, 0.02);
   histMC[19]   = new TH1F ("phoEtaWOPhiWFull5x5MC", "phoEtaWOPhiWFull5x5MC", 150, 0, 1.5);  
   histMC[20]   = new TH1F ("phoHovEMC", "phoHovEMC", 20, 0, 0.05);
   histMC[21]   = new TH1F ("phoHovECorrMC", "phoHovECorrMC", 20, 0, 0.05);
   histMC[22]   = new TH1F ("metMC", "metMC", 10, 150, 1000);
   histMC[23]   = new TH1F ("phoEtaMC","phoEtaMC",10,-1.4442,1.4442);

   /*histgjets[0]   = new TH1F("PFECALClusIsoCorrGj","PFECALClusIsoCorrGj",15,0,30);
   histgjets[1]   = new TH1F("PFHCALClusIsoCorrGj","PFHCALClusIsoCorrGj",20,0,40);
   histgjets[2]   = new TH1F("TrkIsoCorrGj","TrkIsoCorrGj",10,0,20);
   histgjets[3]   = new TH1F("phoBDTGj","phoBDTGj",25,0,1);
   histgjets[4]   = new TH1F("phoEtGj","phoEtGj",30,230,1000);
   histgjets[5]   = new TH1F("phoR9Full5x5Gj", "phoR9Full5x5Gj",100, 0.1, 1.5 );
   histgjets[6]   = new TH1F ("phoS4Full5x5Gj", "phoS4Full5x5Gj", 100, 0.4, 1.1);
   histgjets[7]   = new TH1F ("phoEmaxOESCrFull5x5Gj", "phoEmaxOESCrFull5x5Gj", 150, 0, 1);
   histgjets[8]   = new TH1F ("phoE2ndOESCrFull5x5Gj", "phoE2ndOESCrFull5x5Gj", 150, 0, 0.6);
   histgjets[9]   = new TH1F ("phoE1x3OESCrFull5x5Gj", "phoE1x3OESCrFull5x5Gj", 150, 0.2, 1.1);
   histgjets[10]   = new TH1F ("phoE2x5OESCrFull5x5Gj", "phoE2x5OESCrFull5x5Gj", 150, 0.6, 1.2);
   histgjets[11]   = new TH1F ("phoE5x5OESCrFull5x5Gj", "phoE5x5OESCrFull5x5Gj", 150, 0.6, 1.2);
   histgjets[12]   = new TH1F ("pho2x2OE3x3Full5x5Gj", "pho2x2OE3x3Full5x5Gj", 100, 0.5, 1);
   histgjets[13]   = new TH1F ("phoSigmaIEtaIEtaGj", "phoSigmaIEtaIEtaGj", 100, 0, 0.04);
   histgjets[14]   = new TH1F ("phoSigmaIEtaIPhiGj", "phoSigmaIEtaIPhiGj", 100, -0.001, 0.001);
   histgjets[15]   = new TH1F ("phoSigmaIPhiIPhiGj", "phoSigmaIPhiIPhiGj", 100, 0, 0.07);
   histgjets[16]   = new TH1F ("phoSieieOSipipFull5x5Gj", "phoSieieOSipipFull5x5Gj",150, 0, 1.5);
   histgjets[17]   = new TH1F ("phoPhiWidthGj", "phoPhiWidthGj", 100, 0., 0.06);
   histgjets[18]   = new TH1F ("phoEtaWidthGj", "phoEtaWidthGj", 100, 0, 0.02);
   histgjets[19]   = new TH1F ("phoEtaWOPhiWFull5x5Gj", "phoEtaWOPhiWFull5x5Gj", 150, 0, 1.5);
   histgjets[20]   = new TH1F ("phoHovEGj", "phoHovEGj", 20, 0, 0.05);
   histgjets[21]   = new TH1F ("phoHovECorrGj", "phoHovECorrGj", 20, 0, 0.05);

   histwjets[0]   = new TH1F("PFECALClusIsoCorrWj","PFECALClusIsoCorrWj",15,0,30);
   histwjets[1]   = new TH1F("PFHCALClusIsoCorrWj","PFHCALClusIsoCorrWj",20,0,40);
   histwjets[2]   = new TH1F("TrkIsoCorrWj","TrkIsoCorrWj",10,0,20);
   histwjets[3]   = new TH1F("phoBDTWj","phoBDTWj",25,0,1);
   histwjets[4]   = new TH1F("phoEtWj","phoEtWj",30,230,1000);
   histwjets[5]   = new TH1F("phoR9Full5x5Wj", "phoR9Full5x5Wj",100, 0.1, 1.5 );
   histwjets[6]   = new TH1F ("phoS4Full5x5Wj", "phoS4Full5x5Wj", 100, 0.4, 1.1);
   histwjets[7]   = new TH1F ("phoEmaxOESCrFull5x5Wj", "phoEmaxOESCrFull5x5Wj", 150, 0, 1);
   histwjets[8]   = new TH1F ("phoE2ndOESCrFull5x5Wj", "phoE2ndOESCrFull5x5Wj", 150, 0, 0.6);
   histwjets[9]   = new TH1F ("phoE1x3OESCrFull5x5Wj", "phoE1x3OESCrFull5x5Wj", 150, 0.2, 1.1);
   histwjets[10]   = new TH1F ("phoE2x5OESCrFull5x5Wj", "phoE2x5OESCrFull5x5Wj", 150, 0.6, 1.2);
   histwjets[11]   = new TH1F ("phoE5x5OESCrFull5x5Wj", "phoE5x5OESCrFull5x5Wj", 150, 0.6, 1.2);
   histwjets[12]   = new TH1F ("pho2x2OE3x3Full5x5Wj", "pho2x2OE3x3Full5x5Wj", 100, 0.5, 1);
   histwjets[13]   = new TH1F ("phoSigmaIEtaIEtaWj", "phoSigmaIEtaIEtaWj", 100, 0, 0.04);
   histwjets[14]   = new TH1F ("phoSigmaIEtaIPhiWj", "phoSigmaIEtaIPhiWj", 100, -0.001, 0.001);
   histwjets[15]   = new TH1F ("phoSigmaIPhiIPhiWj", "phoSigmaIPhiIPhiWj", 100, 0, 0.07);
   histwjets[16]   = new TH1F ("phoSieieOSipipFull5x5Wj", "phoSieieOSipipFull5x5Wj",150, 0, 1.5);
   histwjets[17]   = new TH1F ("phoPhiWidthWj", "phoPhiWidthWj", 100, 0., 0.06);
   histwjets[18]   = new TH1F ("phoEtaWidthWj", "phoEtaWidthWj", 100, 0, 0.02);
   histwjets[19]   = new TH1F ("phoEtaWOPhiWFull5x5Wj", "phoEtaWOPhiWFull5x5Wj", 150, 0, 1.5);
   histwjets[20]   = new TH1F ("phoHovEWj", "phoHovEWj", 20, 0, 0.05);
   histwjets[21]   = new TH1F ("phoHovECorrWj", "phoHovECorrWj", 20, 0, 0.05);

   histwlnu[0]   = new TH1F("PFECALClusIsoCorrWlnu","PFECALClusIsoCorrWlnu",15,0,30);
   histwlnu[1]   = new TH1F("PFHCALClusIsoCorrWlnu","PFHCALClusIsoCorrWlnu",20,0,40);
   histwlnu[2]   = new TH1F("TrkIsoCorrWlnu","TrkIsoCorrWlnu",10,0,20);
   histwlnu[3]   = new TH1F("phoBDTWlnu","phoBDTWlnu",25,0,1);
   histwlnu[4]   = new TH1F("phoEtWlnu","phoEtWlnu",30,230,1000);
   histwlnu[5]   = new TH1F("phoR9Full5x5Wlnu", "phoR9Full5x5Wlnu",100, 0.1, 1.5 );
   histwlnu[6]   = new TH1F ("phoS4Full5x5Wlnu", "phoS4Full5x5Wlnu", 100, 0.4, 1.1);
   histwlnu[7]   = new TH1F ("phoEmaxOESCrFull5x5Wlnu", "phoEmaxOESCrFull5x5Wlnu", 150, 0, 1);
   histwlnu[8]   = new TH1F ("phoE2ndOESCrFull5x5Wlnu", "phoE2ndOESCrFull5x5Wlnu", 150, 0, 0.6);
   histwlnu[9]   = new TH1F ("phoE1x3OESCrFull5x5Wlnu", "phoE1x3OESCrFull5x5Wlnu", 150, 0.2, 1.1);
   histwlnu[10]   = new TH1F ("phoE2x5OESCrFull5x5Wlnu", "phoE2x5OESCrFull5x5Wlnu", 150, 0.6, 1.2);
   histwlnu[11]   = new TH1F ("phoE5x5OESCrFull5x5Wlnu", "phoE5x5OESCrFull5x5Wlnu", 150, 0.6, 1.2);
   histwlnu[12]   = new TH1F ("pho2x2OE3x3Full5x5Wlnu", "pho2x2OE3x3Full5x5Wlnu", 100, 0.5, 1);
   histwlnu[13]   = new TH1F ("phoSigmaIEtaIEtaWlnu", "phoSigmaIEtaIEtaWlnu", 100, 0, 0.04);
   histwlnu[14]   = new TH1F ("phoSigmaIEtaIPhiWlnu", "phoSigmaIEtaIPhiWlnu", 100, -0.001, 0.001);
   histwlnu[15]   = new TH1F ("phoSigmaIPhiIPhiWlnu", "phoSigmaIPhiIPhiWlnu", 100, 0, 0.07);
   histwlnu[16]   = new TH1F ("phoSieieOSipipFull5x5Wlnu", "phoSieieOSipipFull5x5Wlnu",150, 0, 1.5);
   histwlnu[17]   = new TH1F ("phoPhiWidthWlnu", "phoPhiWidthWlnu", 100, 0., 0.06);
   histwlnu[18]   = new TH1F ("phoEtaWidthWlnu", "phoEtaWidthWlnu", 100, 0, 0.02);
   histwlnu[19]   = new TH1F ("phoEtaWOPhiWFull5x5Wlnu", "phoEtaWOPhiWFull5x5Wlnu", 150, 0, 1.5);
   histwlnu[20]   = new TH1F ("phoHovEWlnu", "phoHovEWlnu", 20, 0, 0.05);
   histwlnu[21]   = new TH1F ("phoHovECorrWlnu", "phoHovECorrWlnu", 20, 0, 0.05);

   histqcd[0]   = new TH1F("PFECALClusIsoCorrQcd","PFECALClusIsoCorrQcd",15,0,30);
   histqcd[1]   = new TH1F("PFHCALClusIsoCorrQcd","PFHCALClusIsoCorrQcd",20,0,40);
   histqcd[2]   = new TH1F("TrkIsoCorrQcd","TrkIsoCorrQcd",10,0,20);
   histqcd[3]   = new TH1F("phoBDTQcd","phoBDTQcd",25,0,1);
   histqcd[4]   = new TH1F("phoEtQcd","phoEtQcd",30,230,1000);
   histqcd[5]   = new TH1F("phoR9Full5x5Qcd", "phoR9Full5x5Qcd",100, 0.1, 1.5 );
   histqcd[6]   = new TH1F ("phoS4Full5x5Qcd", "phoS4Full5x5Qcd", 100, 0.4, 1.1);
   histqcd[7]   = new TH1F ("phoEmaxOESCrFull5x5Qcd", "phoEmaxOESCrFull5x5Qcd", 150, 0, 1);
   histqcd[8]   = new TH1F ("phoE2ndOESCrFull5x5Qcd", "phoE2ndOESCrFull5x5Qcd", 150, 0, 0.6);
   histqcd[9]   = new TH1F ("phoE1x3OESCrFull5x5Qcd", "phoE1x3OESCrFull5x5Qcd", 150, 0.2, 1.1);
   histqcd[10]   = new TH1F ("phoE2x5OESCrFull5x5Qcd", "phoE2x5OESCrFull5x5Qcd", 150, 0.6, 1.2);
   histqcd[11]   = new TH1F ("phoE5x5OESCrFull5x5Qcd", "phoE5x5OESCrFull5x5Qcd", 150, 0.6, 1.2);
   histqcd[12]   = new TH1F ("pho2x2OE3x3Full5x5Qcd", "pho2x2OE3x3Full5x5Qcd", 100, 0.5, 1);
   histqcd[13]   = new TH1F ("phoSigmaIEtaIEtaQcd", "phoSigmaIEtaIEtaQcd", 100, 0, 0.04);
   histqcd[14]   = new TH1F ("phoSigmaIEtaIPhiQcd", "phoSigmaIEtaIPhiQcd", 100, -0.001, 0.001);
   histqcd[15]   = new TH1F ("phoSigmaIPhiIPhiQcd", "phoSigmaIPhiIPhiQcd", 100, 0, 0.07);
   histqcd[16]   = new TH1F ("phoSieieOSipipFull5x5Qcd", "phoSieieOSipipFull5x5Qcd",150, 0, 1.5);
   histqcd[17]   = new TH1F ("phoPhiWidthQcd", "phoPhiWidthQcd", 100, 0., 0.06);
   histqcd[18]   = new TH1F ("phoEtaWidthQcd", "phoEtaWidthQcd", 100, 0, 0.02);
   histqcd[19]   = new TH1F ("phoEtaWOPhiWFull5x5Qcd", "phoEtaWOPhiWFull5x5Qcd", 150, 0, 1.5);
   histqcd[20]   = new TH1F ("phoHovEQcd", "phoHovEQcd", 20, 0, 0.05);
   histqcd[21]   = new TH1F ("phoHovECorrQcd", "phoHovECorrQcd", 20, 0, 0.05);

   histttgjets[0]   = new TH1F("PFECALClusIsoCorrTTGj","PFECALClusIsoCorrTTGj",15,0,30);
   histttgjets[1]   = new TH1F("PFHCALClusIsoCorrTTGj","PFHCALClusIsoCorrTTGj",20,0,40);
   histttgjets[2]   = new TH1F("TrkIsoCorrTTGj","TrkIsoCorrTTGj",10,0,20);
   histttgjets[3]   = new TH1F("phoBDTTTGj","phoBDTTTGj",25,0,1);
   histttgjets[4]   = new TH1F("phoEtTTGj","phoEtTTGj",30,230,1000);
   histttgjets[5]   = new TH1F("phoR9Full5x5TTGj", "phoR9Full5x5TTGj",100, 0.1, 1.5 );
   histttgjets[6]   = new TH1F ("phoS4Full5x5TTGj", "phoS4Full5x5TTGj", 100, 0.4, 1.1);
   histttgjets[7]   = new TH1F ("phoEmaxOESCrFull5x5TTGj", "phoEmaxOESCrFull5x5TTGj", 150, 0, 1);
   histttgjets[8]   = new TH1F ("phoE2ndOESCrFull5x5TTGj", "phoE2ndOESCrFull5x5TTGj", 150, 0, 0.6);
   histttgjets[9]   = new TH1F ("phoE1x3OESCrFull5x5TTGj", "phoE1x3OESCrFull5x5TTGj", 150, 0.2, 1.1);
   histttgjets[10]   = new TH1F ("phoE2x5OESCrFull5x5TTGj", "phoE2x5OESCrFull5x5TTGj", 150, 0.6, 1.2);
   histttgjets[11]   = new TH1F ("phoE5x5OESCrFull5x5TTGj", "phoE5x5OESCrFull5x5TTGj", 150, 0.6, 1.2);
   histttgjets[12]   = new TH1F ("pho2x2OE3x3Full5x5TTGj", "pho2x2OE3x3Full5x5TTGj", 100, 0.5, 1);
   histttgjets[13]   = new TH1F ("phoSigmaIEtaIEtaTTGj", "phoSigmaIEtaIEtaTTGj", 100, 0, 0.04);
   histttgjets[14]   = new TH1F ("phoSigmaIEtaIPhiTTGj", "phoSigmaIEtaIPhiTTGj", 100, -0.001, 0.001);
   histttgjets[15]   = new TH1F ("phoSigmaIPhiIPhiTTGj", "phoSigmaIPhiIPhiTTGj", 100, 0, 0.07);
   histttgjets[16]   = new TH1F ("phoSieieOSipipFull5x5TTGj", "phoSieieOSipipFull5x5TTGj",150, 0, 1.5);
   histttgjets[17]   = new TH1F ("phoPhiWidthTTGj", "phoPhiWidthTTGj", 100, 0., 0.06);
   histttgjets[18]   = new TH1F ("phoEtaWidthTTGj", "phoEtaWidthTTGj", 100, 0, 0.02);
   histttgjets[19]   = new TH1F ("phoEtaWOPhiWFull5x5TTGj", "phoEtaWOPhiWFull5x5TTGj", 150, 0, 1.5);
   histttgjets[20]   = new TH1F ("phoHovETTGj", "phoHovETTGj", 20, 0, 0.05);
   histttgjets[21]   = new TH1F ("phoHovECorrTTGj", "phoHovECorrTTGj", 20, 0, 0.05);

   histww[0]   = new TH1F("PFECALClusIsoCorrWW","PFECALClusIsoCorrWW",15,0,30);
   histww[1]   = new TH1F("PFHCALClusIsoCorrWW","PFHCALClusIsoCorrWW",20,0,40);
   histww[2]   = new TH1F("TrkIsoCorrWW","TrkIsoCorrWW",10,0,20);
   histww[3]   = new TH1F("phoBDTWW","phoBDTWW",25,0,1);
   histww[4]   = new TH1F("phoEtWW","phoEtWW",30,230,1000);
   histww[5]   = new TH1F("phoR9Full5x5WW", "phoR9Full5x5WW",100, 0.1, 1.5 );
   histww[6]   = new TH1F ("phoS4Full5x5WW", "phoS4Full5x5WW", 100, 0.4, 1.1);
   histww[7]   = new TH1F ("phoEmaxOESCrFull5x5WW", "phoEmaxOESCrFull5x5WW", 150, 0, 1);
   histww[8]   = new TH1F ("phoE2ndOESCrFull5x5WW", "phoE2ndOESCrFull5x5WW", 150, 0, 0.6);
   histww[9]   = new TH1F ("phoE1x3OESCrFull5x5WW", "phoE1x3OESCrFull5x5WW", 150, 0.2, 1.1);
   histww[10]   = new TH1F ("phoE2x5OESCrFull5x5WW", "phoE2x5OESCrFull5x5WW", 150, 0.6, 1.2);
   histww[11]   = new TH1F ("phoE5x5OESCrFull5x5WW", "phoE5x5OESCrFull5x5WW", 150, 0.6, 1.2);
   histww[12]   = new TH1F ("pho2x2OE3x3Full5x5WW", "pho2x2OE3x3Full5x5WW", 100, 0.5, 1);
   histww[13]   = new TH1F ("phoSigmaIEtaIEtaWW", "phoSigmaIEtaIEtaWW", 100, 0, 0.04);
   histww[14]   = new TH1F ("phoSigmaIEtaIPhiWW", "phoSigmaIEtaIPhiWW", 100, -0.001, 0.001);
   histww[15]   = new TH1F ("phoSigmaIPhiIPhiWW", "phoSigmaIPhiIPhiWW", 100, 0, 0.07);
   histww[16]   = new TH1F ("phoSieieOSipipFull5x5WW", "phoSieieOSipipFull5x5WW",150, 0, 1.5);
   histww[17]   = new TH1F ("phoPhiWidthWW", "phoPhiWidthWW", 100, 0., 0.06);
   histww[18]   = new TH1F ("phoEtaWidthWW", "phoEtaWidthWW", 100, 0, 0.02);
   histww[19]   = new TH1F ("phoEtaWOPhiWFull5x5WW", "phoEtaWOPhiWFull5x5WW", 150, 0, 1.5);
   histww[20]   = new TH1F ("phoHovEWW", "phoHovEWW", 20, 0, 0.05);
   histww[21]   = new TH1F ("phoHovECorrWW", "phoHovECorrWW", 20, 0, 0.05);*/

   histwele[0]   = new TH1F("PFECALClusIsoCorrWele","PFECALClusIsoCorrWele",15,0,30);
   histwele[1]   = new TH1F("PFHCALClusIsoCorrWele","PFHCALClusIsoCorrWele",20,0,40);
   histwele[2]   = new TH1F("TrkIsoCorrWele","TrkIsoCorrWele",10,0,20);
   histwele[3]   = new TH1F("phoBDTWele","phoBDTWele",25,0,1);
   histwele[4]   = new TH1F("phoEtWele","phoEtWele",10,200,1000);
   histwele[5]   = new TH1F("phoR9Full5x5Wele", "phoR9Full5x5Wele",100, 0.1, 1.5 );
   histwele[6]   = new TH1F ("phoS4Full5x5Wele", "phoS4Full5x5Wele", 100, 0.4, 1.1);
   histwele[7]   = new TH1F ("phoEmaxOESCrFull5x5Wele", "phoEmaxOESCrFull5x5Wele", 150, 0, 1);
   histwele[8]   = new TH1F ("phoE2ndOESCrFull5x5Wele", "phoE2ndOESCrFull5x5Wele", 150, 0, 0.6);
   histwele[9]   = new TH1F ("phoE1x3OESCrFull5x5Wele", "phoE1x3OESCrFull5x5Wele", 150, 0.2, 1.1);
   histwele[10]   = new TH1F ("phoE2x5OESCrFull5x5Wele", "phoE2x5OESCrFull5x5Wele", 150, 0.6, 1.2);
   histwele[11]   = new TH1F ("phoE5x5OESCrFull5x5Wele", "phoE5x5OESCrFull5x5Wele", 150, 0.6, 1.2);
   histwele[12]   = new TH1F ("pho2x2OE3x3Full5x5Wele", "pho2x2OE3x3Full5x5Wele", 100, 0.5, 1);
   histwele[13]   = new TH1F ("phoSigmaIEtaIEtaWele", "phoSigmaIEtaIEtaWele", 100, 0, 0.04);
   histwele[14]   = new TH1F ("phoSigmaIEtaIPhiWele", "phoSigmaIEtaIPhiWele", 100, -0.001, 0.001);
   histwele[15]   = new TH1F ("phoSigmaIPhiIPhiWele", "phoSigmaIPhiIPhiWele", 100, 0, 0.07);
   histwele[16]   = new TH1F ("phoSieieOSipipFull5x5Wele", "phoSieieOSipipFull5x5Wele",150, 0, 1.5);
   histwele[17]   = new TH1F ("phoPhiWidthWele", "phoPhiWidthWele", 100, 0., 0.06);
   histwele[18]   = new TH1F ("phoEtaWidthWele", "phoEtaWidthWele", 100, 0, 0.02);
   histwele[19]   = new TH1F ("phoEtaWOPhiWFull5x5Wele", "phoEtaWOPhiWFull5x5Wele", 150, 0, 1.5);
   histwele[20]   = new TH1F ("phoHovEWele", "phoHovEWele", 20, 0, 0.05);
   histwele[21]   = new TH1F ("phoHovECorrWele", "phoHovECorrWele", 20, 0, 0.05);
   histwele[22]   = new TH1F ("metWele", "metWele", 10, 150, 1000);
   histwele[23]   = new TH1F("phoEtaWele","phoEtaWele",10,-1.4442,1.4442);
 
   histwz[0]   = new TH1F("PFECALClusIsoCorrWZ","PFECALClusIsoCorrWZ",15,0,30);
   histwz[1]   = new TH1F("PFHCALClusIsoCorrWZ","PFHCALClusIsoCorrWZ",20,0,40);
   histwz[2]   = new TH1F("TrkIsoCorrWZ","TrkIsoCorrWZ",10,0,20);
   histwz[3]   = new TH1F("phoBDTWZ","phoBDTWZ",25,0,1);
   histwz[4]   = new TH1F("phoEtWZ","phoEtWZ",10,200,1000);
   histwz[5]   = new TH1F("phoR9Full5x5WZ", "phoR9Full5x5WZ",100, 0.1, 1.5 );
   histwz[6]   = new TH1F ("phoS4Full5x5WZ", "phoS4Full5x5WZ", 100, 0.4, 1.1);
   histwz[7]   = new TH1F ("phoEmaxOESCrFull5x5WZ", "phoEmaxOESCrFull5x5WZ", 150, 0, 1);
   histwz[8]   = new TH1F ("phoE2ndOESCrFull5x5WZ", "phoE2ndOESCrFull5x5WZ", 150, 0, 0.6);
   histwz[9]   = new TH1F ("phoE1x3OESCrFull5x5WZ", "phoE1x3OESCrFull5x5WZ", 150, 0.2, 1.1);
   histwz[10]   = new TH1F ("phoE2x5OESCrFull5x5WZ", "phoE2x5OESCrFull5x5WZ", 150, 0.6, 1.2);
   histwz[11]   = new TH1F ("phoE5x5OESCrFull5x5WZ", "phoE5x5OESCrFull5x5WZ", 150, 0.6, 1.2);
   histwz[12]   = new TH1F ("pho2x2OE3x3Full5x5WZ", "pho2x2OE3x3Full5x5WZ", 100, 0.5, 1);
   histwz[13]   = new TH1F ("phoSigmaIEtaIEtaWZ", "phoSigmaIEtaIEtaWZ", 100, 0, 0.04);
   histwz[14]   = new TH1F ("phoSigmaIEtaIPhiWZ", "phoSigmaIEtaIPhiWZ", 100, -0.001, 0.001);
   histwz[15]   = new TH1F ("phoSigmaIPhiIPhiWZ", "phoSigmaIPhiIPhiWZ", 100, 0, 0.07);
   histwz[16]   = new TH1F ("phoSieieOSipipFull5x5WZ", "phoSieieOSipipFull5x5WZ",150, 0, 1.5);
   histwz[17]   = new TH1F ("phoPhiWidthWZ", "phoPhiWidthWZ", 100, 0., 0.06);
   histwz[18]   = new TH1F ("phoEtaWidthWZ", "phoEtaWidthWZ", 100, 0, 0.02);
   histwz[19]   = new TH1F ("phoEtaWOPhiWFull5x5WZ", "phoEtaWOPhiWFull5x5WZ", 150, 0, 1.5);
   histwz[20]   = new TH1F ("phoHovEWZ", "phoHovEWZ", 20, 0, 0.05);
   histwz[21]   = new TH1F ("phoHovECorrWZ", "phoHovECorrWZ", 20, 0, 0.05);
   histwz[22]   = new TH1F ("metWZ", "metWZ", 10, 150, 1000);
   histwz[23]   = new TH1F("phoEtaWZ","phoEtaWZ",10,-1.4442,1.4442);

   histzz[0]   = new TH1F("PFECALClusIsoCorrZZ","PFECALClusIsoCorrZZ",15,0,30);
   histzz[1]   = new TH1F("PFHCALClusIsoCorrZZ","PFHCALClusIsoCorrZZ",20,0,40);
   histzz[2]   = new TH1F("TrkIsoCorrZZ","TrkIsoCorrZZ",10,0,20);
   histzz[3]   = new TH1F("phoBDTZZ","phoBDTZZ",25,0,1);
   histzz[4]   = new TH1F("phoEtZZ","phoEtZZ",10,200,1000);
   histzz[5]   = new TH1F("phoR9Full5x5ZZ", "phoR9Full5x5ZZ",100, 0.1, 1.5 );
   histzz[6]   = new TH1F ("phoS4Full5x5ZZ", "phoS4Full5x5ZZ", 100, 0.4, 1.1);
   histzz[7]   = new TH1F ("phoEmaxOESCrFull5x5ZZ", "phoEmaxOESCrFull5x5ZZ", 150, 0, 1);
   histzz[8]   = new TH1F ("phoE2ndOESCrFull5x5ZZ", "phoE2ndOESCrFull5x5ZZ", 150, 0, 0.6);
   histzz[9]   = new TH1F ("phoE1x3OESCrFull5x5ZZ", "phoE1x3OESCrFull5x5ZZ", 150, 0.2, 1.1);
   histzz[10]   = new TH1F ("phoE2x5OESCrFull5x5ZZ", "phoE2x5OESCrFull5x5ZZ", 150, 0.6, 1.2);
   histzz[11]   = new TH1F ("phoE5x5OESCrFull5x5ZZ", "phoE5x5OESCrFull5x5ZZ", 150, 0.6, 1.2);
   histzz[12]   = new TH1F ("pho2x2OE3x3Full5x5ZZ", "pho2x2OE3x3Full5x5ZZ", 100, 0.5, 1);
   histzz[13]   = new TH1F ("phoSigmaIEtaIEtaZZ", "phoSigmaIEtaIEtaZZ", 100, 0, 0.04);
   histzz[14]   = new TH1F ("phoSigmaIEtaIPhiZZ", "phoSigmaIEtaIPhiZZ", 100, -0.001, 0.001);
   histzz[15]   = new TH1F ("phoSigmaIPhiIPhiZZ", "phoSigmaIPhiIPhiZZ", 100, 0, 0.07);
   histzz[16]   = new TH1F ("phoSieieOSipipFull5x5ZZ", "phoSieieOSipipFull5x5ZZ",150, 0, 1.5);
   histzz[17]   = new TH1F ("phoPhiWidthZZ", "phoPhiWidthZZ", 100, 0., 0.06);
   histzz[18]   = new TH1F ("phoEtaWidthZZ", "phoEtaWidthZZ", 100, 0, 0.02);
   histzz[19]   = new TH1F ("phoEtaWOPhiWFull5x5ZZ", "phoEtaWOPhiWFull5x5ZZ", 150, 0, 1.5);
   histzz[20]   = new TH1F ("phoHovEZZ", "phoHovEZZ", 20, 0, 0.05);
   histzz[21]   = new TH1F ("phoHovECorrZZ", "phoHovECorrZZ", 20, 0, 0.05);
   histzz[22]   = new TH1F ("metZZ", "metZZ", 10, 150, 1000);
   histzz[23]   = new TH1F("phoEtZZ","phoEtZZ",10,-1.4442,1.4442);

   Int_t data_nentries = (Int_t)data_tree->GetEntries();
   Int_t signalMC_nentries = (Int_t)signalMC_tree->GetEntries();
   /*Int_t gjets_nentries = (Int_t)gjets_tree->GetEntries();
   Int_t wjets_nentries = (Int_t)wjets_tree->GetEntries();
   Int_t wlnu_nentries = (Int_t)wlnu_tree->GetEntries();
   Int_t qcd_nentries = (Int_t)qcd_tree->GetEntries();
   Int_t ttgjets_nentries = (Int_t)ttgjets_tree->GetEntries();
   Int_t ww_nentries = (Int_t)ww_tree->GetEntries();*/
   Int_t wele_nentries = (Int_t)wele_tree->GetEntries();
   Int_t wz_nentries = (Int_t)wz_tree->GetEntries();
   Int_t zz_nentries = (Int_t)zz_tree->GetEntries();
   
   for (Int_t i=0; i<data_nentries; i++) {
      data_tree->GetEntry(i);
      histdata[0]->Fill(ph_pFECALClusIsoCorr, FRweight);
      histdata[1]->Fill(ph_pFHCALClusIsoCorr, FRweight);
      histdata[2]->Fill(ph_TrkIsoCorr, FRweight);
      histdata[3]->Fill(ph_BDTpred, FRweight);
      histdata[4]->Fill(ph_et, FRweight);
      histdata[5]->Fill(ph_R9Full5x5, FRweight);
      histdata[6]->Fill(ph_S4Full5x5, FRweight);
      histdata[7]->Fill(ph_EmaxOESCrFull5x5, FRweight);
      histdata[8]->Fill(ph_E2ndOESCrFull5x5, FRweight);
      histdata[9]->Fill(ph_E1x3OESCrFull5x5, FRweight);
      histdata[10]->Fill(ph_E2x5OESCrFull5x5, FRweight);
      histdata[11]->Fill(ph_E5x5OESCrFull5x5, FRweight);
      histdata[12]->Fill(ph_2x2OE3x3Full5x5, FRweight);
      histdata[13]->Fill(ph_SigmaIEtaIEta, FRweight);
      histdata[14]->Fill(ph_SigmaIEtaIPhi, FRweight);
      histdata[15]->Fill(ph_SigmaIPhiIPhi, FRweight);
      histdata[16]->Fill(ph_SieieOSipipFull5x5, FRweight);
      histdata[17]->Fill(ph_PhiWidth, FRweight);
      histdata[18]->Fill(ph_EtaWidth, FRweight);
      histdata[19]->Fill(ph_EtaWOPhiWFull5x5, FRweight);
      histdata[20]->Fill(ph_HoverE, FRweight);
      histdata[21]->Fill(ph_HoverECorr, FRweight);
      histdata[22]->Fill(met, FRweight);
      histdata[23]->Fill(ph_eta, FRweight);
   }

   for (Int_t i=0; i<signalMC_nentries; i++) {
      signalMC_tree->GetEntry(i);
      histMC[0]->Fill(ph_pFECALClusIsoCorrMC, puWeightMC);
      histMC[1]->Fill(ph_pFHCALClusIsoCorrMC, puWeightMC);
      histMC[2]->Fill(ph_TrkIsoCorrMC, puWeightMC);
      histMC[3]->Fill(ph_BDTpredMC, puWeightMC);
      histMC[4]->Fill(ph_etMC, puWeightMC);
      histMC[5]->Fill(ph_R9Full5x5MC, puWeightMC);
      histMC[6]->Fill(ph_S4Full5x5MC, puWeightMC);
      histMC[7]->Fill(ph_EmaxOESCrFull5x5MC, puWeightMC);
      histMC[8]->Fill(ph_E2ndOESCrFull5x5MC, puWeightMC); 
      histMC[9]->Fill(ph_E1x3OESCrFull5x5MC, puWeightMC);
      histMC[10]->Fill(ph_E2x5OESCrFull5x5MC, puWeightMC);
      histMC[11]->Fill(ph_E5x5OESCrFull5x5MC, puWeightMC);
      histMC[12]->Fill(ph_2x2OE3x3Full5x5MC, puWeightMC);
      histMC[13]->Fill(ph_SigmaIEtaIEtaMC, puWeightMC);
      histMC[14]->Fill(ph_SigmaIEtaIPhiMC, puWeightMC);
      histMC[15]->Fill(ph_SigmaIPhiIPhiMC, puWeightMC);
      histMC[16]->Fill(ph_SieieOSipipFull5x5MC, puWeightMC);
      histMC[17]->Fill(ph_PhiWidthMC, puWeightMC);
      histMC[18]->Fill(ph_EtaWidthMC, puWeightMC);
      histMC[19]->Fill(ph_EtaWOPhiWFull5x5MC, puWeightMC);
      histMC[20]->Fill(ph_HoverEMC, puWeightMC);
      histMC[21]->Fill(ph_HoverECorrMC, puWeightMC);
      histMC[22]->Fill(metMC, puWeightMC*FRweightMC);
      histMC[23]->Fill(ph_etaMC, puWeightMC*FRweightMC);
   }
   
   /*for (Int_t i=0; i<gjets_nentries; i++) {
      gjets_tree->GetEntry(i);
      histgjets[0]->Fill(ph_pFECALClusIsoCorrGj, puWeightGj);
      histgjets[1]->Fill(ph_pFHCALClusIsoCorrGj, puWeightGj);
      histgjets[2]->Fill(ph_TrkIsoCorrGj, puWeightGj);
      histgjets[3]->Fill(ph_BDTpredGj, puWeightGj);
      histgjets[4]->Fill(ph_etGj, puWeightGj);
      histgjets[5]->Fill(ph_R9Full5x5Gj, puWeightGj);
      histgjets[6]->Fill(ph_S4Full5x5Gj, puWeightGj);
      histgjets[7]->Fill(ph_EmaxOESCrFull5x5Gj, puWeightGj);
      histgjets[8]->Fill(ph_E2ndOESCrFull5x5Gj, puWeightGj);
      histgjets[9]->Fill(ph_E1x3OESCrFull5x5Gj, puWeightGj);
      histgjets[10]->Fill(ph_E2x5OESCrFull5x5Gj, puWeightGj);
      histgjets[11]->Fill(ph_E5x5OESCrFull5x5Gj, puWeightGj);
      histgjets[12]->Fill(ph_2x2OE3x3Full5x5Gj, puWeightGj);
      histgjets[13]->Fill(ph_SigmaIEtaIEtaGj, puWeightGj);
      histgjets[14]->Fill(ph_SigmaIEtaIPhiGj, puWeightGj);
      histgjets[15]->Fill(ph_SigmaIPhiIPhiGj, puWeightGj);
      histgjets[16]->Fill(ph_SieieOSipipFull5x5Gj, puWeightGj);
      histgjets[17]->Fill(ph_PhiWidthGj, puWeightGj);
      histgjets[18]->Fill(ph_EtaWidthGj, puWeightGj);
      histgjets[19]->Fill(ph_EtaWOPhiWFull5x5Gj, puWeightGj);
      histgjets[20]->Fill(ph_HoverEGj, puWeightGj);
      histgjets[21]->Fill(ph_HoverECorrGj, puWeightGj);
   }

   for (Int_t i=0; i<wjets_nentries; i++) {
      wjets_tree->GetEntry(i);
      histwjets[0]->Fill(ph_pFECALClusIsoCorrWj, puWeightWj);
      histwjets[1]->Fill(ph_pFHCALClusIsoCorrWj, puWeightWj);
      histwjets[2]->Fill(ph_TrkIsoCorrWj, puWeightWj);
      histwjets[3]->Fill(ph_BDTpredWj, puWeightWj);
      histwjets[4]->Fill(ph_etWj, puWeightWj);
      histwjets[5]->Fill(ph_R9Full5x5Wj, puWeightWj);
      histwjets[6]->Fill(ph_S4Full5x5Wj, puWeightWj);
      histwjets[7]->Fill(ph_EmaxOESCrFull5x5Wj, puWeightWj);
      histwjets[8]->Fill(ph_E2ndOESCrFull5x5Wj, puWeightWj);
      histwjets[9]->Fill(ph_E1x3OESCrFull5x5Wj, puWeightWj);
      histwjets[10]->Fill(ph_E2x5OESCrFull5x5Wj, puWeightWj);
      histwjets[11]->Fill(ph_E5x5OESCrFull5x5Wj, puWeightWj);
      histwjets[12]->Fill(ph_2x2OE3x3Full5x5Wj, puWeightWj);
      histwjets[13]->Fill(ph_SigmaIEtaIEtaWj, puWeightWj);
      histwjets[14]->Fill(ph_SigmaIEtaIPhiWj, puWeightWj);
      histwjets[15]->Fill(ph_SigmaIPhiIPhiWj, puWeightWj);
      histwjets[16]->Fill(ph_SieieOSipipFull5x5Wj, puWeightWj);
      histwjets[17]->Fill(ph_PhiWidthWj, puWeightWj);
      histwjets[18]->Fill(ph_EtaWidthWj, puWeightWj);
      histwjets[19]->Fill(ph_EtaWOPhiWFull5x5Wj, puWeightWj);
      histwjets[20]->Fill(ph_HoverEWj, puWeightWj);
      histwjets[21]->Fill(ph_HoverECorrWj, puWeightWj);
   }

   for (Int_t i=0; i<wlnu_nentries; i++) {
      wlnu_tree->GetEntry(i);
      histwlnu[0]->Fill(ph_pFECALClusIsoCorrWlnu, puWeightWlnu);
      histwlnu[1]->Fill(ph_pFHCALClusIsoCorrWlnu, puWeightWlnu);
      histwlnu[2]->Fill(ph_TrkIsoCorrWlnu, puWeightWlnu);
      histwlnu[3]->Fill(ph_BDTpredWlnu, puWeightWlnu);
      histwlnu[4]->Fill(ph_etWlnu, puWeightWlnu);
      histwlnu[5]->Fill(ph_R9Full5x5Wlnu, puWeightWlnu);
      histwlnu[6]->Fill(ph_S4Full5x5Wlnu, puWeightWlnu);
      histwlnu[7]->Fill(ph_EmaxOESCrFull5x5Wlnu, puWeightWlnu);
      histwlnu[8]->Fill(ph_E2ndOESCrFull5x5Wlnu, puWeightWlnu);
      histwlnu[9]->Fill(ph_E1x3OESCrFull5x5Wlnu, puWeightWlnu);
      histwlnu[10]->Fill(ph_E2x5OESCrFull5x5Wlnu, puWeightWlnu);
      histwlnu[11]->Fill(ph_E5x5OESCrFull5x5Wlnu, puWeightWlnu);
      histwlnu[12]->Fill(ph_2x2OE3x3Full5x5Wlnu, puWeightWlnu);
      histwlnu[13]->Fill(ph_SigmaIEtaIEtaWlnu, puWeightWlnu);
      histwlnu[14]->Fill(ph_SigmaIEtaIPhiWlnu, puWeightWlnu);
      histwlnu[15]->Fill(ph_SigmaIPhiIPhiWlnu, puWeightWlnu);
      histwlnu[16]->Fill(ph_SieieOSipipFull5x5Wlnu, puWeightWlnu);
      histwlnu[17]->Fill(ph_PhiWidthWlnu, puWeightWlnu);
      histwlnu[18]->Fill(ph_EtaWidthWlnu, puWeightWlnu);
      histwlnu[19]->Fill(ph_EtaWOPhiWFull5x5Wlnu, puWeightWlnu);
      histwlnu[20]->Fill(ph_HoverEWlnu, puWeightWlnu);
      histwlnu[21]->Fill(ph_HoverECorrWlnu, puWeightWlnu);
   }

   for (Int_t i=0; i<qcd_nentries; i++) {
      qcd_tree->GetEntry(i);
      histqcd[0]->Fill(ph_pFECALClusIsoCorrQcd, puWeightQcd);
      histqcd[1]->Fill(ph_pFHCALClusIsoCorrQcd, puWeightQcd);
      histqcd[2]->Fill(ph_TrkIsoCorrQcd, puWeightQcd);
      histqcd[3]->Fill(ph_BDTpredQcd, puWeightQcd);
      histqcd[4]->Fill(ph_etQcd, puWeightQcd);
      histqcd[5]->Fill(ph_R9Full5x5Qcd, puWeightQcd);
      histqcd[6]->Fill(ph_S4Full5x5Qcd, puWeightQcd);
      histqcd[7]->Fill(ph_EmaxOESCrFull5x5Qcd, puWeightQcd);
      histqcd[8]->Fill(ph_E2ndOESCrFull5x5Qcd, puWeightQcd);
      histqcd[9]->Fill(ph_E1x3OESCrFull5x5Qcd, puWeightQcd);
      histqcd[10]->Fill(ph_E2x5OESCrFull5x5Qcd, puWeightQcd);
      histqcd[11]->Fill(ph_E5x5OESCrFull5x5Qcd, puWeightQcd);
      histqcd[12]->Fill(ph_2x2OE3x3Full5x5Qcd, puWeightQcd);
      histqcd[13]->Fill(ph_SigmaIEtaIEtaQcd, puWeightQcd);
      histqcd[14]->Fill(ph_SigmaIEtaIPhiQcd, puWeightQcd);
      histqcd[15]->Fill(ph_SigmaIPhiIPhiQcd, puWeightQcd);
      histqcd[16]->Fill(ph_SieieOSipipFull5x5Qcd, puWeightQcd);
      histqcd[17]->Fill(ph_PhiWidthQcd, puWeightQcd);
      histqcd[18]->Fill(ph_EtaWidthQcd, puWeightQcd);
      histqcd[19]->Fill(ph_EtaWOPhiWFull5x5Qcd, puWeightQcd);
      histqcd[20]->Fill(ph_HoverEQcd, puWeightQcd);
      histqcd[21]->Fill(ph_HoverECorrQcd, puWeightQcd);
   }

   for (Int_t i=0; i<ttgjets_nentries; i++) {
      ttgjets_tree->GetEntry(i);
      histttgjets[0]->Fill(ph_pFECALClusIsoCorrTTGj, puWeightTTGj);
      histttgjets[1]->Fill(ph_pFHCALClusIsoCorrTTGj, puWeightTTGj);
      histttgjets[2]->Fill(ph_TrkIsoCorrTTGj, puWeightTTGj);
      histttgjets[3]->Fill(ph_BDTpredTTGj, puWeightTTGj);
      histttgjets[4]->Fill(ph_etTTGj, puWeightTTGj);
      histttgjets[5]->Fill(ph_R9Full5x5TTGj, puWeightTTGj);
      histttgjets[6]->Fill(ph_S4Full5x5TTGj, puWeightTTGj);
      histttgjets[7]->Fill(ph_EmaxOESCrFull5x5TTGj, puWeightTTGj);
      histttgjets[8]->Fill(ph_E2ndOESCrFull5x5TTGj, puWeightTTGj);
      histttgjets[9]->Fill(ph_E1x3OESCrFull5x5TTGj, puWeightTTGj);
      histttgjets[10]->Fill(ph_E2x5OESCrFull5x5TTGj, puWeightTTGj);
      histttgjets[11]->Fill(ph_E5x5OESCrFull5x5TTGj, puWeightTTGj);
      histttgjets[12]->Fill(ph_2x2OE3x3Full5x5TTGj, puWeightTTGj);
      histttgjets[13]->Fill(ph_SigmaIEtaIEtaTTGj, puWeightTTGj);
      histttgjets[14]->Fill(ph_SigmaIEtaIPhiTTGj, puWeightTTGj);
      histttgjets[15]->Fill(ph_SigmaIPhiIPhiTTGj, puWeightTTGj);
      histttgjets[16]->Fill(ph_SieieOSipipFull5x5TTGj, puWeightTTGj);
      histttgjets[17]->Fill(ph_PhiWidthTTGj, puWeightTTGj);
      histttgjets[18]->Fill(ph_EtaWidthTTGj, puWeightTTGj);
      histttgjets[19]->Fill(ph_EtaWOPhiWFull5x5TTGj, puWeightTTGj);
      histttgjets[20]->Fill(ph_HoverETTGj, puWeightTTGj);
      histttgjets[21]->Fill(ph_HoverECorrTTGj, puWeightTTGj);
   }

   for (Int_t i=0; i<ww_nentries; i++) {
      ww_tree->GetEntry(i);
      histww[0]->Fill(ph_pFECALClusIsoCorrWW, puWeightWW);
      histww[1]->Fill(ph_pFHCALClusIsoCorrWW, puWeightWW);
      histww[2]->Fill(ph_TrkIsoCorrWW, puWeightWW);
      histww[3]->Fill(ph_BDTpredWW, puWeightWW);
      histww[4]->Fill(ph_etWW, puWeightWW);
      histww[5]->Fill(ph_R9Full5x5WW, puWeightWW);
      histww[6]->Fill(ph_S4Full5x5WW, puWeightWW);
      histww[7]->Fill(ph_EmaxOESCrFull5x5WW, puWeightWW);
      histww[8]->Fill(ph_E2ndOESCrFull5x5WW, puWeightWW);
      histww[9]->Fill(ph_E1x3OESCrFull5x5WW, puWeightWW);
      histww[10]->Fill(ph_E2x5OESCrFull5x5WW, puWeightWW);
      histww[11]->Fill(ph_E5x5OESCrFull5x5WW, puWeightWW);
      histww[12]->Fill(ph_2x2OE3x3Full5x5WW, puWeightWW);
      histww[13]->Fill(ph_SigmaIEtaIEtaWW, puWeightWW);
      histww[14]->Fill(ph_SigmaIEtaIPhiWW, puWeightWW);
      histww[15]->Fill(ph_SigmaIPhiIPhiWW, puWeightWW);
      histww[16]->Fill(ph_SieieOSipipFull5x5WW, puWeightWW);
      histww[17]->Fill(ph_PhiWidthWW, puWeightWW);
      histww[18]->Fill(ph_EtaWidthWW, puWeightWW);
      histww[19]->Fill(ph_EtaWOPhiWFull5x5WW, puWeightWW);
      histww[20]->Fill(ph_HoverEWW, puWeightWW);
      histww[21]->Fill(ph_HoverECorrWW, puWeightWW);
   }*/

   for (Int_t i=0; i<wele_nentries; i++) {
      wele_tree->GetEntry(i);
      histwele[0]->Fill(ph_pFECALClusIsoCorrWele, puWeightWele);
      histwele[1]->Fill(ph_pFHCALClusIsoCorrWele, puWeightWele);
      histwele[2]->Fill(ph_TrkIsoCorrWele, puWeightWele);
      histwele[3]->Fill(ph_BDTpredWele, puWeightWele);
      histwele[4]->Fill(ph_etWele, puWeightWele*FRweightWele);
      histwele[5]->Fill(ph_R9Full5x5Wele, puWeightWele);
      histwele[6]->Fill(ph_S4Full5x5Wele, puWeightWele);
      histwele[7]->Fill(ph_EmaxOESCrFull5x5Wele, puWeightWele);
      histwele[8]->Fill(ph_E2ndOESCrFull5x5Wele, puWeightWele);
      histwele[9]->Fill(ph_E1x3OESCrFull5x5Wele, puWeightWele);
      histwele[10]->Fill(ph_E2x5OESCrFull5x5Wele, puWeightWele);
      histwele[11]->Fill(ph_E5x5OESCrFull5x5Wele, puWeightWele);
      histwele[12]->Fill(ph_2x2OE3x3Full5x5Wele, puWeightWele);
      histwele[13]->Fill(ph_SigmaIEtaIEtaWele, puWeightWele);
      histwele[14]->Fill(ph_SigmaIEtaIPhiWele, puWeightWele);
      histwele[15]->Fill(ph_SigmaIPhiIPhiWele, puWeightWele);
      histwele[16]->Fill(ph_SieieOSipipFull5x5Wele, puWeightWele);
      histwele[17]->Fill(ph_PhiWidthWele, puWeightWele);
      histwele[18]->Fill(ph_EtaWidthWele, puWeightWele);
      histwele[19]->Fill(ph_EtaWOPhiWFull5x5Wele, puWeightWele);
      histwele[20]->Fill(ph_HoverEWele, puWeightWele);
      histwele[21]->Fill(ph_HoverECorrWele, puWeightWele);
      histwele[22]->Fill(metWele, puWeightWele*FRweightWele);
      histwele[23]->Fill(ph_etaWele, puWeightWele*FRweightWele);
   }

   for (Int_t i=0; i<wz_nentries; i++) {
      wz_tree->GetEntry(i);
      histwz[0]->Fill(ph_pFECALClusIsoCorrWZ, puWeightWZ);
      histwz[1]->Fill(ph_pFHCALClusIsoCorrWZ, puWeightWZ);
      histwz[2]->Fill(ph_TrkIsoCorrWZ, puWeightWZ);
      histwz[3]->Fill(ph_BDTpredWZ, puWeightWZ);
      histwz[4]->Fill(ph_etWZ, puWeightWZ*FRweightWZ);
      histwz[5]->Fill(ph_R9Full5x5WZ, puWeightWZ);
      histwz[6]->Fill(ph_S4Full5x5WZ, puWeightWZ);
      histwz[7]->Fill(ph_EmaxOESCrFull5x5WZ, puWeightWZ);
      histwz[8]->Fill(ph_E2ndOESCrFull5x5WZ, puWeightWZ);
      histwz[9]->Fill(ph_E1x3OESCrFull5x5WZ, puWeightWZ);
      histwz[10]->Fill(ph_E2x5OESCrFull5x5WZ, puWeightWZ);
      histwz[11]->Fill(ph_E5x5OESCrFull5x5WZ, puWeightWZ);
      histwz[12]->Fill(ph_2x2OE3x3Full5x5WZ, puWeightWZ);
      histwz[13]->Fill(ph_SigmaIEtaIEtaWZ, puWeightWZ);
      histwz[14]->Fill(ph_SigmaIEtaIPhiWZ, puWeightWZ);
      histwz[15]->Fill(ph_SigmaIPhiIPhiWZ, puWeightWZ);
      histwz[16]->Fill(ph_SieieOSipipFull5x5WZ, puWeightWZ);
      histwz[17]->Fill(ph_PhiWidthWZ, puWeightWZ);
      histwz[18]->Fill(ph_EtaWidthWZ, puWeightWZ);
      histwz[19]->Fill(ph_EtaWOPhiWFull5x5WZ, puWeightWZ);
      histwz[20]->Fill(ph_HoverEWZ, puWeightWZ);
      histwz[21]->Fill(ph_HoverECorrWZ, puWeightWZ);
      histwz[22]->Fill(metWZ, puWeightWZ*FRweightWZ);
      histwz[23]->Fill(ph_etaWZ, puWeightWZ*FRweightWZ);
   }

   for (Int_t i=0; i<zz_nentries; i++) {
      zz_tree->GetEntry(i);
      histzz[0]->Fill(ph_pFECALClusIsoCorrZZ, puWeightZZ);
      histzz[1]->Fill(ph_pFHCALClusIsoCorrZZ, puWeightZZ);
      histzz[2]->Fill(ph_TrkIsoCorrZZ, puWeightZZ);
      histzz[3]->Fill(ph_BDTpredZZ, puWeightZZ);
      histzz[4]->Fill(ph_etZZ, puWeightZZ*FRweightZZ);
      histzz[5]->Fill(ph_R9Full5x5ZZ, puWeightZZ);
      histzz[6]->Fill(ph_S4Full5x5ZZ, puWeightZZ);
      histzz[7]->Fill(ph_EmaxOESCrFull5x5ZZ, puWeightZZ);
      histzz[8]->Fill(ph_E2ndOESCrFull5x5ZZ, puWeightZZ);
      histzz[9]->Fill(ph_E1x3OESCrFull5x5ZZ, puWeightZZ);
      histzz[10]->Fill(ph_E2x5OESCrFull5x5ZZ, puWeightZZ);
      histzz[11]->Fill(ph_E5x5OESCrFull5x5ZZ, puWeightZZ);
      histzz[12]->Fill(ph_2x2OE3x3Full5x5ZZ, puWeightZZ);
      histzz[13]->Fill(ph_SigmaIEtaIEtaZZ, puWeightZZ);
      histzz[14]->Fill(ph_SigmaIEtaIPhiZZ, puWeightZZ);
      histzz[15]->Fill(ph_SigmaIPhiIPhiZZ, puWeightZZ);
      histzz[16]->Fill(ph_SieieOSipipFull5x5ZZ, puWeightZZ);
      histzz[17]->Fill(ph_PhiWidthZZ, puWeightZZ);
      histzz[18]->Fill(ph_EtaWidthZZ, puWeightZZ);
      histzz[19]->Fill(ph_EtaWOPhiWFull5x5ZZ, puWeightZZ);
      histzz[20]->Fill(ph_HoverEZZ, puWeightZZ);
      histzz[21]->Fill(ph_HoverECorrZZ, puWeightZZ);
      histzz[22]->Fill(metZZ, puWeightZZ*FRweightZZ);
      histzz[23]->Fill(ph_etaZZ, puWeightZZ*FRweightZZ);
   }

   TString title[nhist] = {"ph_pFECALClusIsoCorr", "ph_pFHCALClusIsoCorr", "ph_TrkIsoCorr", "ph_BDTpred", "ph_et", "ph_R9Full5x5", "ph_S4Full5x5", "ph_EmaxOESCrFull5x5", "ph_E2ndOESCrFull5x5", "ph_E1x3OESCrFull5x5", "ph_E2x5OESCrFull5x5", "ph_E5x5OESCrFull5x5", "ph_2x2OE3x3Full5x5", "ph_sieie", "ph_sieip", "ph_sipip", "ph_sieieOsipip", "ph_sp", "ph_se", "ph_seOsp", "ph_hoe", "ph_hoeCorr", "pfMET", "ph_eta"};
   TString title_axis[nhist] = {"ECAL Cluster IsoCorr", "HCAL Cluster IsoCorr", "Track IsoCorr", "Photon BDT score", "p_{T}^{#gamma}", "R9Full5x5", "S4Full5x5", "E_{max}/E_{SC}^{raw}", "E_{2}/E_{SC}^{raw}", "E_{1x3}/E_{SC}^{raw}", "E_{2x5}/E_{SC}^{raw}", "E_{5x5}/E_{SC}^{raw}", "E_{2x2}/E_{3x3}", "#sigma_{i#etai#eta}", "#sigma_{i#etai#phi}", "#sigma_{i#phii#phi}", "#sigma_{i#etai#eta}/#sigma_{i#phii#phi}", "#sigma_{#phi}", "#sigma_{#eta}", "#sigma_{#eta}/#sigma_{#phi}", "(H/E)_{#gamma}", "(H/E)_{#gamma} Corr", "MET", "#eta"};
 
   for (int i=0; i<nhist; i++) {
      //Double_t scale = histdata[i]->Integral()/histMC[i]->Integral();
      //histMC[i]->Sumw2(1);
      //histMC[i]->Scale(scale);
   }

   TH1F* h_ratio[nhist];
   for (int i=0; i<nhist; i++) {
      h_ratio[i] = (TH1F*)histdata[i]->Clone("h_ratio");
      h_ratio[i]->Sumw2(1);
      //h_ratio[i]->Divide(histwele[i]);
      h_ratio[i]->Divide(histMC[i]);
      
      TCanvas *c1 = new TCanvas("c1", "c1");
      c1->cd();

      TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
      pad1->SetBottomMargin(0);
      pad1->Draw();
      pad1->cd(); 

      histMC[i]->SetFillColor(kOrange);
      histMC[i]->SetLineColor(kOrange);

      /*histMC[i]->SetFillColor(kOrange-4);
      histMC[i]->SetLineColor(kOrange-4);

      histzz[i]->SetFillColor(kSpring+4);
      histzz[i]->SetLineColor(kSpring+4);

      histwz[i]->SetFillColor(kSpring+2);
      histwz[i]->SetLineColor(kSpring+2);

      histww[i]->SetFillColor(kTeal-6);
      histww[i]->SetLineColor(kTeal-6);

      histwlnu[i]->SetFillColor(kCyan-2);
      histwlnu[i]->SetLineColor(kCyan-2);

      histwjets[i]->SetFillColor(kAzure-5); //kBlue-7
      histwjets[i]->SetLineColor(kAzure-5); //kBlue-7

      histttgjets[i]->SetFillColor(kBlue-9);
      histttgjets[i]->SetLineColor(kBlue-9);

      histqcd[i]->SetFillColor(kPink+2);
      histqcd[i]->SetLineColor(kPink+2);

      histgjets[i]->SetFillColor(kMagenta-8);
      histgjets[i]->SetLineColor(kMagenta-8);

      histwz[i]->SetFillColor(kSpring+2);
      histwz[i]->SetLineColor(kSpring+2);

      histzz[i]->SetFillColor(kSpring+4);
      histzz[i]->SetLineColor(kSpring+4);

      histwele[i]->SetFillColor(kCyan-2);
      histwele[i]->SetLineColor(kCyan-2);*/

      histdata[i]->SetMarkerStyle(20);

      histMC[i]->GetYaxis()->SetTitle("Events");
      histMC[i]->GetYaxis()->SetTitleSize(0.05);
      histMC[i]->GetYaxis()->SetTitleOffset(0.95);
      histMC[i]->GetYaxis()->SetLabelSize(0.045);

      TLegend* leg = new TLegend();
      leg = new TLegend(0.68, 0.5, 0.9, 0.891);
      leg->SetBorderSize(0);
      leg->SetEntrySeparation(0.01);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      leg->AddEntry(histMC[i], "Total MC", "f");
      leg->AddEntry(histdata[i], "Data", "ep");

      c1->Update();

      histMC[i]->SetTitle("");
      histdata[i]->SetTitle("");

      histMC[i]->Draw("hist");
      histdata[i]->Draw("ep, SAME");
 
      leg->Draw("SAME");

      //pad1->SetLogy();

      CMS_lumi(pad1, 17, 0);
      pad1->Update();
      c1->Update();
      c1->cd();

      TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
      pad2->SetTopMargin(0.01);
      pad2->SetBottomMargin(0.2);
      pad2->Draw();
      pad2->cd();

      h_ratio[i]->SetTitle("");
      h_ratio[i]->GetXaxis()->SetTitle(title_axis[i]);
      h_ratio[i]->GetXaxis()->SetTitleFont(42);
      h_ratio[i]->GetXaxis()->SetTitleSize(0.11);
      h_ratio[i]->GetXaxis()->SetTitleOffset(0.8);
      h_ratio[i]->GetXaxis()->SetLabelFont(42);
      h_ratio[i]->GetXaxis()->SetLabelSize(0.1);
      h_ratio[i]->GetYaxis()->SetTitle("Data/MC");
      h_ratio[i]->GetYaxis()->SetTitleSize(0.11);
      h_ratio[i]->GetYaxis()->SetTitleOffset(0.43);
      h_ratio[i]->GetYaxis()->SetLabelSize(0.1);
      h_ratio[i]->GetYaxis()->SetLabelOffset(0.01);
      h_ratio[i]->GetYaxis()->SetNdivisions(505);
      h_ratio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
      h_ratio[i]->SetMarkerStyle(20);

      TLine* line = new TLine(h_ratio[i]->GetXaxis()->GetXmin(), 1.0, h_ratio[i]->GetXaxis()->GetXmax(), 1.0); 
      line->SetLineColor(kRed);
      line->SetLineWidth(1);
      h_ratio[i]->Draw();
      line->Draw("SAME");
      c1->Update();
      
      c1->SaveAs(title[i]+".png");

      TCanvas *c2 = new TCanvas("c2", "c2");
      c2->cd();

      TPad* pad11 = new TPad("pad11", "pad11", 0.0, 0.3, 1.0, 1.0);
      pad11->SetBottomMargin(0);
      pad11->Draw();
      pad11->cd();

      hs[i] = new THStack(title[i]+"_stack", "");

      histwele[i]->SetFillColor(kCyan-2);      
      histwele[i]->SetLineColor(kCyan-2);

      histzz[i]->SetFillColor(kSpring+4);
      histzz[i]->SetLineColor(kSpring+4);

      histwz[i]->SetFillColor(kSpring+2);
      histwz[i]->SetLineColor(kSpring+2);
      /*histww[i]->SetFillColor(kTeal-6);
      histwlnu[i]->SetFillColor(kCyan-2);
      histwjets[i]->SetFillColor(kAzure-5); //kBlue-7
      histttgjets[i]->SetLineColor(kBlue-9);
      histqcd[i]->SetFillColor(kPink+2);
      histgjets[i]->SetFillColor(kMagenta-8);*/

      //hs[i]->Add(histMC[i]);
      //hs[i]->Add(histgjets[i]);
      //hs[i]->Add(histqcd[i]);
      //hs[i]->Add(histttgjets[i]);
      //hs[i]->Add(histwjets[i]);
      //hs[i]->Add(histwlnu[i]);
      //hs[i]->Add(histww[i]);
      hs[i]->Add(histzz[i]);
      hs[i]->Add(histwz[i]);
      hs[i]->Add(histwele[i]);
      hs[i]->Draw("hist");
      histdata[i]->Draw("ep, SAME");

      hs[i]->GetYaxis()->SetTitle("Events");

      TLegend* leg_stack = new TLegend();
      leg_stack = new TLegend(0.68, 0.5, 0.9, 0.891);
      leg_stack->SetBorderSize(0);
      leg_stack->SetEntrySeparation(0.01);
      leg_stack->SetFillColor(0);
      leg_stack->SetFillStyle(0);

      leg_stack->AddEntry(histwele[i], "WToENu", "f");
      leg_stack->AddEntry(histwz[i], "WZ", "f");
      leg_stack->AddEntry(histzz[i], "ZZ", "f");
      leg_stack->AddEntry(histdata[i], "Data", "ep");
      leg_stack->Draw("SAME");
      
      //pad11->SetLogy();

      CMS_lumi(pad11, 17, 0);
      pad11->Update();
      c2->Update();
      c2->cd();

      TPad* pad12 = new TPad("pad12", "pad12", 0.0, 0.0, 1.0, 0.3);
      pad12->SetTopMargin(0.01);
      pad12->SetBottomMargin(0.2);
      pad12->Draw();
      pad12->cd();
     
      h_ratio[i]->Draw();
      line->Draw("SAME");
      c2->Update(); 

      c2->SaveAs(title[i]+"stack.png");
   }
}

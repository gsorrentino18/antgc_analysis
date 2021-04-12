#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <THStack.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMarker.h>

using namespace std;

void Roc(){

   TCanvas * c1 = new TCanvas("c1","c1");
   c1->cd(); 
   gPad->SetGrid();

   const Int_t npoints = 6;

   Double_t bkg_rej[npoints] = {1-0.58555, 1-0.61450918, 1-0.64228817, 1-0.69325984, 1-0.73619064, 1-0.78463767};
   Double_t sig_eff[npoints+1] = {0.9381, 0.9583, 0.9719, 0.9913, 0.9963, 0.9984, 1};

   //Double_t bkg_rej[npoints] = {1-0.6976, 1-0.7303, 1-0.7580, 1-0.8011, 1-0.8322, 1-0.8644};
   //Double_t sig_eff[npoints+1] = {0.9381, 0.9583, 0.9719, 0.9913, 0.9963, 0.9984, 1};

   /*Double_t bkg_rej_err[npoints] = {1-0.6973, 1-0.7301, 1-0.7577, 1-0.8011, 1-0.8322, 1-0.8643};
   Double_t sig_eff_err[npoints+1] = {0.9377, 0.9581, 0.9718, 0.9912, 0.9963, 0.9984, 1};
   TGraphErrors* roc = new TGraph(npoints, sig_eff, bkg_rej, sig_eff_err, bkg_rej_err);

   roc->SetTitle("");
   roc->GetXaxis()->SetTitle("Signal Efficiency (%)");
   roc->GetYaxis()->SetTitle("Background Rejection (%)");
   roc->SetMarkerStyle(20);
   roc->Draw("AP");

   Int_t colors[npoints] = {kBlack, kRed, kMagenta, kGreen, kBlue, kOrange};
   TMarker *m[npoints];
   for (Int_t i=0; i<npoints; i++) {
      m[i] = new TMarker(sig_eff[i], bkg_rej[i], 20);
      //m->SetMarkerSize(2);
      m[i]->SetMarkerColor(colors[i]);
      m[i]->Draw();
   }

   TLegend* leg = new TLegend(0.75, 0.7, 0.9, 0.891);
   leg->SetBorderSize(0);
   leg->SetEntrySeparation(0.01);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(m[0],"VVVLoose","p");
   leg->AddEntry(m[1],"VVLoose","p");
   leg->AddEntry(m[2],"VLoose","p");
   leg->AddEntry(m[3],"Loose","p");
   leg->AddEntry(m[4],"Medium","p");
   leg->AddEntry(m[5],"Tight","p");
   leg->Draw("SAME");
   c1->SaveAs("roc.png");*/

   Double_t punzi[npoints+1] = {0};
   Double_t a = 2.;
   Double_t b = 5.;
   Double_t bkg[npoints+1] = {147.975, 155.292, 162.312, 175.193, 186.042, 198.285, 252.709};
   Double_t sig[npoints+1] = {849.179, 867.466, 879.755, 897.301, 901.895, 903.766, 905.211};
   //Double_t bkg[npoints+1] = {447.463, 468.405, 486.168, 513.82, 533.808, 554.446, 641.417};
   //Double_t sig[npoints+1] = {849.179, 867.466, 879.755, 897.301, 901.895, 903.766, 905.211};
   Double_t sigma_bkg[npoints+1] = {0};
   Double_t sigma_sig[npoints+1] = {0};
   Double_t sigma_eff[npoints+1] = {0};
   Double_t sigma_punzi[npoints+1] = {0};
   Double_t sigma_total[npoints+1] = {0};

   Double_t sigma_rej[npoints+1] = {0}; //for TGraphErrors

   for (Int_t j=0; j<npoints+1; j++) {
      sigma_bkg[j] = sqrt(bkg[j]);
      sigma_sig[j] = sqrt(sig[j]);

      punzi[j] = pow(b,2)/2. + a*sqrt(bkg[j]) + (b/2.)*sqrt(b*b +4*a*sqrt(bkg[j]) + 4*bkg[j]);
      //punzi[j] = a*a/8. +9*b*b/13 + a*sqrt(bkg[j]) + (b/2.)*sqrt(b*b +4*a*sqrt(bkg[j]) + 4*bkg[j]);

      //naive efficiency error 
      //sigma_eff[j] = sig_eff[j]*sqrt( pow(sigma_sig[j]/sig[j], 2) + pow(sigma_sig[npoints]/sig[npoints], 2) );
      //sigma_rej[j] = (1-bkg_rej[j])*sqrt( pow(sigma_bkg[j]/bkg[j], 2) + pow(sigma_bkg[npoints]/bkg[npoints], 2) );

      //poissonian efficiency error
      //sigma_eff[j] = sqrt( sig[j]/pow(sig[npoints],2) );
      //sigma_rej[j] = sqrt( bkg[j]/pow(bkg[npoints],2) );;

      //correct efficiency treatment
      sigma_eff[j] = sqrt( ((sig[j]+1)*(sig[j]+2))/((sig[npoints]+2)*(sig[npoints]+3)) - pow( (sig[j]+1)/(sig[npoints]+2), 2) );
      sigma_rej[j] = sqrt( ((bkg[j]+1)*(bkg[j]+2))/((bkg[npoints]+2)*(bkg[npoints]+3)) - pow( (bkg[j]+1)/(bkg[npoints]+2), 2) );

      Double_t p1 = ((2*a*b/sqrt(bkg[j])) + 4*b)/( 4.*sqrt(4*a*sqrt(bkg[j]) + b*b + 4*bkg[j]) );
      Double_t p2 = a/(2*sqrt(bkg[j]));

      sigma_punzi[j] = sigma_bkg[j] * (p1+p2);
      sigma_total[j] = (punzi[j]/sig_eff[j])*sqrt( pow(sigma_punzi[j]/punzi[j], 2) + pow(sigma_eff[j]/sig_eff[j], 2) );
   }
/*   cout << "Punzi Smin/SignalEff " << endl;
   cout << "VVVLoose Punzi: " << punzi[0] << " " << sig_eff[0] << " " << punzi[0]/sig_eff[0] << endl;
   cout << "VVLoose Punzi: " << punzi[1] << " "  << sig_eff[1] << " " << punzi[1]/sig_eff[1] << endl;
   cout << "VLoose Punzi: " << punzi[2] << " "  << sig_eff[2] << " " << punzi[2]/sig_eff[2] << endl;
   cout << "Loose Punzi: " << punzi[3] << " " << sig_eff[3] << " " << punzi[3]/sig_eff[3] << endl;
   cout << "Medium Punzi: " << punzi[4] << " " << sig_eff[4] << " " << punzi[4]/sig_eff[4] << endl;
   cout << "Tight Punzi: " << punzi[5] << " " << sig_eff[5] << " " << punzi[5]/sig_eff[5] << endl; */

   cout << "Punzi Smin/SignalEff " << endl;
   cout << "VVVLoose: " << punzi[0]/sig_eff[0] << endl;
   cout << "VVLoose:  " << punzi[1]/sig_eff[1] << endl;
   cout << "VLoose:  " << punzi[2]/sig_eff[2] << endl;
   cout << "Loose: " << punzi[3]/sig_eff[3] << endl;
   cout << "Medium: "<< punzi[4]/sig_eff[4] << endl;
   cout << "Tight: "<< punzi[5]/sig_eff[5] << endl;
   cout << "No tau veto: "<< punzi[6]/1. << endl;

   cout << "Error on Punzi" << endl;
   cout << "VVVLoose: " << sigma_punzi[0] << "  " << sigma_eff[0] << "  " << sigma_total[0] << endl;
   cout << "VVLoose:  " << sigma_punzi[1] << "  " << sigma_eff[1] << "  " << sigma_total[1] << endl;
   cout << "VLoose:  " << sigma_punzi[2] << "  " << sigma_eff[2] <<  "  " << sigma_total[2] << endl;
   cout << "Loose: " << sigma_punzi[3] << "  " << sigma_eff[3] <<  "  " << sigma_total[3] <<endl;
   cout << "Medium: "<< sigma_punzi[4] << "  " << sigma_eff[4] <<  "  " << sigma_total[4] <<endl;
   cout << "Tight: "<< sigma_punzi[5] << "  " << sigma_eff[5] <<  "  " << sigma_total[5] <<endl;
   cout << "No tau veto: "<< sigma_punzi[6] << "  " << sigma_eff[6] <<  "  " << sigma_total[6] <<endl;
   
   TGraphErrors* roc = new TGraphErrors(npoints, sig_eff, bkg_rej, sigma_eff, sigma_rej);
   //TGraph* roc = new TGraph(npoints, sig_eff, bkg_rej);

   roc->SetTitle("");
   roc->GetXaxis()->SetTitle("Signal Efficiency (%)");
   roc->GetYaxis()->SetTitle("Background Rejection (%)");
   roc->SetMarkerStyle(20);
   roc->Draw("AP");

   Int_t colors[npoints] = {kBlack, kRed, kMagenta, kGreen, kBlue, kOrange};
   TMarker *m[npoints];
   for (Int_t i=0; i<npoints; i++) {
      m[i] = new TMarker(sig_eff[i], bkg_rej[i], 20);
      //m->SetMarkerSize(2);
      m[i]->SetMarkerColor(colors[i]);
      m[i]->Draw();
   }

   TLegend* leg = new TLegend(0.75, 0.7, 0.9, 0.891);
   leg->SetBorderSize(0);
   leg->SetEntrySeparation(0.01);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(m[0],"VVVLoose","p");
   leg->AddEntry(m[1],"VVLoose","p");
   leg->AddEntry(m[2],"VLoose","p");
   leg->AddEntry(m[3],"Loose","p");
   leg->AddEntry(m[4],"Medium","p");
   leg->AddEntry(m[5],"Tight","p");
   leg->Draw("SAME");
   c1->SaveAs("roc.png");
  
}


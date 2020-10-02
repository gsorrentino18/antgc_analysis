#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TArrow.h"
#include "TPaveText.h"
#include "TLine.h"
#include <TF1.h>
#include "TFractionFitter.h"
#include <math.h>
#include <Fit/Fitter.h>
#include <iostream>
#include <fstream>


std::string 													optionsFile 			=	"estimateFakeOptions.txt";
parseOptions 													options;
std::string 													writeString = "";

TH1F* getDataDenominator(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax);
TH1F* getDataNumerator(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax);
TH1F* getDataSideband(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax);
TH1F* getMCSideband(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax);
TH1F* getPromptTemplate(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax);
void fitNumerator(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax);


TH1F* getDataDenominator(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax){

	std::string 												histName				= 		"data_denominator_pT_" + removeTrailingZeros(pTmin) + "to" + removeTrailingZeros(pTmax) + "_eta_" + removeTrailingZeros(etaMin) + "to" + removeTrailingZeros(etaMax);
	histName 																			=		findAndReplaceAll(histName, ".", "p");
	std::string 												writeFile 				=		options.get("saveDir") + histName + ".root";

	TH1F*														denomHist 					= 		(TH1F*) getHistFromFile(histName, writeFile, 1);

	if(denomHist != nullptr) 																		return denomHist;

	denomHist 																				= 		new TH1F(histName.c_str(), ";BDT Score;", options.getInt("bdtNbins"), 0., 1.);
	denomHist->Sumw2();

	TChain*                                     				inTree    				=		openTChain(std::vector<std::string>({options.get("denomDir") + options.get("dataFile")}), options.get("inTreeName"));
	TTreeReader                             					inputTTreeReader(inTree);
	TTreeReaderAnyValue<Float_t>								phoPt_							(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Float_t>								phoSCeta_						(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Float_t>								phoBDTpred_						(inputTTreeReader, "phoBDTpred");	

	while(inputTTreeReader.Next()){

		Float_t 		absEta = std::abs(phoSCeta_);

		if(absEta >= etaMax || absEta < etaMin) continue;
		if(phoPt_ >= pTmax || phoPt_ < pTmin) continue;

		denomHist->Fill(phoBDTpred_);
	}

	closeTChain(inTree);

	writeToFile(denomHist, writeFile, "RECREATE");

	return denomHist;
};


TH1F* getDataNumerator(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax){

	std::string 												histName				= 		"data_numerator_pT_" + removeTrailingZeros(pTmin) + "to" + removeTrailingZeros(pTmax) + "_eta_" + removeTrailingZeros(etaMin) + "to" + removeTrailingZeros(etaMax);
	histName 																			=		findAndReplaceAll(histName, ".", "p");
	std::string 												writeFile 				=		options.get("saveDir") + histName + ".root";

	TH1F*														numHist 					= 		(TH1F*) getHistFromFile(histName, writeFile, 1);

	if(numHist != nullptr) 																		return numHist;

	numHist 																				= 		new TH1F(histName.c_str(), ";BDT Score;", options.getInt("bdtNbins"), 0., 1.);
	numHist->Sumw2();

	TChain*                                     				inTree    				=		openTChain(std::vector<std::string>({options.get("numDir") + options.get("dataFile")}), options.get("inTreeName"));
	TTreeReader                             					inputTTreeReader(inTree);
	TTreeReaderAnyValue<Float_t>								phoPt_							(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Float_t>								phoSCeta_						(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Float_t>								phoBDTpred_						(inputTTreeReader, "phoBDTpred");	

	while(inputTTreeReader.Next()){

		Float_t 		absEta = std::abs(phoSCeta_);

		if(absEta >= etaMax || absEta < etaMin) continue;
		if(phoPt_ >= pTmax || phoPt_ < pTmin) continue;

		numHist->Fill(phoBDTpred_);
	}

	// closeTChain(inTree);

	writeToFile(numHist, writeFile, "RECREATE");

	return numHist;
};


TH1F* getDataSideband(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax){

	std::string 												histName				= 		"data_sideband_pT_" + removeTrailingZeros(pTmin) + "to" + removeTrailingZeros(pTmax) + "_eta_" + removeTrailingZeros(etaMin) + "to" + removeTrailingZeros(etaMax);
	histName 																			=		findAndReplaceAll(histName, ".", "p");
	std::string 												writeFile 				=		options.get("saveDir") + histName + ".root";

	TH1F*														sbHist 					= 		(TH1F*) getHistFromFile(histName, writeFile, 1);

	if(sbHist != nullptr) 																		return sbHist;

	sbHist 																				= 		new TH1F(histName.c_str(), ";BDT Score;", options.getInt("bdtNbins"), 0., 1.);
	sbHist->Sumw2();

	TChain*                                     				inTree    				=		openTChain(std::vector<std::string>({options.get("sidebandDir") + options.get("dataFile")}), options.get("inTreeName"));
	TTreeReader                             					inputTTreeReader(inTree);
	TTreeReaderAnyValue<Float_t>								phoPt_							(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Float_t>								phoSCeta_						(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Float_t>								phoBDTpred_						(inputTTreeReader, "phoBDTpred");	

	while(inputTTreeReader.Next()){

		Float_t 		absEta = std::abs(phoSCeta_);

		if(absEta >= etaMax || absEta < etaMin) continue;
		if(phoPt_ >= pTmax || phoPt_ < pTmin) continue;

		sbHist->Fill(phoBDTpred_);
	}

	// closeTChain(inTree);

	writeToFile(sbHist, writeFile, "RECREATE");

	return sbHist;
};

TH1F* getMCSideband(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax){
	std::string 												histName				= 		"mc_sideband_pT_" + removeTrailingZeros(pTmin) + "to" + removeTrailingZeros(pTmax) + "_eta_" + removeTrailingZeros(etaMin) + "to" + removeTrailingZeros(etaMax);
	histName 																			=		findAndReplaceAll(histName, ".", "p");
	std::string 												writeFile 				=		options.get("saveDir") + histName + ".root";

	TH1F*														sbHist 				= 		(TH1F*) getHistFromFile(histName, writeFile, 1);

	if(sbHist != nullptr) 															return sbHist;

	sbHist 																			= 		new TH1F(histName.c_str(), ";BDT Score;", options.getInt("bdtNbins"), 0., 1.);
	sbHist->Sumw2();

	std::vector<std::string>										samples 			=		options.getList("GJetsBins");
	for(std::string iSample : samples){
		std::string 							iSamplePath 								=	options.get("sidebandDir") + iSample;
		TChain*									tChain 										=	openTChain((std::vector<std::string>){iSamplePath}, options.get("inTreeName"));
		TTreeReader                             inputTTreeReader(tChain);
		TTreeReaderAnyValue<Float_t>			genWeight_										(inputTTreeReader, "genWeight");
		TTreeReaderAnyValue<Float_t>			puWeight_										(inputTTreeReader, "puWeight");
		TTreeReaderAnyValue<Float_t>			phoSCeta_										(inputTTreeReader, "phoSCeta");
		TTreeReaderAnyValue<Float_t>			phoPt_											(inputTTreeReader, "phoPt");		
		TTreeReaderAnyValue<Float_t>			phoBDTpred_										(inputTTreeReader, "phoBDTpred");	

		Float_t 								sumGenWeight 								= 	std::stof(vLookup(findAndReplaceAll(iSample, ".root", ""), options.get("xSectionMap"), 0, 7));
		Float_t 								xSection 									=	std::stof(vLookup(findAndReplaceAll(iSample, ".root", ""), options.get("xSectionMap"), 0, 2));

		while(inputTTreeReader.Next()){
			
			Float_t 		absEta = std::abs(phoSCeta_);

			if(absEta >= etaMax || absEta < etaMin) continue;
			if(phoPt_ >= pTmax || phoPt_ < pTmin) continue;

			Double_t 							weight 										=  	genWeight_ * puWeight_ * xSection * 1000./sumGenWeight;
			
			sbHist->Fill(phoBDTpred_, weight);
		}
		// closeTChain(tChain);
	}

	writeToFile(sbHist, writeFile, "RECREATE");

	return sbHist;
};


TH1F* getPromptTemplate(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax){

	std::string 												histName				= 		"mc_prompt_pT_" + removeTrailingZeros(pTmin) + "to" + removeTrailingZeros(pTmax) + "_eta_" + removeTrailingZeros(etaMin) + "to" + removeTrailingZeros(etaMax);
	histName 																			=		findAndReplaceAll(histName, ".", "p");
	std::string 												writeFile 				=		options.get("saveDir") + histName + ".root";

	TH1F*														promptHist 				= 		(TH1F*) getHistFromFile(histName, writeFile, 1);

	if(promptHist != nullptr) 															return promptHist;

	promptHist 																			= 		new TH1F(histName.c_str(), ";BDT Score;", options.getInt("bdtNbins"), 0., 1.);
	promptHist->Sumw2();

	std::vector<std::string>										samples 			=		options.getList("GJetsBins");
	for(std::string iSample : samples){
		std::string 							iSamplePath 								=	options.get("numDir") + iSample;
		TChain*									tChain 										=	openTChain((std::vector<std::string>){iSamplePath}, options.get("inTreeName"));
		TTreeReader                             inputTTreeReader(tChain);
		TTreeReaderAnyValue<Float_t>			genWeight_										(inputTTreeReader, "genWeight");
		TTreeReaderAnyValue<Float_t>			puWeight_										(inputTTreeReader, "puWeight");
		TTreeReaderAnyValue<Float_t>			phoSCeta_										(inputTTreeReader, "phoSCeta");
		TTreeReaderAnyValue<Float_t>			phoPt_											(inputTTreeReader, "phoPt");		
		TTreeReaderAnyValue<Float_t>			phoBDTpred_										(inputTTreeReader, "phoBDTpred");	

		Float_t 								sumGenWeight 								= 	std::stof(vLookup(findAndReplaceAll(iSample, ".root", ""), options.get("xSectionMap"), 0, 7));
		Float_t 								xSection 									=	std::stof(vLookup(findAndReplaceAll(iSample, ".root", ""), options.get("xSectionMap"), 0, 2));

		while(inputTTreeReader.Next()){
			
			Float_t 		absEta = std::abs(phoSCeta_);

			if(absEta >= etaMax || absEta < etaMin) continue;
			if(phoPt_ >= pTmax || phoPt_ < pTmin) continue;

			Double_t 							weight 										=  	genWeight_ * puWeight_ * xSection * 1000./sumGenWeight;
			
			promptHist->Fill(phoBDTpred_, weight);
		}
		// closeTChain(tChain);
	}

	writeToFile(promptHist, writeFile, "RECREATE");

	return promptHist;
};

void fitNumerator(Float_t pTmin, Float_t pTmax, Float_t etaMin, Float_t etaMax){

	gROOT->SetBatch();
	gStyle->SetOptStat(0);

	std::string 												fitName				= 		"numerator_fit_pT_" + removeTrailingZeros(pTmin) + "to" + removeTrailingZeros(pTmax) + "_eta_" + removeTrailingZeros(etaMin) + "to" + removeTrailingZeros(etaMax);

	TH1F* 										dataNumHist 								=	getDataNumerator(pTmin, pTmax, etaMin, etaMax);
	dataNumHist->SetLineColor(options.getTColFromHex("dataNumCol"));
	dataNumHist->SetLineWidth(options.getFloat("lineWidth"));
	dataNumHist->SetMarkerColor(options.getTColFromHex("dataNumCol"));
	dataNumHist->SetMarkerStyle(options.getInt("dataNumMkrStyle"));
	dataNumHist->SetMarkerSize(options.getFloat("dataNumMkrSize"));

	TH1F* 										mcPrommptHist 								=	getPromptTemplate(pTmin, pTmax, etaMin, etaMax);
	mcPrommptHist->Scale(dataNumHist->Integral()/mcPrommptHist->Integral());
	mcPrommptHist->SetLineColor(options.getTColFromHex("mcNumCol"));
	mcPrommptHist->SetLineWidth(options.getFloat("lineWidth"));

	TH1F* 										dataSidebandHist							=	getDataSideband(pTmin, pTmax, etaMin, etaMax);
	TH1F* 										mcSidebandHist 								=	getMCSideband(pTmin, pTmax, etaMin, etaMax);
	mcSidebandHist->Scale(options.getFloat("luminosity"));
	dataSidebandHist->Add(mcSidebandHist, -1.);
	dataSidebandHist->Scale(dataNumHist->Integral()/dataSidebandHist->Integral());
	dataSidebandHist->SetLineColor(options.getTColFromHex("dataSidebandCol"));
	dataSidebandHist->SetLineWidth(options.getFloat("lineWidth"));

	TObjArray templates;
	templates.Add(mcPrommptHist);
	templates.Add(dataSidebandHist);

	TFractionFitter* ffitter 	=	new TFractionFitter(dataNumHist, &templates, "V");

	std::vector<ROOT::Fit::ParameterSettings> & fitParameters = ffitter->GetFitter()->Config().ParamsSettings();
	fitParameters[0].Set("prompt", options.getFloat("promptInit"), options.getDouble("fitStepSize"), 0., 1.);
	fitParameters[1].Set("fake", 1. - options.getFloat("promptInit"), options.getDouble("fitStepSize"), 0., 1.);

	Int_t fstatus = ffitter->Fit();
	std::cout <<"\tFit status: "<<fstatus<<std::endl;

	TH1F * fitResult = (TH1F*) ffitter->GetPlot()->Clone(fitName.c_str());
	fitResult->SetLineColor(options.getTColFromHex("fitCol"));
	fitResult->SetLineWidth(options.getFloat("lineWidth"));
	fitResult->SetLineStyle(options.getInt("fitLineStyle"));

	

	Double_t promptFraction, promptFractionError;
	Double_t fakeFraction, fakeFractionError;

	ffitter->GetResult(0, promptFraction, promptFractionError);
	ffitter->GetResult(1, fakeFraction, fakeFractionError);

	TLegend 									legend(options.getDouble("legx1"), options.getDouble("legy1"), options.getDouble("legx2"), options.getDouble("legy2"));

	std::string 								legStr 				= 		removeTrailingZeros(etaMin) + "#leq|#eta|<" +removeTrailingZeros(etaMax) + "  " + removeTrailingZeros(pTmin) + "#leqp_{T}<" + removeTrailingZeros(pTmax);
	legend.SetHeader(legStr.c_str(), "C");

	legend.SetTextSize(options.getDouble("legTextSize"));
	legend.SetNColumns(options.getInt("legNcols"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetLineColor(options.getTColFromHex("legFillColor"));
	legend.SetBorderSize(options.getInt("legBorderWidth"));
	// legend.SetMargin(options.getFloat("legMargin"));

	legStr 															= 	"Control numerator (Data, N=" + removeTrailingZeros(dataNumHist->GetEntries()) + ")";
	legend.AddEntry(dataNumHist, legStr.c_str(), "LPE");

	legStr 															= 	"Prompt template (MC, " + to_string_with_precision(promptFraction*100., options.getInt("fracPrecision")) + "#pm" + to_string_with_precision(promptFractionError*100., options.getInt("fracPrecision")) + " %)";
	legend.AddEntry(mcPrommptHist, legStr.c_str(), "L");

	legStr 															= 	"Fake template (Data, " + to_string_with_precision(fakeFraction*100., options.getInt("fracPrecision")) + "#pm" + to_string_with_precision(fakeFraction*100., options.getInt("fracPrecision")) + " %)";
	legend.AddEntry(dataSidebandHist, legStr.c_str(), "L");

	legStr 															= 	"Fit (#chi^{2}/dof=" + to_string_with_precision(ffitter->GetChisquare(), options.getInt("fitLegendPrecision1")) +
	+ "/" + removeTrailingZeros(ffitter->GetNDF()) +	", P=" + to_string_with_precision(ffitter->GetProb() * 100., options.getInt("fitLegendPrecision2")) + "%)";
	legend.AddEntry(fitResult, legStr.c_str(), "L");


	mcPrommptHist->Scale(promptFraction);
	dataSidebandHist->Scale(fakeFraction);

	THStack 								hStack((fitName + "_stack").c_str(), "");
	hStack.Add(mcPrommptHist, "HIST");
	hStack.Add(dataSidebandHist, "HIST");
	hStack.Add(fitResult, "HIST");
	hStack.Add(dataNumHist, "PE");


	TCanvas canvas((fitName+"_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);

	TPad 										pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad0.SetFillStyle(4000);
	pad0.SetFillColor(0);
	pad0.SetFrameFillStyle(4000);
	pad0.SetGrid(1,1);

	canvas.Draw();
	canvas.cd();
	pad0.Draw();
	pad0.cd();

	pad0.cd();

	hStack.Draw("nostack");
	legend.Draw();

	hStack.GetXaxis()->SetTitle("BDT Score");
	hStack.GetYaxis()->SetTitle("Events");
	hStack.GetXaxis()->CenterTitle();
	hStack.GetYaxis()->CenterTitle();
	hStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	hStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	hStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	hStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	hStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	hStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	hStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	hStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	std::string writePath 									= 	options.get("saveDir") + "/" + fitName;
	canvas.SaveAs((writePath+".png").c_str());
	canvas.SaveAs((writePath+".pdf").c_str());


	pad0.Clear();
	// hStack.Clear();
	legend.Clear();

	TH1F* 		denomHist 											=	getDataNumerator(pTmin, pTmax, etaMin, etaMax);
	denomHist->SetLineColor(options.getTColFromHex("dataDenomCol"));
	denomHist->SetLineWidth(options.getFloat("lineWidth"));

	mcSidebandHist->SetLineWidth(options.getFloat("lineWidth"));

	std::string 						ratioPlotName 				=	"fakeRatio_pT_" + removeTrailingZeros(pTmin) + "to" + removeTrailingZeros(pTmax) + "_eta_" + removeTrailingZeros(etaMin) + "to" + removeTrailingZeros(etaMax);

	THStack 								h2Stack((ratioPlotName + "_stack").c_str(), "");
	h2Stack.Add(denomHist, "HIST");
	h2Stack.Add(dataSidebandHist, "HIST");

	Float_t fakeN = dataSidebandHist->Integral(dataSidebandHist->GetXaxis()->FindBin(options.getFloat("bdtThres")), dataSidebandHist->GetNbinsX());
	Float_t fakeNerr = fakeN * fakeFractionError/fakeFraction;

	legStr 															= 	"Control numerator fake =" + to_string_with_precision(fakeN, options.getInt("denomPrecision")) +
	+ "#pm" + to_string_with_precision(fakeNerr, options.getInt("denomPrecision"));
	legend.AddEntry(dataSidebandHist, legStr.c_str(), "L");

	legStr 															= 	"Control denominator =" + removeTrailingZeros(denomHist->GetEntries());
	legend.AddEntry(denomHist, legStr.c_str(), "L");


	Float_t 				fakeRatio								=	fakeN/denomHist->GetEntries();
	Float_t 				fakeRatioerr							=	fakeRatio * std::sqrt(std::pow(fakeFractionError/fakeFraction, 2) + std::pow(1./denomHist->GetEntries(), 2));
	legStr 															= 	"Fake ratio =" + to_string_with_precision(fakeRatio*100., options.getInt("fakeRatioPrecision")) + "#pm" +to_string_with_precision(fakeRatioerr * 100., options.getInt("fakeRatioPrecision")) +" %";

	legend.AddEntry(&canvas, legStr.c_str(), "");

	legStr 				= 		removeTrailingZeros(etaMin) + "#leq|#eta|<" +removeTrailingZeros(etaMax) + "  " + removeTrailingZeros(pTmin) + "#leqp_{T}<" + removeTrailingZeros(pTmax);
	legend.SetHeader(legStr.c_str(), "C");

	pad0.cd();
	h2Stack.Draw("nostack");
	legend.Draw();

	h2Stack.GetXaxis()->SetTitle("BDT Score");
	h2Stack.GetYaxis()->SetTitle("Events");
	h2Stack.GetXaxis()->CenterTitle();
	h2Stack.GetYaxis()->CenterTitle();
	h2Stack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	h2Stack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	h2Stack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	h2Stack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	h2Stack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	h2Stack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	h2Stack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	h2Stack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	writePath 									= 	options.get("saveDir") + "/" + ratioPlotName;
	canvas.SaveAs((writePath+".png").c_str());
	canvas.SaveAs((writePath+".pdf").c_str());

	// clearHeap();


	writeString 			+=		removeTrailingZeros(pTmin) + ", " + removeTrailingZeros(pTmax) + ", " + removeTrailingZeros(etaMin) + ", " + removeTrailingZeros(etaMax) + ", " + removeTrailingZeros(fakeRatio) +  ", " + removeTrailingZeros(fakeRatioerr) + "\n";
};


void init(){

	writeString="pTmin,pTmax,etaMin,etaMax,Fr,FrErr\n";

	options.parseIt(optionsFile, "==",1);

	std::vector<std::string> etaBins 																				=	options.getList("etaBins", ",");
	std::vector<std::string> pTbins																					= 	options.getList("pTbins", ",");


	for(std::string etaBin : etaBins){

		std::vector<Float_t> etaBinF = strToFloatList(etaBin, "-");

		for(std::string ptBin : pTbins){

			std::vector<Float_t> ptBinF = strToFloatList(ptBin, "-");

			fitNumerator(ptBinF[0], ptBinF[1], etaBinF[0], etaBinF[1]);
		}		
	}

	ofstream frOutFile(options.get("saveDir") + options.get("fakeRatioOutFile"));
	frOutFile << writeString;
	frOutFile.close();

};
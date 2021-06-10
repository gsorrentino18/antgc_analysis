#include "TEventList.h"
#include "TF1.h"
#include "TSQLResult.h"
#include <boost/algorithm/string/find.hpp>
#include "/local/cms/user/gsorrent/antgc_analysis/wgCRstudies/ROOTmacros/extra_tools.cc"

// std::vector<std::string> cuts = {	"No Cuts",	"HLT_Photon200_v*",	"goodVertices",	"globalSuperTightHalo2016Filter",	"HBHENoiseFilter",	"HBHENoiseIsoFilter",	"EcalDeadCellTriggerPrimitiveFilter",	"BadPFMuonFilter",	"ecalBadCalibFilter",	"|#eta_{SC}^{#gamma}|<1.4442",	"p_{T}^{#gamma}>200 GeV",	"H/E<0.0421",	"ECALClusIso<5.11 GeV",	"HCALClusIso<17.02 GeV",	"TkrIso<3.50",	"BDT>0.49",	"MIP_{tot}<4.9 GeV",	"Pixel Veto",	">0 #gamma candidate",	"e(tight,p_{T}>30 GeV,|#eta|<2.5)",	"1 e candidate", "-",	"1 #gamma candidate",	"#slash{E}_{T}>50 GeV",	"p_{T}^{#gamma}>220 GeV",	"m_T(e+#slash{E}_{T}})<160 GeV",	"p_{T}^{#gamma}/p_{T}^{(e+#slash{E}_{T})}}<1.4",	"#Delta#phi(#gamma, e+#slash{E}_{T})>0.5",	"#Delta#phi_{min}(#slash{E}_{T},jet)>0.5",	"Veto #mu(p_{T}>10 GeV,loose ID+Iso,|#eta|<2.5)",	"Veto #tau(p_{T}>20 GeV,loose Deep #tau,|#eta|<2.4)"}

std::string 													optionsFile 			=	"plotDataMCoptions.txt";
parseOptions 													options;
ofstream 														yieldsFile;
TH1F* 															getTotalMCBgHist();
void 															createBaysBlocksBinning();
TH1F* 															getJetFakesHist();
TH1F* 															getEleFakesHist();
TH1F* 															getDataHist();
TH1F* 															makePlot();
void 															init();


TH1F* 															getTotalMCBgHist() {

	Float_t 													lumi 					=		options.getFloat("lumi");
	Bool_t 														applyEWkfac 			= 		options.getInt("applyEWkfac");
	Bool_t 														applyQCDkfac 			= 		options.getInt("applyQCDkfac");
	Bool_t 														applyGJetsKfac 			= 		options.getInt("applyGJetsKfac");
	Float_t														GJetsKfac 				= 		options.getFloat("GJetsKfac");

	std::vector<Double_t>										binEdges 				=		strToDoubleList(options.getList("var", ";")[2], ",");
	Bool_t 														useVariableBinning 		=		(binEdges.size() > 3);
	if (binEdges.size() < 3) {
		std::cout << "Error! Invalid binning " << options.get("var") << std::endl;
		return nullptr;
	}

	TH1F*														totalMCBgHist 			= 		useVariableBinning ? new TH1F("totalMCBgHist", "Total MC Background", binEdges.size() - 1, binEdges.data())
	        : new TH1F("totalMCBgHist", "Total MC Background", std::round(binEdges[0]), binEdges[1], binEdges[2]);

	Bool_t 														overFlowInLastXbin 		=		options.getInt("overFlowInLastXbin");
	Double_t 													xLastBinLowEdge			=		totalMCBgHist->GetXaxis()->GetBinLowEdge(totalMCBgHist->GetNbinsX());

	for (const std::string & iSample 	: 		options.getList("bgBins", ";")) {

		for (const std::string & iSampleBin 	: 		split_string(split_string(iSample, ":")[0])) {

			std::string                             				sampleBinName     	    =		findAndReplaceAll(getFileName(iSampleBin), ".root", "");
			std::string 											sampleBinPath 			= 		options.get("SRdir")  + "/" +iSampleBin;

			std::string 											xSecStr 				=		vLookup(sampleBinName, options.get("xSectionMap"), 0, 2);
			trim(xSecStr);
			if (xSecStr.empty()) {
				std::cout<<"Error! Cross section not found for sample "<<sampleBinName<<std::endl;
				return nullptr;
			}
			Float_t 												xSection 				=		std::stod(xSecStr);
			TH1F* 													cutFlowGenWeight  		=		(TH1F*) getHistFromFile(options.get("cutFlowHist"), sampleBinPath);
			Double_t		 										sumGenWeight   			=		cutFlowGenWeight->GetBinContent(1);
			delete cutFlowGenWeight;

			TChain*                                 				sampleBinTree          	=		openTChain((std::vector<std::string>) {sampleBinPath}, options.get("inTree"));
			if ( sampleBinTree->GetEntries() == 0) {
				std::cout<<"No entries in "<<sampleBinName<<"! Skipping..."<<std::endl;
				continue;
			}
			sampleBinTree->Draw(">>preSelection", options.getCSTR("mcCut"), "goff");
			TEventList*												preSelection 			= 		(TEventList*) gDirectory->Get("preSelection");

			sampleBinTree->SetBranchStatus("*",0);
			Float_t 												puWeight;
			Float_t 												genWeight;
			Float_t 												lhePhoPt;
			Char_t 													Wsign;
			Float_t 												theVariable;
			sampleBinTree->SetBranchStatus("puWeight", 1);
			sampleBinTree->SetBranchStatus("genWeight", 1);
			sampleBinTree->SetBranchStatus("lhePhoPt", 1);
			sampleBinTree->SetBranchStatus("Wsign", 1);
			sampleBinTree->SetBranchStatus(options.getList("var", ";")[0].c_str(), 1);
			sampleBinTree->SetBranchAddress("puWeight", &puWeight);
			sampleBinTree->SetBranchAddress("genWeight", &genWeight);
			sampleBinTree->SetBranchAddress("lhePhoPt", &lhePhoPt);
			sampleBinTree->SetBranchAddress("Wsign", &Wsign);
			sampleBinTree->SetBranchAddress(options.getList("var", ";")[0].c_str(), &theVariable);


			kFactorMap 													kFacQCDZNuNuG;
			kFactorMap 													kFacEWZNuNuG;
			kFactorMap 													kFacQCDWpLNuG;
			kFactorMap 													kFacEWWpLNuG;
			kFactorMap 													kFacQCDWmLNuG;
			kFactorMap 													kFacEWWmLNuG;
			Bool_t 													applyZNuNuGKfac			=		match("*ZGTo2NuG*", sampleBinName);
			Bool_t 													applyWLNuGKfac			=		match("*WGToLNuG*", sampleBinName);
			Bool_t 													applyGJetsKfacToSample	=		match("GJets*", sampleBinName);

			if (applyZNuNuGKfac) {
				if (applyQCDkfac) 	kFacQCDZNuNuG.init(options.getList("kFacQCDZNuNuG")[0], options.getList("kFacQCDZNuNuG")[1]);
				if (applyEWkfac) 	kFacEWZNuNuG.init(options.getList("kFacEWZNuNuG")[0], options.getList("kFacEWZNuNuG")[1]);
			} else if (applyWLNuGKfac) {
				if (applyQCDkfac) 	kFacQCDWpLNuG.init(options.getList("kFacQCDWpLNuG")[0], options.getList("kFacQCDWpLNuG")[1]);
				if (applyEWkfac) 	kFacEWWpLNuG.init(options.getList("kFacEWWpLNuG")[0], options.getList("kFacEWWpLNuG")[1]);

				if (applyQCDkfac) 	kFacQCDWmLNuG.init(options.getList("kFacQCDWmLNuG")[0], options.getList("kFacQCDWmLNuG")[1]);
				if (applyEWkfac) 	kFacEWWmLNuG.init(options.getList("kFacEWWmLNuG")[0], options.getList("kFacEWWmLNuG")[1]);
			}

			for (Int_t iEvt = 0; iEvt < preSelection->GetN(); iEvt++) {
				Long64_t 									iEntry				= 	preSelection->GetEntry(iEvt);
				sampleBinTree->GetEntry(iEntry);

				Float_t 									kFac 				=	1.;

				if (applyZNuNuGKfac) {
					kFac 														= 	kFacQCDZNuNuG.getKFac(lhePhoPt) * kFacEWZNuNuG.getKFac(lhePhoPt);
				} else if (applyWLNuGKfac) {
					kFac 														= 	(Wsign > 0) ? kFacQCDWpLNuG.getKFac(lhePhoPt) * kFacEWWpLNuG.getKFac(lhePhoPt) :
					        kFacQCDWmLNuG.getKFac(lhePhoPt) * kFacEWWmLNuG.getKFac(lhePhoPt);
				}

				if (applyGJetsKfacToSample) 					kFac 				*=	GJetsKfac;

				if (theVariable < 20.) cout<<theVariable<<std::endl;

				Float_t	 									weight 				=  	lumi * puWeight * genWeight * xSection*1000. * kFac/sumGenWeight;

				if (overFlowInLastXbin && (theVariable > xLastBinLowEdge)) totalMCBgHist->Fill(xLastBinLowEdge, weight);
				else totalMCBgHist->Fill(theVariable, weight);
			}

			delete preSelection;
			closeTChain(sampleBinTree);
		}

	}

	return totalMCBgHist;
};


void createBaysBlocksBinning() {

	std::vector<Double_t>										binEdges 				=		strToDoubleList(options.getList("var", ";")[2], ",");
	Bool_t 														useVariableBinning 		=		(binEdges.size() > 3);
	Double_t 													xMin 					=		useVariableBinning ? binEdges[0] : binEdges[1];
	Double_t 													xMax 					=		useVariableBinning ? binEdges.back() : binEdges[2];

	//// temporarily update binning options so that other histogramming functions use binning defined above
	//// options.getList("var", ";")[3] = nBins
	std::string 												tmpBinningStr			=		";" + options.getList("var", ";")[3] + ","+ removeTrailingZeros(xMin) + "," + removeTrailingZeros(xMax);
	std::string 												varOptions 				=		options.get("var");
	boost::iterator_range<string::iterator> binningStrStart = boost::find_nth(varOptions, ";", 1);
	boost::iterator_range<string::iterator> binningStrEnd = boost::find_nth(varOptions, ";", 2);
	varOptions.replace(binningStrStart.begin(), binningStrEnd.begin(), tmpBinningStr.c_str(), tmpBinningStr.length());
	options.optMap["var"] = varOptions;

	TH1F* 														rebinSourceHist 		=		nullptr;
	if (options.getInt("rebinWithData") && options.getInt("includeData")) {
		rebinSourceHist																	=		getDataHist();
	} else {
		rebinSourceHist 																=		getTotalMCBgHist();
		TH1F* 														jetFakesHist 		=		getJetFakesHist();
		TH1F* 														eleFakesHist 		=		getEleFakesHist();
		rebinSourceHist->Add(jetFakesHist);
		rebinSourceHist->Add(eleFakesHist);
		delete jetFakesHist;
		delete eleFakesHist;
	}

	TH1F* 														rebinnedHist 			=		(TH1F*) BayesianBlocks::rebin(rebinSourceHist, options.getDouble("baysBlocksP"));
	delete rebinSourceHist;

	//// update binning options with Bays blocks binning
	tmpBinningStr 																		=		";" + removeTrailingZeros(xMin) + ",";
	for (Int_t i = 1; i <= rebinnedHist->GetNbinsX(); i++) {
		if (i == 1 && std::abs(rebinnedHist->GetXaxis()->GetBinLowEdge(i) - xMin)<1.e-20) continue;
		tmpBinningStr 																	+=		removeTrailingZeros(rebinnedHist->GetXaxis()->GetBinLowEdge(i)) + ",";
	}
	tmpBinningStr 																		+=		removeTrailingZeros(xMax);

	varOptions 																			=		options.get("var");
	binningStrStart 																	= 		boost::find_nth(varOptions, ";", 1);
	binningStrEnd 																		= 		boost::find_nth(varOptions, ";", 2);
	varOptions.replace(binningStrStart.begin(), binningStrEnd.begin(), tmpBinningStr.c_str(), tmpBinningStr.length());
	options.optMap["var"] = varOptions;

	std::cout<<"******************************************************************************************************************************************************************************************"<<std::endl;
	std::cout<<"\t\tBinning: "<<tmpBinningStr<<std::endl;
	std::cout<<"******************************************************************************************************************************************************************************************"<<std::endl;

	delete rebinnedHist;
};


TH1F* 			getJetFakesHist() {

	std::string 												histName				= 		"jetFakes_" + options.getList("var", ";")[0];
	std::string                             					sampleName     	     	=		findAndReplaceAll(getFileName(options.getList("data")[0]), ".root", "");
	std::string 												writeFile 				=		options.get("saveDir") + "/" + histName + ".root";
	TH1F*														jetFakesHist 			= 		options.getInt("useExisting") ? (TH1F*) getHistFromFile(histName, writeFile, 1) : nullptr;
	if (jetFakesHist) return jetFakesHist;

	TF1 														jetFakeRatio("jetFakeRatio", options.getList("jetFakesEqn")[0].c_str(), 0., 3000000.);

	std::vector<Double_t>										binEdges 				=		strToDoubleList(options.getList("var", ";")[2], ",");
	Bool_t 														useVariableBinning 		=		(binEdges.size() > 3);
	if (binEdges.size() < 3) {
		std::cout << "Error! Invalid binning " << options.get("var") << std::endl;
		return nullptr;
	}

	std::string 												titles 					=		"Jet Fakes ;" + options.getList("var", ";")[1] + ";Events";
	jetFakesHist 																		= 		useVariableBinning ? new TH1F(histName.c_str(), titles.c_str(), binEdges.size() - 1, binEdges.data())
	        : new TH1F(histName.c_str(), titles.c_str(), std::round(binEdges[0]), binEdges[1], binEdges[2]);
	jetFakesHist->Sumw2();
	jetFakesHist->SetFillColorAlpha(hex2rootColor(options.getList("jetFakesEqn")[1]), options.getFloat("histAlpha"));
	jetFakesHist->SetLineColor(hex2rootColor(options.getList("jetFakesEqn")[1]));

	std::string 												samplePath 				= 		options.get("jetFakes")	+	"/"	+	sampleName	+	".root";
	TChain*                                 					sampleTree         		=		openTChain((std::vector<std::string>) {samplePath}, options.get("inTree"));
	sampleTree->Draw(">>preSelection", options.getCSTR("jetFakesCut"), "goff");
	TEventList*													preSelection			= 		(TEventList*) gDirectory->Get("preSelection");

	sampleTree->SetBranchStatus("*", 0);
	Float_t 													theVariable;
	Float_t 													phoPt;
	sampleTree->SetBranchStatus("phoPt", 1);
	sampleTree->SetBranchAddress("phoPt", &phoPt);
	Bool_t 														varIsPhoPt				=		(options.getList("var", ";")[0] == "phoPt");
	if (!varIsPhoPt) {
		sampleTree->SetBranchStatus(options.getList("var", ";")[0].c_str(), 1);
		sampleTree->SetBranchAddress(options.getList("var", ";")[0].c_str(), &theVariable);
	}

	Bool_t 														overFlowInLastXbin 		=		options.getInt("overFlowInLastXbin");
	Double_t 													xLastBinLowEdge			=		jetFakesHist->GetXaxis()->GetBinLowEdge(jetFakesHist->GetNbinsX());

	for (Int_t iEvt = 0; iEvt < preSelection->GetN(); iEvt++) {
		Long64_t 												iEntry					= 	preSelection->GetEntry(iEvt);
		sampleTree->GetEntry(iEntry);
		Float_t 												jetFakeWeight 			=	jetFakeRatio.Eval(phoPt);
		Float_t 												fillVar 				=	varIsPhoPt ? phoPt : theVariable;

		if (overFlowInLastXbin && (fillVar > xLastBinLowEdge)) jetFakesHist->Fill(xLastBinLowEdge, jetFakeWeight);
		else jetFakesHist->Fill(fillVar, jetFakeWeight);
	}

	delete preSelection;
	closeTChain(sampleTree);

	writeToFile(jetFakesHist, writeFile);

	return jetFakesHist;
};


TH1F* 			getEleFakesHist() {

	std::string 												histName				= 		"eleFakes_" + options.getList("var", ";")[0];
	std::string                             					sampleName     	     	=		findAndReplaceAll(getFileName(options.getList("data")[0]), ".root", "");
	std::string 												writeFile 				=		options.get("saveDir") + "/" + histName + ".root";

	TH1F*														iSampleHist 			= 		options.getInt("useExisting") ? (TH1F*) getHistFromFile(histName, writeFile, 1) : nullptr;

	if (iSampleHist) return iSampleHist;

	TF1 														eleFakeRatio("eleFakeRatio", options.getList("eleFakesEqn")[0].c_str(), -1.6, 1.6);

	std::vector<Double_t>										binEdges 				=		strToDoubleList(options.getList("var", ";")[2], ",");

	Bool_t 														useVariableBinning 		=		(binEdges.size() > 3);
	if (binEdges.size() < 3) {
		std::cout << "Error! Invalid binning " << options.get("var") << std::endl;
		return nullptr;
	}

	std::string 												titles 					=		"Electron Fakes ;" + options.getList("var", ";")[1] + ";Events";
	iSampleHist 																		= 		useVariableBinning ? new TH1F(histName.c_str(), titles.c_str(), binEdges.size() - 1, binEdges.data())
	        : new TH1F(histName.c_str(), titles.c_str(), std::round(binEdges[0]), binEdges[1], binEdges[2]);
	iSampleHist->Sumw2();
	iSampleHist->SetFillColorAlpha(hex2rootColor(options.getList("eleFakesEqn")[1]), options.getFloat("histAlpha"));

	iSampleHist->SetLineColor(hex2rootColor(options.getList("eleFakesEqn")[1]));

	std::string 												samplePath 				= 		options.get("eleFakes")	+	"/"	+	sampleName	+	".root";
	TChain*                                 					sampleTree         		=		openTChain((std::vector<std::string>) {samplePath}, options.get("inTree"));
	sampleTree->Draw(">>preSelection", options.getCSTR("eFakesCut"), "goff");
	TEventList*													preSelection			= 		(TEventList*) gDirectory->Get("preSelection");

	sampleTree->SetBranchStatus("*", 0);
	Float_t 													theVariable;
	Float_t 													phoEta;
	sampleTree->SetBranchStatus("phoEta", 1);
	sampleTree->SetBranchAddress("phoEta", &phoEta);

	Bool_t 														varIsPhoEta				=		(options.getList("var", ";")[0] == "phoEta");
	if (!varIsPhoEta) {
		sampleTree->SetBranchStatus(options.getList("var", ";")[0].c_str(), 1);
		sampleTree->SetBranchAddress(options.getList("var", ";")[0].c_str(), &theVariable);
	}

	Bool_t 														overFlowInLastXbin 		=		options.getInt("overFlowInLastXbin");
	Double_t 													xLastBinLowEdge			=		iSampleHist->GetXaxis()->GetBinLowEdge(iSampleHist->GetNbinsX());

	for (Int_t iEvt = 0; iEvt < preSelection->GetN(); iEvt++) {
		Long64_t 												iEntry					= 	preSelection->GetEntry(iEvt);
		sampleTree->GetEntry(iEntry);
		Float_t 												eleFakeWeight 			=	eleFakeRatio.Eval(phoEta);

		Float_t 												fillVar 				=	varIsPhoEta ? phoEta : theVariable;

		if (overFlowInLastXbin && (fillVar > xLastBinLowEdge)) iSampleHist->Fill(xLastBinLowEdge, eleFakeWeight);
		else iSampleHist->Fill(fillVar, eleFakeWeight);
	}

	delete preSelection;
	closeTChain(sampleTree);

	writeToFile(iSampleHist, writeFile);

	return iSampleHist;

};


TH1F* 															getDataHist() {
	std::string 												histName				= 		"data_" + options.getList("var", ";")[0];
	std::string                             					sampleName     	     	=		findAndReplaceAll(getFileName(options.getList("data")[0]), ".root", "");
	std::string 												writeFile 				=		options.get("saveDir") + "/" + histName + ".root";

	TH1F*														iSampleHist 			= 		options.getInt("useExisting") ? (TH1F*) getHistFromFile(histName, writeFile, 1) : nullptr;
	if (iSampleHist) return iSampleHist;

	std::vector<Double_t>										binEdges 				=		strToDoubleList(options.getList("var", ";")[2], ",");

	Bool_t 														useVariableBinning 		=		(binEdges.size() > 3);
	if (binEdges.size() < 3) {
		std::cout << "Error! Invalid binning " << options.get("var") << std::endl;
		return nullptr;
	}

	std::string 												titles 					=		"Data;" + options.getList("var", ";")[1] + ";Events";
	iSampleHist 																		= 		useVariableBinning ? new TH1F(histName.c_str(), titles.c_str(), binEdges.size() - 1, binEdges.data())
	        : new TH1F(histName.c_str(), titles.c_str(), std::round(binEdges[0]), binEdges[1], binEdges[2]);
	iSampleHist->Sumw2();



	iSampleHist->SetMarkerStyle(std::stoi(options.getList("dataAtts")[0]));
	iSampleHist->SetMarkerSize(std::stod(options.getList("dataAtts")[1]));
	iSampleHist->SetMarkerColor(hex2rootColor(options.getList("dataAtts")[3]));
	iSampleHist->SetLineColor(hex2rootColor(options.getList("dataAtts")[3]));
	iSampleHist->SetLineWidth(std::stod(options.getList("dataAtts")[2]));

	std::string 												samplePath 				= 		options.get("SRdir")	+	"/"	+	sampleName	+	".root";
	TChain*                                 					sampleTree         		=		openTChain((std::vector<std::string>) {samplePath}, options.get("inTree"));
	Float_t 													theVariable;
	sampleTree->SetBranchAddress(options.getList("var", ";")[0].c_str(), &theVariable);

	Bool_t 														overFlowInLastXbin 		=		options.getInt("overFlowInLastXbin");
	Double_t 													xLastBinLowEdge			=		iSampleHist->GetXaxis()->GetBinLowEdge(iSampleHist->GetNbinsX());

	TH1F*         												datCutFlowHist  		=		(TH1F*) getHistFromFile(options.get("cutFlowHistUnWeighted"), samplePath);
	Float_t  													initBinCutFlow  		=		datCutFlowHist->GetBinCenter(datCutFlowHist->FindLastBinAbove(1e-20)) + 1.01;
	std::vector<std::string> 									cutList  				=		options.getList("dataCut", "&&");
	Int_t  														nCuts  		 			=		cutList.size();

	for (Int_t iEvt = 0; iEvt < sampleTree->GetEntries(); iEvt++) {

		Int_t 	iEvtPassCuts 															=		0;
		Float_t iEvtBinCutFlow 															= 		initBinCutFlow;

		for (UInt_t jCut = 0; jCut < cutList.size(); jCut++) {
			sampleTree->Draw(">>iEvtPassJcut", cutList[jCut].c_str(), "goff", 1, iEvt);
			TEventList*													iEvtPassJcut			= 		(TEventList*) gDirectory->Get("iEvtPassJcut");
			if (iEvtPassJcut->GetN() == 0) {
				delete iEvtPassJcut;
				break;
			}
			datCutFlowHist->Fill(iEvtBinCutFlow);
			iEvtBinCutFlow 																+= 		1.;
			iEvtPassCuts 																+=		1;
			delete iEvtPassJcut;
		}

		if (iEvtPassCuts != nCuts) continue;

		sampleTree->GetEntry(iEvt);

		if (overFlowInLastXbin && (theVariable > xLastBinLowEdge)) iSampleHist->Fill(xLastBinLowEdge);
		else iSampleHist->Fill(theVariable);
	}

	closeTChain(sampleTree);

	writeToFile(iSampleHist, writeFile);
	writeToFile(datCutFlowHist, writeFile, "UPDATE");

	return iSampleHist;
};


void drawDataMC() {

        std::cout << "Start drawDataMC" << std::endl;

	std::string 												plotName				=		options.getList("var", ";")[0];
        std::cout << "plotName ok" << std::endl;

	Float_t 													lumi 					=		options.getFloat("lumi");

	Bool_t 														applyEWkfac 			= 		options.getInt("applyEWkfac");
	Bool_t 														applyQCDkfac 			= 		options.getInt("applyQCDkfac");
	Bool_t 														applyGJetsKfac 			= 		options.getInt("applyGJetsKfac");
	Float_t														GJetsKfac 				= 		options.getFloat("GJetsKfac");
	Bool_t 														includeData				=		options.getInt("includeData");

	std::vector<Double_t>										binEdges 				=		strToDoubleList(options.getList("var", ";")[2], ",");
	Bool_t 														useVariableBinning 		=		(binEdges.size() > 3);
	if (binEdges.size() < 3) {
		std::cout << "Error! Invalid binning " << options.get("var") << std::endl;
		return;
	}
	Bool_t 														overFlowInLastXbin 		=		options.getInt("overFlowInLastXbin");

	THStack 													hStack("hStack", "");
	TLegend 													legend(options.getDouble("legx1"), options.getDouble("legy1"), options.getDouble("legx2"), options.getDouble("legy2"));
	legend.SetTextSize(options.getDouble("legTextSize"));
	legend.SetNColumns(options.getInt("legNcols"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetLineColor(options.getTColFromHex("legFillColor"));
	legend.SetBorderSize(0);
	legend.SetHeader(options.getCSTR("legTitle"), "C");
	legend.SetTextAlign(12);
	((TLegendEntry*)legend.GetListOfPrimitives()->First())->SetTextSize(options.getFloat("legHeadSize"));

	std::vector<std::pair<Double_t, TH1F*> >						allHists;

	for (const std::string & iSample 	: 		options.getList("bgBins", ";")) {

		std::string 												sampleName				=		removeNonAlpha(split_string(iSample, ":")[1]);
		std::string 												histName				= 		"MC_" + sampleName + "_" +options.getList("var", ";")[0];
		std::string 												writeFile 				=		options.get("saveDir") + "/" + histName + ".root";

		TH1F*														iSampleHist 			= 		options.getInt("useExisting") ? (TH1F*) getHistFromFile(histName, writeFile, 1) : nullptr;

		if (!iSampleHist) {

			std::string 											titles 						=		split_string(iSample, ":")[1] + ";" + options.getList("var", ";")[1] + ";Events";
			iSampleHist																			= 		useVariableBinning ? new TH1F(sampleName.c_str(), titles.c_str(), binEdges.size() - 1, binEdges.data())
			        : new TH1F(sampleName.c_str(), titles.c_str(), std::round(binEdges[0]), binEdges[1], binEdges[2]);
			iSampleHist->Sumw2();
			iSampleHist->SetLineColor(hex2rootColor(split_string(iSample, ":")[2]));
			iSampleHist->SetFillColorAlpha(hex2rootColor(split_string(iSample, ":")[2]), options.getFloat("histAlpha"));

			Double_t 													xLastBinLowEdge			=		iSampleHist->GetXaxis()->GetBinLowEdge(iSampleHist->GetNbinsX());

			for (const std::string & iSampleBin 	: 		split_string(split_string(iSample, ":")[0])) {

				std::string                             				sampleBinName     	    =		findAndReplaceAll(getFileName(iSampleBin), ".root", "");
				std::string 											sampleBinPath 				= 		options.get("SRdir")  + "/" +iSampleBin;
				std::string 											xSecStr 				=		vLookup(sampleBinName, options.get("xSectionMap"), 0, 2);
				trim(xSecStr);
				if (xSecStr.empty()) {
					std::cout<<"Error! Cross section not found for sample "<<sampleBinName<<std::endl;
					return;
				}
				Float_t 												xSection 				=		std::stod(xSecStr);
				TH1F* 													cutFlowGenWeight  		=		(TH1F*) getHistFromFile(options.get("cutFlowHist"), sampleBinPath);
				Double_t		 										sumGenWeight   			=		cutFlowGenWeight->GetBinContent(1);
				delete cutFlowGenWeight;

				TChain*                                 				sampleBinTree          	=		openTChain((std::vector<std::string>) {sampleBinPath}, options.get("inTree"));
				if ( sampleBinTree->GetEntries() == 0) {
					std::cout<<"No entries in "<<sampleBinName<<"! Skipping..."<<std::endl;
					continue;
				}
				sampleBinTree->Draw(">>preSelection", options.getCSTR("mcCut"), "goff");
				TEventList*												preSelection 			= 		(TEventList*) gDirectory->Get("preSelection");

				sampleBinTree->SetBranchStatus("*",0);
				Float_t 												puWeight;
				Float_t 												genWeight;
				Float_t 												lhePhoPt;
				Char_t 													Wsign;
				Float_t 												theVariable;
				sampleBinTree->SetBranchStatus("puWeight", 1);
				sampleBinTree->SetBranchStatus("genWeight", 1);
				sampleBinTree->SetBranchStatus("lhePhoPt", 1);
				sampleBinTree->SetBranchStatus("Wsign", 1);
				sampleBinTree->SetBranchStatus(options.getList("var", ";")[0].c_str(), 1);
				sampleBinTree->SetBranchAddress("puWeight", &puWeight);
				sampleBinTree->SetBranchAddress("genWeight", &genWeight);
				sampleBinTree->SetBranchAddress("lhePhoPt", &lhePhoPt);
				sampleBinTree->SetBranchAddress("Wsign", &Wsign);
				sampleBinTree->SetBranchAddress(options.getList("var", ";")[0].c_str(), &theVariable);


				kFactorMap 													kFacQCDZNuNuG;
				kFactorMap 													kFacEWZNuNuG;
				kFactorMap 													kFacQCDWpLNuG;
				kFactorMap 													kFacEWWpLNuG;
				kFactorMap 													kFacQCDWmLNuG;
				kFactorMap 													kFacEWWmLNuG;
				Bool_t 													applyZNuNuGKfac			=		match("*ZGTo2NuG*", sampleBinName);
				Bool_t 													applyWLNuGKfac			=		match("*WGToLNuG*", sampleBinName);
				Bool_t 													applyGJetsKfacToSample	=		applyGJetsKfac && match("GJets*", sampleBinName);

				//// k-factors != 1 only if kFactorMap is initialized
				if (applyZNuNuGKfac) {
					if (applyQCDkfac) 	kFacQCDZNuNuG.init(options.getList("kFacQCDZNuNuG")[0], options.getList("kFacQCDZNuNuG")[1]);
					if (applyEWkfac) 	kFacEWZNuNuG.init(options.getList("kFacEWZNuNuG")[0], options.getList("kFacEWZNuNuG")[1]);
				} else if (applyWLNuGKfac) {
					if (applyQCDkfac) 	kFacQCDWpLNuG.init(options.getList("kFacQCDWpLNuG")[0], options.getList("kFacQCDWpLNuG")[1]);
					if (applyEWkfac) 	kFacEWWpLNuG.init(options.getList("kFacEWWpLNuG")[0], options.getList("kFacEWWpLNuG")[1]);

					if (applyQCDkfac) 	kFacQCDWmLNuG.init(options.getList("kFacQCDWmLNuG")[0], options.getList("kFacQCDWmLNuG")[1]);
					if (applyEWkfac) 	kFacEWWmLNuG.init(options.getList("kFacEWWmLNuG")[0], options.getList("kFacEWWmLNuG")[1]);
				}

				if (applyGJetsKfacToSample) {
					std::cout<<"Applying GJets k-fac "<<GJetsKfac<<" to sample "<<sampleBinName<<std::endl;
				}

				for (Int_t iEvt = 0; iEvt < preSelection->GetN(); iEvt++) {
					Long64_t 									iEntry				= 	preSelection->GetEntry(iEvt);
					sampleBinTree->GetEntry(iEntry);

					Float_t 									kFac 				=	1.;

					if (applyZNuNuGKfac) {
						kFac 														= 	kFacQCDZNuNuG.getKFac(lhePhoPt) * kFacEWZNuNuG.getKFac(lhePhoPt);
					} else if (applyWLNuGKfac) {
						kFac 														= 	(Wsign > 0) ? kFacQCDWpLNuG.getKFac(lhePhoPt) * kFacEWWpLNuG.getKFac(lhePhoPt) : kFacQCDWmLNuG.getKFac(lhePhoPt) * kFacEWWmLNuG.getKFac(lhePhoPt);
					}

					if (applyGJetsKfacToSample) 				kFac 				*=	GJetsKfac;

					Float_t	 									weight 				=  	lumi * puWeight * genWeight * xSection*1000. * kFac/sumGenWeight;

					if (overFlowInLastXbin && (theVariable > xLastBinLowEdge)) iSampleHist->Fill(xLastBinLowEdge, weight);
					else iSampleHist->Fill(theVariable, weight);
				}

				delete preSelection;
				closeTChain(sampleBinTree);
			}

			writeToFile(iSampleHist, writeFile);
		}
		allHists.push_back(std::make_pair(iSampleHist->Integral(), iSampleHist));
	}

	TH1F* 														jetFakesHist 			=		getJetFakesHist();
	TH1F* 														eleFakesHist 			=		getEleFakesHist();
	allHists.push_back(std::make_pair(jetFakesHist->Integral(), jetFakesHist));
	allHists.push_back(std::make_pair(eleFakesHist->Integral(), eleFakesHist));

	std::sort(allHists.rbegin(), allHists.rend());


	for (std::vector<std::pair<Double_t,TH1F*>>::reverse_iterator iHist = allHists.rbegin(); iHist!= allHists.rend(); iHist++) {
		if (iHist->second->Integral() < 10.e-30) continue;
		iHist->second->Scale(1., "width");
		hStack.Add(iHist->second, "HIST");
	}

	TCanvas canvas("canvas", "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);

	TPad pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad0.SetFillStyle(4000);
	pad0.SetFillColor(0);
	pad0.SetFrameFillStyle(4000);
	pad0.SetGrid(1,1);


	canvas.Draw();
	canvas.cd();
	pad0.Draw();
	pad0.cd();

	hStack.Draw("HIST");

	hStack.GetXaxis()->SetTitle("");
	hStack.GetXaxis()->CenterTitle();
	hStack.GetXaxis()->SetTitleSize(0);
	hStack.GetXaxis()->SetLabelSize(0);
	hStack.GetXaxis()->SetTitleOffset(0);
	hStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));
	hStack.GetXaxis()->SetNoExponent(1);
	hStack.GetYaxis()->SetTitle(options.getList("var", ";")[4].c_str());
	hStack.GetYaxis()->CenterTitle();
	hStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	hStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	hStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));
	hStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));

	TH1F* 														bgSum 					=		(TH1F*)	((TH1F*)hStack.GetStack()->Last())->Clone("totalBG");
	bgSum->SetLineWidth(0.);
	bgSum->SetFillStyle(std::stoi(options.getList("statUnc")[1]));
	bgSum->SetFillColorAlpha(hex2rootColor(options.getList("statUnc")[0]), std::stod(options.getList("statUnc")[2]));
	writeToFile(bgSum, options.get("saveDir") + "/" +"totalBG.root");

	Double_t 													sumBGIntErr;
	Double_t 													sumBGInt 				=		((TH1F*)bgSum)->IntegralAndError(1, hStack.GetXaxis()->GetNbins(), sumBGIntErr, "width");
	Double_t 													sigInt;
	Double_t 													sigIntErr;
	Double_t 													onlyBGint 				=		0.;
	Double_t 													onlyBGErr 				=		0.;
	Double_t  													smallestBin 			=		std::numeric_limits<Double_t>::max();

	for (std::pair<Double_t,TH1F*> & iHist : allHists) {
		Double_t 												iHistIntErr;
		Double_t 												iHistInt 				=		iHist.second->IntegralAndError(1, iHist.second->GetNbinsX(), iHistIntErr, "width");
		if (iHistInt  < 1e-30) continue;
		Float_t 												iFraction				=		100. * iHistInt/sumBGInt;
		Float_t 												iFractionErr			=		(iFraction == 0.) ? 0. : iFraction * std::sqrt(std::pow(iHistIntErr/iHistInt, 2) + std::pow(sumBGIntErr/sumBGInt, 2) - 2. * (iHistIntErr/iHistInt) * (sumBGIntErr/sumBGInt));

		std::string 											iLegStr 				=		std::string(iHist.second->GetTitle()) + " (" + to_string_with_precision(iFraction, options.getInt("fracPrec")) +
		        "#pm" + to_string_with_precision(iFractionErr, options.getInt("fracPrec")) + ")%";
		legend.AddEntry(iHist.second, iLegStr.c_str(), "F");

		yieldsFile<<iHist.second->GetTitle()<<";"<<removeTrailingZeros(iHistInt)<<";"<<removeTrailingZeros(iHistIntErr)<<";"<<removeTrailingZeros(iFraction)<<";"<<removeTrailingZeros(iFractionErr)<<"\n";
		std::cout<<iHist.second->GetTitle()<<";"<<removeTrailingZeros(iHistInt)<<";"<<removeTrailingZeros(iHistIntErr)<<";"<<removeTrailingZeros(iFraction)<<";"<<removeTrailingZeros(iFractionErr)<<"\n";

		if (smallestBin > iHist.second->GetMinimum(1e-20)) smallestBin = iHist.second->GetMinimum(1e-20);

		if (match(options.get("sig"), iHist.second->GetName())) {
			sigInt = iHistInt;
			sigIntErr = iHistIntErr;
		} else {
			onlyBGint += iHistInt;
			onlyBGErr += iHistIntErr*iHistIntErr;
		}
	}

	onlyBGErr = std::sqrt(onlyBGErr);


	bgSum->Draw("SAME E2");

	yieldsFile<<"TotalPrediction"<<";"<<removeTrailingZeros(sumBGInt)<<";"<<removeTrailingZeros(sumBGIntErr)<<";"<<100<<";"<<0<<"\n";
	std::cout<<"TotalPrediction"<<";"<<removeTrailingZeros(sumBGInt)<<";"<<removeTrailingZeros(sumBGIntErr)<<";"<<100<<";"<<0<<"\n";

	std::string 												statLeg 					=		options.getList("statUnc")[5] + " (" + to_string_with_precision(sumBGInt, options.getInt("fracPrec")) + "#pm"+ to_string_with_precision(sumBGIntErr, options.getInt("fracPrec")) + " events)";
	legend.AddEntry(bgSum, statLeg.c_str(), "F");


	//// Draw Data
	TH1F* 														dataHist 					=		nullptr;

	if (includeData) {
		dataHist 																			=		getDataHist();
		Double_t 													dataNerr;
		Double_t 													dataN 					=		dataHist->IntegralAndError(1, dataHist->GetNbinsX(), dataNerr);
		std::string 												dataLeg 				=		options.getList("dataAtts")[4] + " (" + std::to_string(Int_t(dataHist->Integral())) + "#pm" + to_string_with_precision(dataNerr,options.getInt("fracPrec")) +" events)";
		legend.AddEntry(dataHist, dataLeg.c_str(), "PE");
		dataHist->Scale(1., "width");
		dataHist->Draw("SAME P E");

		if (smallestBin > dataHist->GetMinimum(1e-20)) smallestBin = dataHist->GetMinimum(1e-20);

		yieldsFile<<dataHist->GetTitle()<<";"<<removeTrailingZeros(dataN)<<";"<<removeTrailingZeros(dataNerr)<<";"<<-1<<";"<<-1<<"\n";
		std::cout<<dataHist->GetTitle()<<";"<<removeTrailingZeros(dataN)<<";"<<removeTrailingZeros(dataNerr)<<";"<<-1<<";"<<-1<<"\n";
	}



	if (options.getInt("showSoSqrtB")) {
		Double_t 												sOverSqrtB 					=		sigInt/std::sqrt(onlyBGint);
		Double_t 												sOverSqrtBErr				=		sigIntErr*sigIntErr/onlyBGint + sigInt*sigInt*onlyBGErr*onlyBGErr/(4.*onlyBGint*onlyBGint*onlyBGint);
		sOverSqrtBErr																		=		std::sqrt(sOverSqrtBErr);
		std::string 											sOverSqrtBEntry 			=		"S/#sqrt{B} = " + ToSciString(sOverSqrtB, sOverSqrtBErr, options.getInt("sObprec"), "#pm", "^");
		legend.AddEntry(&legend, sOverSqrtBEntry.c_str(), "");
	}

	legend.Draw();


	Double_t yMin 				=	includeData ? std::min(dataHist->GetMinimum(1e-20), bgSum->GetMinimum(1e-20)) : bgSum->GetMinimum(1e-20);
	yMin 						=	std::min(smallestBin, yMin);
	Double_t yMax 				=	includeData ? std::max(bgSum->GetMaximum(), dataHist->GetMaximum()) : bgSum->GetMaximum();
	Float_t linUnit  			=	(yMax - yMin)/(1. - options.getDouble("yMinPadding") - options.getDouble("yMaxPadding"));

	if (yMin - options.getDouble("yMinPadding")*linUnit < 0) {
		yMin 					= 	0.;
		linUnit 				=	yMax/(1. - options.getDouble("yMaxPadding"));
	} else {
		yMin 					=	yMin - options.getDouble("yMinPadding")*linUnit;
	}
	yMax 						=	yMax + options.getDouble("yMaxPadding")*linUnit;
	hStack.SetMinimum(yMin);
	hStack.SetMaximum(yMax);

	cmsLegend(&pad0, to_string_with_precision(options.getFloat("lumi"), 2) + " fb^{-1} (13 TeV)", "Preliminary 2017", options.getFloat("titleScale"));


	//// Ratio plot
	TPad pad1("pad1", "", options.getDouble("pad1x1"), options.getDouble("pad1y1"), options.getDouble("pad1x2"), options.getDouble("pad1y2"));
	pad1.SetMargin(options.getDouble("pad1marginL"), options.getDouble("pad1marginR"), options.getDouble("pad1marginB"), options.getDouble("pad1marginT"));
	pad1.SetFillStyle(4000);
	pad1.SetFillColor(0);
	pad1.SetFrameFillStyle(4000);
	pad1.SetGrid(1,1);
	canvas.cd();
	pad1.Draw();
	pad1.cd();


	TH1F* 														ratioHist 				=		nullptr;
	if (includeData) {
		ratioHist 																		=		(TH1F*) dataHist->Clone("ratioHist");
		ratioHist->Reset();
		ratioHist->Divide(dataHist, bgSum, 1., 1., "B");
	} else {
		ratioHist 																		=		(TH1F*) bgSum->Clone("ratioHist");
		ratioHist->Reset();
	}

	ratioHist->Draw("PE");
	ratioHist->SetTitle("");
	ratioHist->SetMarkerStyle(std::stoi(options.getList("ratioAtts")[0]));
	ratioHist->SetMarkerSize(std::stod(options.getList("ratioAtts")[1]));
	ratioHist->SetMarkerColor(hex2rootColor(options.getList("ratioAtts")[3]));
	ratioHist->SetLineColor(hex2rootColor(options.getList("ratioAtts")[3]));
	ratioHist->SetLineWidth(std::stod(options.getList("ratioAtts")[2]));
	ratioHist->GetXaxis()->SetTitle(options.getList("var", ";")[1].c_str());
	ratioHist->GetXaxis()->CenterTitle();
	ratioHist->GetXaxis()->SetTitleSize(options.getDouble("pad1axisTitleSize"));
	ratioHist->GetXaxis()->SetLabelSize(options.getDouble("pad1axisLabelSize"));
	ratioHist->GetXaxis()->SetTitleOffset(options.getDouble("pad1XtitleOffset"));
	ratioHist->GetXaxis()->SetNdivisions(options.getInt("pad1xNdivs"));
	ratioHist->GetXaxis()->SetNoExponent(1);
	ratioHist->GetYaxis()->SetTitle(options.getList("ratioAtts")[4].c_str());
	ratioHist->GetYaxis()->CenterTitle();
	ratioHist->GetYaxis()->SetTitleSize(options.getDouble("pad1axisTitleSize"));
	ratioHist->GetYaxis()->SetLabelSize(options.getDouble("pad1axisLabelSize"));
	ratioHist->GetYaxis()->SetTitleOffset(options.getDouble("pad1YtitleOffset"));
	ratioHist->GetYaxis()->SetNdivisions(options.getInt("pad1yNdivs"));
	ratioHist->GetYaxis()->SetRangeUser(options.getFloatList("ratioRange")[0], options.getFloatList("ratioRange")[1]);

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("saveDir")+plotName+".png").c_str());
	canvas.SaveAs((options.get("saveDir")+plotName+".pdf").c_str());

	pad0.cd();
	pad0.SetLogy();

	yMin 														=	includeData ? std::min(dataHist->GetMinimum(1e-20), bgSum->GetMinimum(1e-20)) : bgSum->GetMinimum(1e-20);
	yMin 														=	std::min(smallestBin, yMin);
	yMax 														=	includeData ? std::max(bgSum->GetMaximum(), dataHist->GetMaximum()) : bgSum->GetMaximum();
	Float_t logUnit 											= 	std::log10(yMax/yMin)/(1. - options.getDouble("yMinPadding") - options.getDouble("yMaxPadding"));
	yMin 														= 	yMin/std::pow(10., options.getDouble("yMinPadding")*logUnit);
	yMax 														=	yMax*std::pow(10., options.getDouble("yMaxPadding")*logUnit);

	hStack.SetMaximum(yMax);
	hStack.SetMinimum(yMin);

	if (std::log10(yMax/yMin) < options.getFloat("setMoreLogLabels")) hStack.GetYaxis()->SetMoreLogLabels();

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("saveDir")+ "/" +plotName+"_logY.png").c_str());
	canvas.SaveAs((options.get("saveDir")+ "/"+plotName+"_logY.pdf").c_str());
};


void init() {

	options.parseIt(optionsFile, "==",1);

	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	gStyle->SetHatchesLineWidth(std::stod(options.getList("statUnc")[3]));
	gStyle->SetHatchesSpacing(std::stod(options.getList("statUnc")[4]));
	TGaxis::SetExponentOffset(options.getFloat("yExpOffstX"), options.getFloat("yExpOffstY"),"y");
	TGaxis::SetMaxDigits(options.getInt("yMaxDigits"));

	mkdir(options.get("saveDir"));

	options.optMap["cutFlowMade"] = "0";

	for (std::string iVar : options.getList("vars")) {

		if (!options.keyExists(iVar)) {
			std::cout<<"Binning for var "<<iVar<<" not found! Skipping..."<<std::endl;
			continue;
		}

		options.optMap["var"] = iVar + ";" + options.get(iVar);

		yieldsFile.open(options.get("saveDir") + "/" + options.getList("var", ";")[0] + "_yields.txt");
		yieldsFile<<"Source;Yield;YieldUnc;Frac;FracUnc\n";

		if (options.getInt("useBaysBlocksBinning") && std::stoi(options.getList("var", ";")[3]) > 0) createBaysBlocksBinning();

		drawDataMC();

		yieldsFile.close();

		clearHeap();
	}

	clearHeap();
};

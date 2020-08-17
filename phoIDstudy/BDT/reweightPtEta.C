#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TMap.h"
#include "TObjString.h"
#include "TString.h"
#include <iostream>
#include <fstream>

const Double_t ECAL_EB_ETA_BINS[29] = {-1.4442, -1.3, -1.2, -1.1, -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4442};

void 											makePtEtaHists(TH1D &_PtHist, TH1D &_EtaHist, TH2D & _PtEtaHist, TH1D &_PtNoXsecHist, TH1D &_EtaNoXsecHist, TH2D & _PtEtaNoXsecHist, const std::vector<std::string> & _samplePaths, std::string _directory);
void 											extractPtEtaWeights();
void 											createTrees4BDT();
void 											create1Tree4BDT(std::string filePath, Int_t _sampleIndex, Bool_t _isSignal);
void 											initialize();
std::string 									optionsFile = "reweightPtEtaOptions.txt";
std::string 									weightsFile;
std::vector<Float_t> 							etaBinning;
std::vector<Float_t> 							ptBinning;
parseOptions                                	options;
std::vector<std::string>                    	signalSamples;
std::vector<std::string>						bgSamples;
Float_t 										pTmin;
Float_t 										pTmax;
Float_t 										absEtaMax;

void 											initialize(){
	extractPtEtaWeights();
	createTrees4BDT();
};


void extractPtEtaWeights(){
	
	gROOT->SetBatch();

	options.parseIt(optionsFile, "==");

	//// load eta bins
	std::vector<Float_t> 						absEtaBinning 								=	options.getFloatList("etaBinning");
	absEtaMax 																				=	*std::max_element(absEtaBinning.begin(), absEtaBinning.end());
	for(std::vector<Float_t>::reverse_iterator iEtaBin = absEtaBinning.rbegin(); iEtaBin != absEtaBinning.rend(); ++iEtaBin){
		etaBinning.push_back(-(*iEtaBin));
		etaBinning.push_back(*iEtaBin);
	}
	std::sort (etaBinning.begin(), etaBinning.end());
	etaBinning.erase(std::unique(etaBinning.begin(), etaBinning.end() ), etaBinning.end());

	//// load pT bins
	ptBinning 																				=	options.getFloatList("ptBinning");
	pTmin																					= 	ptBinning[0];
	pTmax																					= 	ptBinning.back();

	mkdir(options.get("outDir"));
	weightsFile																				=	options.get("outDir")+"/weights.root";
	TFile                                       outFile(weightsFile.c_str(), "RECREATE");
	outFile.cd();

	//// Make histograms for signal
	TH1D										signalPt("signalPt", "Signal;p_{T} (GeV)", ptBinning.size()-1, ptBinning.data());
	signalPt.Sumw2();

	TH1D										signalEta("signalEta", "Signal;#eta", etaBinning.size()-1, etaBinning.data());
	signalEta.Sumw2();

	TH2D										signalPtEta("signalPtEta", "Signal;p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	signalPtEta.Sumw2();

	TH1D										signalPtNoXsec("signalPtNoXsec", "Signal (No weights);p_{T} (GeV)", ptBinning.size()-1, ptBinning.data());
	signalPtNoXsec.Sumw2();

	TH1D										signalEtaNoXsec("signalEtaNoXsec", "Signal (No weights);#eta", etaBinning.size()-1, etaBinning.data());
	signalEtaNoXsec.Sumw2();

	TH2D 										signalPtEtaNoXsec("signalPtEtaNoXsec", "Signal (No weights);p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	signalPtEtaNoXsec.Sumw2();	

	signalSamples																			= 	options.getList("signalSamples", ","); std::sort(signalSamples.begin(), signalSamples.end());

	std::cout<<"Getting signal histograms..."<<std::endl;
	makePtEtaHists(signalPt, signalEta, signalPtEta, signalPtNoXsec, signalEtaNoXsec, signalPtEtaNoXsec,  signalSamples, options.get("signalInDir"));

	normalizeHist(signalPt, 1000000.);
	normalizeHist(signalEta, 1000000.);
	normalizeHist(signalPtEta, 1000000.);

	normalizeHist(signalPtNoXsec, 1000000.);
	normalizeHist(signalEtaNoXsec, 1000000.);
	normalizeHist(signalPtEtaNoXsec, 1000000.);

	//// Make histograms for background
	outFile.cd();
	TH1D										bgPt("bgPt", "Background;p_{T} (GeV)", ptBinning.size()-1, ptBinning.data());
	bgPt.Sumw2();

	TH1D										bgEta("bgEta", "Background;#eta", etaBinning.size()-1, etaBinning.data());
	bgEta.Sumw2();

	TH2D 										bgPtEta("bgPtEta", "Background;p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	bgPtEta.Sumw2();

	TH1D										bgPtNoXsec("bgPtNoXsec", "Background;p_{T} (GeV)", ptBinning.size()-1, ptBinning.data());
	bgPtNoXsec.Sumw2();

	TH1D										bgEtaNoXsec("bgEtaNoXsec", "Background (No weights);#eta", etaBinning.size()-1, etaBinning.data());
	bgEtaNoXsec.Sumw2();

	TH2D 										bgPtEtaNoXsec("bgPtEtaNoXsec", "Background (No weights);p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	bgPtEtaNoXsec.Sumw2();
	
	bgSamples																				= 	options.getList("bgSamples", ","); std::sort(signalSamples.begin(), signalSamples.end());
	
	std::cout<<"Getting background histograms..."<<std::endl;
	makePtEtaHists(bgPt, bgEta, bgPtEta, bgPtNoXsec, bgEtaNoXsec, bgPtEtaNoXsec, bgSamples, options.get("backgroundInDir"));

	normalizeHist(bgPt, 1000000.);
	normalizeHist(bgEta, 1000000.);
	normalizeHist(bgPtEta, 1000000.);

	normalizeHist(bgPtNoXsec, 1000000.);
	normalizeHist(bgEtaNoXsec, 1000000.);
	normalizeHist(bgPtEtaNoXsec, 1000000.);

	//// Make weights
	TH2D* 										weightsSoB 									=  (TH2D*) signalPtEta.Clone("weightsSoB");
	weightsSoB->SetDirectory(outFile.GetDirectory(0));
	weightsSoB->Divide(&bgPtEta);

	TH2D* 									weightsBoS 										=  (TH2D*) bgPtEta.Clone("weightsBoS");
	weightsBoS->SetDirectory(outFile.GetDirectory(0));
	weightsBoS->Divide(&signalPtEta);


	//// Make weights with no X-section
	TH2D* 										weightsSoBnoXsec							=  (TH2D*) signalPtEtaNoXsec.Clone("weightsSoBnoXsec");
	weightsSoBnoXsec->SetDirectory(outFile.GetDirectory(0));
	weightsSoBnoXsec->Divide(&bgPtEtaNoXsec);
	
	TH2D* 										weightsBoSnoXsec							=  (TH2D*) bgPtEtaNoXsec.Clone("weightsBoSnoXsec");
	weightsBoS->SetDirectory(outFile.GetDirectory(0));
	weightsBoS->Divide(&signalPtEtaNoXsec);

	outFile.Write();
	outFile.Close();

	std::cout<<"Distributions and weights written to "<<weightsFile<<std::endl;
	clearHeap();
};


void makePtEtaHists(TH1D &_PtHist, TH1D &_EtaHist, TH2D & _PtEtaHist, TH1D &_PtNoXsecHist, TH1D &_EtaNoXsecHist, TH2D & _PtEtaNoXsecHist, const std::vector<std::string> & _samplePaths, std::string _directory){
	
	std::cout<<"Filling Eta-Pt histograms @ "<<getCurrentTime()<<std::endl;

	for(const std::string & iSamplePath : _samplePaths){

		std::string                             sampleName          =	findAndReplaceAll(getFileName(iSamplePath), ".root", "");
		TChain*                                 signalTree          =	openTChain((std::vector<std::string>){_directory+"/"+iSamplePath}, options.get("inTreeName"));
		TTreeReader                             inputTTreeReader(signalTree);
		TTreeReaderAnyValue<Float_t>            puWeight(inputTTreeReader, "puWeight");
		TTreeReaderAnyValue<Float_t>            genWeight(inputTTreeReader, "genWeight");
		TTreeReaderAnyValue<Float_t>            phoPt(inputTTreeReader, "phoPt");
		TTreeReaderAnyValue<Float_t>            phoEta(inputTTreeReader, "phoEta");
		Double_t 								xSection 			=	std::stod(vLookup(sampleName, options.get("xSectionMap"), 0, 2));
		Double_t 								sumGenWeight		=	std::stod(vLookup(sampleName, options.get("xSectionMap"), 0, 7));

		std::cout<<"\t"<<iSamplePath <<"\tx-Section = "<<xSection<<"\tsumGenWeight = "<<sumGenWeight;

		UInt_t 									addedEntries 		=	0;
		while(inputTTreeReader.Next()){
			
			if(phoPt < pTmin) continue;
			if(phoPt > pTmax) continue;
			if(std::abs(phoEta) > absEtaMax) continue;

			Double_t weight =  puWeight * genWeight * xSection*1000000./sumGenWeight;
			
			_PtHist.Fill(phoPt, weight);
			_EtaHist.Fill(phoEta, weight);
			_PtEtaHist.Fill(phoPt, phoEta, weight);

			_PtNoXsecHist.Fill(phoPt);
			_EtaNoXsecHist.Fill(phoEta);
			_PtEtaNoXsecHist.Fill(phoPt, phoEta);

			addedEntries++;
		}

		closeTChain(signalTree);

		std::cout<<"\tN = "<<addedEntries<<std::endl;
	};

	std::cout<<"Histograms filled @ "<<getCurrentTime()<<std::endl;
};


void createTrees4BDT(){

	gROOT->SetBatch();

	std::ofstream sampleMapFile;
	sampleMapFile.open(options.get("outDir") + "/" + options.get("sampleMapFile"));

	std::cout<<"----------------------------------------------------------------------------------------------"<<std::endl;
	std::cout<<"Reweighting signal..."<<std::endl;
	for(UInt_t iBin = 0; iBin < signalSamples.size(); iBin++){
		std::string                             filePath            = 	options.get("signalInDir") + "/" + signalSamples[iBin] ;
		create1Tree4BDT(filePath, iBin, 1);

		sampleMapFile<<std::to_string(iBin)<<","<<filePath<<std::endl;
	}
	std::cout<<"**********************************************************************************************"<<std::endl;



	std::cout<<"----------------------------------------------------------------------------------------------"<<std::endl;
	std::cout<<"Reweighting background..."<<std::endl;
	for(UInt_t iBin = 0; iBin < bgSamples.size(); iBin++){
		Int_t 									sampleIndex         = 	signalSamples.size() + iBin;
		std::string                             filePath            = 	options.get("backgroundInDir") + "/" + bgSamples[iBin] ;
		create1Tree4BDT(filePath, iBin, 0);

		sampleMapFile<<std::to_string(sampleIndex)<<","<<filePath<<std::endl;
	}
	std::cout<<"**********************************************************************************************"<<std::endl;

	sampleMapFile.close();

	std::cout<<getCurrentTime()<<" Done!"<<std::endl;

	clearHeap();
};


void create1Tree4BDT(std::string filePath, Int_t _sampleIndex, Bool_t _isSignal){

	//// Import input tree
	TFile                                       sampleFile(filePath.c_str(), "READ");
	TTree*                                      sampleTree          = 	(TTree*) sampleFile.Get(options.get("inTreeName").c_str());
	ULong64_t                                   nEntries            = 	sampleTree->GetEntries();
	
	Float_t 									puWeight;
	Float_t 									genWeight;
	Float_t 									phoPt;
	Float_t 									phoEta;
	sampleTree->SetBranchStatus("*",1);
	sampleTree->SetBranchAddress("puWeight", &puWeight);
	sampleTree->SetBranchAddress("genWeight", &genWeight);
	sampleTree->SetBranchAddress("phoPt", &phoPt);
	sampleTree->SetBranchAddress("phoEta", &phoEta);

	std::string                             sampleName          =	findAndReplaceAll(getFileName(filePath), ".root", "");
	Double_t 								xSection 			=	std::stod(vLookup(sampleName, options.get("xSectionMap"), 0, 2));
	Double_t 								sumGenWeight		=	std::stod(vLookup(sampleName, options.get("xSectionMap"), 0, 7));
	
	std::cout<<"\t"<<filePath <<"\tx-Section = "<<xSection<<"\tsumGenWeight = "<<sumGenWeight<<std::endl;
	

	//Create output tree
	std::string 								outFilePath									=	options.get("outDir")+ "/" + std::to_string(_isSignal) + "_" + sampleName + ".root";
	TFile                                       outFile(outFilePath.c_str(), "RECREATE");
	outFile.cd();
	TTree*                                      outTree             		= (TTree*) sampleTree->CloneTree(0);
	Float_t 									xSecW;
	Float_t                                     PtEtaRwBG;
	Float_t                                     PtEtaRwSignal;
	Float_t                                     PtEtaNoXsecRwBG;
	Float_t                                     PtEtaNoXsecRwSignal;
	Bool_t                                      isSignal 					= 	_isSignal;
	UChar_t                                     sampleIndex 				= 	_sampleIndex;
	
	outTree->Branch("xSecW", &xSecW);
	outTree->Branch("PtEtaRwBG", &PtEtaRwBG);
	outTree->Branch("PtEtaRwSignal", &PtEtaRwSignal);
	outTree->Branch("PtEtaNoXsecRwBG", &PtEtaNoXsecRwBG);
	outTree->Branch("PtEtaNoXsecRwSignal", &PtEtaNoXsecRwSignal);
	outTree->Branch("isSignal", &isSignal);
	outTree->Branch("sampleIndex", &sampleIndex);
	outTree->SetDirectory(outFile.GetDirectory(0));

	//// Fill tree
	TH2D*                                       weightsBoS     				=   (TH2D*) getHistFromFile("weightsBoS", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       weightsSoB     				=   (TH2D*) getHistFromFile("weightsSoB", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       weightsNoXsecBoS    	 	=   (TH2D*) getHistFromFile("weightsBoSnoXsec", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       weightsNoXsecSoB     		=   (TH2D*) getHistFromFile("weightsSoBnoXsec", weightsFile, 0, outFile.GetDirectory(""));

	outFile.cd();
	TH2D                                        preRwPtEta("preRwPtEta", "Pre Re-weight;p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	TH1D                                        preRwPt("preRwPt", "Pre Re-weight;p_{T} (GeV)", ptBinning.size()-1, ptBinning.data());
	TH1D                                        preRwEta("preRwEta", "Pre Re-weight;#eta", etaBinning.size()-1, etaBinning.data());
	TH2D                                        preRwNoXsecPtEta("preRwNoXsecPtEta", "Pre Re-weight (No #times-section);p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	TH1D                                        preRwNoXsecPt("preRwNoXsecPt", "Pre Re-weight (No #times-section);p_{T} (GeV)", ptBinning.size()-1, ptBinning.data());
	TH1D                                        preRwNoXsecEta("preRwNoXsecEta", "Pre Re-weight (No #times-section);#eta", etaBinning.size()-1, etaBinning.data());

	TH2D                                        postRwPtEta("postRwPtEta", "Reweighted;p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	TH1D                                        postRwPt("postRwPt", "Reweighted;p_{T} (GeV)", ptBinning.size()-1, ptBinning.data());
	TH1D                                        postRwEta("postRwEta", "Reweighted;#eta", etaBinning.size()-1, etaBinning.data());
	TH2D                                        postRwNoXsecPtEta("postRwNoXsecPtEta", "Reweighted (No #times-section);p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	TH1D                                        postRwNoXsecPt("postRwNoXsecPt", "Reweighted (No #times-section);p_{T} (GeV)", ptBinning.size()-1, ptBinning.data());
	TH1D                                        postRwNoXsecEta("postRwNoXsecEta", "Reweighted (No #times-section);#eta", etaBinning.size()-1, etaBinning.data());

	for(ULong64_t iEntry = 0; iEntry < nEntries; iEntry++){

		sampleTree->GetEntry(iEntry);

		if(phoPt < pTmin) continue;
		if(phoPt > pTmax) continue;
		if(std::abs(phoEta) > absEtaMax) continue;
		
		xSecW                             		 							=	puWeight * genWeight * xSection * 1000000. / sumGenWeight;

		preRwPtEta.Fill(phoPt, phoEta, xSecW);
		preRwPt.Fill(phoPt, xSecW);
		preRwEta.Fill(phoEta, xSecW);
		preRwNoXsecPtEta.Fill(phoPt, phoEta);
		preRwNoXsecPt.Fill(phoPt);
		preRwNoXsecEta.Fill(phoEta);

		if(_isSignal){
			PtEtaRwBG 														=	1.;
			PtEtaRwSignal 													=	weightsBoS->GetBinContent(weightsBoS->FindBin(phoPt, phoEta));
			PtEtaNoXsecRwBG 												=	1.;
			PtEtaNoXsecRwSignal												=	weightsNoXsecBoS->GetBinContent(weightsNoXsecBoS->FindBin(phoPt, phoEta));

			postRwPtEta.Fill(phoPt, phoEta, xSecW*PtEtaRwSignal);
			postRwPt.Fill(phoPt, xSecW*PtEtaRwSignal);
			postRwEta.Fill(phoEta, xSecW*PtEtaRwSignal);
			postRwNoXsecPtEta.Fill(phoPt, phoEta, PtEtaNoXsecRwSignal);
			postRwNoXsecPt.Fill(phoPt, PtEtaNoXsecRwSignal);
			postRwNoXsecEta.Fill(phoEta, PtEtaNoXsecRwSignal);

		} else{
			PtEtaRwBG 														=	weightsSoB->GetBinContent(weightsSoB->FindBin(phoPt, phoEta));
			PtEtaRwSignal 													=	1.;
			PtEtaNoXsecRwBG 												=	weightsNoXsecSoB->GetBinContent(weightsNoXsecSoB->FindBin(phoPt, phoEta));
			PtEtaNoXsecRwSignal												=	1.;

			postRwPtEta.Fill(phoPt, phoEta, xSecW*PtEtaRwBG);
			postRwPt.Fill(phoPt, xSecW*PtEtaRwBG);
			postRwEta.Fill(phoEta, xSecW*PtEtaRwBG);
			postRwNoXsecPtEta.Fill(phoPt, phoEta, PtEtaNoXsecRwBG);
			postRwNoXsecPt.Fill(phoPt, PtEtaNoXsecRwBG);
			postRwNoXsecEta.Fill(phoEta, PtEtaNoXsecRwBG);
		}

		outTree->Fill();
	}

	sampleTree->Delete();
	sampleFile.Close();



	TH2D*										signalPtEta 				=   (TH2D*) getHistFromFile("signalPtEta", weightsFile, 0, outFile.GetDirectory(""));
	TH1D*										signalPt 					=   (TH1D*) getHistFromFile("signalPt", weightsFile, 0, outFile.GetDirectory(""));
	TH1D*										signalEta 					=   (TH1D*) getHistFromFile("signalEta", weightsFile, 0, outFile.GetDirectory(""));
	TH2D* 										signalPtEtaNoXsec 			=   (TH2D*) getHistFromFile("signalPtEtaNoXsec", weightsFile, 0, outFile.GetDirectory(""));
	TH1D*										signalPtNoXsec 				=   (TH1D*) getHistFromFile("signalPtNoXsec", weightsFile, 0, outFile.GetDirectory(""));
	TH1D*										signalEtaNoXsec 			=   (TH1D*) getHistFromFile("signalEtaNoXsec", weightsFile, 0, outFile.GetDirectory(""));

	TH2D*										bgPtEta 					=   (TH2D*) getHistFromFile("bgPtEta", weightsFile, 0, outFile.GetDirectory(""));
	TH1D*										bgPt 						=   (TH1D*) getHistFromFile("bgPt", weightsFile, 0, outFile.GetDirectory(""));
	TH1D*										bgEta 						=   (TH1D*) getHistFromFile("bgEta", weightsFile, 0, outFile.GetDirectory(""));
	TH2D* 										bgPtEtaNoXsec 				=   (TH2D*) getHistFromFile("bgPtEtaNoXsec", weightsFile, 0, outFile.GetDirectory(""));
	TH1D*										bgPtNoXsec 					=   (TH1D*) getHistFromFile("bgPtNoXsec", weightsFile, 0, outFile.GetDirectory(""));
	TH1D*										bgEtaNoXsec 				=   (TH1D*) getHistFromFile("bgEtaNoXsec", weightsFile, 0, outFile.GetDirectory(""));

	TH2D*                                       postRwPtEtaOverTarget 		= 	(TH2D*) postRwPtEta.Clone("postRwPtEtaOverTarget");
	TH1D*                                       postRwPtOverTarget 			= 	(TH1D*) postRwPt.Clone("postRwPtOverTarget");
	TH1D*                                       postRwEtaOverTarget 		=	(TH1D*) postRwEta.Clone("postRwEtaOverTarget");
	TH2D*                                       postRwNoXsecPtEtaOverTarget =	(TH2D*) postRwNoXsecPtEta.Clone("postRwNoXsecPtEtaOverTarget");
	TH1D*                                       postRwNoXsecPtOverTarget 	= 	(TH1D*) postRwNoXsecPt.Clone("postRwNoXsecPtOverTarget");
	TH1D*                                       postRwNoXsecEtaOverTarget 	= 	(TH1D*) postRwNoXsecEta.Clone("postRwNoXsecEtaOverTarget");

	normalizeHist(postRwPtEtaOverTarget, 1000000.);
	normalizeHist(postRwPtOverTarget, 1000000.);
	normalizeHist(postRwEtaOverTarget, 1000000.);
	normalizeHist(postRwNoXsecPtEtaOverTarget, 1000000.);
	normalizeHist(postRwNoXsecPtOverTarget, 1000000.);
	normalizeHist(postRwNoXsecEtaOverTarget, 1000000.);

	if(_isSignal){
		postRwPtEtaOverTarget->Divide(bgPtEta);
		postRwPtOverTarget->Divide(bgPt);
		postRwEtaOverTarget->Divide(bgEta);
		postRwNoXsecPtEtaOverTarget->Divide(bgPtEtaNoXsec);
		postRwNoXsecPtOverTarget->Divide(bgPtNoXsec);
		postRwNoXsecEtaOverTarget->Divide(bgEtaNoXsec);
	} else{
		postRwPtEtaOverTarget->Divide(signalPtEta);
		postRwPtOverTarget->Divide(signalPt);
		postRwEtaOverTarget->Divide(signalEta);
		postRwNoXsecPtEtaOverTarget->Divide(signalPtEtaNoXsec);
		postRwNoXsecPtOverTarget->Divide(signalPtNoXsec);
		postRwNoXsecEtaOverTarget->Divide(signalEtaNoXsec);
	}

	outFile.cd();

	postRwPtEtaOverTarget->Write();
	postRwPtOverTarget->Write();
	postRwEtaOverTarget->Write();
	postRwNoXsecPtEtaOverTarget->Write();
	postRwNoXsecPtOverTarget->Write();
	postRwNoXsecEtaOverTarget->Write();

	outFile.Write();
	outFile.Close();

	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"||-- Processed "<< _sampleIndex<<"\t"<<filePath<<std::endl;
	std::cout<<"||-- File written to "<< outFilePath<<std::endl<<std::endl<<std::endl;

	clearHeap();

};
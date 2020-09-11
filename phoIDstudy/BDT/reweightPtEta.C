#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TMap.h"
#include "TObjString.h"
#include "TString.h"
#include <iostream>
#include <fstream>

void 											makePtEtaHists( TH2D & _PtEtaHist, TH2D & _PtEtaNoXsecHist, const std::vector<std::string> & _samplePaths, std::string _directory);
void 											extractPtEtaWeights();
void 											createTrees4BDT();
void 											create1Tree4BDT(std::string filePath, Int_t _sampleIndex, Bool_t _isSignal);
void 											initialize();
std::string 									optionsFile = "reweightPtEtaOptions.txt";
std::string 									weightsFile;
std::vector<Double_t> 							etaBinning;
std::vector<Double_t> 							ptBinning;
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
	std::vector<Double_t> 						absEtaBinning 								=	options.getDoubleList("etaBinning");
	absEtaMax 																				=	*std::max_element(absEtaBinning.begin(), absEtaBinning.end());
	for(std::vector<Double_t>::reverse_iterator iEtaBin = absEtaBinning.rbegin(); iEtaBin != absEtaBinning.rend(); ++iEtaBin){
		etaBinning.push_back(-(*iEtaBin));
		etaBinning.push_back(*iEtaBin);
	}
	std::sort (etaBinning.begin(), etaBinning.end());
	etaBinning.erase(std::unique(etaBinning.begin(), etaBinning.end() ), etaBinning.end());

	//// load pT bins
	ptBinning 																				=	options.getDoubleList("ptBinning");
	pTmin																					= 	ptBinning[0];
	pTmax																					= 	ptBinning.back();

	mkdir(options.get("outDir"));
	weightsFile																				=	options.get("outDir")+"/weights.root";
	TFile                                       outFile(weightsFile.c_str(), "RECREATE");
	outFile.cd();

	//// Make histograms for signal
	TH2D										signalPtEta("signalPtEta", "Signal;p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	signalPtEta.Sumw2();

	TH2D 										signalPtEtaNoXsec("signalPtEtaNoXsec", "Signal (No weights);p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	signalPtEtaNoXsec.Sumw2();	

	signalSamples																			= 	options.getList("signalSamples", ","); std::sort(signalSamples.begin(), signalSamples.end());

	std::cout<<"Getting signal histograms..."<<std::endl;
	makePtEtaHists(signalPtEta,signalPtEtaNoXsec,  signalSamples, options.get("signalInDir"));

	signalPtEta.Scale(1./signalPtEta.Integral());
	signalPtEtaNoXsec.Scale(1./signalPtEtaNoXsec.Integral());

	signalPtEta.Scale(1., "width");
	signalPtEtaNoXsec.Scale(1., "width");


	//// Make histograms for background
	outFile.cd();
	TH2D 										bgPtEta("bgPtEta", "Background;p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	bgPtEta.Sumw2();

	TH2D 										bgPtEtaNoXsec("bgPtEtaNoXsec", "Background (No weights);p_{T} (GeV);#eta", ptBinning.size()-1, ptBinning.data(), etaBinning.size()-1, etaBinning.data());
	bgPtEtaNoXsec.Sumw2();
	
	bgSamples																				= 	options.getList("bgSamples", ","); std::sort(signalSamples.begin(), signalSamples.end());
	
	std::cout<<"Getting background histograms..."<<std::endl;
	makePtEtaHists(bgPtEta, bgPtEtaNoXsec, bgSamples, options.get("backgroundInDir"));

	bgPtEta.Scale(1./bgPtEta.Integral());
	bgPtEtaNoXsec.Scale(1./bgPtEtaNoXsec.Integral());

	bgPtEta.Scale(1., "width");
	bgPtEtaNoXsec.Scale(1., "width");

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


void makePtEtaHists(TH2D & _PtEtaHist, TH2D & _PtEtaNoXsecHist, const std::vector<std::string> & _samplePaths, std::string _directory){
	
	std::cout<<"Filling Eta-Pt histograms @ "<<getCurrentTime()<<std::endl;

	for(const std::string & iSamplePath : _samplePaths){

		std::string                             sampleName          =	findAndReplaceAll(getFileName(iSamplePath), ".root", "");
		TChain*                                 signalTree          =	openTChain((std::vector<std::string>){_directory+"/"+iSamplePath}, options.get("inTreeName"));
		TTreeReader                             inputTTreeReader(signalTree);
		TTreeReaderAnyValue<Float_t>            puWeight(inputTTreeReader, "puWeight");
		TTreeReaderAnyValue<Float_t>            genWeight(inputTTreeReader, "genWeight");
		TTreeReaderAnyValue<Float_t>            phoPt(inputTTreeReader, "phoPt");
		TTreeReaderAnyValue<Float_t>            phoSCeta(inputTTreeReader, "phoSCeta");
		Double_t 								xSection 			=	std::stod(vLookup(sampleName, options.get("xSectionMap"), 0, 2));
		Double_t 								sumGenWeight		=	std::stod(vLookup(sampleName, options.get("xSectionMap"), 0, 7));

		std::cout<<"\t"<<iSamplePath <<"\tx-Section = "<<xSection<<"\tsumGenWeight = "<<sumGenWeight;

		UInt_t 									addedEntries 		=	0;
		while(inputTTreeReader.Next()){
			
			if(phoPt < pTmin) continue;
			if(phoPt > pTmax) continue;
			if(std::abs(phoSCeta) > absEtaMax) continue;

			Double_t weight =  puWeight * genWeight * xSection*100./sumGenWeight;
			
			_PtEtaHist.Fill(phoPt, phoSCeta, weight);
			_PtEtaNoXsecHist.Fill(phoPt, phoSCeta);

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

	std::vector<Double_t> 						absEtaBinning 								=	options.getDoubleList("etaBinning");
	absEtaMax 																				=	*std::max_element(absEtaBinning.begin(), absEtaBinning.end());

	ptBinning 																				=	options.getDoubleList("ptBinning");
	pTmin																					= 	ptBinning[0];
	pTmax																					= 	ptBinning.back();

	//// Import input tree
	TFile                                       sampleFile(filePath.c_str(), "READ");
	TTree*                                      sampleTree          = 	(TTree*) sampleFile.Get(options.get("inTreeName").c_str());
	ULong64_t                                   nEntries            = 	sampleTree->GetEntries();
	
	Float_t 									puWeight;
	Float_t 									genWeight;
	Float_t 									phoPt;
	Float_t 									phoSCeta;
	sampleTree->SetBranchStatus("*",1);
	sampleTree->SetBranchAddress("puWeight", &puWeight);
	sampleTree->SetBranchAddress("genWeight", &genWeight);
	sampleTree->SetBranchAddress("phoPt", &phoPt);
	sampleTree->SetBranchAddress("phoSCeta", &phoSCeta);

	std::string                             sampleName          =	findAndReplaceAll(getFileName(filePath), ".root", "");
	Double_t 								xSection 			=	std::stod(vLookup(sampleName, options.get("xSectionMap"), 0, 2));
	Double_t 								sumGenWeight		=	std::stod(vLookup(sampleName, options.get("xSectionMap"), 0, 7));
	
	std::cout<<"\t"<<filePath <<"\tx-Section = "<<xSection<<"\tsumGenWeight = "<<sumGenWeight<<std::endl;
	

	//Create output tree
	std::string 								outFilePath									=	options.get("outDir")+ "/" + std::to_string(_isSignal) + "_" + sampleName + ".root";
	TFile                                       outFile(outFilePath.c_str(), "RECREATE");
	outFile.cd();
	TTree*                                      outTree             		= (TTree*) sampleTree->CloneTree(0);
	Double_t 									xSecW;
	Double_t 									flatPtEtaRw;
	Double_t 									flatPtEtaRwNoXsec;
	Double_t                                     PtEtaRwBG;
	Double_t                                     PtEtaRwSignal;
	Double_t                                     PtEtaNoXsecRwBG;
	Double_t                                     PtEtaNoXsecRwSignal;
	Bool_t                                      isSignal 					= 	_isSignal;
	UChar_t                                     sampleIndex 				= 	_sampleIndex;
	
	outTree->Branch("xSecW", &xSecW);
	outTree->Branch("flatPtEtaRw", &flatPtEtaRw);
	outTree->Branch("flatPtEtaRwNoXsec", &flatPtEtaRwNoXsec);
	outTree->Branch("PtEtaRwBG", &PtEtaRwBG);
	outTree->Branch("PtEtaRwSignal", &PtEtaRwSignal);
	outTree->Branch("PtEtaNoXsecRwBG", &PtEtaNoXsecRwBG);
	outTree->Branch("PtEtaNoXsecRwSignal", &PtEtaNoXsecRwSignal);
	outTree->Branch("isSignal", &isSignal);
	outTree->Branch("sampleIndex", &sampleIndex);
	outTree->SetDirectory(outFile.GetDirectory(0));

	
	TH2D*                                       signalPtEta     			=   (TH2D*) getHistFromFile("signalPtEta", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       signalPtEtaNoXsec  			=   (TH2D*) getHistFromFile("signalPtEtaNoXsec", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       bgPtEta     				=   (TH2D*) getHistFromFile("bgPtEta", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       bgPtEtaNoXsec     			=   (TH2D*) getHistFromFile("bgPtEtaNoXsec", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       weightsBoS     				=   (TH2D*) getHistFromFile("weightsBoS", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       weightsSoB     				=   (TH2D*) getHistFromFile("weightsSoB", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       weightsNoXsecBoS    	 	=   (TH2D*) getHistFromFile("weightsBoSnoXsec", weightsFile, 0, outFile.GetDirectory(""));
	TH2D*                                       weightsNoXsecSoB     		=   (TH2D*) getHistFromFile("weightsSoBnoXsec", weightsFile, 0, outFile.GetDirectory(""));

	//// Fill tree
	outFile.cd();
	
	for(ULong64_t iEntry = 0; iEntry < nEntries; iEntry++){

		sampleTree->GetEntry(iEntry);

		if(phoPt < pTmin) continue;
		if(phoPt > pTmax) continue;
		if(std::abs(phoSCeta) > absEtaMax) continue;
		
		xSecW                             		 							=	puWeight * genWeight * xSection * 1000. / sumGenWeight;

		if(_isSignal){

			flatPtEtaRw 													=	1./signalPtEta->GetBinContent(signalPtEta->FindBin(phoPt, phoSCeta));
			flatPtEtaRwNoXsec 												=	1./signalPtEtaNoXsec->GetBinContent(signalPtEtaNoXsec->FindBin(phoPt, phoSCeta));
			PtEtaRwBG 														=	1.;
			PtEtaRwSignal 													=	weightsBoS->GetBinContent(weightsBoS->FindBin(phoPt, phoSCeta));
			PtEtaNoXsecRwBG 												=	1.;
			PtEtaNoXsecRwSignal												=	weightsNoXsecBoS->GetBinContent(weightsNoXsecBoS->FindBin(phoPt, phoSCeta));

		} else{
			flatPtEtaRw 													=	1./bgPtEta->GetBinContent(bgPtEta->FindBin(phoPt, phoSCeta));
			flatPtEtaRwNoXsec 												=	1./bgPtEtaNoXsec->GetBinContent(bgPtEtaNoXsec->FindBin(phoPt, phoSCeta));
			PtEtaRwBG 														=	weightsSoB->GetBinContent(weightsSoB->FindBin(phoPt, phoSCeta));
			PtEtaRwSignal 													=	1.;
			PtEtaNoXsecRwBG 												=	weightsNoXsecSoB->GetBinContent(weightsNoXsecSoB->FindBin(phoPt, phoSCeta));
			PtEtaNoXsecRwSignal												=	1.;
		}

		outTree->Fill();
	}

	sampleTree->Delete();
	sampleFile.Close();

	outFile.Write();
	outFile.Close();

	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"||-- Processed "<< _sampleIndex<<"\t"<<filePath<<std::endl;
	std::cout<<"||-- File written to "<< outFilePath<<std::endl<<std::endl<<std::endl;

	clearHeap();

};
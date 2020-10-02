#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include <stdio.h>

void optimizeCuts(){

	TChain*                                     featsTree      							=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/training/dataV2/mergedSamples.root"}), "fullEBTree");
	TChain*                                     bdtTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/training/dataV2/trainingV2/BDTresults_2020_09_12_19_08_19.root"}), "BDTresults");
	featsTree->AddFriend(bdtTree);
	TTreeReader                             	inputTTreeReader(featsTree);
	TTreeReaderAnyValue<long long>				isSignal_										(inputTTreeReader, "isSignalF");
	TTreeReaderAnyValue<long long>				isTest_											(inputTTreeReader, "isTestF");
	TTreeReaderAnyValue<long long>				isValidation_									(inputTTreeReader, "isValidationF");
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Float_t>				phoSigmaIEtaIEta_								(inputTTreeReader, "phoSigmaIEtaIEta");
	TTreeReaderAnyValue<Float_t>				phoHoverE_										(inputTTreeReader, "phoHoverE");
	TTreeReaderAnyValue<Float_t>				phoPFECALClusIsoCorr_							(inputTTreeReader, "phoPFECALClusIsoCorr");
	TTreeReaderAnyValue<Float_t>				phoPFHCALClusIsoCorr_							(inputTTreeReader, "phoPFHCALClusIsoCorr");
	TTreeReaderAnyValue<Float_t>				phoTkrIsoCorr_									(inputTTreeReader, "phoTkrIsoCorr");
	TTreeReaderAnyValue<Double_t>				flatPtEtaRwNoXsecF_								(inputTTreeReader, "flatPtEtaRwNoXsecF");

	TFile outFile("trainingV2/cutOptmizationV2_sietaieta.root", "RECREATE");
	TMVA::Factory factory("optimizeCuts", &outFile, "V:VerboseLevel=Verbose:!Silent:Color:Correlations:DrawProgressBar:AnalysisType=Classification");
	TMVA::DataLoader dataloader("trainingV2/dataset_sietaieta");
	dataloader.AddVariable("phoSigmaIEtaIEta", "#sigma_{i#etai#eta}", "", 'F');
	dataloader.AddVariable("phoHoverE", "H/E", "", 'F');
	dataloader.AddVariable("phoPFECALClusIsoCorr", "I_{ECAL}^{PFClus}", "GeV", 'F');
	dataloader.AddVariable("phoPFHCALClusIsoCorr", "I_{HCAL}^{PFClus}", "GeV", 'F');
	dataloader.AddVariable("phoTkrIsoCorr", "I_{TkrH}", "GeV", 'F');
	
	Double_t 									trainS0 								=		0;
	Double_t 									trainS1 								= 		0;
	Double_t 									testS0  								= 		0;
	Double_t 									testS1  								= 		0;

	Int_t 										currEntry 								= 		0;
	while(inputTTreeReader.Next()){

		if(std::abs(phoSCeta_) > 1.4442) continue;
		if(!isSignal_) continue;
		
		if(isValidation_) trainS0 += flatPtEtaRwNoXsecF_;
		else if(isTest_) testS0 += flatPtEtaRwNoXsecF_;
		
		std::vector<Double_t> vars({phoSigmaIEtaIEta_, phoHoverE_, phoPFECALClusIsoCorr_, phoPFHCALClusIsoCorr_, phoTkrIsoCorr_});

		if(isValidation_) dataloader.AddSignalTrainingEvent(vars, flatPtEtaRwNoXsecF_);
		else if(isTest_) dataloader.AddSignalTestEvent(vars, flatPtEtaRwNoXsecF_);

		if(isValidation_) trainS1 += flatPtEtaRwNoXsecF_;
		else if(isTest_) testS1 += flatPtEtaRwNoXsecF_;

		currEntry++;
	}

	std::cout<<"Signal:"<<std::endl<<"\tTrain S0 = "<<to_string_with_precision(trainS0, 10)<<"\tTrain S1 = "<<to_string_with_precision(trainS1, 10)<<std::endl<<
	"\t\tTrain reduction factor = "<<to_string_with_precision(trainS1/trainS0, 10)<<std::endl<<
	"\tTest S0 = "<<to_string_with_precision(testS0, 10)<<"\tTest S1 = "<<to_string_with_precision(testS1, 10)<<std::endl<<
	"\t\tTest reduction factor = "<<to_string_with_precision(testS1/testS0, 10)<<std::endl;

	trainS0 = 0;
	trainS1 = 0;
	testS0 = 0;
	testS1 = 0;

	inputTTreeReader.SetEntry(-1);
	currEntry = 0;
	while(inputTTreeReader.Next()){
		if(std::abs(phoSCeta_) > 1.4442) continue;

		if(isSignal_) continue;

		if(isValidation_) trainS0 += flatPtEtaRwNoXsecF_;
		else if(isTest_) testS0 += flatPtEtaRwNoXsecF_;

		std::vector<Double_t> vars({phoSigmaIEtaIEta_, phoHoverE_, phoPFECALClusIsoCorr_, phoPFHCALClusIsoCorr_, phoTkrIsoCorr_});
		
		if(isValidation_) dataloader.AddBackgroundTrainingEvent(vars, flatPtEtaRwNoXsecF_);
		else if(isTest_) dataloader.AddBackgroundTestEvent(vars, flatPtEtaRwNoXsecF_);

		if(isValidation_) trainS1 += flatPtEtaRwNoXsecF_;
		else if(isTest_) testS1 += flatPtEtaRwNoXsecF_;

		currEntry++;
	}

	std::cout<<"Background:"<<std::endl<<"\tTrain S0 = "<<to_string_with_precision(trainS0, 10)<<"\tTrain S1 = "<<to_string_with_precision(trainS1, 10)<<std::endl<<
	"\t\tTrain reduction factor = "<<to_string_with_precision(trainS1/trainS0, 10)<<std::endl<<
	"\tTest S0 = "<<to_string_with_precision(testS0, 10)<<"\tTest S1 = "<<to_string_with_precision(testS1, 10)<<std::endl<<
	"\t\tTest reduction factor = "<<to_string_with_precision(testS1/testS0, 10)<<std::endl;

	closeTChain(featsTree);
	closeTChain(bdtTree);

	dataloader.PrepareTrainingAndTestTree("(phoHoverE < 0.05) && (phoPFECALClusIsoCorr<12.) && (phoPFHCALClusIsoCorr < 20.) && (phoTkrIsoCorr < 10)",
		"(phoHoverE < 0.05) && (phoPFECALClusIsoCorr<12.) && (phoPFHCALClusIsoCorr < 20.) && (phoTkrIsoCorr < 10)",
		"NormMode=EqualNumEvents:ScaleWithPreselEff=True:V:VerboseLevel=Verbose:SplitSeed="+std::to_string(rand() % 10000000));
	factory.BookMethod(&dataloader, TMVA::Types::kCuts, "CutsGA", "H:V:VerbosityLevel=Verbose:VarProp[0]=FMin:VarProp[1]=FMin:VarProp[2]=FMin:VarProp[3]=FMin:VarProp[4]=FMin:FitMethod=GA:EffSel:PopSize=600:Steps=50:Cycles=4:SC_steps=10:SC_rate=5:SC_factor=0.95:Seed=0"); 

	factory.TrainAllMethods();
	factory.TestAllMethods();

	outFile.Write();
	outFile.Close();

	std::cout<<"Complete!"<<std::endl;
};
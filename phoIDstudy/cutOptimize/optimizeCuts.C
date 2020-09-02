#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include <stdio.h>

void optimizeCuts(){

	TChain*                                     featsTree      							=		openTChain(std::vector<std::string>({"/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root"}), "fullEB_Tree");
	TChain*                                     isoTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/isoShuffledTree.root"}), "fullEB_isoTree");
	TChain*                                     bdtTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/optimizedV1/BDTresults_0.root"}), "fullEB_BDT_Tree");
	featsTree->AddFriend(isoTree);
	featsTree->AddFriend(bdtTree);
	TTreeReader                             	inputTTreeReader(featsTree);
	TTreeReaderAnyValue<Bool_t>					isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Bool_t>					isTrain_										(inputTTreeReader, "isTrain");
	TTreeReaderAnyValue<Bool_t>					isValidation_									(inputTTreeReader, "isValidation");
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Double_t>				bdtScore_										(inputTTreeReader, "bdtScore");
	TTreeReaderAnyValue<Float_t>				phoHoverE_										(inputTTreeReader, "phoHoverE");
	TTreeReaderAnyValue<Float_t>				phoPFClusEcalIsoCorr_							(inputTTreeReader, "phoPFClusEcalIsoCorr");
	TTreeReaderAnyValue<Float_t>				phoPFClusHcalIsoCorr_							(inputTTreeReader, "phoPFClusHcalIsoCorr");
	TTreeReaderAnyValue<Float_t>				phoTrkSumPtHollowConeDR03Corr_					(inputTTreeReader, "phoTrkSumPtHollowConeDR03Corr");
	TTreeReaderAnyValue<Double_t>				bdtWeight_										(inputTTreeReader, "bdtWeightF");

	TFile outFile("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/cutOptimize/data/cutOptmizationV4.root", "RECREATE");
	TMVA::Factory factory("optimizeCutsEGMHoe", &outFile, "V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
	TMVA::DataLoader dataloader("dataset");
	dataloader.AddVariable("bdtScore", "BDT Score", "", 'F');
	dataloader.AddVariable("phoHoverE", "H/E", "", 'F');
	dataloader.AddVariable("phoPFClusEcalIsoCorr", "I_{ECAL}^{PFClus}", "GeV", 'F');
	dataloader.AddVariable("phoPFClusHcalIsoCor", "I_{HCAL}^{PFClus}", "GeV", 'F');
	dataloader.AddVariable("phoTrkSumPtHollowConeDR03Corr", "I_{TkrH}", "GeV", 'F');
	
	Double_t 									trainS0 								=		0;
	Double_t 									trainS1 								= 		0;
	Double_t 									testS0  								= 		0;
	Double_t 									testS1  								= 		0;

	Int_t 										currEntry 								= 		0;
	while(inputTTreeReader.Next()){

		if(std::abs(phoSCeta_) > 1.4442) continue;
		if(!isSignal_) continue;
		
		if(isValidation_) trainS0 += bdtWeight_;
		else if(!isTrain_) testS0 += bdtWeight_;
		
		//V2
		// if(phoHoverE_ > 0.1) continue;
		// if(phoPFClusEcalIsoCorr_ > 15.) continue;
		// if(phoPFClusHcalIsoCorr_ > 50.) continue;
		// if(phoTrkSumPtHollowConeDR03Corr_ > 10.) continue;

		//V3
		if(phoHoverE_ > 0.05) continue;
		if(phoPFClusEcalIsoCorr_ > 10.) continue;
		if(phoPFClusHcalIsoCorr_ > 20.) continue;
		if(phoTrkSumPtHollowConeDR03Corr_ > 8.) continue;
		
		std::vector<Double_t> vars({bdtScore_, phoHoverE_, phoPFClusEcalIsoCorr_, phoPFClusHcalIsoCorr_, phoTrkSumPtHollowConeDR03Corr_});

		if(isValidation_) dataloader.AddSignalTrainingEvent(vars, bdtWeight_);
		else if(!isTrain_) dataloader.AddSignalTestEvent(vars, bdtWeight_);

		if(isValidation_) trainS1 += bdtWeight_;
		else if(!isTrain_) testS1 += bdtWeight_;

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

		if(isValidation_) trainS0 += bdtWeight_;
		else if(!isTrain_) testS0 += bdtWeight_;

		if(phoHoverE_ > 0.05) continue;
		if(phoPFClusEcalIsoCorr_ > 10.) continue;
		if(phoPFClusHcalIsoCorr_ > 20.) continue;
		if(phoTrkSumPtHollowConeDR03Corr_ > 8.) continue;

		std::vector<Double_t> vars({bdtScore_, phoHoverE_, phoPFClusEcalIsoCorr_, phoPFClusHcalIsoCorr_, phoTrkSumPtHollowConeDR03Corr_});
		
		if(isValidation_) dataloader.AddBackgroundTrainingEvent(vars, bdtWeight_);
		else if(!isTrain_) dataloader.AddBackgroundTestEvent(vars, bdtWeight_);

		if(isValidation_) trainS1 += bdtWeight_;
		else if(!isTrain_) testS1 += bdtWeight_;

		currEntry++;
	}

	std::cout<<"Background:"<<std::endl<<"\tTrain S0 = "<<to_string_with_precision(trainS0, 10)<<"\tTrain S1 = "<<to_string_with_precision(trainS1, 10)<<std::endl<<
	"\t\tTrain reduction factor = "<<to_string_with_precision(trainS1/trainS0, 10)<<std::endl<<
	"\tTest S0 = "<<to_string_with_precision(testS0, 10)<<"\tTest S1 = "<<to_string_with_precision(testS1, 10)<<std::endl<<
	"\t\tTest reduction factor = "<<to_string_with_precision(testS1/testS0, 10)<<std::endl;

	closeTChain(featsTree);
	closeTChain(isoTree);
	closeTChain(bdtTree);

	dataloader.PrepareTrainingAndTestTree( "", "", "NormMode=EqualNumEvents:V" );
	factory.BookMethod(&dataloader, TMVA::Types::kCuts, "CutsGA", "H:!V:VarProp[0]=FMax:VarProp[1]=FMin:VarProp[2]=FMin:VarProp[3]=FMin:VarProp[4]=FMin:FitMethod=GA:EffSel:PopSize=1000:Steps=60:Cycles=3:SC_steps=10:SC_rate=5:SC_factor=0.95:SaveBestCycle=15");

	factory.TrainAllMethods();
	factory.TestAllMethods();

	outFile.Write();
	outFile.Close();

	std::cout<<"Complete!"<<std::endl;
};
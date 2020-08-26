#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"

void makeIdTree(){


	//// Import input tree
	TChain*                                     isoTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/isoShuffledTree.root"}), "fullEB_isoTree");
	TChain*                                     bdtTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/optimizedV1/BDTresults_0.root"}), "fullEB_BDT_Tree");
	TChain*                                     featsTree      							=		openTChain(std::vector<std::string>({"/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root"}), "fullEB_Tree");
	isoTree->AddFriend(bdtTree);
	isoTree->AddFriend(featsTree);
	TTreeReader                             	inputTTreeReader(isoTree);
	TTreeReaderAnyValue<Double_t>				bdtScore_										(inputTTreeReader, "bdtScore");
	TTreeReaderAnyValue<Float_t>				phoHoverE_										(inputTTreeReader, "phoHoverE");
	TTreeReaderAnyValue<Float_t>				phoPFClusEcalIsoCorr_							(inputTTreeReader, "phoPFClusEcalIsoCorr");
	TTreeReaderAnyValue<Float_t>				phoPFClusHcalIsoCorr_							(inputTTreeReader, "phoPFClusHcalIsoCorr");
	TTreeReaderAnyValue<Float_t>				phoTrkSumPtHollowConeDR03Corr_					(inputTTreeReader, "phoTrkSumPtHollowConeDR03Corr");

	//// Create output tree
	TFile                                       outFile("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/efficiency/data/idShuffledTree.root", "RECREATE");
	outFile.cd();
	TTree                                      outTree											("fullEB_IdTree", "Barrel Isolation Tree");
	Bool_t  									pass95;
	Bool_t 	 									pass90;
	Bool_t 	 									pass80;
	Bool_t 	 									pass70;
	
	outTree.Branch("pass95", &pass95);
	outTree.Branch("pass90", &pass90);
	outTree.Branch("pass80", &pass80);
	outTree.Branch("pass70", &pass70);
	outTree.SetDirectory(outFile.GetDirectory(0));


	Int_t 										currEntry = 0;
	while(inputTTreeReader.Next()){
		if(currEntry % 1000000 == 0) std::cout<<getCurrentTime()<<"\tReading entry\t"<<currEntry<<std::endl;

		pass95 = (bdtScore_ >= 8.49e-02) && (phoHoverE_ < 4.50e-02) && (phoPFClusEcalIsoCorr_ < 3.41) && (phoPFClusHcalIsoCorr_ < 9.28) && (phoTrkSumPtHollowConeDR03Corr_ < 4.21);
		pass90 = (bdtScore_ >= 2.05e-01) && (phoHoverE_ < 3.66e-02) && (phoPFClusEcalIsoCorr_ < 4.02) && (phoPFClusHcalIsoCorr_ < 5.21) && (phoTrkSumPtHollowConeDR03Corr_ < 3.33);
		pass80 = (bdtScore_ >= 3.29e-01) && (phoHoverE_ < 4.87e-02) && (phoPFClusEcalIsoCorr_ < 4.47) && (phoPFClusHcalIsoCorr_ < 3.03) && (phoTrkSumPtHollowConeDR03Corr_ < 1.43);
		pass70 = (bdtScore_ >= 6.47e-01) && (phoHoverE_ < 3.16e-02) && (phoPFClusEcalIsoCorr_ < 1.53) && (phoPFClusHcalIsoCorr_ < 9.26) && (phoTrkSumPtHollowConeDR03Corr_ < 1.75);

		outTree.Fill();

		currEntry++;
	}

	closeTChain(bdtTree);
	closeTChain(isoTree);
	closeTChain(featsTree);

	outFile.Write();
	outFile.Close();

	std::cout<<"Complete!"<<std::endl;

};
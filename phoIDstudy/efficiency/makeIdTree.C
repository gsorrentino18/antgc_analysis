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
	TTreeReaderAnyValue<Float_t>				phoSigmaIEtaIEta_								(inputTTreeReader, "phoSigmaIEtaIEta");	
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

	Bool_t  									cutbased95;
	Bool_t 	 									cutbased90;
	Bool_t 	 									cutbased80;
	Bool_t 	 									cutbased70;
	
	outTree.Branch("pass95", &pass95);
	outTree.Branch("pass90", &pass90);
	outTree.Branch("pass80", &pass80);
	outTree.Branch("pass70", &pass70);


	outTree.Branch("cutbased95", &cutbased95);
	outTree.Branch("cutbased90", &cutbased90);
	outTree.Branch("cutbased80", &cutbased80);
	outTree.Branch("cutbased70", &cutbased70);

	outTree.SetDirectory(outFile.GetDirectory(0));


	Int_t 										currEntry = 0;
	while(inputTTreeReader.Next()){
		if(currEntry % 1000000 == 0) std::cout<<getCurrentTime()<<"\tReading entry\t"<<currEntry<<std::endl;

		// pass95 = (bdtScore_ >= 8.49e-02) && (phoHoverE_ < 4.50e-02) && (phoPFClusEcalIsoCorr_ < 3.41) && (phoPFClusHcalIsoCorr_ < 9.28) && (phoTrkSumPtHollowConeDR03Corr_ < 4.21);
		// pass90 = (bdtScore_ >= 2.05e-01) && (phoHoverE_ < 3.66e-02) && (phoPFClusEcalIsoCorr_ < 4.02) && (phoPFClusHcalIsoCorr_ < 5.21) && (phoTrkSumPtHollowConeDR03Corr_ < 3.33);
		// pass80 = (bdtScore_ >= 3.29e-01) && (phoHoverE_ < 4.87e-02) && (phoPFClusEcalIsoCorr_ < 4.47) && (phoPFClusHcalIsoCorr_ < 3.03) && (phoTrkSumPtHollowConeDR03Corr_ < 1.43);
		// pass70 = (bdtScore_ >= 6.47e-01) && (phoHoverE_ < 3.16e-02) && (phoPFClusEcalIsoCorr_ < 1.53) && (phoPFClusHcalIsoCorr_ < 9.26) && (phoTrkSumPtHollowConeDR03Corr_ < 1.75);

		pass95 = (bdtScore_ >= 8.4858950660326768e-02) && (phoHoverE_ < 4.5018283831009365e-02) && (phoPFECALClusIsoCorr_ < 3.4103181645009730e+00) && (phoPFHCALClusIsoCorr_ < 9.2786995496122131e+00) && (phoTkrIsoCorr_ < 4.2117599454665431e+00);
		pass90 = (bdtScore_ >= 2.0507127927211427e-01) && (phoHoverE_ < 3.6568641101583359e-02) && (phoPFECALClusIsoCorr_ < 4.0163389219598731e+00) && (phoPFHCALClusIsoCorr_ < 5.2111989884892580e+00) && (phoTkrIsoCorr_ < 3.3294343175354939e+00);
		pass80 = (bdtScore_ >= 3.2887671669257112e-01) && (phoHoverE_ < 4.8717872871108611e-02) && (phoPFECALClusIsoCorr_ < 4.4686313648949589e+00) && (phoPFHCALClusIsoCorr_ < 3.0344907163887740e+00) && (phoTkrIsoCorr_ < 1.4285177877985955e+00);
		pass70 = (bdtScore_ >= 6.4674869718575889e-01) && (phoHoverE_ < 3.1618057666772573e-02) && (phoPFECALClusIsoCorr_ < 1.5315069640623058e-01) && (phoPFHCALClusIsoCorr_ < 9.2591138752480902e+00) && (phoTkrIsoCorr_ < 1.7539121705839822e+00);


		cutbased95 = (phoSigmaIEtaIEta_ < 1.0367983332487710e-02) && (phoHoverE_ < 4.9562589785941503e-02) && (phoPFClusEcalIsoCorr_ < 2.5519119822574599e+00) && (phoPFClusHcalIsoCorr_ < 5.7111424647345075e+00) && (phoTrkSumPtHollowConeDR03Corr_ < 3.4514573046845207e+00);
		cutbased90 = (phoSigmaIEtaIEta_ <1.0367983332487710e-02) && (phoHoverE_ < 3.6599900598484954e-02) && (phoPFClusEcalIsoCorr_ < 2.3816563604751853e+00) && (phoPFClusHcalIsoCorr_ < 4.3056564294352881e+00) && (phoTrkSumPtHollowConeDR03Corr_ < 1.9607086152426874e+00);
		cutbased80 = (phoSigmaIEtaIEta_ <1.0367983332487710e-02) && (phoHoverE_ < 1.8096978991278500e-02) && (phoPFClusEcalIsoCorr_ < -6.8794942900652245e-02) && (phoPFClusHcalIsoCorr_ < 3.2546532037354154e+00) && (phoTrkSumPtHollowConeDR03Corr_ < 2.0260429003153302e+00);
		cutbased70 = (phoSigmaIEtaIEta_ <1.0634545004917687e-02) && (phoHoverE_ < 3.4559647710403776e-02) && (phoPFClusEcalIsoCorr_ < 1.9707257306003800e+00) && (phoPFClusHcalIsoCorr_ < 4.3810222291263301e+00) && (phoTrkSumPtHollowConeDR03Corr_ < 2.2012193177264550e-01);

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
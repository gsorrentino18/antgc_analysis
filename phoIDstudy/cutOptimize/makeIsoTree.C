#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include <algorithm>
#include <random>

void makeIsoTree(){

	effectiveAreaMap							phoPFClusEcalIsoEA("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/90pc/phoPFClusEcalIso.txt", 1, ",", 0);
	isoPtScalingMap								phoPFClusEcalIsoPtScaling("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/isoPtScaling/90pc//phoPFClusEcalIso.txt", 0, 1, ",", 0);
	effectiveAreaMap							phoPFClusHcalIsoEA("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/90pc/phoPFClusHcalIso.txt", 1, ",", 0);
	isoPtScalingMap								phoPFClusHcalIsoPtScaling("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/isoPtScaling/90pc//phoPFClusHcalIso.txt", 1, 1, ",", 0);
	effectiveAreaMap							phoTrkSumPtHollowConeDR03EA("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/90pc/phoTrkSumPtHollowConeDR03.txt", 1, ",", 0);

	//// Import input tree
	TChain*                                     inTree         							=		openTChain(std::vector<std::string>({"/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root"}), "fullEB_Tree");
	TTreeReader                             	inputTTreeReader(inTree);
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Float_t>				phoPt_											(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Bool_t>					isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<UChar_t>				phoQualityBits_									(inputTTreeReader, "phoQualityBits");
	TTreeReaderAnyValue<Float_t>				rho_											(inputTTreeReader, "rho");
	TTreeReaderAnyValue<Float_t>				phoPFClusEcalIso_								(inputTTreeReader, "phoPFClusEcalIso");
	TTreeReaderAnyValue<Float_t>				phoPFClusHcalIso_								(inputTTreeReader, "phoPFClusHcalIso");
	TTreeReaderAnyValue<Float_t>				phoTrkSumPtHollowConeDR03_						(inputTTreeReader, "phoTrkSumPtHollowConeDR03");
	

	//// Create output tree
	TFile                                       outFile("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/isoShuffledTree.root", "RECREATE");
	outFile.cd();
	TTree*                                      outTree											("fullEB_isoTree", "Barrel Isolation Tree");
	Float_t 									phoPFClusEcalIsoCorr;
	Float_t 									phoPFClusHcalIsoCorr;
	Float_t 									phoTrkSumPtHollowConeDR03Corr;
	
	outTree.Branch("phoPFClusEcalIsoCorr", &phoPFClusEcalIsoCorr);
	outTree.Branch("phoPFClusHcalIsoCorr", &phoPFClusHcalIsoCorr);
	outTree.Branch("phoTrkSumPtHollowConeDR03Corr", &phoTrkSumPtHollowConeDR03Corr);
	outTree.SetDirectory(outFile.GetDirectory(0));


	Int_t 										currEntry = 0;
	while(inputTTreeReader.Next()){
		if(currEntry % 100000 == 0) std::cout<<getCurrentTime()<<"\tReading entry\t"<<currEntry<<std::endl;

		Float_t 								absSCEta								=		std::abs(phoSCeta_);
		Bool_t 									isEB 									=		(absSCEta <= 1.4442);
		phoPFClusEcalIsoCorr 															=		phoPFClusEcalIso_ - (isEB ? (rho_*phoPFClusEcalIsoEA.getEffectiveArea(absSCEta) +  phoPFClusEcalIsoPtScaling.getPtScaling(absSCEta, phoPt_)) : 0.);
		phoPFClusHcalIsoCorr															=		phoPFClusHcalIso_ - (isEB ? rho_*phoPFClusHcalIsoEA.getEffectiveArea(absSCEta) +  phoPFClusHcalIsoPtScaling.getPtScaling(absSCEta, phoPt_) : 0.);
		phoTrkSumPtHollowConeDR03Corr													=		phoTrkSumPtHollowConeDR03_ - (isEB ? rho_*phoTrkSumPtHollowConeDR03EA.getEffectiveArea(absSCEta) : 0.);

		outTree->Fill();

		currEntry++;
	}

	closeTChain(inTree);

	outFile.Write();
	outFile.Close();

	std::cout<<"Complete!"<<std::endl;

};
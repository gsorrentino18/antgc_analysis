#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include <algorithm>
#include <random>

void randomizeRows(){

	//// Import input tree
	TFile                                       inFile("/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamples.root", "READ");
	TTree*                                      inTree          = 	(TTree*) inFile.Get("fullEB_Tree");
	ULong64_t                                   nEntries            = 	inTree->GetEntries();
	Int_t iBaskets = inTree->LoadBaskets(31000000000);
	std::cout<<"iBaskets="<<iBaskets<<std::endl;

	//// Randomize entry list
	std::vector<UInt_t> 						randomizedEntries(nEntries);
	std::iota(randomizedEntries.begin(), randomizedEntries.end(), 0);
	std::default_random_engine generator;
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);

	//// Create output tree
	TFile                                       outFile("/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/mergedSamplesShuffled.root", "RECREATE");
	outFile.cd();
	TTree*                                      outTree             		= (TTree*) inTree->CloneTree(0);
	Float_t 									splitRand;
	Bool_t 										isTrain;
	Bool_t 										isValidation;
	outTree->Branch("splitRand", &splitRand);
	outTree->Branch("isTrain", &isTrain);
	outTree->Branch("isValidation", &isValidation);
	outTree->SetDirectory(outFile.GetDirectory(0));


	Int_t 										currEntry = 0;
	for(UInt_t iRandEntry : randomizedEntries){

		inTree->GetEntry(iRandEntry);

		splitRand = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		(splitRand < 0.7) ? isTrain = 1 : isTrain = 0;
		(splitRand < 0.29) ? isValidation = 1 : isValidation = 0;

		outTree->Fill();

		if(currEntry % 1000 == 0) std::cout<<getCurrentTime()<<"\tReading entry\t"<<currEntry<<"\t"<<iRandEntry <<std::endl;

		currEntry++;
	}

	inTree->Delete();
	inFile.Close();


	outFile.Write();
	outFile.Close();

};
// #include "extra_tools.cc"
// R__LOAD_LIBRARY(/panfs/roc/msisoft/mpfr/3.1.2-p11_gcc4.9.2/lib/libmpfr.so.4.1.2)
// #include </panfs/roc/msisoft/mpfr/3.1.2-p11_gcc4.9.2/include/mpfr.h>
// #include </panfs/roc/msisoft/mpfr/3.1.2-p11_gcc4.9.2/include/mpf2mpfr.h>

#include "TFile.h"
#include "TTree.h"
#include "iostream"
#include "sstream"
#include "fstream"
#include <algorithm>
#include <random>



void randomizeRows(){

	
	//// Import input tree
	TFile                                       inFile("dataV2/mergedSamples.root", "READ");
	TTree*                                      inTree          = 	(TTree*) inFile.Get("fullEBTree");
	ULong64_t                                   nEntries            = 	inTree->GetEntries();
	Int_t iBaskets = inTree->LoadBaskets(31000000000);
	std::cout<<"iBaskets="<<iBaskets<<std::endl;
	
	//// Randomize entry list
	std::vector<UInt_t> 						randomizedEntries(nEntries);
	std::iota(randomizedEntries.begin(), randomizedEntries.end(), 0);
	std::default_random_engine generator;

	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	generator.seed(time(0));
	std::shuffle(std::begin(randomizedEntries), std::end(randomizedEntries), generator);
	
	//// Create output tree
	TFile                                       outFile("dataV2/mergedSamplesShuffled.root", "RECREATE");
	outFile.cd();
	TTree*                                      outTree             		= (TTree*) inTree->CloneTree(0);
	Float_t 									splitRand;
	Bool_t 										isTrain;
	Bool_t 										isValidation;
	Bool_t 										isTest;
	
	outTree->Branch("splitRand", &splitRand);
	outTree->Branch("isTrain", &isTrain);
	outTree->Branch("isValidation", &isValidation);
	outTree->Branch("isTest", &isTest);
	outTree->SetDirectory(outFile.GetDirectory(0));
	

	Int_t 										currEntry = 0;
	for(UInt_t iRandEntry : randomizedEntries){

		inTree->GetEntry(iRandEntry);

		splitRand = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

		isTrain = 0;
		isValidation = 0;
		isTest = 0;

		if(splitRand < 0.4) isTrain = 1;
		else if (splitRand < 0.7) isValidation = 1;
		else isTest = 1;

		outTree->Fill();

		if(currEntry % 500000 == 0) std::cout<<"\tReading entry\t"<<currEntry<<"\t"<<iRandEntry <<std::endl;

		currEntry++;
	}
	
	inTree->Delete();
	inFile.Close();


	outFile.Write();
	outFile.Close();

};
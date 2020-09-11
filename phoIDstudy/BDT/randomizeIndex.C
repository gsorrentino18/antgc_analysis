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
	inTree->Delete();
	inFile.Close();

	//// Randomize entry list
	//// Create output tree
	TFile                                       outFile("dataV2/mergedSamplesShuffledIndices.root", "RECREATE");
	outFile.cd();
	TTree*                                      outTree             		= new TTree("randomized_split", "");
	Float_t 									splitRand;
	Bool_t 										isTrain;
	Bool_t 										isValidation;
	Bool_t 										isTest;
	
	outTree->Branch("splitRand", &splitRand);
	outTree->Branch("isTrain", &isTrain);
	outTree->Branch("isValidation", &isValidation);
	outTree->Branch("isTest", &isTest);
	outTree->SetDirectory(outFile.GetDirectory(0));
	
	srand(time(0));

	for(UInt_t iRandEntry = 0; iRandEntry < nEntries; iRandEntry++){

		splitRand = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

		isTrain = 0;
		isValidation = 0;
		isTest = 0;

		if(splitRand < 0.4) isTrain = 1;
		else if (splitRand < 0.7) isValidation = 1;
		else isTest = 1;

		outTree->Fill();

		if(iRandEntry % 500000 == 0) std::cout<<"\tEntry\t"<<iRandEntry <<std::endl;
	}

	outFile.Write();
	outFile.Close();

};
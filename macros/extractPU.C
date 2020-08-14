#include "extra_tools.cc"

std::string INFILES="/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5NtuplesSummary.txt";
std::string OUTDIR="/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/";
std::string HISTNAME="ggNtuplizer/hPUTruew";
std::string GENWHISTNAME="ggNtuplizer/hSumGenWeight";
std::string TREENAME="ggNtuplizer/EventTree";


void extractPU(std::string _inFileList, std::string _outFile, std::string _histName = HISTNAME, Bool_t _verbose=0);
void batchExtractPU(std::string _inLists=INFILES, std::string _outDir=OUTDIR, std::string _histName = HISTNAME, Bool_t _verbose=0);

void extractEffEntries(std::string _inFileList, std::string _treeName = TREENAME, Bool_t _verbose=0);
void batchExtractEffEntries(std::string _inLists=INFILES, std::string _treeName = TREENAME, Bool_t _verbose=0);


void extractPU(std::string _inFileList, std::string _outFile, std::string _histName, Bool_t _verbose){

	if(_verbose) std::cout<<"Extracting pileup histograms from ntuples listed in "<<_inFileList<<std::endl;

	std::vector<std::string> ntuples = getNonemptyLines(_inFileList);

	TH1F *puHist = (TH1F*) getHistFromFile(_histName, ntuples[0], 0);
	TH1F *genWeightHist = (TH1F*) getHistFromFile(GENWHISTNAME, ntuples[0], 0);

	if(!puHist) {
		std::cout<<"Error! Cannot read "<<_histName<<" from "<<ntuples[0]<<std::endl;
		return;
	}

	UInt_t nEntries = puHist->GetEntries();

	for(UInt_t i = 1; i < ntuples.size(); i++){
		TH1F *puHist_i = (TH1F*) getHistFromFile(_histName, ntuples[i], 0);
		TH1F *genWeightHist_i = (TH1F*) getHistFromFile(GENWHISTNAME, ntuples[i], 0);

		if(!puHist_i) {
			std::cout<<"Error! Cannot read "<<_histName<<" from "<<ntuples[i]<<std::endl;
			continue;
		}

		nEntries += puHist_i->GetEntries();
		puHist->Add(puHist_i);
		genWeightHist->Add(genWeightHist_i);
		delete puHist_i;
		delete genWeightHist_i;

		
	}

	writeToFile(puHist, _outFile, "RECREATE", _verbose);

	if(_verbose) std::cout<<" \t Entries = "<<puHist->GetEntries()<<std::endl;

	std::cout<<getFileName(_outFile)<<","<<nEntries<<","<<to_string_with_precision(genWeightHist->GetBinContent(1),20)<<std::endl;

	delete  puHist;
	delete genWeightHist;
};

void extractEffEntries(std::string _inFileList, std::string _treeName, Bool_t _verbose){
	TChain *tChain = openTChain(_inFileList, _treeName);
	TTreeReader                             tReader(tChain);
	TTreeReaderAnyValue<Float_t>			genWeight_										(tReader, "genWeight");

	Double_t 								sumW =0;
	Double_t 								sumW2 =0;
	ULong64_t								Nevts=0;
	Double_t 								wMax = - std::numeric_limits<Double_t>::max();
	Double_t 								wMin = std::numeric_limits<Double_t>::max();
	while(tReader.Next()){
		// Double_t errorSearch = (Double_t) std::abs(genWeight_) / (Double_t) std::numeric_limits<Float_t>::max();
		// if(errorSearch >0.9) std::cout<<genWeight_<<std::endl;
		// if(std::abs(genWeight_) < std::numeric_limits<Float_t>::min()*20) std::cout<<genWeight_<<std::endl;

		if(genWeight_ > wMax) wMax = genWeight_;
		if(genWeight_ < wMin) wMin = genWeight_;

		sumW+=genWeight_;
		sumW2+=genWeight_ * genWeight_;
		Nevts++;
	}

	closeTChain(tChain);

	std::cout<<getFileName(_inFileList)<<","<<Nevts<<","<<to_string_with_precision(sumW,20)<<","<<to_string_with_precision(sumW2,20)<<","<<to_string_with_precision(sumW*sumW/sumW2,20)<<std::endl;
	std::cout<<"Min="<<to_string_with_precision(wMin,10)<<"\tMax="<<to_string_with_precision(wMax,10)<<std::endl;
};

void batchExtractEffEntries(std::string _inLists, std::string _treeName, Bool_t _verbose){
	std::vector<std::string> ntupleLists = getLinesRegex(_inLists, "^((?!SinglePhoton).)*$");

	for(std::string it : ntupleLists){
		extractEffEntries(it, _treeName, _verbose);
	}
};


void batchExtractPU(std::string _inLists, std::string _outDir, std::string _histName, Bool_t _verbose){

	std::vector<std::string> ntupleLists = getLinesRegex(_inLists, "^((?!SinglePhoton).)*$");

	mkdir(_outDir);

	for(std::string it : ntupleLists){
		std::string outFile = _outDir + "/" + getFileName(it);
		outFile = findAndReplaceAll(outFile, ".txt", ".root");
		extractPU(it, outFile, _histName);
	}
};

void readEntries(std::string _inFileList, std::string _treeName=TREENAME, Bool_t _verbose=0){
	TChain *tChain = openTChain(_inFileList, _treeName);
	tChain->SetBranchStatus("*",0);
	// TTreeReader                             tReader(tChain);
	// // TTreeReaderAnyValue<Int_t>			run_										(tReader, "run");
	ULong64_t								Nevts=tChain->GetEntries();
	// while(tReader.Next()){
	// 	Nevts++;
	// }
	closeTChain(tChain);

	std::cout<<getFileName(_inFileList)<<","<<Nevts<<std::endl;
};

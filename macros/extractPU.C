#include "extra_tools.cc"

std::string INFILES="/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5NtuplesSummary.txt";
std::string OUTDIR="/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/";
std::string HISTNAME="ggNtuplizer/hPUTruew";


void extractPU(std::string _inFileList, std::string _outFile, std::string _histName = HISTNAME);
void batchExtractPU(std::string _inLists, std::string _outDir, std::string _histName = HISTNAME);



void extractPU(std::string _inFileList, std::string _outFile, std::string _histName){

	std::cout<<"Extracting pileup histograms from ntuples listed in "<<_inFileList<<std::endl;

	std::vector<std::string> ntuples = getNonemptyLines(_inFileList);

	TH1F *puHist = (TH1F*) getHistFromFile(_histName, ntuples[0], 0);

	if(!puHist) {
		std::cout<<"Error! Cannot read "<<_histName<<" from "<<ntuples[0]<<std::endl;
		return;
	}

	for(UInt_t i = 1; i < ntuples.size(); i++){
		TH1F *puHist_i = (TH1F*) getHistFromFile(_histName, ntuples[i], 0);

		if(!puHist_i) {
			std::cout<<"Error! Cannot read "<<_histName<<" from "<<ntuples[i]<<std::endl;
			continue;
		}

		puHist->Add(puHist_i);
		puHist_i->Delete();
	}

	writeToFile(puHist, _outFile);

	std::cout<<" \t Entries = "<<puHist->GetEntries()<<std::endl;

	std::cout<<getFileName(_outFile)<<","<<puHist->GetEntries()<<","<<to_string_with_precision(puHist->GetEffectiveEntries())<<std::endl;

	puHist->Delete();
};




void batchExtractPU(std::string _inLists=INFILES, std::string _outDir=OUTDIR, std::string _histName){

	std::vector<std::string> ntupleLists = getLinesRegex(_inLists, "^((?!SinglePhoton).)*$");

	mkdir(_outDir);

	for(std::string it : ntupleLists){
		std::string outFile = _outDir + "/" + getFileName(it);
		outFile = findAndReplaceAll(outFile, ".txt", ".root");
		extractPU(it, outFile, _histName);
	}
};

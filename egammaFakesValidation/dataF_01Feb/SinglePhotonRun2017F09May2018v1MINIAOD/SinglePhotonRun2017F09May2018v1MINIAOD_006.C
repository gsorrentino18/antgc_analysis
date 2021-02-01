#include "/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataF_01Feb//egammafakeValidation.cc"

void SinglePhotonRun2017F09May2018v1MINIAOD_006(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	egammafakeValidation("/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataF_01Feb//SinglePhotonRun2017F09May2018v1MINIAOD/SinglePhotonRun2017F09May2018v1MINIAOD_006", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataF_01Feb//SinglePhotonRun2017F09May2018v1MINIAOD//SinglePhotonRun2017F09May2018v1MINIAOD_006.root",
				-1.,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

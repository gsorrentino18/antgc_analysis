#include "/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataF_01FebV2//egammafakeValidationV2.cc"

void SinglePhotonRun2017F09May2018v1MINIAOD_017(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	egammafakeValidation("/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataF_01FebV2//SinglePhotonRun2017F09May2018v1MINIAOD/SinglePhotonRun2017F09May2018v1MINIAOD_017", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataF_01FebV2//SinglePhotonRun2017F09May2018v1MINIAOD//SinglePhotonRun2017F09May2018v1MINIAOD_017.root",
				-1.,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

#include "/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataB_01FebV2//egammafakeValidationV2.cc"

void SinglePhotonRun2017B31Mar2018v1MINIAOD_013(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	egammafakeValidation("/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataB_01FebV2//SinglePhotonRun2017B31Mar2018v1MINIAOD/SinglePhotonRun2017B31Mar2018v1MINIAOD_013", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/dataB_01FebV2//SinglePhotonRun2017B31Mar2018v1MINIAOD//SinglePhotonRun2017B31Mar2018v1MINIAOD_013.root",
				-1.,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

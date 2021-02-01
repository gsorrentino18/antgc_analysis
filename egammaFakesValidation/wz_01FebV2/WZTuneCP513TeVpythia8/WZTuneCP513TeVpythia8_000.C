#include "/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/wz_01FebV2//egammafakeValidationV2.cc"

void WZTuneCP513TeVpythia8_000(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	egammafakeValidation("/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/wz_01FebV2//WZTuneCP513TeVpythia8/WZTuneCP513TeVpythia8_000", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/wz_01FebV2//WZTuneCP513TeVpythia8//WZTuneCP513TeVpythia8_000.root",
				27.57,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/WZTuneCP513TeVpythia8.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

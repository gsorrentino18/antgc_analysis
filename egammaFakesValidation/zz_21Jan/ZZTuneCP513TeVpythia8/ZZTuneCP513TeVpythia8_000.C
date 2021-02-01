#include "/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/zz_21Jan//egammafakeValidation.cc"

void ZZTuneCP513TeVpythia8_000(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	egammafakeValidation("/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/zz_21Jan//ZZTuneCP513TeVpythia8/ZZTuneCP513TeVpythia8_000", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/zz_21Jan//ZZTuneCP513TeVpythia8//ZZTuneCP513TeVpythia8_000.root",
				12.14,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/ZZTuneCP513TeVpythia8.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

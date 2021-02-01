#include "/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/wele_28Jan//egammafakeValidation.cc"

void WToENuM100TuneCP513TeVpythia8_000(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	egammafakeValidation("/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/wele_28Jan//WToENuM100TuneCP513TeVpythia8/WToENuM100TuneCP513TeVpythia8_000", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/egammaFakesValidation/wele_28Jan//WToENuM100TuneCP513TeVpythia8//WToENuM100TuneCP513TeVpythia8_000.root",
				174.,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/WToENuM100TuneCP513TeVpythia8.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

#include "/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis_skim/phoID/2020_08_13_Prompt//genPhoMatcher.cc"

void WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8_000(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	genPhoMatcher("/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/batch/jobs/2020_08_13_Propmt//WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8/WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8_000", 
				"/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis_skim/phoID/2020_08_13_Prompt//WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8//WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8_000.root",
				.009687,
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/WGToLNuG01J5fPtG490TuneCP513TeVamcatnloFXFXpythia8.root",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/pileup_2017_data.root",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/EGM/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt,",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/EGM/effAreaPhotons_cone03_pfWorstChargedHadrons_70percentBased.txt",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/EGM/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/EGM/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

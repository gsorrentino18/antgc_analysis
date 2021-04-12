#include "/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis_skim/phoID/2020_08_07/promptDR0p01DPt0p1//GJetsDR0p4HT400To600TuneCP513TeVmadgraphMLMpythia8v2//genPhoMatcher.cc"

void GJetsDR0p4HT400To600TuneCP513TeVmadgraphMLMpythia8v2_013(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	genPhoMatcher("/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/batch/jobs/2020_08_07/promptDR0p01DPt0p1//GJetsDR0p4HT400To600TuneCP513TeVmadgraphMLMpythia8v2/GJetsDR0p4HT400To600TuneCP513TeVmadgraphMLMpythia8v2_013", 
				"/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis_skim/phoID/2020_08_07/promptDR0p01DPt0p1//GJetsDR0p4HT400To600TuneCP513TeVmadgraphMLMpythia8v2//GJetsDR0p4HT400To600TuneCP513TeVmadgraphMLMpythia8v2_013.root",
				124.6,
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/GJetsDR0p4HT400To600TuneCP513TeVmadgraphMLMpythia8v2.root",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/pileup_2017_data.root",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt,",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfWorstChargedHadrons_70percentBased.txt,",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt,",
				"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

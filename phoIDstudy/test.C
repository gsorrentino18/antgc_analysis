#include "genPhoMatcher.cc"

void test(){
	genPhoMatcher t(
		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.txt",
		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
		40.66,
		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/GJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/pileup_2017_data.root",
		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt",
		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfWorstChargedHadrons_70percentBased.txt",
		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt",
		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt");
}; 

// #include "fakePhoFinder.cc"

// void test(){
// 	fakePhoFinder t(
// 		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.txt",
// 		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
// 		40.66,
// 		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/GJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
// 		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/METv5_pileup/pileup_2017_data.root",
// 		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt",
// 		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfWorstChargedHadrons_70percentBased.txt",
// 		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt",
// 		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/effAreas/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt");
// }; 
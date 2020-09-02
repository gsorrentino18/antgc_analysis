// #include "fakeControlFinder.cc"

// void test(){
// 	fakeControlFinder t(
// 		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.txt",
// 		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
// 		40.66,
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/GJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfWorstChargedHadrons_70percentBased.txt",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt");
// }; 

// #include "fakePhoFinderV2.cc"

// void test(){
// 	fakePhoFinder t(
// 		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.txt",
// 		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
// 		40.66,
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/GJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfWorstChargedHadrons_70percentBased.txt",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt",
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt",
// 		"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/optimizedV0/aNTGC_photon_BDT.model");
// }; 

#include "fakePhoFinderV2.cc"

void test(){
	fakePhoFinder t(
		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.txt",
		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
		40.66,
		"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/GJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root");
}; 


// #include "genPhoMatcherV3.cc"

// void test(){
// 	genPhoMatcher t(
// 		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.txt",
// 		"testGJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root",
// 		40.66,
// 		"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/GJetsDR0p4HT600ToInfTuneCP513TeVmadgraphMLMpythia8.root");
// }; 
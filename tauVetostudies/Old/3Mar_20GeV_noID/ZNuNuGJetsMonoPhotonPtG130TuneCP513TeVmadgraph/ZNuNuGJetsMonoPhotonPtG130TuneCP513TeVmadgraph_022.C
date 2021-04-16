#include "/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_20GeV_noID//tauVetoV3.cc"

void ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_022(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	tauVeto("/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_20GeV_noID//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph/ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_022", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_20GeV_noID//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_022.root",
				.2038,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

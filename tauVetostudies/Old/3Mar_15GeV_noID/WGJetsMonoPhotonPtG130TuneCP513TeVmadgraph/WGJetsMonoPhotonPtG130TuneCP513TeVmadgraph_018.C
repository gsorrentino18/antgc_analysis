#include "/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_15GeV_noID//tauVetoV3.cc"

void WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_018(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	tauVeto("/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_15GeV_noID//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph/WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_018", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/3Mar_15GeV_noID//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_018.root",
				0.7158,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

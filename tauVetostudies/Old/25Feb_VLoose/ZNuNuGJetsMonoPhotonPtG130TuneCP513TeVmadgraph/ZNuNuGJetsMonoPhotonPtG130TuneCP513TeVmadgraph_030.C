#include "/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/25Feb_VLoose//tauVeto.cc"

void ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_030(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	tauVeto("/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/25Feb_VLoose//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph/ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_030", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/25Feb_VLoose//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph//ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_030.root",
				.2038,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

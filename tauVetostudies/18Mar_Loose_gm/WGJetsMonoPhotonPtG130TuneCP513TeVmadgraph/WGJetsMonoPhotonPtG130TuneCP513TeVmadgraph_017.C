#include "/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/18Mar_Loose_gm//tauVeto.cc"

void WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_017(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	tauVeto("/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/18Mar_Loose_gm//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph/WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_017", 
				"/data/cmszfs1/user/gsorrent/antgc_analysis/tauVetostudies/18Mar_Loose_gm//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph//WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph_017.root",
				0.7158,
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/WGJetsMonoPhotonPtG130TuneCP513TeVmadgraph.root",
				"/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

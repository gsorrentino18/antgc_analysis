#include "#ccfilepath"

void #macroname(){
	std::cout<<getCurrentTime()<<std::endl;
	std::cout<<"Begin root macro..."<<std::endl;
	
	#className("#fileList", 
				"#outfilepath",
				#xSec,
				"#mcPU",
				"#dataPU",
				"#CHEffArea",
				"#WCHEffArea",
				"#PhoEffArea",
				"#NHEffArea",
				"#BDTpath");

	std::cout<<"End root macro!"<<std::endl;
	std::cout<<getCurrentTime()<<std::endl;
};

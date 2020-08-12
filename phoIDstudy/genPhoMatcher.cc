//////////////////////////////////////////////////////////////////////
//  Mohammad Abrar Wadud, Univeristy of Minnesota           		//
//	August/04/2020                                            		//
//  Photon ID Study   												//
//////////////////////////////////////////////////////////////////////

#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"

R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/py2-xgboost/0.82/lib/python3.6/site-packages/xgboost/lib/libxgboost.so)

#include <xgboost/c_api.h>

#ifndef GENPHOMATCHER
#define GENPHOMATCHER

// Barrel-Endcap transition region
#define BETRetaMin 1.4442
#define BETRetaMax 1.566
#define HBetaMax 1.3920                     // ref. Josh H.
#define ZMASS 91.1876
#define pi7 3.1415927
#define REPORT_EVERY 50000
#define CUTFLOWSTEPS 30


const Double_t ECAL_ETA_BINS[57] = {-5., -3., -2.9, -2.7, -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -BETRetaMax, -BETRetaMin, -1.3, -1.2, -1.1, -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0., 0.1, 0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, BETRetaMin, BETRetaMax, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 2.9, 3., 5.};
const Double_t ECAL_ABS_ETA_BINS[13] = { 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, BETRetaMin, BETRetaMax, 1.8, 2.0, 2.2, 2.5};
const Double_t ECAL_EB_ETA_BINS[57] = {-BETRetaMin, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05, -1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, BETRetaMin};
const Double_t ECAL_FINE_ETA_BINS[85] = {-5., -3., -2.9, -2.7, -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -BETRetaMax, -BETRetaMin, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05, -1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, BETRetaMin, BETRetaMax, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 2.9, 3., 5.};
const Double_t ECAL_EB_ABS_ETA_BINS[8] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, BETRetaMin};


////////////////////////////////////////////////// Container for categories //////////////////////////////////////////////////////////////////////////////////////
struct eventType{
	TTree *                             tree = nullptr;

	// Cut efficiency tracking
	Float_t                             lastCutStep = 0.;
	TH1F *                              cutFlowCount = nullptr;
	TH1F *                              cutFlowGenWeight = nullptr;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class genPhoMatcher{
public:
	genPhoMatcher(std::string FILELIST, std::string OUTFILE, Float_t XSECTION=-1., std::string MCPILEUPHIST="", std::string DATAPILEUPHIST="",
		std::string CHARGED_HADRONIC_EFFECTIVE_AREAS="", std::string WORST_CHARGED_HADRONIC_EFFECTIVE_AREAS="", std::string PHOTONIC_EFFECTIVE_AREAS="", std::string NEUTRAL_HADRONIC_EFFECTIVE_AREAS="",
		std::string BDT_PATH="/local/cms/user/wadud/aNTGCmet/antgcpreselector/macros/photonIDmva/BDT/data/GJets_Central_NewPromptDef/resultsKDE/tuned/removeExtras/trained/aNTGC_photon_BDT.model");

	~genPhoMatcher(){
		XGBoosterFree(phoBDT_h);
		std::cout<<"END @ "<<getCurrentTime()<<std::endl;
		std::cout<<"*************************************************************************************************************************************************"<<std::endl;
	};

private:
	Bool_t              isMC = false;
	Bool_t              doPUreweight = false;
	Float_t             xSec = -1.;

	effectiveAreaMap    ChHadEffAreas;
	effectiveAreaMap    WChHadEffAreas;
	effectiveAreaMap    PhoEffAreas;
	effectiveAreaMap    NeuHadEffAreas;

	void                analyze();
	Bool_t              selectEvent();

	Short_t				photonIsPrompt(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	Bool_t				photonIsFake(Short_t _phoIndex, Float_t _deltaRmax = 0.3);
	Short_t matchWithRecoPho(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);

	TFile *             outFile = nullptr;

	/////////////////////////////////////////// Pileup Reweighting /////////////////////////////////////////////////////////
	PileupReWeighting   puReweighter;
	TH1F                pileupPreweight{"pileupUnweighted", "Unweighted Pileup; True # of Interactions", 200, 0., 200.};
	TH1F                pileupPostweight{"pileupWeighted", "Weighted Pileup; True # of Interactions", 200, 0., 200.};
	TH1F                rhoPreweight{"rhoUnweighted", "Unweighted #rho; #rho", 200, 0., 200.};
	TH1F                rhoPostweight{"rhoWeighted", "Weighted #rho; #rho", 200, 0., 200.};
	TH1F                nvtxPreweight{"nvtxUnweighted", "Unweighted # of Vertices; # of Vertices", 200, 0., 200.};
	TH1F                nvtxPostweight{"nvtxWeighted", "Weighted # of Vertices; # of Vertices", 200, 0., 200.};
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////// Input TTree ///////////////////////////////////////////////////////////
	Bool_t                                  initNtuples(std::string FILELIST);
	TTreeReader                             inputTTreeReader;
	TChain *                                inputTree = nullptr;

	TTreeReaderAnyValue<Int_t>              _run;
	TTreeReaderAnyValue<Long64_t>           _event;
	TTreeReaderAnyValue<UShort_t>           _lumis;
	TTreeReaderAnyValue<UChar_t>            _nVtx;
	TTreeReaderAnyValue<Float_t>            _rho;
	TTreeReaderAnyValue<ULong64_t>          _HLTPho;
	TTreeReaderAnyValue<UShort_t>           _beamHaloSummary;
	TTreeReaderAnyValue<UShort_t>           _metFilters;

	TTreeReaderAnyValue<UChar_t>            _puTrue;
	TTreeReaderAnyValue<Float_t>            _genWeight;
	TTreeReaderAnyValue<UShort_t>			_nMC;
	TTreeReaderVectorValue<Int_t>           _mcPID;
	TTreeReaderVectorValue<Float_t>			_mcPt;
	TTreeReaderVectorValue<Float_t>         _mcEta;
	TTreeReaderVectorValue<Float_t>         _mcPhi;
	TTreeReaderVectorValue<UShort_t>        _mcStatusFlag;
	TTreeReaderVectorValue<Short_t>         _mcStatus;

	TTreeReaderAnyValue<UShort_t>           _nPho;
	TTreeReaderVectorValue<Float_t>         _phoCalibEt;
	TTreeReaderVectorValue<Float_t>         _phoEt;
	TTreeReaderVectorValue<Float_t>         _phoEta;
	TTreeReaderVectorValue<Float_t>         _phoPhi;	
	TTreeReaderVectorValue<Float_t>         _phoSeedTime;
	TTreeReaderVectorValue<UChar_t>         _phoFiducialRegion;

	TTreeReaderVectorValue<UChar_t>         _phoQualityBits;
	TTreeReaderVectorValue<Float_t>         _phoR9Full5x5;
	TTreeReaderVectorValue<Float_t>         _phoSigmaIEtaIEtaFull5x5;
	TTreeReaderVectorValue<Float_t>         _phoSigmaIEtaIPhiFull5x5;
	TTreeReaderVectorValue<Float_t>         _phoSigmaIPhiIPhiFull5x5;	
	TTreeReaderVectorValue<Float_t>         _phoE2x2Full5x5;
	TTreeReaderVectorValue<Float_t>			_phoE5x5Full5x5;

	TTreeReaderVectorValue<Float_t>			_phoMaxEnergyXtal;
	TTreeReaderVectorValue<Float_t>			_phoE2ndFull5x5;
	TTreeReaderVectorValue<Float_t>			_phoE1x3Full5x5;
	TTreeReaderVectorValue<Float_t>			_phoE1x5Full5x5;
	TTreeReaderVectorValue<Float_t>			_phoE2x5Full5x5;

	TTreeReaderVectorValue<Float_t>			_phoPFClusEcalIso;
	TTreeReaderVectorValue<Float_t>			_phoPFClusHcalIso;
	TTreeReaderVectorValue<Float_t>			_phoTrkSumPtSolidConeDR04;
	TTreeReaderVectorValue<Float_t>			_phoTrkSumPtHollowConeDR04;
	TTreeReaderVectorValue<Float_t>			_phoTrkSumPtSolidConeDR03;
	TTreeReaderVectorValue<Float_t>			_phoTrkSumPtHollowConeDR03;
	TTreeReaderVectorValue<Float_t>			_phoECALIso;
	TTreeReaderVectorValue<Float_t>			_phoHCALIso;

	TTreeReaderVectorValue<Float_t>         _phoHoverE;
	TTreeReaderVectorValue<Float_t>         _phoPFChIso;
	TTreeReaderVectorValue<Float_t>         _phoPFPhoIso;
	TTreeReaderVectorValue<Float_t>         _phoPFNeuIso;
	TTreeReaderVectorValue<Float_t>         _phoPFChWorstIso;
	TTreeReaderVectorValue<Float_t>         _phoIDMVA;
	TTreeReaderVectorValue<UChar_t>         _phoIDbit;
	TTreeReaderVectorValue<Float_t>         _phoMIPTotEnergy;
	
	TTreeReaderVectorValue<Short_t>         _phoDirectEcalSCindex;
	TTreeReaderVectorValue<Float_t>         _ecalSCeta;
	TTreeReaderVectorValue<Float_t>         _ecalSCphi;
	TTreeReaderVectorValue<Float_t>         _ecalSCEn;
	TTreeReaderVectorValue<Float_t>         _ecalSCRawEn;
	TTreeReaderVectorValue<Float_t>         _ecalSCetaWidth;
	TTreeReaderVectorValue<Float_t>         _ecalSCphiWidth;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////// Output variables //////////////////////////////////////////////////////////
	Char_t              fillPhoVars(Short_t _phoIndex);

	Int_t               run_;
	Long64_t            event_;
	UShort_t            lumis_;	
	Float_t 			rho_;
	UChar_t				nVtx_;

	Float_t             genWeight_ = 1.;
	Float_t             puWeight_ = 1.;
	UChar_t				puTrue_;
	UShort_t			genStatusFlag_;
	Short_t				genStatus_;
	Float_t 			deltaRgenPho_;
	Float_t 			relDeltaPtGenPho_;
	Float_t 			deltaRPt_;
	Int_t				genPDGid_ = 0;

	Float_t 			phoPt_;
	Float_t 			phoEta_;
	Float_t 			phoPhi_;
	Float_t 			phoSeedTime_;
	
	UChar_t				phoQualityBits_;
	Float_t 			phoR9Full5x5_;
	Float_t 			phoS4Full5x5_;

	Float_t 			phoS1Full5x5_;
	Float_t 			phoS2ndFull5x5_;
	Float_t 			phoU21Full5x5_;
	Float_t 			phoS1x3Full5x5_;
	Float_t 			phoS2x5Full5x5_;
	Float_t 			phoS25Full5x5_;

	Float_t 			phoSigmaIEtaIEta_;
	Float_t 			phoSigmaIPhiIPhi_;
	Float_t 			phoSigmaIEtaIPhi_;
	Float_t 			phoE2x2Full5x5_;
	Float_t 			phoE5x5Full5x5_;

	Float_t				phoMaxEnergyXtal_;
	Float_t				phoE2ndFull5x5_;
	Float_t				phoE1x3Full5x5_;
	Float_t				phoE1x5Full5x5_;
	Float_t				phoE2x5Full5x5_;

	Float_t				phoPFClusEcalIso_;
	Float_t				phoPFClusHcalIso_;
	Float_t				phoTrkSumPtSolidConeDR04_;
	Float_t				phoTrkSumPtHollowConeDR04_;
	Float_t				phoTrkSumPtSolidConeDR03_;
	Float_t				phoTrkSumPtHollowConeDR03_;
	Float_t				phoECALIso_;
	Float_t				phoHCALIso_;

	Float_t 			phoHoverE_;
	Float_t 			phoPFChIsoRaw_;
	Float_t 			phoPFPhoIsoRaw_;
	Float_t 			phoPFNeuIsoRaw_;
	Float_t 			phoPFChWorstIsoRaw_;
	Float_t 			phoPFChIsoCorr_;
	Float_t 			phoPFChWorstIsoCorr_;
	Float_t 			phoPFPhoIsoCorr_;
	Float_t 			phoPFNeuIsoCorr_;	
	Float_t 			phoIDMVA_;
	UChar_t             phoIDbit_;
	Float_t 			phoBDTprediction_;
	Float_t 			phoMIP_;	

	Float_t 			phoRelTightPFChIsoCorr_;
	Float_t 			phoRelTightPFPhoIsoCorr_;
	Float_t 			phoRelTightPFNeuIsoCorr_;
	Float_t 			phoRelMedPFChIsoCorr_;
	Float_t 			phoRelMedPFPhoIsoCorr_;
	Float_t 			phoRelMedPFNeuIsoCorr_;
	Float_t 			phoRelLoosePFChIsoCorr_;
	Float_t 			phoRelLoosePFPhoIsoCorr_;
	Float_t 			phoRelLoosePFNeuIsoCorr_;

	Float_t 			phoSCet_;
	Float_t 			phoSCrawet_;
	Float_t 			phoSCeta_;
	Float_t 			phoSCphi_;
	Float_t 			phoSCEn_;
	Float_t 			phoSCRawEn_;
	Float_t 			phoSCetaWidth_;
	Float_t 			phoSCphiWidth_;	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////Buffer variables///////////////////////////////////////////////////////
	Short_t matchedGenPhoIndex;
	UChar_t nPromptPho_;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////// XGBoost ////////////////////////////////////////////////////////////
	BoosterHandle phoBDT_h;
	Bool_t predictBDT = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////Categories////////////////////////////////////////////////////////////
	Bool_t          initEventTypes();
	void            initEventType(eventType & evType, std::string typeName, std::string typeTitle);
	void            fillEventType(eventType & evType);
	void            registerCutFlow(eventType & evType);
	void            registerAllCutFlow();

	eventType 		fullEB;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
genPhoMatcher::genPhoMatcher(std::string FILELIST, std::string OUTFILE, Float_t XSECTION, std::string MCPILEUPHIST, std::string DATAPILEUPHIST,
	std::string CHARGED_HADRONIC_EFFECTIVE_AREAS, std::string WORST_CHARGED_HADRONIC_EFFECTIVE_AREAS, std::string PHOTONIC_EFFECTIVE_AREAS, std::string NEUTRAL_HADRONIC_EFFECTIVE_AREAS, std::string BDT_PATH){

	std::cout<<"*************************************************************************************************************************************************"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"Running genPhoMatcher"<<std::endl<<
	"\n\nInput parameters:"<<std::endl<<
	"\t\tFile list = "<<FILELIST<<std::endl<<
	"\t\tOutput file = "<<OUTFILE<<std::endl<<
	"\t\tCross section = "<<XSECTION<<std::endl<<
	"\t\tMC pileup histogram = "<<MCPILEUPHIST<<std::endl<<
	"\t\tData pileup histogram = "<<DATAPILEUPHIST<<std::endl<<
	"\t\tCharged Hadronic Effective Areas = "<<CHARGED_HADRONIC_EFFECTIVE_AREAS<<std::endl<<
	"\t\tPhotonic Effective Areas = "<<PHOTONIC_EFFECTIVE_AREAS<<std::endl<<
	"\t\tNeutral Hadronic Effective Areas = "<<NEUTRAL_HADRONIC_EFFECTIVE_AREAS<<std::endl;
	

	xSec = XSECTION;
	if(XSECTION > 0.) isMC = true;
	if(isMC && file_exists(MCPILEUPHIST) && file_exists(DATAPILEUPHIST)) doPUreweight = true;

	std::cout<<"\t\tSample is simulation = "<<std::boolalpha<<isMC<<std::endl<<
	"\t\tDo pileup reweight = "<<doPUreweight<<"\n\n"<<std::endl;

	// doPUreweight = 0;

	if(doPUreweight) {
		std::cout<<"Pileup reweighting:"<<std::endl;
		puReweighter.init(MCPILEUPHIST, DATAPILEUPHIST, "hPUTruew", "pileup");
	}

	if(file_exists(CHARGED_HADRONIC_EFFECTIVE_AREAS)) {
		std::cout<<"Charged hadronic effective area:"<<std::endl;
		ChHadEffAreas.init(CHARGED_HADRONIC_EFFECTIVE_AREAS);
	}

	if(file_exists(WORST_CHARGED_HADRONIC_EFFECTIVE_AREAS)) {
		std::cout<<"Worst Charged hadronic effective area:"<<std::endl;
		WChHadEffAreas.init(WORST_CHARGED_HADRONIC_EFFECTIVE_AREAS);
	}
	
	if(file_exists(PHOTONIC_EFFECTIVE_AREAS)){
		std::cout<<"Photonic effective area:"<<std::endl;
		PhoEffAreas.init(PHOTONIC_EFFECTIVE_AREAS);
	}
	if(file_exists(NEUTRAL_HADRONIC_EFFECTIVE_AREAS)){
		std::cout<<"Neutral hadronic effective area:"<<std::endl;
		NeuHadEffAreas.init(NEUTRAL_HADRONIC_EFFECTIVE_AREAS);
	}

	if(file_exists(BDT_PATH)){
		std::cout<<"Loading BDT model from "<<BDT_PATH <<std::endl;
		XGBoosterCreate(NULL, 0, &phoBDT_h);
		XGBoosterSetParam(phoBDT_h, "seed", "0");
		Int_t mLdSuccess = XGBoosterLoadModel(phoBDT_h, BDT_PATH.c_str());
		// XGBoosterSetParam(phoBDT_h, "nthread", "1");
		if(mLdSuccess == 0) predictBDT=1;
		else{
			std::cout<<"Failed to load BDT model!"<<std::endl;
		}
	}

	initNtuples(FILELIST);

	outFile = new TFile(OUTFILE.c_str(), "RECREATE");

	initEventTypes();

	analyze();

	outFile->Write();
	outFile->Close();

	closeTChain(inputTree);

	std::cout<<"\n\nOutput written to file\t"<<OUTFILE <<std::endl<<"Complete!"<<std::endl<<getCurrentTime()<<std::endl;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	genPhoMatcher::photonIsPrompt(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){

	Short_t matchedPromptGenPho = -999;
	Float_t minDeltaR = 999.;

	for(UShort_t iGenP=0; iGenP < _nMC; iGenP++){
		if(_mcPID[iGenP] != 22) continue;

		UShort_t iGenPStFl = _mcStatusFlag[iGenP];
		if(!getBit(iGenPStFl,1)) continue;

		Float_t dRiGenPho = deltaR(_phoEta[_phoIndex], _phoPhi[_phoIndex], _mcEta[iGenP], _mcPhi[iGenP]);
		if(dRiGenPho > _deltaRmax) continue;

		// Float_t relDeltaPtiGenPho = std::abs(_mcPt[iGenP] - _phoCalibEt[_phoIndex])/_mcPt[iGenP];
		// if(relDeltaPtiGenPho > _relDeltaPtMax) continue;
		// if(relDeltaPtiGenPho < _relDeltaPtMin) continue;

		if(dRiGenPho < minDeltaR){
			minDeltaR = dRiGenPho;
			matchedPromptGenPho = iGenP;
		}
	}

	return matchedPromptGenPho;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	genPhoMatcher::matchWithRecoPho(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){

	Short_t matchedRecoPho = -999;
	Float_t minDeltaRPt = 999.;

	for(UShort_t iPho = 0; iPho < _nPho; iPho++){

		UChar_t iPhoFidReg = _phoFiducialRegion[iPho];
		if(!getBit(iPhoFidReg, 0)) continue; 			// skip if not EB (0 = EB, 1 = EE, 2 = EB-EE gap)

		// if(_phoCalibEt[iPho] < 200.) continue;
		
		Float_t relDeltaPtiGenPho = (_phoCalibEt[iPho] - _mcPt[_genIndex])/_mcPt[_genIndex];
		if(relDeltaPtiGenPho > _relDeltaPtMax) continue;
		if(relDeltaPtiGenPho < _relDeltaPtMin) continue;

		Float_t dRiGenPho = deltaR(_phoEta[iPho], _phoPhi[iPho], _mcEta[_genIndex], _mcPhi[_genIndex]);
		if(dRiGenPho > _deltaRmax) continue;

		// Float_t iPhoGenDeltaRPt = std::sqrt(relDeltaPtiGenPho*relDeltaPtiGenPho + dRiGenPho*dRiGenPho);
		Float_t iPhoGenDeltaRPt = dRiGenPho;

		if(iPhoGenDeltaRPt < minDeltaRPt){
			matchedRecoPho = iPho;
			minDeltaRPt = iPhoGenDeltaRPt;
		}
	}

	return matchedRecoPho;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t genPhoMatcher::photonIsFake(Short_t _phoIndex, Float_t _deltaRmax){
	// Checks for the existence of a prompt gen photon within a given deltaR cone
	Bool_t photonIsFake = 1;
	Float_t photonEta = _phoEta[_phoIndex];
	Float_t photonPhi = _phoPhi[_phoIndex];

	for(UShort_t iGenP=0; iGenP<_nMC; iGenP++){
		if(_mcPID[iGenP] != 22) continue;

		UShort_t iGenPStFl = _mcStatusFlag[iGenP];
		if(!getBit(iGenPStFl,1)) continue;

		Float_t dRiGenPho = deltaR(photonEta, photonPhi, _mcEta[iGenP], _mcPhi[iGenP]);
		if(dRiGenPho > _deltaRmax) continue;

		photonIsFake = 0;

		break;
	}

	return photonIsFake;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t genPhoMatcher::selectEvent(){
	//// reset event cut flow
	fullEB.lastCutStep = 0.;
	registerAllCutFlow();

	/////////////////////////////////////////////////Global Event Cuts /////////////////////////////////////////////////////
	// ULong64_t HLTPho    = (_HLTPho);                //// 9 = HLT_Photon200_v, 11 = HLT_Photon300_NoHE_v, 19 = HLT_Photon135_PFMET100_v, no trigger on 2016 MC
	// if(!getBit(HLTPho, 9))      return 0;           //// 200 GeV photon trigger
	// registerAllCutFlow();

	///////////////////////////////////////////// Select highest pT photon /////////////////////////////////////////////////
	Short_t highetsPtGenIndex 	= -9999;
	Float_t highestPt 			= -999.;
	Bool_t	passedPrompt		= 0;
	Bool_t 	passedGenFidCut		= 0;
	Bool_t 	passedGenPtCut 		= 0;
	Bool_t 	passedRecoMatch		= 0;
	Short_t matchedRecoPhoton 	= -999;
	
	nPromptPho_ = 0;

	for(UShort_t iGenP=0; iGenP<_nMC; iGenP++){
		if(_mcPID[iGenP] != 22) continue;

		UShort_t iGenPStFl = _mcStatusFlag[iGenP];
		if(!getBit(iGenPStFl,1)) continue;
		passedPrompt = 1;

		nPromptPho_++;

		if(std::abs(_mcEta[iGenP]) > 1.4442) continue;
		passedGenFidCut = 1;

		if(std::abs(_mcPt[iGenP]) < 200.) continue;
		passedGenPtCut = 1;

		Short_t iGenRecoMatch = matchWithRecoPho(iGenP, 0.05, -0.1, 0.1);
		if(iGenRecoMatch < 0) continue;
		passedRecoMatch = 1;

		if(_mcPt[iGenP] > highestPt){
			highetsPtGenIndex = iGenP;
			highestPt = _mcPt[iGenP];
			matchedRecoPhoton = iGenRecoMatch;
		}
	}
	
	if(passedPrompt) registerAllCutFlow();
	if(passedGenFidCut) registerAllCutFlow();
	if(passedGenPtCut) registerAllCutFlow();
	if(passedRecoMatch) registerAllCutFlow();

	if(highetsPtGenIndex < 0) return 0;

	matchedGenPhoIndex = highetsPtGenIndex;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fillPhoVars(matchedRecoPhoton);
	fillEventType(fullEB);

	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t genPhoMatcher::initEventTypes(){

	pileupPreweight.SetDirectory(outFile->GetDirectory(""));
	pileupPostweight.SetDirectory(outFile->GetDirectory(""));
	rhoPreweight.SetDirectory(outFile->GetDirectory(""));
	rhoPostweight.SetDirectory(outFile->GetDirectory(""));
	nvtxPreweight.SetDirectory(outFile->GetDirectory(""));
	nvtxPostweight.SetDirectory(outFile->GetDirectory(""));

	initEventType(fullEB, "fullEB", "Full ECAL Barrel");
	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genPhoMatcher::analyze(){
	std::cout<<"--------------------------------------------------------------------------------------------------"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"Analyzing events.."<<std::endl;
	ULong64_t current_entry =0;
	while(inputTTreeReader.Next()){

		if(current_entry % REPORT_EVERY == 0){
			std::cout<<"\t"<< getCurrentTime()<<"\tAnalyzing entry\t"<<current_entry<<
			",\t\tevent\t"<<(_event)<<"\t\tFile " <<inputTree->GetCurrentFile()->GetName() <<std::endl;
		}

		// if(current_entry > 1000) break;
		rhoPreweight.Fill(_rho, genWeight_);
		nvtxPreweight.Fill(_nVtx, genWeight_);

		if(isMC) {  
			if(doPUreweight){
				genWeight_ = _genWeight;
				puWeight_ = puReweighter.weight(_puTrue);
				Float_t genPUweight = genWeight_ * puWeight_;
				pileupPostweight.Fill(_puTrue, genPUweight);
				rhoPostweight.Fill(_rho, genPUweight);
				nvtxPostweight.Fill(_nVtx, genPUweight);
			} else{
				genWeight_ = _genWeight;
				puWeight_ = genWeight_;    
			}
			pileupPreweight.Fill(_puTrue, genWeight_);
		}

		selectEvent();

		current_entry++; 
	};

	std::cout<<"Done analyzing!"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"--------------------------------------------------------------------------------------------------"<<std::endl;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Char_t genPhoMatcher::fillPhoVars(Short_t _phoIndex){

	run_ 					= _run;
	event_ 					= _event;
	lumis_ 					= _lumis;
	rho_ 					= _rho;
	nVtx_					= _nVtx;

	if(isMC){
		puTrue_					= _puTrue;

		Short_t iPhoGenIndex = matchedGenPhoIndex; //photonIsPrompt(_phoIndex, 0.3, 0.5);
		if(iPhoGenIndex > -1) {
			genStatusFlag_ 			= 	_mcStatusFlag[iPhoGenIndex];
			genStatus_ 				=	_mcStatus[iPhoGenIndex];
			deltaRgenPho_ 			=	deltaR(_mcEta[iPhoGenIndex], _mcPhi[iPhoGenIndex], _phoEta[_phoIndex], _phoPhi[_phoIndex]);
			relDeltaPtGenPho_ 		=	(_phoCalibEt[_phoIndex] - _mcPt[iPhoGenIndex])/_mcPt[iPhoGenIndex];
			deltaRPt_ 				= std::sqrt(deltaRgenPho_*deltaRgenPho_ + relDeltaPtGenPho_*relDeltaPtGenPho_);
			genPDGid_ 				= 	_mcPID[iPhoGenIndex];
		} else {
			genStatusFlag_ 			= 	9999;
			genStatus_ 				= 	-9999;
			deltaRgenPho_ 			=	-9999;
			relDeltaPtGenPho_ 		=	-9999;
			deltaRPt_ 				= 	-9999;
			genPDGid_				= -9999;
		}
	}

	Short_t phoSCindex 		= _phoDirectEcalSCindex[_phoIndex];
	Float_t phoAbsSCEta 	= std::abs(phoSCeta_);

	phoPt_ 					= _phoCalibEt[_phoIndex];
	phoEta_ 				= _phoEta[_phoIndex];
	phoPhi_ 				= _phoPhi[_phoIndex];
	phoSeedTime_ 			= _phoSeedTime[_phoIndex];

	phoQualityBits_			= _phoQualityBits[_phoIndex];
	phoR9Full5x5_ 			= _phoR9Full5x5[_phoIndex];
	phoS4Full5x5_			= _phoE2x2Full5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoS1Full5x5_ 			= _phoMaxEnergyXtal[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoS2ndFull5x5_ 		= _phoE2ndFull5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoU21Full5x5_ 			= _phoE2ndFull5x5[_phoIndex]/_phoMaxEnergyXtal[_phoIndex];
	phoS1x3Full5x5_ 		= _phoE1x3Full5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoS2x5Full5x5_ 		= _phoE2x5Full5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoS25Full5x5_ 			= _phoE5x5Full5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];

	phoSigmaIEtaIEta_ 		= _phoSigmaIEtaIEtaFull5x5[_phoIndex];
	phoSigmaIEtaIPhi_ 		= _phoSigmaIEtaIPhiFull5x5[_phoIndex];
	phoSigmaIPhiIPhi_ 		= _phoSigmaIPhiIPhiFull5x5[_phoIndex];	
	phoE2x2Full5x5_			= _phoE2x2Full5x5[_phoIndex];
	phoE5x5Full5x5_			= _phoE5x5Full5x5[_phoIndex];

	phoMaxEnergyXtal_		= _phoMaxEnergyXtal[_phoIndex];
	phoE2ndFull5x5_			= _phoE2x2Full5x5[_phoIndex];
	phoE1x3Full5x5_			= _phoE1x3Full5x5[_phoIndex];
	phoE1x5Full5x5_			= _phoE1x5Full5x5[_phoIndex];
	phoE2x5Full5x5_			= _phoE2x2Full5x5[_phoIndex];

	phoPFClusEcalIso_		= _phoPFClusEcalIso[_phoIndex];
	phoPFClusHcalIso_ 		= _phoPFClusHcalIso[_phoIndex];
	phoTrkSumPtSolidConeDR04_	= _phoTrkSumPtSolidConeDR04[_phoIndex];
	phoTrkSumPtHollowConeDR04_ 	= _phoTrkSumPtHollowConeDR04[_phoIndex];
	phoTrkSumPtSolidConeDR03_	= _phoTrkSumPtSolidConeDR03[_phoIndex];
	phoTrkSumPtHollowConeDR03_ 	= _phoTrkSumPtHollowConeDR03[_phoIndex];
	phoECALIso_ 			= _phoECALIso[_phoIndex];
	phoHCALIso_ 			= _phoHCALIso[_phoIndex];

	phoHoverE_ 				= _phoHoverE[_phoIndex];
	phoPFChIsoRaw_ 			= _phoPFChIso[_phoIndex];
	phoPFPhoIsoRaw_ 		= _phoPFPhoIso[_phoIndex];
	phoPFNeuIsoRaw_ 		= _phoPFNeuIso[_phoIndex];
	phoPFChWorstIsoRaw_ 	= _phoPFChWorstIso[_phoIndex];
	phoPFChIsoCorr_ 		= std::max(_phoPFChIso[_phoIndex] - _rho * ChHadEffAreas.getEffectiveArea(phoAbsSCEta), (Float_t)0.);
	phoPFChWorstIsoCorr_ 	= std::max(_phoPFChWorstIso[_phoIndex] - _rho * WChHadEffAreas.getEffectiveArea(phoAbsSCEta), (Float_t)0.);
	phoPFPhoIsoCorr_ 		= std::max(_phoPFPhoIso[_phoIndex] - _rho * PhoEffAreas.getEffectiveArea(phoAbsSCEta), (Float_t)0.);
	phoPFNeuIsoCorr_ 		= std::max(_phoPFNeuIso[_phoIndex] - _rho * NeuHadEffAreas.getEffectiveArea(phoAbsSCEta), (Float_t)0.);
	phoIDMVA_ 				= _phoIDMVA[_phoIndex];
	phoIDbit_ 				= _phoIDbit[_phoIndex];
	phoMIP_ 				= _phoMIPTotEnergy[_phoIndex];

	phoRelTightPFChIsoCorr_		= phoPFChIsoCorr_/0.65;
	phoRelTightPFPhoIsoCorr_	= phoPFPhoIsoCorr_/(2.044 + 0.004017*_phoEt[_phoIndex]);
	phoRelTightPFNeuIsoCorr_	= phoPFNeuIsoCorr_/(0.317 + 0.01512*_phoEt[_phoIndex] + 2.259e-05*_phoEt[_phoIndex]*_phoEt[_phoIndex]);
	phoRelMedPFChIsoCorr_ 		= phoPFChIsoCorr_/1.141;
	phoRelMedPFPhoIsoCorr_ 		= phoPFPhoIsoCorr_/(2.08 + 0.004017*_phoEt[_phoIndex]);
	phoRelMedPFNeuIsoCorr_ 		= phoPFNeuIsoCorr_/(1.189 + 0.01512*_phoEt[_phoIndex] + 2.259e-05*_phoEt[_phoIndex]*_phoEt[_phoIndex]);
	phoRelLoosePFChIsoCorr_ 	= phoPFChIsoCorr_/1.694;
	phoRelLoosePFPhoIsoCorr_ 	= phoPFPhoIsoCorr_/(2.876 + 0.004017*_phoEt[_phoIndex]);
	phoRelLoosePFNeuIsoCorr_ 	= phoPFNeuIsoCorr_/(24.032 + 0.01512*_phoEt[_phoIndex] + 2.259e-05*_phoEt[_phoIndex]*_phoEt[_phoIndex]);
	
	phoSCet_ 				= (_ecalSCEn[phoSCindex]) / std::cosh(phoSCeta_);
	phoSCrawet_				= (_ecalSCRawEn[phoSCindex]) / std::cosh(phoSCeta_);
	phoSCeta_ 				= _ecalSCeta[phoSCindex];
	phoSCphi_ 				= _ecalSCphi[phoSCindex];
	phoSCEn_ 				= _ecalSCEn[phoSCindex];
	phoSCRawEn_			= _ecalSCRawEn[phoSCindex];
	phoSCetaWidth_ 		= _ecalSCetaWidth[phoSCindex];
	phoSCphiWidth_ 		= _ecalSCphiWidth[phoSCindex];
	
	// BDT prediction
	// Features in order:	'ecalSCetaWidth', 'ecalSCphiWidth', 'phoHoverE', 'phoR9Full5x5', 'phoS4realFull5x5', 'phoSigmaIEtaIEta', 'phoSigmaIEtaIPhi', 'phoSigmaIPhiIPhi'
	if(predictBDT){
		float feats[1][8];
		feats[0][0] = phoSCetaWidth_;
		feats[0][1] = phoSCphiWidth_;
		feats[0][2] = phoHoverE_;
		feats[0][3] = phoR9Full5x5_;
		feats[0][4] = phoS4Full5x5_;
		feats[0][5] = phoSigmaIEtaIEta_;
		feats[0][6] = phoSigmaIEtaIPhi_;
		feats[0][7] = phoSigmaIPhiIPhi_;
		DMatrixHandle dTest;
		XGDMatrixCreateFromMat((float*)feats, 1, 8, -1, &dTest);
		bst_ulong out_len;
		const float *prediction;
		XGBoosterPredict(phoBDT_h, dTest, 0, 0, &out_len, &prediction);
		assert(out_len == 1);
		XGDMatrixFree(dTest);
		phoBDTprediction_ = prediction[0];
	}

	return 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t genPhoMatcher::initNtuples(std::string FILELIST){

	inputTree = openTChain(FILELIST, "ggNtuplizer/EventTree");
	inputTTreeReader.SetTree(inputTree);

	std::cout<<"**************************************************************************************************************************************************************"<<std::endl<<
	"Initializing branches in input ntuples..."<<std::endl;

	_run.set(inputTTreeReader, "run");
	_event.set(inputTTreeReader, "event");
	_lumis.set(inputTTreeReader, "lumis");
	_HLTPho.set(inputTTreeReader, "HLTPho");
	_beamHaloSummary.set(inputTTreeReader, "beamHaloSummary");
	_metFilters.set(inputTTreeReader, "metFilters");
	_nVtx.set(inputTTreeReader, "nVtx");
	_rho.set(inputTTreeReader, "rho");

	if(isMC){
		_puTrue.set(inputTTreeReader, "puTrue");
		_genWeight.set(inputTTreeReader, "genWeight");
		_nMC.set(inputTTreeReader, "nMC");
		_mcPID.set(inputTTreeReader, "mcPID");
		_mcPt.set(inputTTreeReader, "mcPt");
		_mcEta.set(inputTTreeReader, "mcEta");
		_mcPhi.set(inputTTreeReader, "mcPhi");
		_mcStatusFlag.set(inputTTreeReader, "mcStatusFlag");
		_mcStatus.set(inputTTreeReader, "mcStatus");
	}

	_nPho.set(inputTTreeReader, "nPho");
	_phoEt.set(inputTTreeReader, "phoEt");
	_phoCalibEt.set(inputTTreeReader, "phoCalibEt");
	_phoEta.set(inputTTreeReader, "phoEta");
	_phoPhi.set(inputTTreeReader, "phoPhi");
	_phoSeedTime.set(inputTTreeReader, "phoSeedTime");
	_phoFiducialRegion.set(inputTTreeReader, "phoFiducialRegion");

	_phoQualityBits.set(inputTTreeReader, "phoQualityBits");
	_phoR9Full5x5.set(inputTTreeReader, "phoR9Full5x5");
	_phoSigmaIEtaIEtaFull5x5.set(inputTTreeReader, "phoSigmaIEtaIEtaFull5x5");
	_phoSigmaIEtaIPhiFull5x5.set(inputTTreeReader, "phoSigmaIEtaIPhiFull5x5");
	_phoSigmaIPhiIPhiFull5x5.set(inputTTreeReader, "phoSigmaIPhiIPhiFull5x5");
	_phoE2x2Full5x5.set(inputTTreeReader, "phoE2x2Full5x5");
	_phoE5x5Full5x5.set(inputTTreeReader, "phoE5x5Full5x5");

	_phoMaxEnergyXtal.set(inputTTreeReader, "phoMaxEnergyXtal");
	_phoE2ndFull5x5.set(inputTTreeReader, "phoE2ndFull5x5");
	_phoE1x3Full5x5.set(inputTTreeReader, "phoE1x3Full5x5");
	_phoE1x5Full5x5.set(inputTTreeReader, "phoE1x5Full5x5");
	_phoE2x5Full5x5.set(inputTTreeReader, "phoE2x5Full5x5");

	_phoPFClusEcalIso.set(inputTTreeReader, "phoPFClusEcalIso");
	_phoPFClusHcalIso.set(inputTTreeReader, "phoPFClusHcalIso");
	_phoTrkSumPtSolidConeDR04.set(inputTTreeReader, "phoTrkSumPtSolidConeDR04");
	_phoTrkSumPtHollowConeDR04.set(inputTTreeReader, "phoTrkSumPtHollowConeDR04");
	_phoTrkSumPtSolidConeDR03.set(inputTTreeReader, "phoTrkSumPtSolidConeDR03");
	_phoTrkSumPtHollowConeDR03.set(inputTTreeReader, "phoTrkSumPtHollowConeDR03");
	_phoECALIso.set(inputTTreeReader, "phoECALIso");
	_phoHCALIso.set(inputTTreeReader, "phoHCALIso");
	
	_phoHoverE.set(inputTTreeReader, "phoHoverE");
	_phoPFChIso.set(inputTTreeReader, "phoPFChIso");
	_phoPFPhoIso.set(inputTTreeReader, "phoPFPhoIso");
	_phoPFNeuIso.set(inputTTreeReader, "phoPFNeuIso");
	_phoPFChWorstIso.set(inputTTreeReader, "phoPFChWorstIso");
	_phoIDMVA.set(inputTTreeReader, "phoIDMVA");
	_phoIDbit.set(inputTTreeReader, "phoIDbit");
	_phoMIPTotEnergy.set(inputTTreeReader, "phoMIPTotEnergy");
	
	_phoDirectEcalSCindex.set(inputTTreeReader, "phoDirectEcalSCindex");
	_ecalSCeta.set(inputTTreeReader, "ecalSCeta");
	_ecalSCphi.set(inputTTreeReader, "ecalSCphi");
	_ecalSCEn.set(inputTTreeReader, "ecalSCEn");
	_ecalSCRawEn.set(inputTTreeReader, "ecalSCRawEn");
	_ecalSCetaWidth.set(inputTTreeReader, "ecalSCetaWidth");
	_ecalSCphiWidth.set(inputTTreeReader, "ecalSCphiWidth");

	std::cout<<"Branches initialized!"<<std::endl<<
	"**************************************************************************************************************************************************************"<<std::endl;
	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genPhoMatcher::initEventType(eventType & evType, std::string typeName, std::string typeTitle){

	mkTFileDir(outFile, typeName);

	evType.cutFlowCount = new TH1F((typeName+"_cutFlowCount").c_str(), "Cut Flow (Unweighted)", CUTFLOWSTEPS, 0., (Float_t)CUTFLOWSTEPS);
	evType.cutFlowGenWeight = new TH1F((typeName+"_cutFlowGenWeight").c_str(), "Cut Flow (Gen Weighted)", CUTFLOWSTEPS, 0., (Float_t)CUTFLOWSTEPS);
	evType.cutFlowCount->SetDirectory(outFile->GetDirectory(typeName.c_str()));
	evType.cutFlowGenWeight->SetDirectory(outFile->GetDirectory(typeName.c_str()));

	//////// initialize tree
	evType.tree = new TTree((typeName+"_Tree").c_str(), typeTitle.c_str());
	evType.tree->SetDirectory(outFile->GetDirectory(typeName.c_str()));

	evType.tree->Branch("run", &run_);
	evType.tree->Branch("event", &event_);
	evType.tree->Branch("lumis", &lumis_);
	evType.tree->Branch("rho", &rho_);
	evType.tree->Branch("nVtx", &nVtx_);
	evType.tree->Branch("puTrue", &puTrue_);

	evType.tree->Branch("puWeight", &puWeight_);
	evType.tree->Branch("genWeight", &genWeight_);
	evType.tree->Branch("genStatusFlag", &genStatusFlag_);
	evType.tree->Branch("genStatus", &genStatus_);
	evType.tree->Branch("deltaRgenPho", &deltaRgenPho_);
	evType.tree->Branch("relDeltaPtGenPho", &relDeltaPtGenPho_);
	evType.tree->Branch("deltaRPt", &deltaRPt_);
	evType.tree->Branch("genPDGid", &genPDGid_);
	evType.tree->Branch("nPromptPho", &nPromptPho_);	

	evType.tree->Branch("phoPt", &phoPt_);
	evType.tree->Branch("phoEta", &phoEta_);
	evType.tree->Branch("phoPhi", &phoPhi_);
	evType.tree->Branch("phoSeedTime", &phoSeedTime_);

	evType.tree->Branch("phoQualityBits", &phoQualityBits_);
	evType.tree->Branch("phoR9Full5x5", &phoR9Full5x5_);
	evType.tree->Branch("phoS4Full5x5", &phoS4Full5x5_);

	evType.tree->Branch("phoS1Full5x5", &phoS1Full5x5_);
	evType.tree->Branch("phoS2ndFull5x5", &phoS2ndFull5x5_);
	evType.tree->Branch("phoU21Full5x5", &phoU21Full5x5_);	
	evType.tree->Branch("phoS1x3Full5x5", &phoS1x3Full5x5_);
	evType.tree->Branch("phoS2x5Full5x5", &phoS2x5Full5x5_);
	evType.tree->Branch("phoS25Full5x5", &phoS25Full5x5_);

	evType.tree->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta_);
	evType.tree->Branch("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi_);
	evType.tree->Branch("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi_);
	evType.tree->Branch("phoE2x2Full5x5", &phoE2x2Full5x5_);
	evType.tree->Branch("phoE5x5Full5x5", &phoE5x5Full5x5_);

	evType.tree->Branch("phoMaxEnergyXtal", &phoMaxEnergyXtal_);
	evType.tree->Branch("phoE2ndFull5x5", &phoE2ndFull5x5_);
	evType.tree->Branch("phoE1x3Full5x5", &phoE1x3Full5x5_);
	evType.tree->Branch("phoE1x5Full5x5", &phoE1x5Full5x5_);
	evType.tree->Branch("phoE2x5Full5x5", &phoE2x5Full5x5_);

	evType.tree->Branch("phoPFClusEcalIso", &phoPFClusEcalIso_);
	evType.tree->Branch("phoPFClusHcalIso", &phoPFClusHcalIso_);
	evType.tree->Branch("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04_);
	evType.tree->Branch("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04_);
	evType.tree->Branch("phoTrkSumPtSolidConeDR03", &phoTrkSumPtSolidConeDR03_);
	evType.tree->Branch("phoTrkSumPtHollowConeDR03", &phoTrkSumPtHollowConeDR03_);
	evType.tree->Branch("phoECALIso", &phoECALIso_);
	evType.tree->Branch("phoHCALIso", &phoHCALIso_);

	evType.tree->Branch("phoHoverE", &phoHoverE_);
	evType.tree->Branch("phoPFChIsoRaw", &phoPFChIsoRaw_);
	evType.tree->Branch("phoPFPhoIsoRaw", &phoPFPhoIsoRaw_);
	evType.tree->Branch("phoPFNeuIsoRaw", &phoPFNeuIsoRaw_);
	evType.tree->Branch("phoPFChWorstIsoRaw", &phoPFChWorstIsoRaw_);
	evType.tree->Branch("phoPFChIsoCorr", &phoPFChIsoCorr_);
	evType.tree->Branch("phoPFChWorstIsoCorr", &phoPFChWorstIsoCorr_);
	evType.tree->Branch("phoPFPhoIsoCorr", &phoPFPhoIsoCorr_);
	evType.tree->Branch("phoPFNeuIsoCorr", &phoPFNeuIsoCorr_);
	evType.tree->Branch("phoEGMidMVA", &phoIDMVA_);
	if(predictBDT) evType.tree->Branch("phoBDTprediction", &phoBDTprediction_);
	evType.tree->Branch("phoIDbit", &phoIDbit_);
	evType.tree->Branch("phoMIP", &phoMIP_);
	
	evType.tree->Branch("phoRelTightPFChIso", &phoRelTightPFChIsoCorr_);
	evType.tree->Branch("phoRelTightPFPhoIso", &phoRelTightPFPhoIsoCorr_);
	evType.tree->Branch("phoRelTightPFNeuIso", &phoRelTightPFNeuIsoCorr_);
	evType.tree->Branch("phoRelMedPFChIso", &phoRelMedPFChIsoCorr_);
	evType.tree->Branch("phoRelMedPFPhoIso", &phoRelMedPFPhoIsoCorr_);
	evType.tree->Branch("phoRelMedPFNeuIso", &phoRelMedPFNeuIsoCorr_);
	evType.tree->Branch("phoRelLoosePFChIso", &phoRelLoosePFChIsoCorr_);
	evType.tree->Branch("phoRelLoosePFPhoIso", &phoRelLoosePFPhoIsoCorr_);
	evType.tree->Branch("phoRelLoosePFNeuIso", &phoRelLoosePFNeuIsoCorr_);
	
	evType.tree->Branch("phoSCet", &phoSCet_);
	evType.tree->Branch("phoSCrawet", &phoSCrawet_);
	evType.tree->Branch("phoSCeta", &phoSCeta_);
	evType.tree->Branch("phoSCphi", &phoSCphi_);
	evType.tree->Branch("phoSCEn", &phoSCEn_);
	evType.tree->Branch("phoSCRawEn", &phoSCRawEn_);
	evType.tree->Branch("phoSCetaWidth", &phoSCetaWidth_);
	evType.tree->Branch("phoSCphiWidth", &phoSCphiWidth_);	
	
	std::cout<<"Created output tree:\t"<<typeName<<"\t"<<typeTitle<<std::endl;
	// evType.tree->Print();
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genPhoMatcher::fillEventType(eventType & evType){
	evType.tree->Fill();

	registerCutFlow(evType);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genPhoMatcher::registerAllCutFlow(){
	registerCutFlow(fullEB);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reset lastCutStep before each event
void genPhoMatcher::registerCutFlow(eventType & evType){
	evType.cutFlowCount->Fill(evType.lastCutStep);
	evType.cutFlowGenWeight->Fill(evType.lastCutStep, genWeight_);
	evType.lastCutStep = evType.lastCutStep + 1.;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

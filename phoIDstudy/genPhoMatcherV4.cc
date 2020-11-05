//////////////////////////////////////////////////////////////////////
//  Mohammad Abrar Wadud, Univeristy of Minnesota           		//
//	August/04/2020                                            		//
//  Photon ID Study   												//
//////////////////////////////////////////////////////////////////////

#include "/local/cms/user/gsorrent/antgc_analysis/macros/extra_tools.cc"

R__ADD_INCLUDE_PATH(/local/cms/user/wadud/aNTGCmet/xgboost/include/xgboost/)
R__LOAD_LIBRARY(/local/cms/user/wadud/aNTGCmet/xgboost/lib/libxgboost.so)
#include </local/cms/user/wadud/aNTGCmet/xgboost/include/xgboost/c_api.h>


#ifndef GENPHOMATCHER
#define GENPHOMATCHER

// Barrel-Endcap transition region
#define BETRetaMin 1.4442
#define BETRetaMax 1.566
#define HBetaMax 1.3920                     // ref. Josh H.
#define ZMASS 91.1876
#define pi7 3.1415927
#define REPORT_EVERY 10000
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
	genPhoMatcher(std::string FILELIST, std::string OUTFILE, Float_t XSECTION=-1., std::string MCPILEUPHIST="", std::string DATAPILEUPHIST="/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root",
		std::string PFECALCLUS_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoPFClusEcalIso.txt",
		std::string PFHCALCLUS_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoPFClusHcalIso.txt",
		std::string TKRISO_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoTrkSumPtHollowConeDR03.txt",
		std::string PFECALCLUS_PTSCALING="/hdfs/cms/user/wadud/anTGC/analysis_data/isoPtScaling/90pc/phoPFClusEcalIso.txt",
		std::string PFHCALCLUS_PTSCALING="/hdfs/cms/user/wadud/anTGC/analysis_data/isoPtScaling/90pc/phoPFClusHcalIso.txt",
		std::string BDT_PATH="/home/wadud/Public/trainingV2/aNTGC_photon_BDT_2020_09_12_19_08_19.model");

	~genPhoMatcher(){
		XGBoosterFree(phoBDT_h);
		std::cout<<"END @ "<<getCurrentTime()<<std::endl;
		std::cout<<"*************************************************************************************************************************************************"<<std::endl;
	};

private:
	Bool_t              isMC = false;
	Bool_t              doPUreweight = false;
	Float_t             xSec = -1.;

	effectiveAreaMap    PFHCALClusEffAreas;
	effectiveAreaMap    PFECALClusEffAreas;
	effectiveAreaMap    TkrEffAreas;

	isoPtScalingMap 	PFHCALClusPtScaling;
	isoPtScalingMap 	PFECALClusPtScaling;

	void                analyze();
	Bool_t              selectEvent();

	Short_t				photonIsPrompt(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	Bool_t				photonIsFake(Short_t _phoIndex, Float_t _deltaRmax = 0.3);
	Short_t 			matchWithRecoPho(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
        Short_t                         matchWithRecoEle(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	Short_t 			matchWithTrigPho(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
        Short_t                         electronIsTag(); 
        Short_t				photonIsProbe();

	TFile *             outFile = nullptr;

	/////////////////////////////////////////// Pileup Reweighting /////////////////////////////////////////////////////////
	PileupReWeighting   puReweighter;
	TH1F                pileupPreweight{"pileupUnweighted", "Unweighted Pileup; True # of Interactions", 200, 0., 200.};
	TH1F                pileupPostweight{"pileupWeighted", "Weighted Pileup; True # of Interactions", 200, 0., 200.};
	TH1F                rhoPreweight{"rhoUnweighted", "Unweighted #rho; #rho", 200, 0., 200.};
	TH1F                rhoPostweight{"rhoWeighted", "Weighted #rho; #rho", 200, 0., 200.};
	TH1F                nvtxPreweight{"nvtxUnweighted", "Unweighted # of Vertices; # of Vertices", 200, 0., 200.};
	TH1F                nvtxPostweight{"nvtxWeighted", "Weighted # of Vertices; # of Vertices", 200, 0., 200.};
        TH1D                invmass{"invmass", "Invariant egamma mass", 120, 70, 120};
        TH1D                invmass2{"invmass w/o selection", "Invariant egamma mass", 500, 0, 500};
        TH1D                deltaPhi{"deltaPhi", "deltaPhi", 100, -5, 5};
        TH1D                deltaEta{"deltaEta", "deltaEta", 100, -5, 5};
        TH1D                deltaRs{"deltaRs", "deltaR", 100, 0, 1};
        TH1D                deltaPt{"deltaPt", "Invariant egamma mass", 100, -100, 100};
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
	TTreeReaderAnyValue<UShort_t>		_nMC;
	TTreeReaderVectorValue<Int_t>           _mcPID;
	TTreeReaderVectorValue<Float_t>		_mcPt;
	TTreeReaderVectorValue<Float_t>         _mcEta;
	TTreeReaderVectorValue<Float_t>         _mcPhi;
	TTreeReaderVectorValue<UShort_t>        _mcStatusFlag;
	TTreeReaderVectorValue<Short_t>         _mcStatus;
	TTreeReaderVectorValue<Short_t>     	_mcIndex;
        //TTreeReaderVectorValue<Short_t>         _pho_gen_index;

	TTreeReaderAnyValue<UShort_t>           _nPho;
	TTreeReaderVectorValue<Float_t>         _phoCalibEt;
        TTreeReaderVectorValue<Float_t>         _phoCalibE;
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
	TTreeReaderVectorValue<Float_t>		_phoE5x5Full5x5;

	TTreeReaderVectorValue<Float_t>		_phoMaxEnergyXtal;
	TTreeReaderVectorValue<Float_t>		_phoE2ndFull5x5;
	TTreeReaderVectorValue<Float_t>		_phoE1x3Full5x5;
	TTreeReaderVectorValue<Float_t>		_phoE1x5Full5x5;
	TTreeReaderVectorValue<Float_t>		_phoE2x5Full5x5;

	TTreeReaderVectorValue<Float_t>		_phoPFClusEcalIso;
	TTreeReaderVectorValue<Float_t>		_phoPFClusHcalIso;
	TTreeReaderVectorValue<Float_t>		_phoTrkSumPtSolidConeDR04;
	TTreeReaderVectorValue<Float_t>		_phoTrkSumPtHollowConeDR04;
	TTreeReaderVectorValue<Float_t>		_phoTrkSumPtSolidConeDR03;
	TTreeReaderVectorValue<Float_t>		_phoTrkSumPtHollowConeDR03;
	TTreeReaderVectorValue<Float_t>		_phoECALIso;
	TTreeReaderVectorValue<Float_t>		_phoHCALIso;

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

	TTreeReaderAnyValue<Float_t>            _pfMET;
	TTreeReaderAnyValue<Float_t>            _pfMETPhi;
	TTreeReaderAnyValue<Float_t>            _pfMET_metSig;
	TTreeReaderAnyValue<Float_t>            _pfMET_EtSig;

	TTreeReaderAnyValue<UShort_t>           _nEle;
        TTreeReaderVectorValue<Float_t>         _elePt;
        TTreeReaderVectorValue<Float_t>         _eleEta;
        TTreeReaderVectorValue<Float_t>         _elePhi;
	TTreeReaderVectorValue<Float_t>         _eleCalibPt;
        TTreeReaderVectorValue<Float_t>         _eleCalibEn;
	TTreeReaderVectorValue<UChar_t>         _eleIDbit;

	TTreeReaderAnyValue<UShort_t>           _nMu;
	TTreeReaderVectorValue<Float_t>         _muPt;
	TTreeReaderVectorValue<Int_t>           _muIDbit;

	TTreeReaderAnyValue<UShort_t>           _nAK4CHSJet;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Pt;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Eta;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Phi;
	TTreeReaderVectorValue<Char_t>          _AK4CHSJet_PUFullID;
	TTreeReaderVectorValue<Char_t>          _AK4CHSJet_ID;


	TTreeReaderAnyValue<UChar_t>		_ntrgObjPho;
	TTreeReaderVectorValue<UChar_t>		_trgObjPhoBits;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoPt;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoEta;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoPhi;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////// Output variables //////////////////////////////////////////////////////////
	Char_t              fillPhoVars(Short_t _phoIndex, Short_t _eleIndex);

	Int_t		 run_;
	Long64_t	 event_;
	UShort_t	 lumis_;	
	Float_t	 	 rho_;
	UChar_t	 	 nVtx_;

	UShort_t		metFilters_;
	UShort_t		beamHaloSummary_;

	Float_t         genWeight_ = 1.;
	Float_t         puWeight_ = 1.;
	UChar_t		puTrue_;
	UShort_t	genStatusFlag_;
	Short_t		genStatus_;
	Float_t 	deltaRgenPho_;
	Float_t 	relDeltaPtGenPho_;
	Float_t 	deltaRPt_;
	Int_t		genPDGid_ = 0;

	Float_t 		deltaRtrg_;
	Float_t 		deltaPttrg_;

	Float_t 		phoPt_;
	Float_t 		phoEta_;
	Float_t 		phoPhi_;
        Float_t                 phoEn_;
	Float_t 		phoSeedTime_;

        Float_t                 elePt_;
        Float_t                 eleEta_;
        Float_t                 elePhi_;
        Float_t                 eleEn_;
	
	UChar_t			phoGenBits_;
	UChar_t			phoQualityBits_;
	Float_t 		phoR9Full5x5_;
	Float_t 		phoS4Full5x5_;

	Float_t 		phoEmaxOESCrFull5x5_;
	Float_t 		phoE2ndOESCrFull5x5_;
	Float_t 		phoE2ndOEmaxFull5x5_;
	Float_t 		phoE1x3OESCrFull5x5_;
	Float_t 		phoE2x5OESCrFull5x5_;
	Float_t 		phoE5x5OESCrFull5x5_;

	Float_t 		phoEmaxOE3x3Full5x5_;
	Float_t 		phoE2ndOE3x3Full5x5_;
	Float_t 		pho2x2OE3x3Full5x5_;
	Float_t 		phoSieieOSipipFull5x5_;
	Float_t			phoEtaWOPhiWFull5x5_;

	Float_t 		phoSigmaIEtaIEta_;
	Float_t 		phoSigmaIPhiIPhi_;
	Float_t 		phoSigmaIEtaIPhi_;

	Float_t 		phoE2x2Full5x5_;
	Float_t 		phoE3x3Full5x5_;
	Float_t 		phoE5x5Full5x5_;
	Float_t			phoMaxEnergyXtal_;
	Float_t			phoE2ndFull5x5_;
	Float_t			phoE1x3Full5x5_;
	Float_t			phoE1x5Full5x5_;
	Float_t			phoE2x5Full5x5_;

	Float_t			phoPFClusEcalIso_;
	Float_t			phoPFClusHcalIso_;
	Float_t			phoTrkSumPtSolidConeDR04_;
	Float_t			phoTrkSumPtHollowConeDR04_;
	Float_t			phoTrkSumPtSolidConeDR03_;
	Float_t			phoTrkSumPtHollowConeDR03_;
	Float_t			phoECALIso_;
	Float_t			phoHCALIso_;
	Float_t 		phoPFECALClusIsoCorr_;
	Float_t 		phoPFHCALClusIsoCorr_;
	Float_t 		phoTkrIsoCorr_;	

	Float_t 		phoHoverE_;
	Float_t 		phoPFChIsoRaw_;
	Float_t 		phoPFPhoIsoRaw_;
	Float_t 		phoPFNeuIsoRaw_;
	Float_t 		phoPFChWorstIsoRaw_;
	Float_t 		phoIDMVA_;
	UChar_t                 phoIDbit_;
	Float_t 		phoBDTpred_;
	UChar_t 		phoPFClusIDbits_;

	Bool_t 			pass95_ = 0;
	Bool_t 			pass90_ = 0;
	Bool_t 			pass80_ = 0;
	Bool_t			pass70_ = 0;

	Float_t 		phoMIP_;	

	Float_t 		phoSCet_;
	Float_t 		phoSCrawet_;
	Float_t 		phoSCeta_;
	Float_t 		phoSCphi_;
	Float_t 		phoSCEn_;
	Float_t 		phoSCRawEn_;
	Float_t 		phoEtaWidth_;
	Float_t 		phoPhiWidth_;

	UChar_t 		lepVeto_;

	Float_t 		met_;
	Float_t			metPhi_;
	Float_t 		metSig_;
	Float_t 		EtSig_;

	Float_t 		mT_;
	Float_t 		deltaPhiMetPho_;
	Float_t 		phoPtOverMet_;


	UChar_t 		nJet30_;
	Float_t 		Ht30_;
	Float_t 		minDeltaPhiMetJet30_;
	Float_t 		minDeltaRPhoJet30_;
        Double_t                InvMass_;
        Double_t                deltaPhi_;
        Double_t                deltaEta_;
        Double_t                deltaR_;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////Buffer variables///////////////////////////////////////////////////////
	Short_t			phoTrigMatch;
	Short_t			matchedGenPhoIndex;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////// XGBoost ////////////////////////////////////////////////////////////
	DMatrixHandle 		dTest;
	BoosterHandle 		phoBDT_h;
	Bool_t 			predictBDT = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////Categories////////////////////////////////////////////////////////////
	Bool_t          	initEventTypes();
	void            	initEventType(eventType & evType, std::string typeName, std::string typeTitle);
	void            	fillEventType(eventType & evType);
	void            	registerCutFlow(eventType & evType);
	void            	registerAllCutFlow();
	eventType		fullEB;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
genPhoMatcher::genPhoMatcher(std::string FILELIST, std::string OUTFILE, Float_t XSECTION, std::string MCPILEUPHIST, std::string DATAPILEUPHIST,
	std::string PFECALCLUS_EFFECTIVE_AREAS, std::string PFHCALCLUS_EFFECTIVE_AREAS, std::string TKRISO_EFFECTIVE_AREAS,
	std::string PFECALCLUS_PTSCALING,	std::string PFHCALCLUS_PTSCALING,
	std::string BDT_PATH){

	std::cout<<"*************************************************************************************************************************************************"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"Running genPhoMatcher"<<std::endl<<
	"\n\nInput parameters:"<<std::endl<<
	"\t\tFile list = "<<FILELIST<<std::endl<<
	"\t\tOutput file = "<<OUTFILE<<std::endl<<
	"\t\tCross section = "<<XSECTION<<std::endl<<
	"\t\tMC pileup histogram = "<<MCPILEUPHIST<<std::endl<<
	"\t\tData pileup histogram = "<<DATAPILEUPHIST<<std::endl<<
	"\t\tECAL PFCLuster Isolation Effective Areas = "<<PFECALCLUS_EFFECTIVE_AREAS<<std::endl<<
	"\t\tHCAL PFCLuster Isolation Effective Areas = "<<PFHCALCLUS_EFFECTIVE_AREAS<<std::endl<<
	"\t\tTracker Isolation Effective Areas = "<<TKRISO_EFFECTIVE_AREAS<<std::endl<<
	"\t\tECAL PFCLuster Isolation pT Scaling = "<<PFECALCLUS_PTSCALING<<std::endl<<
	"\t\tHCAL PFCLuster Isolation pT Scaling = "<<PFHCALCLUS_PTSCALING<<std::endl<<
	"\t\tBDT model file = "<<BDT_PATH<<std::endl;

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

	if(file_exists(PFECALCLUS_EFFECTIVE_AREAS)) {
		std::cout<<"PF ECAL Cluster effective areas:"<<std::endl;
		PFECALClusEffAreas.init(PFECALCLUS_EFFECTIVE_AREAS, 1, ",", 0);
	}

	if(file_exists(PFHCALCLUS_EFFECTIVE_AREAS)) {
		std::cout<<"PF HCAL Cluster effective areas:"<<std::endl;
		PFHCALClusEffAreas.init(PFHCALCLUS_EFFECTIVE_AREAS, 1, ",", 0);
	}
	
	if(file_exists(TKRISO_EFFECTIVE_AREAS)){
		std::cout<<"Tracker isolation effective areas:"<<std::endl;
		TkrEffAreas.init(TKRISO_EFFECTIVE_AREAS, 1, ",", 0);
	}

	if(file_exists(PFECALCLUS_PTSCALING)){
		std::cout<<"ECAL pT scaling:"<<std::endl;
		PFECALClusPtScaling.init(PFECALCLUS_PTSCALING, 0, 1, ",", 0);
	}

	if(file_exists(PFHCALCLUS_PTSCALING)){
		std::cout<<"HCAL pT scaling:"<<std::endl;
		PFHCALClusPtScaling.init(PFHCALCLUS_PTSCALING, 1, 1, ",", 0);
	}

	if(file_exists(BDT_PATH)){
		std::cout<<"\nLoading BDT model from "<<BDT_PATH <<std::endl;
		XGBoosterCreate(NULL, 0, &phoBDT_h);
		XGBoosterSetParam(phoBDT_h, "seed", "0");
		Int_t mLdSuccess = XGBoosterLoadModel(phoBDT_h, BDT_PATH.c_str());
		if(mLdSuccess == 0) predictBDT=1;
		else{
			std::cout<<"Failed to load BDT model!"<<std::endl;
		}
	}

	std::cout<<"\nCreating TChain... " <<std::endl;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t genPhoMatcher::electronIsTag(){

        Short_t eleTag = -999;
        Short_t highestPtEleIndex = -999;
        Short_t highestPt = -999;

        for(UShort_t iEle=0; iEle < _nEle; iEle++){

           if (_eleCalibPt[iEle] < 20) continue;

           UShort_t iHEEP = _eleIDbit[iEle];
           if(!getBit(iHEEP,4)) continue;

           if(_eleCalibPt[iEle] > highestPt){
             highestPtEleIndex = iEle;
             highestPt = _eleCalibPt[iEle];
           }
           eleTag = highestPtEleIndex;  
        }
        //std::cout << "eleTag: "<< eleTag << std::endl;
        return eleTag;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Short_t genPhoMatcher::photonIsProbe(){

        Short_t phoProbe = -999;
        Short_t highestPtPhoIndex = -999;
        Short_t highestPt = -999;
        
        for(UShort_t iPho=0; iPho<_nPho; iPho++){

           if(_phoCalibEt[iPho] < 200 ) continue;
           if(_phoHoverE[iPho] > 0.05 ) continue;

              if(_phoCalibEt[iPho] > highestPt){
              highestPtPhoIndex = iPho;
              highestPt = _phoCalibEt[iPho];
           }
           phoProbe = highestPtPhoIndex;
        }
        //std::cout << "phoProbe: "<< phoProbe << std::endl;
        return phoProbe;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	genPhoMatcher::matchWithRecoPho(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){

	Short_t matchedRecoPho = -999;
	Float_t minDeltaR = 999.;

	for(UShort_t iPho = 0; iPho < _nPho; iPho++){

		UChar_t iPhoFidReg = _phoFiducialRegion[iPho];
		if(!getBit(iPhoFidReg, 0)) continue; 			// skip if not EB (0 = EB, 1 = EE, 2 = EB-EE gap)

                if(_phoCalibEt[iPho] < 200.) continue;
                if(_phoHoverE[iPho] > 0.05 ) continue;

		Float_t relDeltaPtiGenPho = (_phoCalibEt[iPho] - _mcPt[_genIndex])/_mcPt[_genIndex];
		if(relDeltaPtiGenPho > _relDeltaPtMax) continue;
		if(relDeltaPtiGenPho < _relDeltaPtMin) continue;

		Float_t dRiGenPho = deltaR(_phoEta[iPho], _phoPhi[iPho], _mcEta[_genIndex], _mcPhi[_genIndex]);
		if(dRiGenPho > _deltaRmax) continue;

		if(dRiGenPho < minDeltaR){
			matchedRecoPho = iPho;
			minDeltaR = dRiGenPho;
		}
	}

	return matchedRecoPho;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Short_t genPhoMatcher::matchWithRecoEle(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){

        Short_t matchedRecoEle = -999;
        Float_t minDeltaR = 999.;

        for(UShort_t iEle = 0; iEle < _nEle; iEle++){

                UShort_t iHEEP = _eleIDbit[iEle];
                if(!getBit(iHEEP,4)) continue;

                if (_elePt[iEle] < 60) continue;
                 
                Float_t relDeltaPtiGenEle = (_elePt[iEle] - _mcPt[_genIndex])/_mcPt[_genIndex];
                if(relDeltaPtiGenEle > _relDeltaPtMax) continue;
                if(relDeltaPtiGenEle < _relDeltaPtMin) continue;
                Float_t dRiGenEle = deltaR(_eleEta[iEle], _elePhi[iEle], _mcEta[_genIndex], _mcPhi[_genIndex]);
                if(dRiGenEle > _deltaRmax) continue;

                if(dRiGenEle < minDeltaR){
                        matchedRecoEle = iEle;
                        minDeltaR = dRiGenEle;
                //std::cout << "electron index: " << iEle << "    over all electrons in the event: " << _nEle << endl;
                }
        }

        return matchedRecoEle;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	genPhoMatcher::matchWithTrigPho(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){
	Float_t minDeltaR = 999.;
	Short_t matchedTrigPho = -999;

	for(UShort_t iTrgPho = 0; iTrgPho < _ntrgObjPho; iTrgPho++){

		Float_t relDeltaPt = (_phoCalibEt[_phoIndex] - _trgObjPhoPt[iTrgPho])/_trgObjPhoPt[iTrgPho];
		if(relDeltaPt > _relDeltaPtMax) continue;
		if(relDeltaPt < _relDeltaPtMin) continue;

		Float_t dR = deltaR(_phoEta[_phoIndex], _phoPhi[_phoIndex], _trgObjPhoEta[iTrgPho], _trgObjPhoPhi[iTrgPho]);
		if(dR > _deltaRmax) continue;

		if(dR < minDeltaR){
			matchedTrigPho = iTrgPho;
			minDeltaR = dR;
		}
	}

	return matchedTrigPho;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t genPhoMatcher::selectEvent(){
	//// reset event cut flow
	fullEB.lastCutStep = 0.;
	registerAllCutFlow();

	// ULong64_t HLTPho    = (_HLTPho);                //// 9 = HLT_Photon200_v, 11 = HLT_Photon300_NoHE_v, 19 = HLT_Photon135_PFMET100_v, no trigger on 2016 MC
	// if(!getBit(HLTPho, 9))      return 0;           //// 200 GeV photon trigger
	// registerAllCutFlow();

        /*if (isMC) {
           Short_t highestPtGenIndexPho    = -9999;
           Float_t highestPtPho            = -999.;
           Short_t highestPtGenIndexEle    = -9999;
           Float_t highestPtEle            = -999.;
           Bool_t	passedPrompt		= 0;
           Bool_t 	passedGenFidCut		= 0;
           Bool_t 	passedGenPtCut		= 0;
           Bool_t 	passedRecoMatchPho	= 0;
           Bool_t       passedRecoMatchEle      = 0;
           Bool_t 	passedRecoPtCutPho 	= 0;
           Bool_t       passedRecoPtCutEle      = 0;
           Bool_t 	passedRecoFidCut	= 0;
           Bool_t       passedRecoHoverECut     = 0;
           // Bool_t 	passedPixelCut 		= 0;
           Short_t matchedRecoPhoton 	= -999;
           Short_t matchedRecoElectron    = -999;
           
           for(UShort_t iGenP=0; iGenP<_nMC; iGenP++){

              UShort_t iGenPStFl = _mcStatusFlag[iGenP];
              if(!getBit(iGenPStFl,1)) continue;
              passedPrompt = 1;

              if(std::abs(_mcEta[iGenP]) > 1.6) continue;
              passedGenFidCut = 1;

              if(fabs(_mcPID[iGenP]) == 11) {

                 Short_t iGenRecoMatchEle = matchWithRecoEle(iGenP, 0.05, -0.1, 0.1);
                 if(iGenRecoMatchEle < 0) continue;
                 passedRecoMatchEle = 1;

                 //if (_mcPt[iGenRecoMatchEle] < 60) continue;
                 passedRecoPtCutEle = 1;

                 //UShort_t iHEEP = _eleIDbit[iGenRecoMatchEle];
                 //if(!getBit(iHEEP,4)) continue;

                 if(_mcPt[iGenP] > highestPtEle){
                    highestPtGenIndexEle = iGenP;
                    highestPtEle = _mcPt[iGenP];
                    matchedRecoElectron = iGenRecoMatchEle;
                 }
                 //std::cout << "matchedRecoElectron: " << matchedRecoElectron << std::endl;
              }

              if(_mcPID[iGenP] == 22) {

                 Short_t iGenRecoMatchPho = matchWithRecoPho(iGenP, 0.05, -0.1, 0.1);
                 if(iGenRecoMatchPho < 0) continue;
                 passedRecoMatchPho = 1;
           
                 //if(std::abs(_ecalSCeta[_phoDirectEcalSCindex[iGenRecoMatchPho]]) > BETRetaMin) continue;
                 passedRecoFidCut = 1;
           
                 //if(_phoCalibEt[iGenRecoMatchPho] < 200.) continue;
                 passedRecoPtCutPho = 1;
              
                 //if(_phoHoverE[iGenRecoMatchPho] > 0.05 ) continue;
                 passedRecoHoverECut= 1;
   
                 // UChar_t iQualityBits = _phoQualityBits[iGenRecoMatch];
                 // if(getBit(iQualityBits,0)) continue;                             	// 0=has pixel seed, 1=electron veto
                 // passedPixelCut = 1;
           
                 if(_mcPt[iGenP] > highestPtPho){
                    highestPtGenIndexPho = iGenP;
                    highestPtPho = _mcPt[iGenP];
                    matchedRecoPhoton = iGenRecoMatchPho;
                 }
                 //std::cout << "matchedRecoPhoton: " << matchedRecoPhoton << std::endl;
              }

           }

           if(passedPrompt) registerAllCutFlow();
           if(passedGenFidCut) registerAllCutFlow();
           if(passedRecoMatchPho) registerAllCutFlow();
           if(passedRecoMatchEle) registerAllCutFlow();
           if(passedRecoFidCut) registerAllCutFlow();
           //if(passedGenPtCut) registerAllCutFlow();
           if(passedRecoPtCutPho) registerAllCutFlow();
           if(passedRecoPtCutEle) registerAllCutFlow();
           if(passedRecoHoverECut)registerAllCutFlow();
           // if(passedPixelCut) registerAllCutFlow();
           
           if((matchedRecoElectron < 0) || (matchedRecoPhoton < 0)) return 0;        

           if((matchedRecoElectron >= 0) && (matchedRecoPhoton >= 0)) {
              matchedGenPhoIndex = highestPtGenIndexPho;

              TLorentzVector v_ele(0.,0.,0.,0.);
              TLorentzVector v_pho(0.,0.,0.,0.);

              v_ele.SetPtEtaPhiE(_eleCalibPt[matchedRecoElectron], _eleEta[matchedRecoElectron], _elePhi[matchedRecoElectron], _eleCalibEn[matchedRecoElectron]);
              v_pho.SetPtEtaPhiE(_phoCalibEt[matchedRecoPhoton], _phoEta[matchedRecoPhoton], _phoPhi[matchedRecoPhoton], _phoCalibE[matchedRecoPhoton]);

              Double_t eg_mass = (v_ele+v_pho).M();
              if ((eg_mass > 80.) && (eg_mass < 110.)) {
                 std::cout << "invariant mass: " << eg_mass << std::endl;
                 fillPhoVars(matchedRecoPhoton, matchedRecoElectron); ////attenzione a questo indice, matchedRecoPhoton in origine
                 fillEventType(fullEB);
                 return 1;
              }
           }
        }*/

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         
        /*Short_t tagElectron = electronIsTag();
        Short_t probePhoton = photonIsProbe();

        if((tagElectron < 0) || (probePhoton < 0)) return 0;

        if((tagElectron >= 0) && (probePhoton >= 0)) {
           cout << "tagEle" << tagElectron << endl;
           std::cout << "Electron pt: "<<_eleCalibPt[tagElectron] <<   std::endl;
           cout << "probePho" << probePhoton << endl;
           std::cout << "Photon pt: "<<_phoCalibEt[probePhoton] << "   phoHoverE: " << _phoHoverE[probePhoton] <<   std::endl;
           TLorentzVector v_ele(0.,0.,0.,0.);
           TLorentzVector v_pho(0.,0.,0.,0.);
           cout << "_phoEta" << _phoEta[probePhoton] << "    _phoCalibE" << _phoCalibE[probePhoton] <<  std::endl;
           v_ele.SetPtEtaPhiE(_eleCalibPt[tagElectron], _eleEta[tagElectron], _elePhi[tagElectron], _eleCalibEn[tagElectron]);
           v_pho.SetPtEtaPhiE(_phoCalibEt[probePhoton], _phoEta[probePhoton], _phoPhi[probePhoton], _phoCalibE[probePhoton]);*/

           /*Double_t deltaPhi_ = fabs(_elePhi[tagElectron] - _phoPhi[probePhoton]);
           if (deltaPhi_ > acos(-1)) {
              deltaPhi_ = 2*acos(-1) - deltaPhi_;
           }
           Double_t deltaEta_ = _eleEta[tagElectron] - _phoEta[probePhoton];
           Double_t deltaR_ = sqrt((deltaPhi_*deltaPhi_) + (deltaEta_*deltaEta_));

           if (deltaR_ < 0.3) return 0;*/
           //Double_t eg_mass = (v_ele+v_pho).M();
           //if ((eg_mass > 80.) && (eg_mass < 110.)) {
              //fillPhoVars(probePhoton, tagElectron);
              //fillEventType(fullEB);
              return 1;
          // } 
           //else {
              //return 0;
           //}
        //}
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
        invmass.SetDirectory(outFile->GetDirectory(""));
        invmass2.SetDirectory(outFile->GetDirectory(""));
        deltaPhi.SetDirectory(outFile->GetDirectory(""));
        deltaEta.SetDirectory(outFile->GetDirectory(""));
        deltaRs.SetDirectory(outFile->GetDirectory(""));
        deltaPt.SetDirectory(outFile->GetDirectory(""));

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
				puWeight_ = 1.;    
			}
			pileupPreweight.Fill(_puTrue, genWeight_);
		}

		selectEvent();

               Short_t tagElectron = -999;
               Short_t probePhoton = -999;
               Short_t tagEle = -999;
               Short_t probePho = -999;
               Short_t highestPtEleIndex = -999;
               Short_t highestPt = -999;
               Double_t eg_mass = 0;

               for(UShort_t iEle=0; iEle < _nEle; iEle++){

                  UShort_t iHEEP = _eleIDbit[iEle];
                  if( (_eleCalibPt[iEle] > 60) && (getBit(iHEEP,4)) ) {
                     if(_eleCalibPt[iEle] > highestPt){
                        highestPt = _eleCalibPt[iEle];
                        tagElectron = iEle;
                     }
                  }
                  if (tagElectron < 0) continue;
                  TLorentzVector v_ele(0.,0.,0.,0.);
                  v_ele.SetPtEtaPhiM(_eleCalibPt[tagElectron], _eleEta[tagElectron], _elePhi[tagElectron], 0.511/1000.);  

                  //Short_t probePhoton = -999;
                  Short_t highestPtPhoIndex = -999;
                  Short_t highestPt2 = -999; 
 
                  for(UShort_t iPho=0; iPho<_nPho; iPho++){
                     if( (_phoCalibEt[iPho] > 200 ) && (_phoHoverE[iPho] < 0.05 ) ) {
                        if(_phoCalibEt[iPho] > highestPt2){
                           highestPt2 = _phoCalibEt[iPho];
                           probePhoton = iPho;
                        }
                     }
                     if (probePhoton < 0) continue;
                     TLorentzVector v_pho(0.,0.,0.,0.);
                     v_pho.SetPtEtaPhiM(_phoCalibEt[probePhoton], _phoEta[probePhoton], _phoPhi[probePhoton], 0);
                     
                     eg_mass = (v_ele+v_pho).M();
                     if ((eg_mass > 80.) && (eg_mass < 110.)) {
                        tagEle = tagElectron;
                        probePho = probePhoton;        
                     }
                  }
               }

               if((tagEle >= 0) && (probePho >= 0)) {

                  Double_t deltaPhi_ = fabs(_elePhi[tagEle] - _phoPhi[probePho]);
                  Double_t deltaEta_ = _eleEta[tagEle] - _phoEta[probePho];
                  Double_t deltaR_ = sqrt((deltaPhi_*deltaPhi_) + (deltaEta_*deltaEta_));
                  Double_t deltaPt_ = _eleCalibPt[tagEle] - _phoCalibEt[probePho];
                  if (deltaPhi_ > acos(-1)) {
                     deltaPhi_ = 2*acos(-1) - deltaPhi_;
                  }
                  cout << "Tag found: " << tagElectron << std::endl;
                  std::cout << "Electron pt: "<<_eleCalibPt[tagEle] <<   std::endl;
                  cout << "Probe found: " << probePhoton << std::endl;
                  std::cout << "Photon pt: "<<_phoCalibEt[probePho] << std::endl;

                  invmass.Fill(eg_mass, genWeight_);
                  deltaPhi.Fill(deltaPhi_, genWeight_);
                  deltaEta.Fill(deltaEta_, genWeight_);
                  deltaRs.Fill(deltaR_, genWeight_);
                  deltaPt.Fill(deltaPt_, genWeight_);
                  
                  fillPhoVars(probePhoton, tagElectron);
                  fillEventType(fullEB);
               } //tag and probe > 0 cycle 
     
		current_entry++; 
	};

	std::cout<<"Done analyzing!"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"--------------------------------------------------------------------------------------------------"<<std::endl;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Char_t genPhoMatcher::fillPhoVars(Short_t _phoIndex, Short_t _eleIndex){

	run_ 					= _run;
	event_ 					= _event;
	lumis_ 					= _lumis;
	rho_ 					= _rho;
	nVtx_					= _nVtx;

	metFilters_ 		= _metFilters;
	beamHaloSummary_	= _beamHaloSummary;

	Short_t phoSCindex 	= _phoDirectEcalSCindex[_phoIndex];
	phoSCeta_ 		= _ecalSCeta[phoSCindex];
	Float_t phoAbsSCEta 	= std::abs(phoSCeta_);

	phoPt_ 			= _phoCalibEt[_phoIndex];
	phoEta_			= _phoEta[_phoIndex];
	phoPhi_			= _phoPhi[_phoIndex];
        phoEn_                  = _phoCalibE[_phoIndex];
	phoSeedTime_ 		= _phoSeedTime[_phoIndex];

        elePt_                  = _eleCalibPt[_eleIndex];
        eleEta_                 = _eleEta[_eleIndex];
        elePhi_                 = _elePhi[_eleIndex];
        eleEn_                  = _eleCalibEn[_eleIndex];

	phoQualityBits_		= _phoQualityBits[_phoIndex];
	phoR9Full5x5_ 		= _phoR9Full5x5[_phoIndex];
	phoS4Full5x5_		= _phoE2x2Full5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoEmaxOESCrFull5x5_ 	= _phoMaxEnergyXtal[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoE2ndOESCrFull5x5_ 	= _phoE2ndFull5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoE2ndOEmaxFull5x5_ 	= _phoE2ndFull5x5[_phoIndex]/_phoMaxEnergyXtal[_phoIndex];
	phoE1x3OESCrFull5x5_ 	= _phoE1x3Full5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoE2x5OESCrFull5x5_ 	= _phoE2x5Full5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];
	phoE5x5OESCrFull5x5_ 	= _phoE5x5Full5x5[_phoIndex]/_ecalSCRawEn[phoSCindex];


	phoSigmaIEtaIEta_ 	= _phoSigmaIEtaIEtaFull5x5[_phoIndex];
	phoSigmaIEtaIPhi_ 	= _phoSigmaIEtaIPhiFull5x5[_phoIndex];
	phoSigmaIPhiIPhi_	= _phoSigmaIPhiIPhiFull5x5[_phoIndex];

	phoE2x2Full5x5_		= _phoE2x2Full5x5[_phoIndex];
	phoE3x3Full5x5_		= phoR9Full5x5_ * _ecalSCRawEn[phoSCindex];
	phoE5x5Full5x5_		= _phoE5x5Full5x5[_phoIndex];
	phoMaxEnergyXtal_	= _phoMaxEnergyXtal[_phoIndex];
	phoE2ndFull5x5_		= _phoE2ndFull5x5[_phoIndex];
	phoE1x3Full5x5_		= _phoE1x3Full5x5[_phoIndex];
	phoE1x5Full5x5_		= _phoE1x5Full5x5[_phoIndex];
	phoE2x5Full5x5_		= _phoE2x2Full5x5[_phoIndex];

	phoEmaxOE3x3Full5x5_ 	= phoMaxEnergyXtal_/phoE3x3Full5x5_;
	phoE2ndOE3x3Full5x5_ 	= phoE2ndFull5x5_/phoE3x3Full5x5_;
	pho2x2OE3x3Full5x5_	= phoE2x2Full5x5_/phoE3x3Full5x5_;
	phoSieieOSipipFull5x5_	= phoSigmaIEtaIEta_/phoSigmaIPhiIPhi_;

	phoPFClusEcalIso_		= _phoPFClusEcalIso[_phoIndex];
	phoPFClusHcalIso_ 		= _phoPFClusHcalIso[_phoIndex];
	phoTrkSumPtSolidConeDR04_	= _phoTrkSumPtSolidConeDR04[_phoIndex];
	phoTrkSumPtHollowConeDR04_ 	= _phoTrkSumPtHollowConeDR04[_phoIndex];
	phoTrkSumPtSolidConeDR03_	= _phoTrkSumPtSolidConeDR03[_phoIndex];
	phoTrkSumPtHollowConeDR03_ 	= _phoTrkSumPtHollowConeDR03[_phoIndex];
	phoECALIso_ 			= _phoECALIso[_phoIndex];
	phoHCALIso_ 			= _phoHCALIso[_phoIndex];

	phoPFECALClusIsoCorr_ 		= 	phoPFClusEcalIso_ - rho_ * PFECALClusEffAreas.getEffectiveArea(phoAbsSCEta) - PFECALClusPtScaling.getPtScaling(phoAbsSCEta, phoPt_);
	phoPFHCALClusIsoCorr_ 		= 	phoPFClusHcalIso_ - rho_ * PFHCALClusEffAreas.getEffectiveArea(phoAbsSCEta) - PFHCALClusPtScaling.getPtScaling(phoAbsSCEta, phoPt_);
	phoTkrIsoCorr_ 			= 	phoTrkSumPtHollowConeDR03_ - rho_ * TkrEffAreas.getEffectiveArea(phoAbsSCEta);

	phoHoverE_ 		= _phoHoverE[_phoIndex];
	phoPFChIsoRaw_ 		= _phoPFChIso[_phoIndex];
	phoPFPhoIsoRaw_		= _phoPFPhoIso[_phoIndex];
	phoPFNeuIsoRaw_		= _phoPFNeuIso[_phoIndex];
	phoPFChWorstIsoRaw_ 	= _phoPFChWorstIso[_phoIndex];
	phoIDMVA_		= _phoIDMVA[_phoIndex];
	phoIDbit_		= _phoIDbit[_phoIndex];
	phoMIP_ 		= _phoMIPTotEnergy[_phoIndex];

	phoSCet_			= (_ecalSCEn[phoSCindex]) / std::cosh(phoSCeta_);
	phoSCrawet_			= (_ecalSCRawEn[phoSCindex]) / std::cosh(phoSCeta_);
	phoSCphi_ 			= _ecalSCphi[phoSCindex];
	phoSCEn_ 			= _ecalSCEn[phoSCindex];
	phoSCRawEn_			= _ecalSCRawEn[phoSCindex];
	phoEtaWidth_ 			= _ecalSCetaWidth[phoSCindex];
	phoPhiWidth_ 			= _ecalSCphiWidth[phoSCindex];

	phoEtaWOPhiWFull5x5_	= phoEtaWidth_/phoPhiWidth_;

	if(isMC){
		puTrue_			= 	_puTrue;
		phoGenBits_ 		= 	0;
		Short_t iPhoGenIndex 	= 	matchedGenPhoIndex;

		if(iPhoGenIndex > -1) {
			setBit(phoGenBits_,0,1);	
			genStatusFlag_ 			= 	_mcStatusFlag[iPhoGenIndex];
			genStatus_ 			=	_mcStatus[iPhoGenIndex];
			deltaRgenPho_ 			=	deltaR(_mcEta[iPhoGenIndex], _mcPhi[iPhoGenIndex], phoEta_, phoPhi_);
			relDeltaPtGenPho_ 		=	(_phoCalibEt[_phoIndex] - _mcPt[iPhoGenIndex])/_mcPt[iPhoGenIndex];
			deltaRPt_ 			= std::sqrt(deltaRgenPho_*deltaRgenPho_ + relDeltaPtGenPho_*relDeltaPtGenPho_);
			genPDGid_ 			= 	_mcPID[iPhoGenIndex];
		} else {
			genStatusFlag_ 			=  9999;
			genStatus_ 			= -9999;
			deltaRgenPho_ 			= -9999;
			relDeltaPtGenPho_ 		= -9999;
			deltaRPt_ 			= -9999;
			genPDGid_			= -9999;
		}

		setBit(phoGenBits_, 1, 1);
	}

	phoTrigMatch = matchWithTrigPho(_phoIndex, 0.3, 0.5, 1.5);

	if(phoTrigMatch > -1){
		deltaRtrg_  = deltaR(phoEta_, phoPhi_, _trgObjPhoEta[phoTrigMatch], _trgObjPhoPhi[phoTrigMatch]);
		deltaPttrg_ = (_trgObjPhoPt[phoTrigMatch] - phoPt_)/_trgObjPhoPt[phoTrigMatch];
	} else{
		deltaRtrg_  = - 999;
		deltaPttrg_ = - 999;
	}

	phoPFClusIDbits_ 	= 0;

	if(predictBDT){

		std::vector<Float_t> feats{pho2x2OE3x3Full5x5_, phoE1x3OESCrFull5x5_, phoE2ndOESCrFull5x5_, phoE2x5OESCrFull5x5_, phoE5x5OESCrFull5x5_, phoEmaxOESCrFull5x5_, phoEtaWOPhiWFull5x5_, phoEtaWidth_, phoPhiWidth_, phoR9Full5x5_, phoS4Full5x5_, phoSieieOSipipFull5x5_, phoSigmaIEtaIEta_, phoSigmaIEtaIPhi_, phoSigmaIPhiIPhi_};

		XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
		bst_ulong out_len;
		const float *prediction;
		XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
		assert(out_len == 1);
		XGDMatrixFree(dTest);
		phoBDTpred_ = prediction[0];

		pass95_ = (phoBDTpred_ >= 2.1347871547393873e-01) && (phoHoverE_ < 3.1972081872819871e-02) && (phoPFECALClusIsoCorr_ < 5.5702003247919532e+00) && (phoPFHCALClusIsoCorr_ < 1.5870381409324537e+01) && (phoTkrIsoCorr_ < 3.9988330612434098e+00);
		pass90_ = (phoBDTpred_ >= 4.2242258443976971e-01) && (phoHoverE_ < 4.0873656933608248e-02) && (phoPFECALClusIsoCorr_ < 3.7277726987096962e+00) && (phoPFHCALClusIsoCorr_ < 1.1052528830800725e+01) && (phoTkrIsoCorr_ < 2.8711557036145372e+00);
		pass80_ = (phoBDTpred_ >= 8.7244571931769122e-01) && (phoHoverE_ < 4.4450510008151485e-02) && (phoPFECALClusIsoCorr_ < 4.7979935544225469e+00) && (phoPFHCALClusIsoCorr_ < 1.0128719656374308e+01) && (phoTkrIsoCorr_ < 4.1483811040263188e+00);
		pass70_ = (phoBDTpred_ >= 9.5084052145385245e-01) && (phoHoverE_ < 4.8348822338560769e-02) && (phoPFECALClusIsoCorr_ < 3.4110871397501707e+00) && (phoPFHCALClusIsoCorr_ < 1.5238652601828349e+01) && (phoTkrIsoCorr_ < 4.0837764862499419e+00);

		if(pass70_) setBit(phoPFClusIDbits_, 0, 1);
		if(pass80_) setBit(phoPFClusIDbits_, 1, 1);
		if(pass90_) setBit(phoPFClusIDbits_, 2, 1);
		if(pass95_) setBit(phoPFClusIDbits_, 3, 1);
	}

        TLorentzVector v_ele(0.,0.,0.,0.);
        TLorentzVector v_pho(0.,0.,0.,0.);

        v_ele.SetPtEtaPhiE(_eleCalibPt[_eleIndex], _eleEta[_eleIndex], _elePhi[_eleIndex], _eleCalibEn[_eleIndex]);
        v_pho.SetPtEtaPhiE(_phoCalibEt[_phoIndex], _phoEta[_phoIndex], _phoPhi[_phoIndex], _phoCalibE[_phoIndex]);

        Double_t eg_mass = (v_ele+v_pho).M();
        InvMass_ = eg_mass;

        deltaPhi_ = fabs(_elePhi[_eleIndex] - _phoPhi[_phoIndex]);
        if (deltaPhi_ > acos(-1)) {
           deltaPhi_ = 2*acos(-1) - deltaPhi_;
        }
        deltaEta_ = _eleEta[_eleIndex] - _phoEta[_phoIndex];
        deltaR_ = sqrt((deltaPhi_*deltaPhi_) + (deltaEta_*deltaEta_));

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
	_beamHaloSummary.set(inputTTreeReader, "beamHaloSummary");
	_metFilters.set(inputTTreeReader, "metFilters");
	_nVtx.set(inputTTreeReader, "nVtx");
	_rho.set(inputTTreeReader, "rho");
	_HLTPho.set(inputTTreeReader, "HLTPho");
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
		_mcIndex.set(inputTTreeReader, "mcIndex");
                //_pho_gen_index.set(inputTTreeReader, "pho_gen_index");
	}

	_nPho.set(inputTTreeReader, "nPho");
	_phoEt.set(inputTTreeReader, "phoEt");
	_phoCalibEt.set(inputTTreeReader, "phoCalibEt");
        _phoCalibE.set(inputTTreeReader, "phoCalibE");
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

	_metFilters.set(inputTTreeReader, "metFilters");
	_pfMET.set(inputTTreeReader, "pfMET");
	_pfMETPhi.set(inputTTreeReader, "pfMETPhi");
	_pfMET_metSig.set(inputTTreeReader, "pfMET_metSig");
	_pfMET_EtSig.set(inputTTreeReader, "pfMET_EtSig");

	_nEle.set(inputTTreeReader, "nEle");
        _elePt.set(inputTTreeReader, "elePt");
        _eleEta.set(inputTTreeReader, "eleEta");
        _elePhi.set(inputTTreeReader, "elePhi");
	_eleCalibPt.set(inputTTreeReader, "eleCalibPt");
        _eleCalibEn.set(inputTTreeReader, "eleCalibEn");
	_eleIDbit.set(inputTTreeReader, "eleIDbit");

	_nMu.set(inputTTreeReader, "nMu");
	_muPt.set(inputTTreeReader, "muPt");
	_muIDbit.set(inputTTreeReader, "muIDbit");

	_nAK4CHSJet.set(inputTTreeReader, "nAK4CHSJet");
	_AK4CHSJet_Pt.set(inputTTreeReader, "AK4CHSJet_Pt");
	_AK4CHSJet_Eta.set(inputTTreeReader, "AK4CHSJet_Eta");
	_AK4CHSJet_Phi.set(inputTTreeReader, "AK4CHSJet_Phi");
	_AK4CHSJet_PUFullID.set(inputTTreeReader, "AK4CHSJet_PUFullID");
	_AK4CHSJet_ID.set(inputTTreeReader, "AK4CHSJet_ID");

	_ntrgObjPho.set(inputTTreeReader, "ntrgObjPho");
	_trgObjPhoBits.set(inputTTreeReader, "trgObjPhoBits");
	_trgObjPhoPt.set(inputTTreeReader, "trgObjPhoPt");
	_trgObjPhoEta.set(inputTTreeReader, "trgObjPhoEta");
	_trgObjPhoPhi.set(inputTTreeReader, "trgObjPhoPhi");

	std::cout<<"Branches initialized!"<<std::endl<<
	"**************************************************************************************************************************************************************"<<std::endl;
	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genPhoMatcher::initEventType(eventType & evType, std::string typeName, std::string typeTitle){

	mkTFileDir(outFile, typeName);

	evType.cutFlowCount = new TH1F((typeName+"cutFlowCount").c_str(), "Cut Flow (Unweighted)", CUTFLOWSTEPS, 0., (Float_t)CUTFLOWSTEPS);
	evType.cutFlowGenWeight = new TH1F((typeName+"cutFlowGenWeight").c_str(), "Cut Flow (Gen Weighted)", CUTFLOWSTEPS, 0., (Float_t)CUTFLOWSTEPS);
	evType.cutFlowCount->SetDirectory(outFile->GetDirectory(typeName.c_str()));
	evType.cutFlowGenWeight->SetDirectory(outFile->GetDirectory(typeName.c_str()));

	//////// initialize tree
	evType.tree = new TTree((typeName+"Tree").c_str(), typeTitle.c_str());
	evType.tree->SetDirectory(outFile->GetDirectory(typeName.c_str()));

	evType.tree->Branch("run", &run_);
	evType.tree->Branch("event", &event_);
	evType.tree->Branch("lumis", &lumis_);
	evType.tree->Branch("rho", &rho_);
	evType.tree->Branch("nVtx", &nVtx_);
	evType.tree->Branch("puTrue", &puTrue_);
	evType.tree->Branch("metFilters", &metFilters_);
	evType.tree->Branch("beamHaloSummary", &beamHaloSummary_);

	evType.tree->Branch("puWeight", &puWeight_);
	evType.tree->Branch("genWeight", &genWeight_);
	evType.tree->Branch("genStatusFlag", &genStatusFlag_);
	evType.tree->Branch("genStatus", &genStatus_);
	evType.tree->Branch("deltaRgenPho", &deltaRgenPho_);
	evType.tree->Branch("relDeltaPtGenPho", &relDeltaPtGenPho_);
	evType.tree->Branch("deltaRPt", &deltaRPt_);
	evType.tree->Branch("genPDGid", &genPDGid_);	

	evType.tree->Branch("deltaRtrg", &deltaRtrg_);
	evType.tree->Branch("deltaPttrg", &deltaPttrg_);

	evType.tree->Branch("phoPt", &phoPt_);
	evType.tree->Branch("phoEta", &phoEta_);
	evType.tree->Branch("phoPhi", &phoPhi_);
        evType.tree->Branch("phoEn", &phoEn_);
	evType.tree->Branch("phoSeedTime", &phoSeedTime_);

        evType.tree->Branch("elePt", &elePt_);
        evType.tree->Branch("eleEta", &eleEta_);
        evType.tree->Branch("elePhi", &elePhi_);
        evType.tree->Branch("eleEn", &eleEn_);

	evType.tree->Branch("phoGenBits", &phoGenBits_);
	evType.tree->Branch("phoQualityBits", &phoQualityBits_);

	evType.tree->Branch("phoR9Full5x5", &phoR9Full5x5_);
	evType.tree->Branch("phoS4Full5x5", &phoS4Full5x5_);
	evType.tree->Branch("phoEmaxOESCrFull5x5", &phoEmaxOESCrFull5x5_);
	evType.tree->Branch("phoE2ndOESCrFull5x5", &phoE2ndOESCrFull5x5_);
	evType.tree->Branch("phoE2ndOEmaxFull5x5", &phoE2ndOEmaxFull5x5_);	
	evType.tree->Branch("phoE1x3OESCrFull5x5", &phoE1x3OESCrFull5x5_);
	evType.tree->Branch("phoE2x5OESCrFull5x5", &phoE2x5OESCrFull5x5_);
	evType.tree->Branch("phoE5x5OESCrFull5x5", &phoE5x5OESCrFull5x5_);
	evType.tree->Branch("phoEmaxOE3x3Full5x5", &phoEmaxOE3x3Full5x5_);
	evType.tree->Branch("phoE2ndOE3x3Full5x5", &phoE2ndOE3x3Full5x5_);
	evType.tree->Branch("pho2x2OE3x3Full5x5", &pho2x2OE3x3Full5x5_);
	evType.tree->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta_);
	evType.tree->Branch("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi_);
	evType.tree->Branch("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi_);
	evType.tree->Branch("phoSieieOSipipFull5x5", &phoSieieOSipipFull5x5_);
	evType.tree->Branch("phoEtaWidth", &phoEtaWidth_);
	evType.tree->Branch("phoPhiWidth", &phoPhiWidth_);
	evType.tree->Branch("phoEtaWOPhiWFull5x5", &phoEtaWOPhiWFull5x5_);

	evType.tree->Branch("phoMaxEnergyXtal", &phoMaxEnergyXtal_);
	evType.tree->Branch("phoE2ndFull5x5", &phoE2ndFull5x5_);
	evType.tree->Branch("phoE2x2Full5x5", &phoE2x2Full5x5_);
	evType.tree->Branch("phoE3x3Full5x5", &phoE3x3Full5x5_);
	evType.tree->Branch("phoE5x5Full5x5", &phoE5x5Full5x5_);
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
	evType.tree->Branch("phoPFECALClusIsoCorr", &phoPFECALClusIsoCorr_);
	evType.tree->Branch("phoPFHCALClusIsoCorr", &phoPFHCALClusIsoCorr_);
	evType.tree->Branch("phoTkrIsoCorr", &phoTkrIsoCorr_);

	evType.tree->Branch("phoHoverE", &phoHoverE_);
	evType.tree->Branch("phoPFChIsoRaw", &phoPFChIsoRaw_);
	evType.tree->Branch("phoPFPhoIsoRaw", &phoPFPhoIsoRaw_);
	evType.tree->Branch("phoPFNeuIsoRaw", &phoPFNeuIsoRaw_);
	evType.tree->Branch("phoPFChWorstIsoRaw", &phoPFChWorstIsoRaw_);
	evType.tree->Branch("phoEGMidMVA", &phoIDMVA_);
	evType.tree->Branch("phoBDTpred", &phoBDTpred_);
	evType.tree->Branch("phoPFClusIDbits", &phoPFClusIDbits_);	
	evType.tree->Branch("phoIDbit", &phoIDbit_);

	evType.tree->Branch("pass95", &pass95_);	
	evType.tree->Branch("pass90", &pass90_);
	evType.tree->Branch("pass80", &pass80_);
	evType.tree->Branch("pass70", &pass70_);

	evType.tree->Branch("phoMIP", &phoMIP_);
	
	evType.tree->Branch("phoSCet", &phoSCet_);
	evType.tree->Branch("phoSCrawet", &phoSCrawet_);
	evType.tree->Branch("phoSCeta", &phoSCeta_);
	evType.tree->Branch("phoSCphi", &phoSCphi_);
	evType.tree->Branch("phoSCEn", &phoSCEn_);
	evType.tree->Branch("phoSCRawEn", &phoSCRawEn_);

	evType.tree->Branch("lepVeto", &lepVeto_);
	
	evType.tree->Branch("met", &met_);
	evType.tree->Branch("metPhi", &metPhi_);
	evType.tree->Branch("metSig", &metSig_);
	evType.tree->Branch("EtSig", &EtSig_);

	evType.tree->Branch("mT", &mT_);
	evType.tree->Branch("deltaPhiMetPho", &deltaPhiMetPho_);
	evType.tree->Branch("phoPtOverMet", &phoPtOverMet_);

	evType.tree->Branch("nJet30", &nJet30_);
	evType.tree->Branch("Ht30", &Ht30_);
	evType.tree->Branch("minDeltaPhiMetJet30", &minDeltaPhiMetJet30_);
	evType.tree->Branch("minDeltaRPhoJet30", &minDeltaRPhoJet30_);

        evType.tree->Branch("InvMass", &InvMass_);
        evType.tree->Branch("deltaPhi", &deltaPhi_);
        evType.tree->Branch("deltaEta", &deltaEta_);
        	
	std::cout<<"Created output tree:\t"<<typeName<<"\t"<<typeTitle<<std::endl<<std::endl;
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

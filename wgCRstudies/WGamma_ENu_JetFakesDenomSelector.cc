//////////////////////////////////////////////////////////////////////
//  Mohammad Abrar Wadud, Univeristy of Minnesota           		//
//	August/04/2020                                            		//
//  Photon ID Study   												//
//////////////////////////////////////////////////////////////////////

#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"

R__ADD_INCLUDE_PATH(/local/cms/user/wadud/aNTGCmet/xgboost/include/xgboost/)
R__LOAD_LIBRARY(/local/cms/user/wadud/aNTGCmet/xgboost/lib/libxgboost.so)
#include </local/cms/user/wadud/aNTGCmet/xgboost/include/xgboost/c_api.h>


#ifndef WGAMMAENUJETFAKESDENOMSELECTOR
#define WGAMMAENUJETFAKESDENOMSELECTOR

// Barrel-Endcap transition region
#define BETRetaMin 1.4442
#define BETRetaMax 1.566
#define HBetaMax 1.3920                     // ref. Josh H.
#define ZMASS 91.1876
#define pi7 3.1415927
#define REPORT_EVERY 100000
#define CUTFLOWSTEPS 50

const Double_t ECAL_ETA_BINS[57] = { -5., -3., -2.9, -2.7, -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -BETRetaMax, -BETRetaMin, -1.3, -1.2, -1.1, -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, BETRetaMin, BETRetaMax, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 2.9, 3., 5.};
const Double_t ECAL_ABS_ETA_BINS[13] = { 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, BETRetaMin, BETRetaMax, 1.8, 2.0, 2.2, 2.5};
const Double_t ECAL_EB_ETA_BINS[57] = { -BETRetaMin, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05, -1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, BETRetaMin};
const Double_t ECAL_FINE_ETA_BINS[85] = { -5., -3., -2.9, -2.7, -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -BETRetaMax, -BETRetaMin, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05, -1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, BETRetaMin, BETRetaMax, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 2.9, 3., 5.};
const Double_t ECAL_EB_ABS_ETA_BINS[8] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, BETRetaMin};


////////////////////////////////////////////////// Container for categories //////////////////////////////////////////////////////////////////////////////////////
struct eventType {
	TTree *                             tree = nullptr;

	// Cut efficiency tracking
	Float_t                             lastCutStep = 0.;
	TH1F *                              cutFlowCount = nullptr;
	TH1F *                              cutFlowGenWeight = nullptr;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class WGamma_ENu_JetFakesDenomSelector {
 public:
	WGamma_ENu_JetFakesDenomSelector(std::string FILELIST, std::string OUTFILE, Float_t XSECTION = -1., std::string MCPILEUPHIST = "", std::string DATAPILEUPHIST = "/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root",
	                                 std::string PFECALCLUS_EFFECTIVE_AREAS = "/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoPFClusEcalIso.txt",
	                                 std::string PFHCALCLUS_EFFECTIVE_AREAS = "/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoPFClusHcalIso.txt",
	                                 std::string TKRISO_EFFECTIVE_AREAS = "/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoTrkSumPtHollowConeDR03.txt",
	                                 std::string PFECALCLUS_PTSCALING = "/hdfs/cms/user/wadud/anTGC/analysis_data/isoPtScaling/90pc/phoPFClusEcalIso.txt",
	                                 std::string PFHCALCLUS_PTSCALING = "/hdfs/cms/user/wadud/anTGC/analysis_data/isoPtScaling/90pc/phoPFClusHcalIso.txt",
	                                 std::string BDT_PATH = "/hdfs/cms/user/wadud/anTGC/analysis_data/phoBDT/V3/aNTGC_photon_BDT_2020_12_14_20_11_56.model",
                           std::string ELE_SF_PATH_CAND="",
                           std::string ELE_SF_PATH_VETO="",
                           std::string MU_ID_SF_PATH_VETO="/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/scaleFactors/rereco/muon/RunBCDEF_SF_ID_syst.root,NUM_TightID_DEN_genTracks_pt_abseta",
                           std::string MU_ISO_SF_PATH_VETO="/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/scaleFactors/rereco/muon/RunBCDEF_SF_ID_syst.root,NUM_LooseID_DEN_genTracks_pt_abseta",
                           std::string PFPHOISO_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt",
                std::string PFNEUISO_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt",
                std::string PFCHISO_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/EGM/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt");

	~WGamma_ENu_JetFakesDenomSelector() {
		XGBoosterFree(phoBDT_h);
		std::cout << "END @ " << getCurrentTime() << std::endl;
		std::cout << "*************************************************************************************************************************************************" << std::endl;
	};

 private:
	Bool_t              isMC = false;
	Bool_t              doPUreweight = false;
	Float_t             xSec = -1.;

	effectiveAreaMap    PFHCALClusEffAreas;
	effectiveAreaMap    PFECALClusEffAreas;
	effectiveAreaMap    TkrEffAreas;

        effectiveAreaMap    PFNeuIsoEffAreas;
        effectiveAreaMap    PFPhoIsoEffAreas;
        effectiveAreaMap    PFChIsoEffAreas;

	isoPtScalingMap 	PFHCALClusPtScaling;
	isoPtScalingMap 	PFECALClusPtScaling;

	void                analyze();
	Bool_t              selectEvent();

	Short_t				nearestFinalGen(Short_t _phoIndex, Float_t _deltaRmax);
	Short_t				nearestFinalGen(Float_t _eta, Float_t _phi, Float_t _deltaRmax);
	Short_t				photonIsPrompt(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	Bool_t				photonIsFake(Short_t _phoIndex, Float_t _deltaRmax = 0.3);
	Short_t 			matchWithRecoPho(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	Short_t 			matchWithTrigPho(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	Bool_t 				phoIsSigCandidate(Short_t _phoIndex);

	Float_t 			getPhoBDTScore(Short_t iPho);

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
	TTreeReaderVectorValue<Short_t>     	_mcIndex;
	TTreeReaderVectorValue<Char_t>     		_mcPromptStatusType;
	TTreeReaderAnyValue<Float_t>        	_genPho1;
	TTreeReaderAnyValue<Float_t>         	_genPho2;

	TTreeReaderAnyValue<UShort_t>           _nPho;
	TTreeReaderVectorValue<Short_t>			_pho_gen_index;
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

	TTreeReaderAnyValue<Float_t>            _pfMET;
	TTreeReaderAnyValue<Float_t>            _pfMETPhi;
	TTreeReaderAnyValue<Float_t>            _pfMET_metSig;
	TTreeReaderAnyValue<Float_t>            _pfMET_EtSig;

	TTreeReaderAnyValue<UShort_t>           _nEle;
	TTreeReaderVectorValue<Float_t>         _eleCalibPt;
	TTreeReaderVectorValue<Float_t>         _eleEta;
	TTreeReaderVectorValue<Float_t>         _elePhi;
	TTreeReaderVectorValue<UChar_t>         _eleIDbit;
	TTreeReaderVectorValue<UShort_t>        _eleQualityBits;
	TTreeReaderVectorValue<Float_t>         _eleIDMVAIso;
	TTreeReaderVectorValue<Short_t>         _eleDirectEcalSCindex;


	TTreeReaderAnyValue<UShort_t>           _nMu;
	TTreeReaderVectorValue<Float_t>         _muEta;
	TTreeReaderVectorValue<Float_t>         _muPt;
	TTreeReaderVectorValue<Int_t>           _muIDbit;

	TTreeReaderAnyValue<UShort_t> 			_nTau;
	TTreeReaderVectorValue<Float_t>         _tauPt;
	TTreeReaderVectorValue<Float_t>			_tauEta;
	TTreeReaderVectorValue<UInt_t>			_tauIDbitsDeepTau2017v2p1;
	TTreeReaderVectorValue<std::vector<Char_t>>  _tauDMs;

	TTreeReaderAnyValue<UShort_t>           _nAK4CHSJet;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Pt;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Eta;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Phi;
	TTreeReaderVectorValue<Char_t>          _AK4CHSJet_PUFullID;
	TTreeReaderVectorValue<Char_t>          _AK4CHSJet_ID;


	TTreeReaderAnyValue<UChar_t>			_ntrgObjPho;
	TTreeReaderVectorValue<UChar_t>			_trgObjPhoBits;
	TTreeReaderVectorValue<Float_t>			_trgObjPhoPt;
	TTreeReaderVectorValue<Float_t>			_trgObjPhoEta;
	TTreeReaderVectorValue<Float_t>			_trgObjPhoPhi;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////// Output variables //////////////////////////////////////////////////////////
	Char_t              fillVars(Short_t _phoIndex, Short_t _eleIndex);

	Int_t               run_;
	Long64_t            event_;
	UShort_t            lumis_;
	Float_t 			rho_;
	UChar_t				nVtx_;

	UShort_t			metFilters_;
	UShort_t 			beamHaloSummary_;

	Float_t             genWeight_ = 1.;
	Float_t             puWeight_ = 1.;
	UChar_t				puTrue_;
	UShort_t			genStatusFlag_;
	Short_t				genStatus_;
	Float_t 			deltaRgenPho_;
	Float_t 			relDeltaPtGenPho_;
	Float_t 			deltaRPt_;
	Int_t				genPDGid_ = 0;
	Char_t				genPromptStatusType_ = 0;

	Float_t 			lhePhoPt_;
	Char_t 				Wsign_;
	Char_t 				lepSign_;

	Float_t 			deltaRtrg_;
	Float_t 			deltaPttrg_;

	Float_t 			nPhoCand_;
	Float_t 			phoPt_;
	Float_t 			phoEta_;
	Float_t 			phoPhi_;
	Float_t 			phoSeedTime_;

	UChar_t				phoGenBits_;
	UChar_t				phoQualityBits_;
	Float_t 			phoR9Full5x5_;
	Float_t 			phoS4Full5x5_;

	Float_t 			phoEmaxOESCrFull5x5_;
	Float_t 			phoE2ndOESCrFull5x5_;
	Float_t 			phoE2ndOEmaxFull5x5_;
	Float_t 			phoE1x3OESCrFull5x5_;
	Float_t 			phoE2x5OESCrFull5x5_;
	Float_t 			phoE5x5OESCrFull5x5_;

	Float_t 			phoEmaxOE3x3Full5x5_;
	Float_t 			phoE2ndOE3x3Full5x5_;
	Float_t 			pho2x2OE3x3Full5x5_;
	Float_t 			phoSieieOSipipFull5x5_;
	Float_t				phoEtaWOPhiWFull5x5_;

	Float_t 			phoSigmaIEtaIEta_;
	Float_t 			phoSigmaIPhiIPhi_;
	Float_t 			phoSigmaIEtaIPhi_;

	Float_t 			phoE2x2Full5x5_;
	Float_t 			phoE3x3Full5x5_;
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
	Float_t 			phoPFECALClusIsoCorr_;
	Float_t 			phoPFHCALClusIsoCorr_;
	Float_t 			phoTkrIsoCorr_;

	Float_t 			phoHoverE_;
	Float_t 			phoPFChIsoRaw_;
	Float_t 			phoPFPhoIsoRaw_;
	Float_t 			phoPFNeuIsoRaw_;
	Float_t 			phoPFChWorstIsoRaw_;
	Float_t 			phoIDMVA_;
	UChar_t             phoIDbit_;
	Float_t 			phoBDTpred_;
	UChar_t 			phoPFClusIDbits_;

	Bool_t 					pass95_ = 0;
	Bool_t 					pass90_ = 0;
	Bool_t 					pass80_ = 0;
	Bool_t 					pass70_ = 0;

	Float_t 			phoMIP_;

	Float_t 			phoSCet_;
	Float_t 			phoSCrawet_;
	Float_t 			phoSCeta_;
	Float_t 			phoSCphi_;
	Float_t 			phoSCEn_;
	Float_t 			phoSCRawEn_;
	Float_t 			phoEtaWidth_;
	Float_t 			phoPhiWidth_;

	UChar_t 			lepVeto_;

	Float_t 			lepPt_;
	Float_t 			lepEta_;
	Float_t 			lepPhi_;
	Float_t 			mTLepMet_;
	Int_t 				lepGenPID_;
	Float_t 			lepMetRecoilPt_;
	Float_t 			lepMetRecoilPhi_;
	Float_t 			deltaPhiRecoilPho_;
	Float_t 			minDeltaPhiRecoilJet30_;

	Float_t 			met_;
	Float_t				metPhi_;
	Float_t 			metSig_;
	Float_t 			EtSig_;

	Float_t 			mT_;
	Float_t 			deltaPhiMetPho_;
	Float_t 			phoPtOverMet_;


	Float_t 			nJet30_;
	Float_t 			Ht30_;
	Float_t 			minDeltaPhiMetJet30_;
	Float_t 			minDeltaRPhoJet30_;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////Buffer variables///////////////////////////////////////////////////////
	Short_t 			phoTrigMatch;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////// XGBoost ////////////////////////////////////////////////////////////
	DMatrixHandle 		dTest;
	BoosterHandle 		phoBDT_h;
	Bool_t 				predictBDT = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////Categories////////////////////////////////////////////////////////////
	Bool_t          	initEventTypes();
	void            	initEventType(eventType & evType, std::string typeName, std::string typeTitle);
	void            	fillEventType(eventType & evType);
	void            	registerCutFlow(eventType & evType);
	void            	registerAllCutFlow();
	eventType 			fullEB;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
WGamma_ENu_JetFakesDenomSelector::WGamma_ENu_JetFakesDenomSelector(std::string FILELIST, std::string OUTFILE, Float_t XSECTION, std::string MCPILEUPHIST, std::string DATAPILEUPHIST,
    std::string PFECALCLUS_EFFECTIVE_AREAS, std::string PFHCALCLUS_EFFECTIVE_AREAS, std::string TKRISO_EFFECTIVE_AREAS,
    std::string PFECALCLUS_PTSCALING,	std::string PFHCALCLUS_PTSCALING,
    std::string BDT_PATH, std::string ELE_SF_PATH_CAND, std::string ELE_SF_PATH_VETO,
                                       std::string MU_ID_SF_PATH_VETO, std::string MU_ISO_SF_PATH_VETO, std::string PFPHOISO_EFFECTIVE_AREAS, std::string PFNEUISO_EFFECTIVE_AREAS, std::string PFCHISO_EFFECTIVE_AREAS ) {

	std::cout << "*************************************************************************************************************************************************" << std::endl <<
	          getCurrentTime() << std::endl <<
	          "Running WGamma_ENu_JetFakesDenomSelector" << std::endl <<
	          "\n\nInput parameters:" << std::endl <<
	          "\t\tFile list = " << FILELIST << std::endl <<
	          "\t\tOutput file = " << OUTFILE << std::endl <<
	          "\t\tCross section = " << XSECTION << std::endl <<
	          "\t\tMC pileup histogram = " << MCPILEUPHIST << std::endl <<
	          "\t\tData pileup histogram = " << DATAPILEUPHIST << std::endl <<
	          "\t\tECAL PFCLuster Isolation Effective Areas = " << PFECALCLUS_EFFECTIVE_AREAS << std::endl <<
	          "\t\tHCAL PFCLuster Isolation Effective Areas = " << PFHCALCLUS_EFFECTIVE_AREAS << std::endl <<
	          "\t\tTracker Isolation Effective Areas = " << TKRISO_EFFECTIVE_AREAS << std::endl <<
	          "\t\tECAL PFCLuster Isolation pT Scaling = " << PFECALCLUS_PTSCALING << std::endl <<
	          "\t\tHCAL PFCLuster Isolation pT Scaling = " << PFHCALCLUS_PTSCALING << std::endl <<
                  "\t\tBDT model file = " << BDT_PATH <<
                  "\t\tPF PHOISO Effective Areas = "<<PFPHOISO_EFFECTIVE_AREAS<<std::endl<<
                  "\t\tPF NEUISO Effective Areas = "<<PFNEUISO_EFFECTIVE_AREAS<<std::endl<<
                  "\t\tPF CHISO Effective Areas = "<<PFCHISO_EFFECTIVE_AREAS<<std::endl;

	xSec = XSECTION;
	if (XSECTION > 0.) isMC = true;
	if (isMC && file_exists(MCPILEUPHIST) && file_exists(DATAPILEUPHIST)) doPUreweight = true;

	std::cout << "\t\tSample is simulation = " << std::boolalpha << isMC << std::endl <<
	          "\t\tDo pileup reweight = " << doPUreweight << "\n\n" << std::endl;

	if (doPUreweight) {
		std::cout << "Pileup reweighting:" << std::endl;
		puReweighter.init(MCPILEUPHIST, DATAPILEUPHIST, "hPUTruew", "pileup");
	}

	if (file_exists(PFECALCLUS_EFFECTIVE_AREAS)) {
		std::cout << "PF ECAL Cluster effective areas:" << std::endl;
		PFECALClusEffAreas.init(PFECALCLUS_EFFECTIVE_AREAS, 1, ",", 0);
	}

	if (file_exists(PFHCALCLUS_EFFECTIVE_AREAS)) {
		std::cout << "PF HCAL Cluster effective areas:" << std::endl;
		PFHCALClusEffAreas.init(PFHCALCLUS_EFFECTIVE_AREAS, 1, ",", 0);
	}

	if (file_exists(TKRISO_EFFECTIVE_AREAS)) {
		std::cout << "Tracker isolation effective areas:" << std::endl;
		TkrEffAreas.init(TKRISO_EFFECTIVE_AREAS, 1, ",", 0);
	}

	if (file_exists(PFECALCLUS_PTSCALING)) {
		std::cout << "ECAL pT scaling:" << std::endl;
		PFECALClusPtScaling.init(PFECALCLUS_PTSCALING, 0, 1, ",", 0);
	}

	if (file_exists(PFHCALCLUS_PTSCALING)) {
		std::cout << "HCAL pT scaling:" << std::endl;
		PFHCALClusPtScaling.init(PFHCALCLUS_PTSCALING, 1, 1, ",", 0);
	}

	if (file_exists(BDT_PATH)) {
		std::cout << "\nLoading BDT model from " << BDT_PATH << std::endl;
		XGBoosterCreate(NULL, 0, &phoBDT_h);
		XGBoosterSetParam(phoBDT_h, "seed", "0");
		Int_t mLdSuccess = XGBoosterLoadModel(phoBDT_h, BDT_PATH.c_str());
		if (mLdSuccess == 0) predictBDT = 1;
		else {
			std::cout << "Failed to load BDT model!" << std::endl;
		}
	}
        if(file_exists(PFPHOISO_EFFECTIVE_AREAS)) {
          std::cout<<"PF Pho effective areas:"<<std::endl;
          PFPhoIsoEffAreas.init(PFPHOISO_EFFECTIVE_AREAS, 1, "        ", 2);
        }

        if(file_exists(PFNEUISO_EFFECTIVE_AREAS)) {
          std::cout<<"PF Neu effective areas:"<<std::endl;
          PFNeuIsoEffAreas.init(PFNEUISO_EFFECTIVE_AREAS, 1, "        ", 2);
        }

        if(file_exists(PFCHISO_EFFECTIVE_AREAS)){
          std::cout<<"PF CH effective areas:"<<std::endl;
          PFChIsoEffAreas.init(PFCHISO_EFFECTIVE_AREAS, 1, "        ", 2);
        }

	std::cout << "\nCreating TChain... " << std::endl;
	if (!initNtuples(FILELIST)) exit(EXIT_FAILURE);

	outFile = new TFile(OUTFILE.c_str(), "RECREATE");

	initEventTypes();

	analyze();

	outFile->Write();
	outFile->Close();

	closeTChain(inputTree);

	std::cout << "\n\nOutput written to file\t" << OUTFILE << std::endl << "Complete!" << std::endl << getCurrentTime() << std::endl;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t				WGamma_ENu_JetFakesDenomSelector::nearestFinalGen(Float_t _eta, Float_t _phi, Float_t _deltaRmax) {

	Short_t matchedGen = -999;
	Float_t minDeltaR = 999.;

	for (UShort_t iGenP = 0; iGenP < _nMC; iGenP++) {

		if (_mcStatus[iGenP] != 1) continue;

		Float_t dRiGenReco = deltaR(_eta, _phi, _mcEta[iGenP], _mcPhi[iGenP]);
		if (dRiGenReco > _deltaRmax) continue;

		if (dRiGenReco < minDeltaR) {
			minDeltaR = dRiGenReco;
			matchedGen = iGenP;
		}
	}

	return matchedGen;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t				WGamma_ENu_JetFakesDenomSelector::nearestFinalGen(Short_t _phoIndex, Float_t _deltaRmax) {

	Short_t matchedGen = -999;
	Float_t minDeltaR = 999.;

	for (UShort_t iGenP = 0; iGenP < _nMC; iGenP++) {

		if (_mcStatus[iGenP] != 1) continue;

		Float_t dRiGenPho = deltaR(_phoEta[_phoIndex], _phoPhi[_phoIndex], _mcEta[iGenP], _mcPhi[iGenP]);
		if (dRiGenPho > _deltaRmax) continue;

		if (dRiGenPho < minDeltaR) {
			minDeltaR = dRiGenPho;
			matchedGen = iGenP;
		}
	}

	return matchedGen;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	WGamma_ENu_JetFakesDenomSelector::photonIsPrompt(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax) {

	Short_t matchedPromptGenPho = -999;
	Float_t minDeltaR = 999.;

	for (UShort_t iGenP = 0; iGenP < _nMC; iGenP++) {
		if (_mcPID[iGenP] != 22) continue;

		UShort_t iGenPStFl = _mcStatusFlag[iGenP];
		if (!getBit(iGenPStFl, 1)) continue;

		Float_t dRiGenPho = deltaR(_phoEta[_phoIndex], _phoPhi[_phoIndex], _mcEta[iGenP], _mcPhi[iGenP]);
		if (dRiGenPho > _deltaRmax) continue;

		// Float_t relDeltaPtiGenPho = std::abs(_mcPt[iGenP] - _phoCalibEt[_phoIndex])/_mcPt[iGenP];
		// if(relDeltaPtiGenPho > _relDeltaPtMax) continue;
		// if(relDeltaPtiGenPho < _relDeltaPtMin) continue;

		if (dRiGenPho < minDeltaR) {
			minDeltaR = dRiGenPho;
			matchedPromptGenPho = iGenP;
		}
	}

	return matchedPromptGenPho;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	WGamma_ENu_JetFakesDenomSelector::matchWithRecoPho(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax) {

	Short_t matchedRecoPho = -999;
	Float_t minDeltaR = 999.;

	for (UShort_t iPho = 0; iPho < _nPho; iPho++) {

		UChar_t iPhoFidReg = _phoFiducialRegion[iPho];
		if (!getBit(iPhoFidReg, 0)) continue; 			// skip if not EB (0 = EB, 1 = EE, 2 = EB-EE gap)

		// if(_phoCalibEt[iPho] < 200.) continue;

		Float_t relDeltaPtiGenPho = (_phoCalibEt[iPho] - _mcPt[_genIndex]) / _mcPt[_genIndex];
		if (relDeltaPtiGenPho > _relDeltaPtMax) continue;
		if (relDeltaPtiGenPho < _relDeltaPtMin) continue;

		Float_t dRiGenPho = deltaR(_phoEta[iPho], _phoPhi[iPho], _mcEta[_genIndex], _mcPhi[_genIndex]);
		if (dRiGenPho > _deltaRmax) continue;

		if (dRiGenPho < minDeltaR) {
			matchedRecoPho = iPho;
			minDeltaR = dRiGenPho;
		}
	}

	return matchedRecoPho;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	WGamma_ENu_JetFakesDenomSelector::matchWithTrigPho(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax) {
	Float_t minDeltaR = 999.;
	Short_t matchedTrigPho = -999;

	for (UShort_t iTrgPho = 0; iTrgPho < _ntrgObjPho; iTrgPho++) {

		Float_t relDeltaPt = (_phoCalibEt[_phoIndex] - _trgObjPhoPt[iTrgPho]) / _trgObjPhoPt[iTrgPho];
		if (relDeltaPt > _relDeltaPtMax) continue;
		if (relDeltaPt < _relDeltaPtMin) continue;

		Float_t dR = deltaR(_phoEta[_phoIndex], _phoPhi[_phoIndex], _trgObjPhoEta[iTrgPho], _trgObjPhoPhi[iTrgPho]);
		if (dR > _deltaRmax) continue;

		if (dR < minDeltaR) {
			matchedTrigPho = iTrgPho;
			minDeltaR = dR;
		}
	}

	return matchedTrigPho;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t WGamma_ENu_JetFakesDenomSelector::photonIsFake(Short_t _phoIndex, Float_t _deltaRmax) {
	// Checks for the existence of a prompt gen photon within a given deltaR cone
	Bool_t photonIsFake = 1;
	Float_t photonEta = _phoEta[_phoIndex];
	Float_t photonPhi = _phoPhi[_phoIndex];

	for (UShort_t iGenP = 0; iGenP < _nMC; iGenP++) {
		if (_mcPID[iGenP] != 22) continue;

		UShort_t iGenPStFl = _mcStatusFlag[iGenP];
		if (!getBit(iGenPStFl, 1)) continue;

		Float_t dRiGenPho = deltaR(photonEta, photonPhi, _mcEta[iGenP], _mcPhi[iGenP]);
		if (dRiGenPho > _deltaRmax) continue;

		photonIsFake = 0;

		break;
	}

	return photonIsFake;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t 			WGamma_ENu_JetFakesDenomSelector::phoIsSigCandidate(Short_t _phoIndex) {

	Float_t _phoIndexAbsSCEta = std::abs(_ecalSCeta[_phoDirectEcalSCindex[_phoIndex]]);

	if ( _phoIndexAbsSCEta > BETRetaMin) return 0;

	if (_phoCalibEt[_phoIndex] < 225.) return 0;

	if (_phoHoverE[_phoIndex] > 4.2132449148158196e-02) return 0;

	Float_t _phoIndexPFECALClusIsoCorr 		= 	_phoPFClusEcalIso[_phoIndex] - _rho * PFECALClusEffAreas.getEffectiveArea(_phoIndexAbsSCEta) - PFECALClusPtScaling.getPtScaling(_phoIndexAbsSCEta, _phoCalibEt[_phoIndex]);
	if (_phoIndexPFECALClusIsoCorr > 5.1143879813442599e+00) return 0;

	Float_t _phoIndexPFHCALClusIsoCorr 		= 	_phoPFClusHcalIso[_phoIndex] - _rho * PFHCALClusEffAreas.getEffectiveArea(_phoIndexAbsSCEta) - PFHCALClusPtScaling.getPtScaling(_phoIndexAbsSCEta, _phoCalibEt[_phoIndex]);
	if (_phoIndexPFHCALClusIsoCorr > 1.7020221456277284e+01) return 0;

	Float_t _phoIndexTkrIsoCorr 				= 	_phoTrkSumPtHollowConeDR03[_phoIndex] - _rho * TkrEffAreas.getEffectiveArea(_phoIndexAbsSCEta);
	if (_phoIndexTkrIsoCorr > 3.5001887785563284e+00) return 0;

	if (getPhoBDTScore(_phoIndex) < 4.9013550283797025e-01) return 0;

	if (getBit((_phoQualityBits[_phoIndex]), 0)) return 0;

	if (_phoMIPTotEnergy[_phoIndex] > 4.9) return 0;

	return 1;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t WGamma_ENu_JetFakesDenomSelector::selectEvent() {

	//// reset event cut flow
	fullEB.lastCutStep = 0.1;
	registerAllCutFlow();

	nPhoCand_ = 0.0001;

	ULong64_t HLTPho    = (_HLTPho);                //// 9 = HLT_Photon200_v, 11 = HLT_Photon300_NoHE_v, 19 = HLT_Photon135_PFMET100_v, no trigger on 2016 MC
	if (!getBit(HLTPho, 9))      return 0;          //// 200 GeV photon trigger
	registerAllCutFlow();

	UShort_t metFilters = (_metFilters);            //// MET filters

	if (getBit(metFilters, 0))   return 0;          //// Flag_goodVertices
	registerAllCutFlow();

	if (getBit(metFilters, 1))   return 0;          //// Flag_globalSuperTightHalo2016Filter
	registerAllCutFlow();

	if (getBit(metFilters, 2))   return 0;          //// Flag_HBHENoiseFilter
	registerAllCutFlow();

	if (getBit(metFilters, 3))   return 0;          //// Flag_HBHENoiseIsoFilter
	registerAllCutFlow();

	if (getBit(metFilters, 4))   return 0;          //// Flag_EcalDeadCellTriggerPrimitiveFilter
	registerAllCutFlow();

	if (getBit(metFilters, 5))   return 0;          //// Flag_BadPFMuonFilter
	registerAllCutFlow();

	if (getBit(metFilters, 9))   return 0;          //// Updated ecalBadCalibFilter
	registerAllCutFlow();

	Short_t highetsPtIndex 	= -9999;
	Bool_t 	passedFidCut 	= 0;
	Float_t highestPhoPt 	= -999.;
	Bool_t 	passedPtCut 	= 0;
	Bool_t 	passedHoE = 0;
	Bool_t 	passedIso = 0;
	Bool_t 	passedBDT = 0;
	Bool_t passedMipCut = 0;
	Bool_t passedPixelCut = 0;

        Short_t iPhoton = -999;
        Float_t highestPtPho = -999;
        Int_t selPho = 0;
        Bool_t foundPho = false;

        Short_t iElectron = -999;
        Int_t TightEle = 0;
        Bool_t foundEle = false;

	for (UShort_t iPho = 0; iPho < _nPho; iPho++) {

           Float_t rho                    = _rho;
           Short_t phoSCindex      = _phoDirectEcalSCindex[iPho];
           Float_t phoSCeta        = _ecalSCeta[phoSCindex];
           Float_t phoAbsSCEta     = std::abs(phoSCeta);

              if (phoIsSigCandidate(iPho))       nPhoCand_ += 1.;
              
              if(_phoHoverE[iPho] > 0.0590) continue;
              if(_phoSigmaIEtaIEtaFull5x5[iPho] < 0.03) continue;

              float iPhoPFNeuIsoCorr =    _phoPFNeuIso[iPho] - PFNeuIsoEffAreas.getEffectiveArea(phoAbsSCEta) * rho - 0.0117*_phoCalibEt[iPho] - 2.3e-05*std::pow(_phoCalibEt[iPho],2.);
              bool passLooseNeuIso = (iPhoPFNeuIsoCorr < 19.722);

              float iPhoPFPhoIsoCorr = _phoPFPhoIso[iPho] - PFPhoIsoEffAreas.getEffectiveArea(phoAbsSCEta) * rho - 0.0037*_phoCalibEt[iPho];
              bool passLoosePhoIso = (iPhoPFPhoIsoCorr < 4.162);

              float iPhoPFChIsoCorr = _phoPFChIso[iPho] - PFChIsoEffAreas.getEffectiveArea(phoAbsSCEta) * rho;
              bool passLooseChiIso = (iPhoPFChIsoCorr < 2.089);

              if(passLoosePhoIso && passLooseNeuIso && passLooseChiIso) continue;
              if((iPhoPFNeuIsoCorr > 5.*19.722) || (iPhoPFPhoIsoCorr > 5.*4.162)  || (iPhoPFChIsoCorr > 5.*2.089)) continue;

              UChar_t phoPixVeto = _phoQualityBits[iPho];
              UChar_t phoCutBasedID = _phoIDbit[iPho];

              if (phoPixVeto == 0 ) continue;  // if phoPixVeto == 0 -> photon hasPixelSeed
              if (_phoCalibEt[iPho] < 220) continue;
              if (std::abs(_phoEta[iPho]) < 1.566 || std::abs(_phoEta[iPho]) > 2.5) continue;

              iPhoton = iPho;
              selPho++;
	}

        if (selPho != 1) return 0;
        if (iPhoton >=0) foundPho = true;
        if (!foundPho) return 0;

        for(UShort_t iEle=0; iEle<_nEle; iEle++) {

                     if (_eleCalibPt[iEle] < 30) continue;
                     if (std::abs(_eleEta[iEle]) > 2.5) continue;
                     if (!getBit((_eleIDbit[iEle]), 3)) continue;  // ele has tight ID 
                     iElectron = iEle;
                     TightEle++;
        }

        if (TightEle != 1) return 0; // veto evts with 0 or more than 1 tight ele              
        if (iElectron >=0) foundEle = true;
        if (!foundEle) return 0;

	if (passedFidCut) registerAllCutFlow();
	if (passedPtCut) registerAllCutFlow();
	if (passedHoE) registerAllCutFlow();
	if (passedIso) registerAllCutFlow();
	if (passedBDT) registerAllCutFlow();
	if (passedMipCut) registerAllCutFlow();
	if (passedPixelCut) registerAllCutFlow();

	//if (highetsPtIndex < 0) return 0;

	// Data e-->pho FR
	// Median fit eqn : 0.0470368 + 0.00999998|x| + -0.0225408x*x
	// HighSys fit eqn : 0.0563197 + 0.00208314|x| + -0.0168909x*x
	// LowSys fit eqn : 0.0386988 + 0.0142861|x| + -0.0259301x*x

	// Must have 1 medium electron
	/*Short_t highetsPtEleIndex 	= -999;
	Float_t highestElePt 		= -999.;
	Short_t nEle 				=	0;
	for (Int_t iEle = 0; iEle < (_nEle); iEle++) {
		if (_eleCalibPt[iEle] < 30.) continue;
		if (std::abs(_eleEta[iEle]) > 2.5) continue;
		if (!getBit((_eleIDbit[iEle]), 3)) continue; 		// 3=tight ID
		nEle++;
		if (_eleCalibPt[iEle] > highestElePt) {
			highestElePt = _eleCalibPt[iEle];
			highetsPtEleIndex = iEle;
		}
	}

	if (highetsPtEleIndex < 0) return 0;
	registerAllCutFlow();

	if (nEle > 1) return 0;*/

	fillVars(iPhoton, iElectron);

	fillEventType(fullEB);


	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t WGamma_ENu_JetFakesDenomSelector::initEventTypes() {

	pileupPreweight.SetDirectory(outFile->GetDirectory(""));
	pileupPostweight.SetDirectory(outFile->GetDirectory(""));
	rhoPreweight.SetDirectory(outFile->GetDirectory(""));
	rhoPostweight.SetDirectory(outFile->GetDirectory(""));
	nvtxPreweight.SetDirectory(outFile->GetDirectory(""));
	nvtxPostweight.SetDirectory(outFile->GetDirectory(""));

	initEventType(fullEB, "fullEB", "ECAL Barrel");
	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void WGamma_ENu_JetFakesDenomSelector::analyze() {
	std::cout << "--------------------------------------------------------------------------------------------------" << std::endl <<
	          getCurrentTime() << std::endl <<
	          "Analyzing events.." << std::endl;
	ULong64_t current_entry = 0;
	while (inputTTreeReader.Next()) {

		if (current_entry % REPORT_EVERY == 0) {
			std::cout << "\t" << getCurrentTime() << "\tAnalyzing entry\t" << current_entry <<
			          ",\t\tevent\t" << (_event) << "\t\tFile " << inputTree->GetCurrentFile()->GetName() << std::endl;
		}

		// if(current_entry > 1000) break;
		rhoPreweight.Fill(_rho, genWeight_);
		nvtxPreweight.Fill(_nVtx, genWeight_);

		if (isMC) {
			if (doPUreweight) {
				genWeight_ = _genWeight;
				puWeight_ = puReweighter.weight(_puTrue);
				Float_t genPUweight = genWeight_ * puWeight_;
				pileupPostweight.Fill(_puTrue, genPUweight);
				rhoPostweight.Fill(_rho, genPUweight);
				nvtxPostweight.Fill(_nVtx, genPUweight);
			} else {
				genWeight_ = _genWeight;
				puWeight_ = 1.;
			}
			pileupPreweight.Fill(_puTrue, genWeight_);
		}

		selectEvent();

		current_entry++;
	};

	std::cout << "Done analyzing!" << std::endl <<
	          getCurrentTime() << std::endl <<
	          "--------------------------------------------------------------------------------------------------" << std::endl;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t 	WGamma_ENu_JetFakesDenomSelector::getPhoBDTScore(Short_t iPho) {

	Short_t phoSCindex 		= _phoDirectEcalSCindex[iPho];

	std::vector<Float_t> feats{_phoE2x2Full5x5[iPho] / (_phoR9Full5x5[iPho] * _ecalSCRawEn[phoSCindex]),
	                           _phoE1x3Full5x5[iPho] / _ecalSCRawEn[phoSCindex],
	                           _phoE2ndFull5x5[iPho] / _ecalSCRawEn[phoSCindex],
	                           _phoE2x5Full5x5[iPho] / _ecalSCRawEn[phoSCindex],
	                           _phoMaxEnergyXtal[iPho] / _ecalSCRawEn[phoSCindex],
	                           _ecalSCetaWidth[phoSCindex] / _ecalSCphiWidth[phoSCindex],
	                           _ecalSCetaWidth[phoSCindex],
	                           _ecalSCphiWidth[phoSCindex],
	                           _phoR9Full5x5[iPho],
	                           _phoE2x2Full5x5[iPho] / _ecalSCRawEn[phoSCindex],
	                           _phoSigmaIEtaIEtaFull5x5[iPho] / _phoSigmaIPhiIPhiFull5x5[iPho],
	                           _phoSigmaIEtaIEtaFull5x5[iPho],
	                           _phoSigmaIEtaIPhiFull5x5[iPho],
	                           _phoSigmaIPhiIPhiFull5x5[iPho]};

	XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
	bst_ulong out_len;
	const float *prediction;
	XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
	assert(out_len == 1);
	Float_t iPhoBDTpred = prediction[0];
	XGDMatrixFree(dTest);

	return iPhoBDTpred;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Char_t WGamma_ENu_JetFakesDenomSelector::fillVars(Short_t _phoIndex, Short_t _eleIndex) {


	run_ 					= _run;
	event_ 					= _event;
	lumis_ 					= _lumis;
	rho_ 					= _rho;
	nVtx_					= _nVtx;

	metFilters_ 			= _metFilters;
	beamHaloSummary_		= _beamHaloSummary;

	Short_t phoSCindex 		= _phoDirectEcalSCindex[_phoIndex];
	phoSCeta_ 				= _ecalSCeta[phoSCindex];
	Float_t phoAbsSCEta 	= std::abs(phoSCeta_);

	phoPt_ 					= _phoCalibEt[_phoIndex];
	phoEta_ 				= _phoEta[_phoIndex];
	phoPhi_ 				= _phoPhi[_phoIndex];
	phoSeedTime_ 			= _phoSeedTime[_phoIndex];

	phoQualityBits_			= _phoQualityBits[_phoIndex];
	phoR9Full5x5_ 			= _phoR9Full5x5[_phoIndex];
	phoS4Full5x5_			= _phoE2x2Full5x5[_phoIndex] / _ecalSCRawEn[phoSCindex];
	phoEmaxOESCrFull5x5_ 	= _phoMaxEnergyXtal[_phoIndex] / _ecalSCRawEn[phoSCindex];
	phoE2ndOESCrFull5x5_ 	= _phoE2ndFull5x5[_phoIndex] / _ecalSCRawEn[phoSCindex];
	phoE2ndOEmaxFull5x5_ 	= _phoE2ndFull5x5[_phoIndex] / _phoMaxEnergyXtal[_phoIndex];
	phoE1x3OESCrFull5x5_ 	= _phoE1x3Full5x5[_phoIndex] / _ecalSCRawEn[phoSCindex];
	phoE2x5OESCrFull5x5_ 	= _phoE2x5Full5x5[_phoIndex] / _ecalSCRawEn[phoSCindex];
	phoE5x5OESCrFull5x5_ 	= _phoE5x5Full5x5[_phoIndex] / _ecalSCRawEn[phoSCindex];


	phoSigmaIEtaIEta_ 		= _phoSigmaIEtaIEtaFull5x5[_phoIndex];
	phoSigmaIEtaIPhi_ 		= _phoSigmaIEtaIPhiFull5x5[_phoIndex];
	phoSigmaIPhiIPhi_ 		= _phoSigmaIPhiIPhiFull5x5[_phoIndex];

	phoE2x2Full5x5_			= _phoE2x2Full5x5[_phoIndex];
	phoE3x3Full5x5_			= phoR9Full5x5_ * _ecalSCRawEn[phoSCindex];
	phoE5x5Full5x5_			= _phoE5x5Full5x5[_phoIndex];
	phoMaxEnergyXtal_		= _phoMaxEnergyXtal[_phoIndex];
	phoE2ndFull5x5_			= _phoE2ndFull5x5[_phoIndex];
	phoE1x3Full5x5_			= _phoE1x3Full5x5[_phoIndex];
	phoE1x5Full5x5_			= _phoE1x5Full5x5[_phoIndex];
	phoE2x5Full5x5_			= _phoE2x2Full5x5[_phoIndex];

	phoEmaxOE3x3Full5x5_ 	= phoMaxEnergyXtal_ / phoE3x3Full5x5_;
	phoE2ndOE3x3Full5x5_ 	= phoE2ndFull5x5_ / phoE3x3Full5x5_;
	pho2x2OE3x3Full5x5_		= _phoE2x2Full5x5[_phoIndex] / (_phoR9Full5x5[_phoIndex] * _ecalSCRawEn[phoSCindex]);
	phoSieieOSipipFull5x5_	= _phoSigmaIEtaIEtaFull5x5[_phoIndex] / _phoSigmaIPhiIPhiFull5x5[_phoIndex];

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
	phoTkrIsoCorr_ 				= 	phoTrkSumPtHollowConeDR03_ - rho_ * TkrEffAreas.getEffectiveArea(phoAbsSCEta);

	phoHoverE_ 				= _phoHoverE[_phoIndex];
	phoPFChIsoRaw_ 			= _phoPFChIso[_phoIndex];
	phoPFPhoIsoRaw_ 		= _phoPFPhoIso[_phoIndex];
	phoPFNeuIsoRaw_ 		= _phoPFNeuIso[_phoIndex];
	phoPFChWorstIsoRaw_ 	= _phoPFChWorstIso[_phoIndex];
	phoIDMVA_ 				= _phoIDMVA[_phoIndex];
	phoIDbit_ 				= _phoIDbit[_phoIndex];
	phoMIP_ 				= _phoMIPTotEnergy[_phoIndex];

	phoSCet_ 				= (_ecalSCEn[phoSCindex]) / std::cosh(phoSCeta_);
	phoSCrawet_				= (_ecalSCRawEn[phoSCindex]) / std::cosh(phoSCeta_);
	phoSCphi_ 				= _ecalSCphi[phoSCindex];
	phoSCEn_ 				= _ecalSCEn[phoSCindex];
	phoSCRawEn_				= _ecalSCRawEn[phoSCindex];
	phoEtaWidth_ 			= _ecalSCetaWidth[phoSCindex];
	phoPhiWidth_ 			= _ecalSCphiWidth[phoSCindex];

	phoEtaWOPhiWFull5x5_	= _ecalSCetaWidth[phoSCindex] / _ecalSCphiWidth[phoSCindex];

	if (isMC) {
		puTrue_					= _puTrue;

		phoGenBits_ 			= 0;

		Short_t iPhoGenIndex = nearestFinalGen(_phoIndex, 0.1);

		if (iPhoGenIndex > -1) {
			setBit(phoGenBits_, 0, 1);
			genStatusFlag_ 			= 	_mcStatusFlag[iPhoGenIndex];
			genStatus_ 				=	_mcStatus[iPhoGenIndex];
			deltaRgenPho_ 			=	deltaR(_mcEta[iPhoGenIndex], _mcPhi[iPhoGenIndex], phoEta_, phoPhi_);
			relDeltaPtGenPho_ 		=	(_phoCalibEt[_phoIndex] - _mcPt[iPhoGenIndex]) / _mcPt[iPhoGenIndex];
			deltaRPt_ 				= std::sqrt(deltaRgenPho_ * deltaRgenPho_ + relDeltaPtGenPho_ * relDeltaPtGenPho_);
			genPDGid_ 				= 	_mcPID[iPhoGenIndex];
			genPromptStatusType_ 	=	_mcPromptStatusType[iPhoGenIndex];
		} else {
			genStatusFlag_ 			= 	9999;
			genStatus_ 				= 	-9999;
			deltaRgenPho_ 			=	-9999;
			relDeltaPtGenPho_ 		=	-9999;
			deltaRPt_ 				= 	-9999;
			genPDGid_				= -9999;
			genPromptStatusType_ 	=	-99;
		}

		setBit(phoGenBits_, 1, 1);
	}

	if (isMC) {

		Wsign_ = 0;
		lepSign_ = 0;

		Float_t highestLepPt = -999;
		Short_t highetsPtLepIndex = -999;
		for (UShort_t iGenP = 0; iGenP < _nMC; iGenP++) {
			if (abs(_mcPID[iGenP]) == 24) Wsign_ = _mcPID[iGenP];

			if (_mcPt[iGenP] > 	highestLepPt && (std::abs(_mcPID[iGenP]) == 11 || std::abs(_mcPID[iGenP]) == 13 || std::abs(_mcPID[iGenP]) == 15)) {
				highetsPtLepIndex = iGenP;
				highestLepPt = _mcPt[iGenP];
			}

		}

		if (highetsPtLepIndex > -1) lepSign_ = _mcPID[highetsPtLepIndex];

		lhePhoPt_ = std::max(_genPho1, _genPho2);

	}

	phoTrigMatch = matchWithTrigPho(_phoIndex, 0.3, -10., 10.);
	if (phoTrigMatch > -1) {
		deltaRtrg_ = deltaR(phoEta_, phoPhi_, _trgObjPhoEta[phoTrigMatch], _trgObjPhoPhi[phoTrigMatch]);
		deltaPttrg_ = (_trgObjPhoPt[phoTrigMatch] - phoPt_) / _trgObjPhoPt[phoTrigMatch];
	} else {
		deltaRtrg_ = - 999;
		deltaPttrg_ = - 999;
	}

	phoPFClusIDbits_ 		= 0;

	/*if (predictBDT) {
		phoBDTpred_ = getPhoBDTScore(_phoIndex);
		pass95_ = (phoBDTpred_ >= 4.9013550283797025e-01) && (phoHoverE_ < 4.2132449148158196e-02) && (phoPFECALClusIsoCorr_ < 5.1143879813442599e+00) && (phoPFHCALClusIsoCorr_ < 1.7020221456277284e+01) && (phoTkrIsoCorr_ < 3.5001887785563284e+00);
		if (pass95_) setBit(phoPFClusIDbits_, 3, 1);
	}*/

	lepVeto_ = 0;
	for (Int_t i = 0; i < _nEle; i++) {
		if (i == _eleIndex) continue;
		if ((_eleCalibPt[i]) < 10.) continue;
		if (std::abs(_eleEta[i]) > 2.5) continue;

		// if (!getBit((_eleIDbit[i]), 0)) continue;// veto ID
		if (!getBit((_eleIDbit[i]), 1)) continue; //loose ID
		setBit(lepVeto_, 0, 1);
		break;
	}

	for (Int_t i = 0; i < (_nMu); i++) {
		if ((_muPt[i]) < 10.) continue;
		if (std::abs(_muEta[i]) > 2.4) continue;
		// if (!getBit((_muIDbit[i]), 0) || !getBit((_muIDbit[i]), 11)) continue;// 0=loose ID,11=loose tracker iso ID
		if (!getBit((_muIDbit[i]), 0) || !getBit((_muIDbit[i]), 7)) continue;// 0=loose ID,7=loose Iso
		setBit(lepVeto_, 1, 1);
		break;
	}

	for (Int_t i = 0; i < (_nTau); i++) {
		if ((_tauPt[i]) < 20.) continue;
		if (std::abs(_tauEta[i]) > 2.3) continue;
		if ((_tauDMs[i][0] == 5) || (_tauDMs[i][0] == 6)) continue; //veto DM 5 & 6 with decayModeFindingNEW
		UInt_t tmpTauIDbit = _tauIDbitsDeepTau2017v2p1[i];
		// 3 = byLooseDeepTau2017v2p1VSjet, 11 = byLooseDeepTau2017v2p1VSe, 17 = byLooseDeepTau2017v2p1VSmu
		if (!getBit(tmpTauIDbit, 3) || !getBit((tmpTauIDbit), 11) || !getBit((tmpTauIDbit), 17)) continue;
		setBit(lepVeto_, 2, 1);
		break;
	}

	lepPt_			=	_eleCalibPt[_eleIndex];
	lepEta_ 		=	_eleEta[_eleIndex];
	lepPhi_			=	_elePhi[_eleIndex];
	mTLepMet_			=	std::sqrt(2. * _eleCalibPt[_eleIndex] * _pfMET * (1. - std::cos(deltaPhi(_elePhi[_eleIndex], _pfMETPhi))));

	lepGenPID_ 		=	-999;
	if (isMC) {
		Short_t 	iGenLep 		= nearestFinalGen(_eleEta[_eleIndex], _eleEta[_eleIndex], 0.1);
		if (iGenLep > -1) lepGenPID_ 					=  _mcPID[iGenLep];
	}

	Float_t 		lepMetRecoilPx 			=	_pfMET * std::cos(_pfMETPhi) + _eleCalibPt[_eleIndex] * std::cos(_elePhi[_eleIndex]);
	Float_t 		lepMetRecoilPy 			=	_pfMET * std::sin(_pfMETPhi) + _eleCalibPt[_eleIndex] * std::sin(_elePhi[_eleIndex]);

	lepMetRecoilPt_ 						=	std::sqrt(lepMetRecoilPx * lepMetRecoilPx + lepMetRecoilPy * lepMetRecoilPy);
	lepMetRecoilPhi_ 						=	atan2(lepMetRecoilPy, lepMetRecoilPx);
	deltaPhiRecoilPho_ 						=	deltaPhi(_phoPhi[_phoIndex], lepMetRecoilPhi_);

	met_ 			= 	_pfMET;
	metPhi_ 		= 	_pfMETPhi;
	metSig_			=	_pfMET_metSig;
	EtSig_ 			=	_pfMET_EtSig;
	deltaPhiMetPho_	=	deltaPhi(metPhi_, phoPhi_);

	mT_ 			=	std::sqrt(2. * phoPt_ * met_ * (1. - std::cos(deltaPhiMetPho_)));
	phoPtOverMet_	= 	phoPt_ / met_;

	nJet30_ 		=	0.001;
	Ht30_ 			=	0.;

	std::vector<std::pair<UShort_t, Float_t>> jetsOrdered;
	for (UShort_t iJet = 0; iJet < _nAK4CHSJet; iJet++) {

		if (_AK4CHSJet_Pt[iJet] < 30.) continue;

		if (std::abs(_AK4CHSJet_Eta[iJet]) > 5.) continue;

		Char_t iJetPUID = _AK4CHSJet_PUFullID[iJet];
		if (!getBit(iJetPUID, 0) && _AK4CHSJet_Pt[iJet] < 50) continue; 			// loose PU ID

		Char_t iJetID = _AK4CHSJet_ID[iJet];
		if (!getBit(iJetID, 1)) continue; 			// tight ID

		nJet30_++;
		Ht30_ += _AK4CHSJet_Pt[iJet];

		if (deltaR(phoEta_, phoPhi_,_AK4CHSJet_Eta[iJet], _AK4CHSJet_Phi[iJet]) < 0.4) continue;
		if (deltaR(lepEta_, lepPhi_, _AK4CHSJet_Eta[iJet], _AK4CHSJet_Phi[iJet]) < 0.4) continue;

		jetsOrdered.push_back(std::make_pair(iJet, _AK4CHSJet_Pt[iJet]));

	};

	std::sort(jetsOrdered.rbegin(), jetsOrdered.rend());

	minDeltaPhiMetJet30_	=	99999999.;
	minDeltaRPhoJet30_		=	99999999.;
	minDeltaPhiRecoilJet30_ =	99999999.;

	for (UShort_t iJet = 0; iJet < jetsOrdered.size(); iJet++) {

		Float_t iJetDeltaPhiMETJet 	= deltaPhi(_pfMETPhi, _AK4CHSJet_Phi[jetsOrdered[iJet].first]);
		Float_t iJetDeltaRPhoJet 	= deltaR(phoEta_, phoPhi_, _AK4CHSJet_Eta[jetsOrdered[iJet].first], _AK4CHSJet_Phi[jetsOrdered[iJet].first]);
		Float_t iJetDeltaPhiRecoilJet 	= deltaPhi(lepMetRecoilPhi_, _AK4CHSJet_Phi[jetsOrdered[iJet].first]);

		if (iJetDeltaPhiMETJet < minDeltaPhiMetJet30_) {
			minDeltaPhiMetJet30_ = iJetDeltaPhiMETJet;
		}

		if (iJetDeltaRPhoJet < minDeltaRPhoJet30_) {
			minDeltaRPhoJet30_ = iJetDeltaRPhoJet;
		}

		if (iJetDeltaPhiRecoilJet < minDeltaPhiRecoilJet30_) {
			minDeltaPhiRecoilJet30_ = iJetDeltaPhiRecoilJet;
		}
	}

	return 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t WGamma_ENu_JetFakesDenomSelector::initNtuples(std::string FILELIST) {

	inputTree = openTChain(FILELIST, "ggNtuplizer/EventTree");
	if (!inputTree) return 0;
	inputTTreeReader.SetTree(inputTree);

	std::cout << "**************************************************************************************************************************************************************" << std::endl <<
	          "Initializing branches in input ntuples..." << std::endl;

	_run.set(inputTTreeReader, "run");
	_event.set(inputTTreeReader, "event");
	_lumis.set(inputTTreeReader, "lumis");
	_beamHaloSummary.set(inputTTreeReader, "beamHaloSummary");
	_metFilters.set(inputTTreeReader, "metFilters");
	_nVtx.set(inputTTreeReader, "nVtx");
	_rho.set(inputTTreeReader, "rho");
	_HLTPho.set(inputTTreeReader, "HLTPho");
	if (isMC) {
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
		_mcPromptStatusType.set(inputTTreeReader, "mcPromptStatusType");
		_pho_gen_index.set(inputTTreeReader, "pho_gen_index");
		_genPho1.set(inputTTreeReader, "genPho1");
		_genPho2.set(inputTTreeReader, "genPho2");
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

	_metFilters.set(inputTTreeReader, "metFilters");
	_pfMET.set(inputTTreeReader, "pfMET");
	_pfMETPhi.set(inputTTreeReader, "pfMETPhi");
	_pfMET_metSig.set(inputTTreeReader, "pfMET_metSig");
	_pfMET_EtSig.set(inputTTreeReader, "pfMET_EtSig");

	_nEle.set(inputTTreeReader, "nEle");
	_eleCalibPt.set(inputTTreeReader, "eleCalibPt");
	_eleEta.set(inputTTreeReader, "eleEta");
	_elePhi.set(inputTTreeReader, "elePhi");
	_eleIDbit.set(inputTTreeReader, "eleIDbit");
	_eleQualityBits.set(inputTTreeReader, "eleQualityBits");
	_eleIDMVAIso.set(inputTTreeReader, "eleIDMVAIso");
	_eleDirectEcalSCindex.set(inputTTreeReader, "eleDirectEcalSCindex");

	_nMu.set(inputTTreeReader, "nMu");
	_muPt.set(inputTTreeReader, "muPt");
	_muEta.set(inputTTreeReader, "muEta");
	_muIDbit.set(inputTTreeReader, "muIDbit");

	_nTau.set(inputTTreeReader, "nTau");
	_tauPt.set(inputTTreeReader, "tauPt");
	_tauEta.set(inputTTreeReader, "tauEta");
	_tauIDbitsDeepTau2017v2p1.set(inputTTreeReader, "tauIDbitsDeepTau2017v2p1");
	_tauDMs.set(inputTTreeReader, "tauDMs");

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

	std::cout << "Branches initialized!" << std::endl <<
	          "**************************************************************************************************************************************************************" << std::endl;
	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void WGamma_ENu_JetFakesDenomSelector::initEventType(eventType & evType, std::string typeName, std::string typeTitle) {

	mkTFileDir(outFile, typeName);

	evType.cutFlowCount = new TH1F((typeName + "cutFlowCount").c_str(), "Cut Flow (Unweighted)", CUTFLOWSTEPS, 0., (Float_t)CUTFLOWSTEPS);
	evType.cutFlowGenWeight = new TH1F((typeName + "cutFlowGenWeight").c_str(), "Cut Flow (Gen Weighted)", CUTFLOWSTEPS, 0., (Float_t)CUTFLOWSTEPS);
	evType.cutFlowCount->SetDirectory(outFile->GetDirectory(typeName.c_str()));
	evType.cutFlowGenWeight->SetDirectory(outFile->GetDirectory(typeName.c_str()));

	//////// initialize tree
	evType.tree = new TTree((typeName + "Tree").c_str(), typeTitle.c_str());
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
	evType.tree->Branch("mcPromptStatusType", &genPromptStatusType_);

	evType.tree->Branch("lhePhoPt", &lhePhoPt_);
	evType.tree->Branch("Wsign", &Wsign_);
	evType.tree->Branch("lepSign", &lepSign_);

	evType.tree->Branch("deltaRtrg", &deltaRtrg_);
	evType.tree->Branch("deltaPttrg", &deltaPttrg_);


	evType.tree->Branch("nPhoCand", &nPhoCand_);
	evType.tree->Branch("phoPt", &phoPt_);
	evType.tree->Branch("phoEta", &phoEta_);
	evType.tree->Branch("phoPhi", &phoPhi_);
	evType.tree->Branch("phoSeedTime", &phoSeedTime_);

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

	evType.tree->Branch("pass95", &pass95_);
	evType.tree->Branch("pass90", &pass90_);
	evType.tree->Branch("pass80", &pass80_);
	evType.tree->Branch("pass70", &pass70_);

	evType.tree->Branch("phoIDbit", &phoIDbit_);
	evType.tree->Branch("phoMIP", &phoMIP_);

	evType.tree->Branch("phoSCet", &phoSCet_);
	evType.tree->Branch("phoSCrawet", &phoSCrawet_);
	evType.tree->Branch("phoSCeta", &phoSCeta_);
	evType.tree->Branch("phoSCphi", &phoSCphi_);
	evType.tree->Branch("phoSCEn", &phoSCEn_);
	evType.tree->Branch("phoSCRawEn", &phoSCRawEn_);

	evType.tree->Branch("lepVeto", &lepVeto_);

	evType.tree->Branch("lepPt", &lepPt_);
	evType.tree->Branch("lepEta", &lepEta_);
	evType.tree->Branch("lepPhi", &lepPhi_);
	evType.tree->Branch("mTLepMet", &mTLepMet_);
	evType.tree->Branch("lepGenPID", &lepGenPID_);
	evType.tree->Branch("lepMetRecoilPt", &lepMetRecoilPt_);
	evType.tree->Branch("lepMetRecoilPhi", &lepMetRecoilPhi_);
	evType.tree->Branch("deltaPhiRecoilPho", &deltaPhiRecoilPho_);
	evType.tree->Branch("minDeltaPhiRecoilJet30", &minDeltaPhiRecoilJet30_);

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

	std::cout << "Created output tree:\t" << typeName << "\t" << typeTitle << std::endl << std::endl;
	// evType.tree->Print();
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void WGamma_ENu_JetFakesDenomSelector::fillEventType(eventType & evType) {
	evType.tree->Fill();

	registerCutFlow(evType);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void WGamma_ENu_JetFakesDenomSelector::registerAllCutFlow() {
	registerCutFlow(fullEB);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reset lastCutStep before each event
void WGamma_ENu_JetFakesDenomSelector::registerCutFlow(eventType & evType) {
	evType.cutFlowCount->Fill(evType.lastCutStep);
	evType.cutFlowGenWeight->Fill(evType.lastCutStep, genWeight_);
	evType.lastCutStep = evType.lastCutStep + 1.;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

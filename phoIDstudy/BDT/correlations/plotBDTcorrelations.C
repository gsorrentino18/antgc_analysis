#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TArrow.h"
#include "TPaveText.h"
#include "TLine.h"

std::string 													optionsFile 			=	"plotBDTcorrelationsOptions.txt";
parseOptions 													options;
TH2F*															get2DHistogram(std::string _bdtVar, std::string _isoVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax, effectiveAreaMap* _EAmap, isoPtScalingMap* _ptSmap);
void 															getCorrelation(std::string _bdtVar, std::string _isoVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax, ofstream* _isoCorrelationOutFilePath, effectiveAreaMap* _EAmap, isoPtScalingMap* _ptSmap);
TH2F*															get2DHistogramCombIso(std::string _bdtVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax);
void 															getCombIsoCorrelation(std::string _bdtVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax, ofstream* _isoCorrelationOutFilePath);
void 															initialize();
std::map<std::string,std::vector<std::string>>					varPlotInfo;
Bool_t 															truncateCorr;
Bool_t 															plotSignal;
Bool_t 															useSaved;

TH2F*															get2DHistogram(std::string _bdtVar, std::string _isoVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax, effectiveAreaMap* _EAmap, isoPtScalingMap* _ptSmap){

	std::string 								histName 									=	_bdtVar + "_vs_" + _isoVar + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p") + findAndReplaceAll(removeTrailingZeros(_pTmin) + "to" + removeTrailingZeros(_pTmax), ".", "p");
	std::string 								yTitle 										=	varPlotInfo[_isoVar][0];
	yTitle 																					= 	(_EAmap ? (findAndReplaceAll(yTitle, "(", "- #rhoA_{eff} (")) : yTitle); 
	yTitle 																					= 	(_ptSmap ? (findAndReplaceAll(yTitle, "(", "- Scale(p_{T}) (")) : yTitle); 
	std::string 								histTitle 									=	";"	+ yTitle + ";" +  varPlotInfo[_bdtVar][0];
	TH2F* 										histBDTvsIso								=	new TH2F(histName.c_str(), histTitle.c_str(), std::stoi(varPlotInfo[_isoVar][1]), std::stof(varPlotInfo[_isoVar][2]), std::stof(varPlotInfo[_isoVar][3]),
		std::stoi(varPlotInfo[_bdtVar][1]), std::stof(varPlotInfo[_bdtVar][2]), std::stof(varPlotInfo[_bdtVar][3]));
	
	histBDTvsIso->SetContour(options.getInt("2dNcontours"));
	histBDTvsIso->Sumw2();
	
	TChain*										fTree 										=	openTChain((std::vector<std::string>){options.get("featsFile")}, options.get("featsTree"));
	TChain*										bTree 										=	openTChain((std::vector<std::string>){options.get("bdtFile")}, options.get("bdtTree"));
	fTree->AddFriend(bTree);
	TTreeReader                             	inputTTreeReader(fTree);
	TTreeReaderAnyValue<Float_t>				xSecW_											(inputTTreeReader, "xSecW");
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Float_t>				phoPt_											(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Bool_t>					isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Bool_t>					isTrain_										(inputTTreeReader, "isTrain");
	TTreeReaderAnyValue<UChar_t>				phoQualityBits_									(inputTTreeReader, "phoQualityBits");
	TTreeReaderAnyValue<Float_t>				rho_;
	if(_EAmap)									rho_.											set(inputTTreeReader, "rho");
	TTreeReaderAnyValue<Float_t>				isolation_										(inputTTreeReader, _isoVar);
	TTreeReaderAnyValue<Double_t>				bdtVar_											(inputTTreeReader, _bdtVar);

	while(inputTTreeReader.Next()){
		if(isTrain_) continue;
		if(plotSignal != isSignal_) continue;
		if(phoPt_ < _pTmin) continue;
		if(phoPt_ >= _pTmax) continue;
		Float_t 								absSCEta									=	std::abs(phoSCeta_);
		if (absSCEta <= _etaMin) continue;
		if (absSCEta > _etaMax) continue;
		UChar_t 								iPhoQalBits 								=	phoQualityBits_;
		if(getBit(iPhoQalBits,0)) continue;
		Float_t 								corrections 								=	(_EAmap ? rho_ * _EAmap->getEffectiveArea(absSCEta) : 0.) + (_ptSmap ?  _ptSmap->getPtScaling(absSCEta, phoPt_) : 0.);
		Float_t 								correctedIso 								= 	isolation_ - corrections;
		histBDTvsIso->Fill(correctedIso, bdtVar_, xSecW_);
	}

	closeTChain(fTree);
	closeTChain(bTree);
	
	histBDTvsIso->GetXaxis()->CenterTitle();
	histBDTvsIso->GetYaxis()->CenterTitle();
	histBDTvsIso->GetZaxis()->CenterTitle();
	histBDTvsIso->GetZaxis()->RotateTitle();
	histBDTvsIso->GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histBDTvsIso->GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histBDTvsIso->GetZaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histBDTvsIso->GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	histBDTvsIso->GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	histBDTvsIso->GetZaxis()->SetLabelSize(options.getDouble("pad0ZaxisLabelSize"));
	histBDTvsIso->GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	histBDTvsIso->GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	histBDTvsIso->GetZaxis()->SetTitleOffset(options.getDouble("pad0ZtitleOffset"));
	histBDTvsIso->GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	histBDTvsIso->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	return histBDTvsIso;
};


void getCorrelation(std::string _bdtVar, std::string _isoVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax, ofstream* _isoCorrelationOutFilePath, effectiveAreaMap* _EAmap, isoPtScalingMap* _ptSmap){
	
	std::string 								plotName 									=	_bdtVar + "_vs_" + _isoVar + "_" + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p") + "_" +findAndReplaceAll(removeTrailingZeros(_pTmin) + "to" + removeTrailingZeros(_pTmax), ".", "p");
	std::string 								writePath 									= 	options.get("writeDir") + "/" + std::to_string(plotSignal) + "_" + plotName;

	TH2F* 										histBDTvsIso;
	if(useSaved && file_exists(writePath+".root")){
		histBDTvsIso 																		=	(TH2F*)	getHistFromFile(plotName, writePath+".root");
	} else{
		histBDTvsIso																		= 	get2DHistogram(_bdtVar, _isoVar, _etaMin, _etaMax, _pTmin, _pTmax, _EAmap, _ptSmap);
		histBDTvsIso->SetName(plotName.c_str());
		writeToFile(histBDTvsIso, writePath+".root", "RECREATE", 1);
	}
	
	if(options.getInt("autoZrange")) histBDTvsIso->GetZaxis()->SetRangeUser(histBDTvsIso->GetMinimum(std::numeric_limits<Double_t>::min()) / options.getFloat("zMinRatio"), histBDTvsIso->GetMaximum() * options.getFloat("zMaxRatio"));

	TCanvas canvas((plotName+"_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);
	canvas.Draw();
	canvas.cd();

	TPad 										pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad0.SetFillStyle(4000);
	pad0.SetFillColor(0);
	pad0.SetFrameFillStyle(4000);
	pad0.SetGrid(1,1);
	canvas.cd();
	pad0.Draw();
	pad0.cd();
	pad0.SetLogz();
	histBDTvsIso->Draw("colz");

	std::string 								legTitle 									=	(plotSignal ? "Prompt #bullet " : "Fake  ") + removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax) + "  " + removeTrailingZeros(_pTmin) + "<p_{T}#leq" + removeTrailingZeros(_pTmax);
	TLegend 									legend(options.getDouble("legx1"), options.getDouble("legy1"), options.getDouble("legx2"), options.getDouble("legy2"));
	legend.SetHeader(legTitle.c_str(), "C");
	legend.SetTextSize(options.getDouble("legTextSize"));
	legend.SetNColumns(options.getInt("legNcols"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetLineColor(options.getTColFromHex("legFillColor"));
	legend.SetBorderSize(options.getInt("legBorderWidth"));
	legend.SetMargin(options.getFloat("legMargin"));
	legend.SetTextAlign(options.getInt("legTextAlign"));

	TProfile* 									meanBDTvsIso 								= 	histBDTvsIso->ProfileX((plotName+"_pfx").c_str());
	meanBDTvsIso->SetLineColor(options.getTColFromHex("profCol"));
	meanBDTvsIso->SetLineWidth(options.getFloat("profWidth"));
	meanBDTvsIso->SetMarkerColor(options.getTColFromHex("profCol"));
	meanBDTvsIso->SetMarkerSize(options.getFloat("profMkrSz"));
	meanBDTvsIso->SetMarkerStyle(options.getInt("profMkrStl"));
	meanBDTvsIso->Draw("A SAME PE");

	

	Float_t 								pR 											=	histBDTvsIso->GetCorrelationFactor();
	Float_t 								sR 											=	spearmanR(histBDTvsIso, options.getInt("spearmanrNbins"));
	std::string 							profLegStr		= "Mean (r_{P}=" + to_string_with_precision(pR, options.getInt("corrPrecison")) + ", " + "r_{S}=" + to_string_with_precision(sR, options.getInt("corrPrecison")) + ")";
	legend.AddEntry(meanBDTvsIso, profLegStr.c_str(), "LPE");

	std::string outFileString = removeTrailingZeros(_etaMin) + ", " + removeTrailingZeros(_etaMax) + ", " +
	removeTrailingZeros(_pTmin) + ", " + removeTrailingZeros(_pTmax) + ", " +
	to_string_with_precision(pR, 10) + ", " +
	to_string_with_precision(sR, 10);

	(*_isoCorrelationOutFilePath) << outFileString << endl;

	pad0.cd();

	legend.Draw();

	histBDTvsIso->GetYaxis()->SetRangeUser(options.getFloatList("yRange")[0], options.getFloatList("yRange")[1]);
	histBDTvsIso->GetXaxis()->SetRangeUser(options.getFloatList("xRange")[0], options.getFloatList("xRange")[1]);
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	
	canvas.SaveAs((writePath+".png").c_str());
	canvas.SaveAs((writePath+".pdf").c_str());

	clearHeap();
};



void 															initialize(){
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	options.parseIt(optionsFile, "==");
	gStyle->SetPalette(options.getInt("colorPalette"));

	CSVReader 									varFile(options.get("prettNamesFile"), "==");
	std::vector<std::vector<std::string>>		varDat 										=	varFile.getData();
	for(const std::vector<std::string> & iRow : varDat){
		if(iRow.size()==3) {
			std::vector<std::string> 			iVarPlotInfo 								=	{iRow[1]};
			std::vector<std::string> 			binInfo 									=	split_string(iRow[2]);
			iVarPlotInfo.insert(iVarPlotInfo.end(),binInfo.begin(), binInfo.end());
			varPlotInfo[iRow[0]]															=	iVarPlotInfo;
		}
	};

	mkdir(options.get("writeDir"));
	
	std::string 								bdtVar 										=	options.get("bdtVar");
	std::vector<std::string> 					isoVars 									=	options.getList("isoVars", ",");
	std::vector<std::string> 					etaBins 									=	options.getList("etaBins", ";");
	std::vector<std::string> 					pTbins										= 	options.getList("pTbins", ";");
	truncateCorr																			= 	options.getInt("truncateCorr");
	plotSignal																				= 	options.getInt("plotSignal");
	useSaved																				= 	options.getInt("useSaved");
	

	for(std::string isoVar : isoVars){

		std::cout<<"\n\nIsolation: "<<isoVar<<std::endl;

		if(!options.keyExists(isoVar+"EA")){
			std::cout<<"Key "<<isoVar<<"EA"<<" not found! Skipping..."<<std::endl;
			continue;
		}

		Bool_t									applypTcorr									= 	!match("*Trk*", isoVar) && !match("*Ch*", isoVar);
		if(applypTcorr && !options.keyExists(isoVar+"PtScaling")){
			std::cout<<"Key "<<isoVar<<"PtScaling"<<" not found! Skipping..."<<std::endl;
			continue;
		}

		std::string 							isoCorrelationOutFilePath	 				= 	options.get("writeDir") + "/" +std::to_string(plotSignal) + "_" + isoVar + ".txt";
		ofstream* 								isoCorrelationOutFile 						=	new ofstream(isoCorrelationOutFilePath, std::ofstream::out);

		effectiveAreaMap*						isoEA 										=	new effectiveAreaMap(options.get(isoVar+"EA"), 1, ",", 0);
		std::string 							isoPtScalingPath							=	applypTcorr? options.getList(isoVar+"PtScaling")[1] : "";
		Bool_t 		 							isoPtScalingIsQuadratic						=	applypTcorr? (std::stoi(options.getList(isoVar+"PtScaling")[0])==2) : 0;
		isoPtScalingMap*						isoPtScaling								=	applypTcorr ? new isoPtScalingMap(isoPtScalingPath, isoPtScalingIsQuadratic, 1, ",", 0) : nullptr;
		
		for(std::string iPtBinStr : pTbins){
			std::vector<Float_t> iPtMinMax = strToFloatList(iPtBinStr, ",");
			std::cout<<"Pt: "<<iPtMinMax[0]<<"-"<<iPtMinMax[1]<<std::endl;
			for(std::string iEtaBinStr : etaBins){
				std::vector<Float_t> iEtaMinMax = strToFloatList(iEtaBinStr, ",");
				std::cout<<"Eta: "<<iEtaMinMax[0]<<"-"<<iEtaMinMax[1]<<std::endl;
				getCorrelation(bdtVar, isoVar, iEtaMinMax[0], iEtaMinMax[1], iPtMinMax[0], iPtMinMax[1], isoCorrelationOutFile, isoEA, isoPtScaling);
			}
		}
		
		isoCorrelationOutFile->close();
		delete isoCorrelationOutFile;
		delete isoEA;
		delete isoPtScaling;
	}

	// std::string 							isoCorrelationOutFilePath	 				= 	options.get("writeDir") + "/" +std::to_string(plotSignal) + "_CombinedIso.txt";
	// ofstream 								isoCorrelationOutFile 							(isoCorrelationOutFilePath, std::ofstream::out);
	// for(std::string iPtBinStr : pTbins){
	// 	std::vector<Float_t> iPtMinMax = strToFloatList(iPtBinStr, ",");
	// 	std::cout<<"Pt: "<<iPtMinMax[0]<<"-"<<iPtMinMax[1]<<std::endl;
	// 	for(std::string iEtaBinStr : etaBins){
	// 		std::vector<Float_t> iEtaMinMax = strToFloatList(iEtaBinStr, ",");
	// 		std::cout<<"Eta: "<<iEtaMinMax[0]<<"-"<<iEtaMinMax[1]<<std::endl;
	// 		getCombIsoCorrelation(bdtVar, iEtaMinMax[0], iEtaMinMax[1], iPtMinMax[0], iPtMinMax[1], &isoCorrelationOutFile);
	// 	}
	// }
	// isoCorrelationOutFile.close();
};


TH2F*															get2DHistogramCombIso(std::string _bdtVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax){

	std::string 								histName 									=	_bdtVar + "_vs_combIso_" + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p") + findAndReplaceAll(removeTrailingZeros(_pTmin) + "to" + removeTrailingZeros(_pTmax), ".", "p");
	std::string 								yTitle 										=	options.get("combIsoTitle");
	std::string 								histTitle 									=	";"	+ yTitle + ";" +  varPlotInfo[_bdtVar][0];
	std::vector<Float_t> 						combIsoBinning 								= 	options.getFloatList("combIsoBinning");
	TH2F* 										histBDTvsIso								=	new TH2F(histName.c_str(), histTitle.c_str(),
		Int_t(combIsoBinning[0]), combIsoBinning[1], combIsoBinning[2],
		std::stoi(varPlotInfo[_bdtVar][1]), std::stof(varPlotInfo[_bdtVar][2]), std::stof(varPlotInfo[_bdtVar][3]));
	
	histBDTvsIso->SetContour(options.getInt("2dNcontours"));
	histBDTvsIso->Sumw2();
	
	TChain*										fTree 										=	openTChain((std::vector<std::string>){options.get("featsFile")}, options.get("featsTree"));
	TChain*										bTree 										=	openTChain((std::vector<std::string>){options.get("bdtFile")}, options.get("bdtTree"));
	fTree->AddFriend(bTree);
	TTreeReader                             	inputTTreeReader(fTree);
	TTreeReaderAnyValue<Float_t>				xSecW_											(inputTTreeReader, "xSecW");
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Float_t>				phoPt_											(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Bool_t>					isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Bool_t>					isTrain_										(inputTTreeReader, "isTrain");
	TTreeReaderAnyValue<UChar_t>				phoQualityBits_									(inputTTreeReader, "phoQualityBits");
	TTreeReaderAnyValue<Float_t>				rho_ 											(inputTTreeReader, "rho");
	TTreeReaderAnyValue<Float_t>				phoTrkSumPtHollowConeDR03_						(inputTTreeReader, "phoTrkSumPtHollowConeDR03");
	TTreeReaderAnyValue<Float_t>				phoPFClusEcalIso_								(inputTTreeReader, "phoPFClusEcalIso");
	TTreeReaderAnyValue<Float_t>				phoPFClusHcalIso_								(inputTTreeReader, "phoPFClusHcalIso");
	TTreeReaderAnyValue<Double_t>				bdtVar_											(inputTTreeReader, _bdtVar);

	effectiveAreaMap							phoPFClusEcalIsoEA								(options.get("phoPFClusEcalIsoEA"), 1, ",", 0);
	effectiveAreaMap							phoPFClusHcalIsoEA								(options.get("phoPFClusHcalIsoEA"), 1, ",", 0);
	effectiveAreaMap							phoTrkSumPtHollowConeDR03EA						(options.get("phoTrkSumPtHollowConeDR03EA"), 1, ",", 0);

	isoPtScalingMap								phoPFClusEcalIsoPtScaling						(options.getList("phoPFClusEcalIsoPtScaling")[1], 0, 1, ",", 0);
	isoPtScalingMap								phoPFClusHcalIsoPtScaling						(options.getList("phoPFClusHcalIsoPtScaling")[1], 1, 1, ",", 0);

	while(inputTTreeReader.Next()){
		if(isTrain_) continue;
		if(plotSignal != isSignal_) continue;
		if(phoPt_ < _pTmin) continue;
		if(phoPt_ >= _pTmax) continue;
		Float_t 								absSCEta									=	std::abs(phoSCeta_);
		if (absSCEta <= _etaMin) continue;
		if (absSCEta > _etaMax) continue;
		UChar_t 								iPhoQalBits 								=	phoQualityBits_;
		if(getBit(iPhoQalBits,0)) continue;
		Float_t 								corrections 								=	rho_ * (phoPFClusEcalIsoEA.getEffectiveArea(absSCEta) + phoPFClusHcalIsoEA.getEffectiveArea(absSCEta) + phoTrkSumPtHollowConeDR03EA.getEffectiveArea(absSCEta));
		corrections 																		+=	phoPFClusEcalIsoPtScaling.getPtScaling(absSCEta, phoPt_) + phoPFClusHcalIsoPtScaling.getPtScaling(absSCEta, phoPt_);
		Float_t 								correctedIso 								= 	phoPFClusEcalIso_ + phoPFClusHcalIso_ + phoTrkSumPtHollowConeDR03_ - corrections;
		histBDTvsIso->Fill(correctedIso, bdtVar_, xSecW_);
	}

	closeTChain(fTree);
	closeTChain(bTree);
	
	histBDTvsIso->GetXaxis()->CenterTitle();
	histBDTvsIso->GetYaxis()->CenterTitle();
	histBDTvsIso->GetZaxis()->CenterTitle();
	histBDTvsIso->GetZaxis()->RotateTitle();
	histBDTvsIso->GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histBDTvsIso->GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histBDTvsIso->GetZaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histBDTvsIso->GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	histBDTvsIso->GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	histBDTvsIso->GetZaxis()->SetLabelSize(options.getDouble("pad0ZaxisLabelSize"));
	histBDTvsIso->GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	histBDTvsIso->GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	histBDTvsIso->GetZaxis()->SetTitleOffset(options.getDouble("pad0ZtitleOffset"));
	histBDTvsIso->GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	histBDTvsIso->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	return histBDTvsIso;
};




void getCombIsoCorrelation(std::string _bdtVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax, ofstream* _isoCorrelationOutFilePath){
	
	std::string 								plotName 									=	_bdtVar + "_vs_combIso_" + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p") + "_" +findAndReplaceAll(removeTrailingZeros(_pTmin) + "to" + removeTrailingZeros(_pTmax), ".", "p");
	std::string 								writePath 									= 	options.get("writeDir") + "/" + std::to_string(plotSignal) + "_" + plotName;

	TH2F* 										histBDTvsIso;
	if(useSaved && file_exists(writePath+".root")){
		histBDTvsIso 																		=	(TH2F*)	getHistFromFile(plotName, writePath+".root");
	} else{
		histBDTvsIso																		= 	get2DHistogramCombIso(_bdtVar,_etaMin, _etaMax, _pTmin, _pTmax);
		histBDTvsIso->SetName(plotName.c_str());
		writeToFile(histBDTvsIso, writePath+".root", "RECREATE", 1);
	}
	
	if(options.getInt("autoZrange")) histBDTvsIso->GetZaxis()->SetRangeUser(histBDTvsIso->GetMinimum(std::numeric_limits<Double_t>::min()) / options.getFloat("zMinRatio"), histBDTvsIso->GetMaximum() * options.getFloat("zMaxRatio"));

	TCanvas canvas((plotName+"_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);
	canvas.Draw();
	canvas.cd();

	TPad 										pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad0.SetFillStyle(4000);
	pad0.SetFillColor(0);
	pad0.SetFrameFillStyle(4000);
	pad0.SetGrid(1,1);
	canvas.cd();
	pad0.Draw();
	pad0.cd();
	pad0.SetLogz();
	histBDTvsIso->Draw("colz");

	std::string 								legTitle 									=	(plotSignal ? "Prompt #bullet " : "Fake  ") + removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax) + "  " + removeTrailingZeros(_pTmin) + "<p_{T}#leq" + removeTrailingZeros(_pTmax);
	TLegend 									legend(options.getDouble("legx1"), options.getDouble("legy1"), options.getDouble("legx2"), options.getDouble("legy2"));
	legend.SetHeader(legTitle.c_str(), "C");
	legend.SetTextSize(options.getDouble("legTextSize"));
	legend.SetNColumns(options.getInt("legNcols"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetLineColor(options.getTColFromHex("legFillColor"));
	legend.SetBorderSize(options.getInt("legBorderWidth"));
	legend.SetMargin(options.getFloat("legMargin"));
	legend.SetTextAlign(options.getInt("legTextAlign"));

	TProfile* 									meanBDTvsIso 								= 	histBDTvsIso->ProfileX((plotName+"_pfx").c_str());
	meanBDTvsIso->SetLineColor(options.getTColFromHex("profCol"));
	meanBDTvsIso->SetLineWidth(options.getFloat("profWidth"));
	meanBDTvsIso->SetMarkerColor(options.getTColFromHex("profCol"));
	meanBDTvsIso->SetMarkerSize(options.getFloat("profMkrSz"));
	meanBDTvsIso->SetMarkerStyle(options.getInt("profMkrStl"));
	meanBDTvsIso->Draw("A SAME PE");

	

	Float_t 								pR 											=	histBDTvsIso->GetCorrelationFactor();
	Float_t 								sR 											=	spearmanR(histBDTvsIso, options.getInt("spearmanrNbins"));
	std::string 							profLegStr		= "Mean (r_{P}=" + to_string_with_precision(pR, options.getInt("corrPrecison")) + ", " + "r_{S}=" + to_string_with_precision(sR, options.getInt("corrPrecison")) + ")";
	legend.AddEntry(meanBDTvsIso, profLegStr.c_str(), "LPE");

	std::string outFileString = removeTrailingZeros(_etaMin) + ", " + removeTrailingZeros(_etaMax) + ", " +
	 removeTrailingZeros(_pTmin) + ", " + removeTrailingZeros(_pTmax) + ", " +
	to_string_with_precision(pR, 10) + ", " +
	to_string_with_precision(sR, 10);

	(*_isoCorrelationOutFilePath) << outFileString << endl;

	pad0.cd();

	legend.Draw();

	histBDTvsIso->GetYaxis()->SetRangeUser(options.getFloatList("yRange")[0], options.getFloatList("yRange")[1]);
	histBDTvsIso->GetXaxis()->SetRangeUser(options.getFloatList("xRange")[0], options.getFloatList("xRange")[1]);
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	
	canvas.SaveAs((writePath+".png").c_str());
	canvas.SaveAs((writePath+".pdf").c_str());

	clearHeap();
};
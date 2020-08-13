#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TArrow.h"
#include "TPaveText.h"
#include "TLine.h"

std::map<std::string,std::vector<std::string>>					varPlotInfo;
std::string 													optionsFile 			=	"fitPtScalingVSetaOptions.txt";
parseOptions 													options;
std::vector<std::string>										samples;
std::vector<Double_t>											pTbinning;
std::vector<std::string> 										etaBins;
Bool_t 															applyEAcorr;
Bool_t 															truncateCorr;
Float_t 														percentile;
std::vector<Double_t> 											fitRange;
TH2F*															get2DHistogram(std::string _isoVar, Float_t _etaMin, Float_t _etaMax, effectiveAreaMap* _EAmap);
void 															initialize();
void 															getEffectiveArea(std::string _isoVar, ofstream* _pTscalingOutFile, effectiveAreaMap* _EAmap);


TH2F*															get2DHistogram(std::string _isoVar, Float_t _etaMin, Float_t _etaMax, effectiveAreaMap* _EAmap){

	std::string 								histName 									=	_isoVar + "Corr_vs_pT_eta" + findAndReplaceAll(std::to_string(_etaMin) + "to" + std::to_string(_etaMax), ".", "p");
	std::string 								histTitle 									=	";"	+ varPlotInfo["phoPt"][0] + ";" + (applyEAcorr ? (findAndReplaceAll(varPlotInfo[_isoVar][0], "(", " - #rhoA_{eff}^{" + options.get("percentile")+ "%} (")) : varPlotInfo[_isoVar][0]) ;
	TH2F* 										histIsoVsPt								=	new TH2F(histName.c_str(), histTitle.c_str(), pTbinning.size()-1, pTbinning.data(),
		std::stoi(varPlotInfo[_isoVar][1]), std::stof(varPlotInfo[_isoVar][2]), std::stof(varPlotInfo[_isoVar][3]));
	histIsoVsPt->Sumw2();
	
	for(std::string iSample : samples){
		std::string 							iSamplePath 								=	options.get("fileDir") + "/" + iSample + ".root" ;
		TChain*									tChain 										=	openTChain((std::vector<std::string>){iSamplePath}, options.get("inTreeName"));
		TTreeReader                             inputTTreeReader(tChain);
		TTreeReaderAnyValue<Float_t>			genWeight_										(inputTTreeReader, "genWeight");
		TTreeReaderAnyValue<Float_t>			puWeight_										(inputTTreeReader, "puWeight");
		TTreeReaderAnyValue<Float_t>			phoSCeta_										(inputTTreeReader, "phoSCeta");
		TTreeReaderAnyValue<Float_t>			phoPt_											(inputTTreeReader, "phoPt");		
		TTreeReaderAnyValue<Float_t>			rho_											(inputTTreeReader, "rho");	
		TTreeReaderAnyValue<Float_t>			isolation_										(inputTTreeReader, _isoVar);
		TH1F* 									cutFlowGenWeight 							= 	(TH1F*) getHistFromFile(options.get("cutFlowHist"), iSamplePath);
		Double_t 								sumGenWeight 								= 	cutFlowGenWeight->GetBinContent(1);
		delete cutFlowGenWeight;
		Float_t 								xSection 									=	std::stof(vLookup(iSample, options.get("xSectionMap"), 0, 2));

		while(inputTTreeReader.Next()){
			if(phoPt_ < 200.) continue;
			Float_t 							absSCEta									=	std::abs(phoSCeta_);
			if (absSCEta <= _etaMin) continue;
			if (absSCEta > _etaMax) continue;
			Double_t 							weight 										=  	genWeight_ * puWeight_ * xSection * 1000000./sumGenWeight;
			Float_t 							untruncatedEACorrectedIso					=	applyEAcorr ? (isolation_ - rho_ * _EAmap->getEffectiveArea(absSCEta)) : isolation_;
			Float_t 							EAcorrectedIso 								= 	truncateCorr ? std::max((Float_t)0., untruncatedEACorrectedIso) : untruncatedEACorrectedIso;
			histIsoVsPt->Fill(phoPt_, EAcorrectedIso, weight);
		}
		closeTChain(tChain);
	}

	return histIsoVsPt;
};


void 															initialize(){
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	options.parseIt(optionsFile, "==");

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

	samples																					=	options.getList("samples", ",");
	std::sort(samples.begin(), samples.end());
	mkdir(options.get("writeDir"));
	
	std::vector<std::string> 					isoVarStrs 									=	options.getList("isoVar", ";");
	etaBins 																				=	options.getList("etaBins", ";");
	pTbinning																				= 	options.getDoubleList("pTbinning");
	applyEAcorr																				= 	options.getInt("applyEAcorr");
	truncateCorr 																			=	options.getInt("truncateCorr");
	percentile											 									= 	options.getFloat("percentile")/100.;
	fitRange 																				=	{pTbinning[0], pTbinning.back()};

	for(std::string isoVarStr : isoVarStrs){

		std::string 							isoVar 										= 	split_string(isoVarStr, ",")[0];

		if(!options.keyExists(isoVar+"EA")){
			std::cout<<"Key "<<isoVar<<"EA"<<" not found! Skipping..."<<std::endl;
			continue;
		}

		std::string 							pTScalingOutFilePath 						= 	options.get("writeDir") + "/" +isoVar + ".txt";
		
		ofstream* 								pTScalingOutFile							=	new ofstream(pTScalingOutFilePath, std::ofstream::out);
		
		effectiveAreaMap*						isoEA 										=	applyEAcorr ? new effectiveAreaMap(options.get(isoVar+"EA"), 1, ",", 0) : nullptr;
		
		getEffectiveArea(isoVarStr, pTScalingOutFile, isoEA);

		
		pTScalingOutFile->close();
		delete pTScalingOutFile;
		delete isoEA;
	}
	
};


void getEffectiveArea(std::string _isoVarStr, ofstream* _pTscalingOutFile, effectiveAreaMap* _EAmap){

	std::string									isoVar										= 	split_string(_isoVarStr, ",")[0];
	Int_t 										polOder										= 	std::stoi(split_string(_isoVarStr, ",")[1]);

	std::string 								plotName 									=	isoVar + "RhoCorr_vs_pT_eta_" + options.get("percentile") + "pc";
	TCanvas canvas((plotName+"_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);

	TPad 										pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad0.SetFillStyle(4000);
	pad0.SetFillColor(0);
	pad0.SetFrameFillStyle(4000);
	pad0.SetGrid(1,1);

	TLegend 									legend(options.getDouble("legx1"), options.getDouble("legy1"), options.getDouble("legx2"), options.getDouble("legy2"));
	legend.SetTextSize(options.getDouble("legTextSize"));
	legend.SetNColumns(options.getInt("legNcols"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetLineColor(options.getTColFromHex("legFillColor"));
	legend.SetBorderSize(options.getInt("legBorderWidth"));
	legend.SetMargin(options.getFloat("legMargin"));

	canvas.Draw();
	canvas.cd();
	pad0.Draw();
	pad0.cd();

	TH1D* 											firstHist 								= 	nullptr;
	Double_t 										yMax 									= 	-999;
	Double_t 										yMin 									= 	10000;

	for(std::string etaStr : etaBins){

		Float_t 									etaMin 									= 	std::stof(split_string(etaStr, ",")[0]);
		Float_t 									etaMax 									= 	std::stof(split_string(etaStr, ",")[1]);

		TH2F* 										histIsoVsPt								=	get2DHistogram(isoVar, etaMin, etaMax, _EAmap);
		std::string 								etaHistName								=	plotName + "_eta" + split_string(etaStr, ",")[0] +"to" + split_string(etaStr, ",")[1];
		TH1D* 										isoPcVSpT 								= 	histIsoVsPt->QuantilesX(percentile, etaHistName.c_str());

		isoPcVSpT->SetLineColor(hex2rootColor(split_string(etaStr, ",")[2]));
		isoPcVSpT->SetLineWidth(options.getFloat("pcWidth"));
		isoPcVSpT->SetMarkerColor(hex2rootColor(split_string(etaStr, ",")[2]));
		isoPcVSpT->SetMarkerSize(options.getFloat("pcMkrSz"));
		isoPcVSpT->SetMarkerStyle(options.getInt("pcMkrStl"));
		isoPcVSpT->Draw("SAME PE");

		if(!firstHist) {
			firstHist = isoPcVSpT;
			std::string 								histTitle 									=	";"	+ varPlotInfo["phoPt"][0] + ";" + options.get("percentile")+ "%^{ile} " + (applyEAcorr ? (findAndReplaceAll(varPlotInfo[isoVar][0], "(", " - #rhoA_{eff}^{" + options.get("percentile")+ "%} (")) : varPlotInfo[isoVar][0]) ;
			firstHist->SetTitle(histTitle.c_str());
			firstHist->GetXaxis()->CenterTitle();
			firstHist->GetYaxis()->CenterTitle();
			firstHist->GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
			firstHist->GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
			firstHist->GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
			firstHist->GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
			firstHist->GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
			firstHist->GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
			firstHist->GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
			firstHist->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));
		}
		
		std::string outFileString;

		if(polOder == 0){

			TF1* 										linFitLine 								= 	new TF1((etaHistName+"_fitLinear").c_str(), "[0]", fitRange[0], fitRange[1]);
			TFitResultPtr 								linFitResult							=	isoPcVSpT->Fit(linFitLine, options.getCSTR("fitOpt"));
			linFitResult->Print("V");	
			std::string 							linFitText 									= 	"#bf{"+split_string(etaStr, ",")[0] + "#leq|#eta|<" + split_string(etaStr, ",")[1] + "}: ("
			+ to_string_with_precision(linFitLine->GetParameter(0), options.getInt("constantPrecision")) + "#pm"
			+ to_string_with_precision(linFitLine->GetParError(0), options.getInt("constantPrecision")) + ")";
			// + " [#chi^{2}/ndf=" + to_string_with_precision(linFitResult->Chi2 (), options.getInt("chi2precision")) + "/" + std::to_string(linFitResult->Ndf()) + "]";

			std::string 							linFitDef 									=	to_string_with_precision(linFitLine->GetParameter(0), 10);
			TF1* 									linFitCopy 									=  new TF1((etaHistName+"_fitLinearCopy").c_str(), linFitDef.c_str(), histIsoVsPt->GetXaxis()->GetXmin(), histIsoVsPt->GetXaxis()->GetXmax());
			linFitCopy->SetLineColor(hex2rootColor(split_string(etaStr, ",")[2]));
			linFitCopy->SetLineWidth(options.getFloat("fitLnWidth"));
			linFitCopy->SetLineStyle(options.getInt("fitLnStl"));
			linFitCopy->Draw("A SAME L");
			legend.AddEntry(isoPcVSpT, linFitText.c_str() , "LPE");

			outFileString = split_string(etaStr, ",")[0] + ", " + split_string(etaStr, ",")[1] + ", " + to_string_with_precision(linFitLine->GetParameter(0), 10);

			yMax 																				= 	std::max(yMax, linFitLine->Eval(histIsoVsPt->GetXaxis()->GetBinUpEdge(histIsoVsPt->GetNbinsX())));
			yMin 																				= 	std::min(yMin, linFitLine->Eval(histIsoVsPt->GetXaxis()->GetBinLowEdge(1)));

		} else if(polOder == 1){
			TF1* 										linFitLine 								= 	new TF1((etaHistName+"_fitLinear").c_str(), "[0]+[1]*x", fitRange[0], fitRange[1]);
			TFitResultPtr 								linFitResult							=	isoPcVSpT->Fit(linFitLine, options.getCSTR("fitOpt"));
			linFitResult->Print("V");	
			std::string 							linFitText 									= 	"#bf{"+split_string(etaStr, ",")[0] + "#leq|#eta|<" + split_string(etaStr, ",")[1] + "}: ("
			+ to_string_with_precision(linFitLine->GetParameter(1), options.getInt("slopePrecision")) + "#pm" 
			+ to_string_with_precision(linFitLine->GetParError(1), options.getInt("slopePrecision")) + ")p_{T}+("
			+ to_string_with_precision(linFitLine->GetParameter(0), options.getInt("interceptPrecision")) + "#pm"
			+ to_string_with_precision(linFitLine->GetParError(0), options.getInt("interceptPrecision")) + ")";
			// + " [#chi^{2}/ndf=" + to_string_with_precision(linFitResult->Chi2 (), options.getInt("chi2precision")) + "/" + std::to_string(linFitResult->Ndf()) + "]";

			std::string 							linFitDef 									=	to_string_with_precision(linFitLine->GetParameter(1), 10) + "*x + " + to_string_with_precision(linFitLine->GetParameter(0), 10);
			TF1* 									linFitCopy 									=  new TF1((etaHistName+"_fitLinearCopy").c_str(), linFitDef.c_str(), histIsoVsPt->GetXaxis()->GetXmin(), histIsoVsPt->GetXaxis()->GetXmax());
			linFitCopy->SetLineColor(hex2rootColor(split_string(etaStr, ",")[2]));
			linFitCopy->SetLineWidth(options.getFloat("fitLnWidth"));
			linFitCopy->SetLineStyle(options.getInt("fitLnStl"));
			linFitCopy->Draw("A SAME L");
			legend.AddEntry(isoPcVSpT, linFitText.c_str() , "LPE");

			outFileString = split_string(etaStr, ",")[0] + ", " + split_string(etaStr, ",")[1] + ", " + to_string_with_precision(linFitLine->GetParameter(0), 10) + ", " + to_string_with_precision(linFitLine->GetParameter(1), 10);

			yMax 																				= 	std::max(yMax, linFitLine->Eval(histIsoVsPt->GetXaxis()->GetBinUpEdge(histIsoVsPt->GetNbinsX())));
			yMin 																				= 	std::min(yMin, linFitLine->Eval(histIsoVsPt->GetXaxis()->GetBinLowEdge(1)));
		} else if(polOder == 2){
			TF1* 									quadFitLine 								= 	new TF1((etaHistName+"_fitQuad").c_str(), "[0]+[1]*x+[2]*x*x", fitRange[0], fitRange[1]);
			TFitResultPtr 							quadFitResult 								=	isoPcVSpT->Fit(quadFitLine, options.getCSTR("fitOpt"));
			quadFitResult->Print("V");

			std::string 							quadFitText 								= 	"#bf{"+split_string(etaStr, ",")[0] + "#leq|#eta|<" + split_string(etaStr, ",")[1] + "}: ("
			+ to_string_with_precision(quadFitLine->GetParameter(2), options.getInt("c2Precision")) + "#pm" 
			+ to_string_with_precision(quadFitLine->GetParError(2), options.getInt("c2Precision")) + ")p_{T}^{2}+("
			+ to_string_with_precision(quadFitLine->GetParameter(1), options.getInt("c1Precision")) + "#pm"
			+ to_string_with_precision(quadFitLine->GetParError(1), options.getInt("c1Precision")) + ")p_{T}+("
			+ to_string_with_precision(quadFitLine->GetParameter(0), options.getInt("c0Precision")) + "#pm"
			+ to_string_with_precision(quadFitLine->GetParError(0), options.getInt("c0Precision")) + ")";
			// + " [#chi^{2}/ndf=" + to_string_with_precision(quadFitResult->Chi2 (), options.getInt("chi2precision")) + "/" + std::to_string(quadFitResult->Ndf()) + "]";

			std::string 							quadFitDef 									=	to_string_with_precision(quadFitLine->GetParameter(2), 10) + "*x*x +" + to_string_with_precision(quadFitLine->GetParameter(1), 10) + "*x + " + to_string_with_precision(quadFitLine->GetParameter(0), 10);
			TF1* 									quadFitCopy 									=  new TF1((plotName+"_fitQuadCopy").c_str(), quadFitDef.c_str(), histIsoVsPt->GetXaxis()->GetXmin(), histIsoVsPt->GetXaxis()->GetXmax());
			quadFitCopy->SetLineColor(hex2rootColor(split_string(etaStr, ",")[2]));
			quadFitCopy->SetLineWidth(options.getFloat("fitLnWidth"));
			quadFitCopy->SetLineStyle(options.getInt("fitLnStl"));
			quadFitCopy->Draw("A SAME L");
			legend.AddEntry(isoPcVSpT, quadFitText.c_str() , "LPE");

			outFileString 																		= split_string(etaStr, ",")[0] + ", " + split_string(etaStr, ",")[1] +
			", " + to_string_with_precision(quadFitLine->GetParameter(0), 10) +
			", "+ to_string_with_precision(quadFitLine->GetParameter(1), 10) +
			", "+ to_string_with_precision(quadFitLine->GetParameter(2), 10);

			yMax 																				= std::max(yMax, quadFitLine->Eval(histIsoVsPt->GetXaxis()->GetBinUpEdge(histIsoVsPt->GetNbinsX())));
			yMin 																				= 	std::min(yMin, quadFitLine->Eval(histIsoVsPt->GetXaxis()->GetBinLowEdge(1)));
		}

		(*_pTscalingOutFile) << outFileString << endl;
		
	}

	if(polOder == 0)firstHist->GetYaxis()->SetRangeUser(0., options.getFloat("yMaxRatio0")*yMax);
	else firstHist->GetYaxis()->SetRangeUser(std::min(0.,yMin), options.getFloat("yMaxRatio12")*yMax);
	pad0.cd();

	legend.Draw();
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	std::string writePath 									= 	options.get("writeDir") + "/" + plotName;
	canvas.SaveAs((writePath+".png").c_str());
	canvas.SaveAs((writePath+".pdf").c_str());

	clearHeap();
};
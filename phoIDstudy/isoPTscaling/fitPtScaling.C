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
std::string 													optionsFile 			=	"fitPtScalingOptions.txt";
parseOptions 													options;
std::vector<std::string>										samples;
std::vector<Double_t>											pTbinning;
TH2F*															get2DHistogram(std::string _isoVar, Float_t _etaMin, Float_t _etaMax, effectiveAreaMap* _EAmap);
void 															initialize();
void 															getEffectiveArea(std::string _isoVar, Float_t _etaMin, Float_t _etaMax, ofstream* _pTscalingOutFile, effectiveAreaMap* _EAmap);

TH2F*															get2DHistogram(std::string _isoVar, Float_t _etaMin, Float_t _etaMax, effectiveAreaMap* _EAmap){

	std::string 								histName 									=	_isoVar + "Corr_vs_pT_eta" + findAndReplaceAll(std::to_string(_etaMin) + "to" + std::to_string(_etaMax), ".", "p");
	std::string 								histTitle 									=	";"	+ varPlotInfo["phoPt"][0] + ";" + findAndReplaceAll(varPlotInfo[_isoVar][0], "(", " - #rhoA_{eff}^{" + options.get("percentile")+ "%} (");
	TH2F* 										histIsoVsPt								=	new TH2F(histName.c_str(), histTitle.c_str(), pTbinning.size()-1, pTbinning.data(),
		std::stoi(varPlotInfo[_isoVar][1]), std::stof(varPlotInfo[_isoVar][2]), std::stof(varPlotInfo[_isoVar][3]));
	histIsoVsPt->SetContour(options.getInt("2dNcontours"));

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
			// Float_t 							EAcorrectedIso 								= 	std::max((Float_t)0., isolation_ - rho_ * _EAmap->getEffectiveArea(absSCEta));
			// Float_t 							EAcorrectedIso 								= 	isolation_ - rho_ * _EAmap->getEffectiveArea(absSCEta);
			Float_t 							EAcorrectedIso 								= 	isolation_;
			histIsoVsPt->Fill(phoPt_, EAcorrectedIso, weight);
		}
		closeTChain(tChain);
	}
	
	histIsoVsPt->GetXaxis()->CenterTitle();
	histIsoVsPt->GetYaxis()->CenterTitle();
	histIsoVsPt->GetZaxis()->CenterTitle();
	histIsoVsPt->GetZaxis()->RotateTitle();
	histIsoVsPt->GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histIsoVsPt->GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histIsoVsPt->GetZaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histIsoVsPt->GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	histIsoVsPt->GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	histIsoVsPt->GetZaxis()->SetLabelSize(options.getDouble("pad0ZaxisLabelSize"));
	histIsoVsPt->GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	histIsoVsPt->GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	histIsoVsPt->GetZaxis()->SetTitleOffset(options.getDouble("pad0ZtitleOffset"));
	histIsoVsPt->GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	histIsoVsPt->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	return histIsoVsPt;
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

	samples																					=	options.getList("samples", ",");
	std::sort(samples.begin(), samples.end());
	mkdir(options.get("writeDir"));
	
	std::vector<std::string> 					isoVars 									=	options.getList("isoVar", ",");
	std::vector<std::string> 					etaBins 									=	options.getList("etaBins", ";");
	pTbinning																				= options.getDoubleList("pTbinning");

	for(std::string isoVar : isoVars){

		if(!options.keyExists(isoVar+"EA")){
			std::cout<<"Key "<<isoVar<<"EA"<<" not found! Skipping..."<<std::endl;
			continue;
		}

		std::string 							pTScalingOutFilePath 						= 	options.get("writeDir") + "/" +isoVar + ".txt";
		
		ofstream* 								pTScalingOutFile							=	new ofstream(pTScalingOutFilePath, std::ofstream::out);
		
		effectiveAreaMap*						isoEA 										=	new effectiveAreaMap(options.get(isoVar+"EA"), 1, ",", 0);
		
		for(std::string iEtaBinStr : etaBins){
			
			std::vector<Float_t> iEtaMinMax = strToFloatList(iEtaBinStr, ",");
			
			getEffectiveArea(isoVar, iEtaMinMax[0], iEtaMinMax[1], pTScalingOutFile, isoEA);
			
		}
		
		pTScalingOutFile->close();
		delete pTScalingOutFile;
		delete isoEA;
	}
	
};


void getEffectiveArea(std::string _isoVar, Float_t _etaMin, Float_t _etaMax, ofstream* _pTscalingOutFile, effectiveAreaMap* _EAmap){

	TH2F* 										histIsoVsPt								=	get2DHistogram(_isoVar, _etaMin, _etaMax, _EAmap);

	std::string 								plotName 									=	_isoVar + "Corr_vs_pT_eta" + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p") + "_" + options.get("percentile") + "pc";
	TCanvas canvas((plotName+"_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);

	canvas.Draw();
	canvas.cd();

	//// X projection
	TPad padUwing("padUwing", "", options.getFloat("padUwingx1"), options.getFloat("padUwingy1"), options.getFloat("padUwingx2"), options.getFloat("padUwingy2"));
	padUwing.SetMargin(options.getFloat("padUwingmarginL"), options.getFloat("padUwingmarginR"), options.getFloat("padUwingmarginB"), options.getFloat("padUwingmarginT"));
	padUwing.SetFillStyle(0);
	padUwing.SetFillColor(0);
	padUwing.SetFrameFillStyle(0);
	canvas.cd();
	padUwing.Draw();
	padUwing.SetGrid(1,0);
	padUwing.cd();

	TH1D* 										ptProj 									= 	histIsoVsPt->ProjectionX((plotName+"_px").c_str());
	ptProj->SetFillStyle(1001);
	ptProj->SetFillColor(options.getTColFromHex("projXcol"));
	ptProj->SetMinimum(0.);
	ptProj->SetLineWidth(0.);
	ptProj->GetXaxis()->SetTitle("");
	ptProj->GetYaxis()->SetTitle("");
	ptProj->GetXaxis()->SetLabelSize(0);
	ptProj->GetYaxis()->SetLabelSize(0);
	ptProj->SetStats(0);
	ptProj->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));
	ptProj->GetYaxis()->SetNdivisions(options.getInt("padUwingyNdivs"));
	ptProj->Draw("HIST");

	std::vector<Double_t> 						fitRange(pTbinning[0], pTbinning.back());

	padUwing.RedrawAxis();
	padUwing.RedrawAxis("G");
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");

	canvas.cd();
	TLatex 										texWriter;
	texWriter.SetNDC(0);
	texWriter.SetTextAlign(11);
	texWriter.SetTextSize(options.getFloat("categoryInfoTextSize"));
	texWriter.SetTextColor(options.getTColFromHex("categoryInfoTextColor"));
	std::string 								etaTitle 									=	removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax);
	texWriter.DrawLatex(options.getFloat("categoryInfoX"), options.getFloat("categoryInfoY"), etaTitle.c_str());

	//// Y projection
	TPad 										padRwing("padRwing", "", options.getFloat("padRwingx1"), options.getFloat("padRwingy1"), options.getFloat("padRwingx2"), options.getFloat("padRwingy2"));
	padRwing.SetMargin(options.getFloat("padRwingmarginL"), options.getFloat("padRwingmarginR"), options.getFloat("padRwingmarginB"), options.getFloat("padRwingmarginT"));
	padRwing.SetFillStyle(0);
	padRwing.SetFillColor(0);
	padRwing.SetFrameFillStyle(0);
	canvas.cd();
	padRwing.Draw();
	padRwing.SetGrid(0,1);
	padRwing.cd();

	if(options.getInt("padRwingLogY")) {
		padRwing.SetLogx();
	}

	TH1D* 										isoProj 									=	histIsoVsPt->ProjectionY((plotName+"_py").c_str());
	removeNegativeBins(isoProj);
	isoProj->SetFillStyle(1001);
	isoProj->SetFillColor(options.getTColFromHex ("projYcol"));
	isoProj->SetLineWidth(0.);
	isoProj->SetMinimum(0.);
	isoProj->GetXaxis()->SetTitle("");
	isoProj->GetYaxis()->SetTitle("");
	isoProj->GetXaxis()->SetLabelSize(0);
	isoProj->GetYaxis()->SetLabelSize(0);
	isoProj->SetStats(0);
	// isoProj->GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	isoProj->GetXaxis()->SetNdivisions(options.getInt("padRwingyNdivs"));
	isoProj->GetYaxis()->SetNdivisions(options.getInt("padRwingxNdivs"));	
	isoProj->GetYaxis()->SetTickLength(0.);
	isoProj->GetYaxis()->SetTickSize(0.);
	isoProj->SetMaximum(isoProj->GetMaximum()*options.getFloat("yMaxRatio"));
	isoProj->Draw("hbar");

	padRwing.RedrawAxis();
	padRwing.RedrawAxis("G");
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");

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
	histIsoVsPt->Draw("colz");

	TLegend 									legend(options.getDouble("legx1"), options.getDouble("legy1"), options.getDouble("legx2"), options.getDouble("legy2"));
	legend.SetTextSize(options.getDouble("legTextSize"));
	legend.SetNColumns(options.getInt("legNcols"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetLineColor(options.getTColFromHex("legFillColor"));
	legend.SetBorderSize(options.getInt("legBorderWidth"));

	if(options.getInt("autoZrange")) histIsoVsPt->GetZaxis()->SetRangeUser(histIsoVsPt->GetMinimum(std::numeric_limits<Double_t>::min()) / options.getFloat("zMinRatio"), histIsoVsPt->GetMaximum() * options.getFloat("zMaxRatio"));

	Float_t 								percentile 									= 	options.getFloat("percentile")/100.;

	TH1D* 									isoPcVSpT 										= 	histIsoVsPt->QuantilesX(percentile, plotName.c_str());
	isoPcVSpT->SetLineColor(options.getTColFromHex("pcCol"));
	isoPcVSpT->SetLineWidth(options.getFloat("pcWidth"));
	isoPcVSpT->SetMarkerColor(options.getTColFromHex("pcCol"));
	isoPcVSpT->SetMarkerSize(options.getFloat("pcMkrSz"));
	isoPcVSpT->SetMarkerStyle(options.getInt("pcMkrStl"));
	isoPcVSpT->Draw("A SAME PE");
	legend.AddEntry(isoPcVSpT, (options.get("percentile")+"%").c_str() , "LPE");

	TF1* 									linFitLine 									= 	new TF1((plotName+"_fitLinear").c_str(), "[0]+[1]*x", fitRange[0], fitRange[1]);
	TFitResultPtr 							linFitResult 								=	isoPcVSpT->Fit(linFitLine, options.getCSTR("fitOpt"));
	linFitResult->Print("V");	
	std::string 							linFitText 									= 	"("
	+ to_string_with_precision(linFitLine->GetParameter(1), options.getInt("slopePrecision")) + "#pm" 
	+ to_string_with_precision(linFitLine->GetParError(1), options.getInt("slopePrecision")) + ")p_{T}+("
	+ to_string_with_precision(linFitLine->GetParameter(0), options.getInt("interceptPrecision")) + "#pm"
	+ to_string_with_precision(linFitLine->GetParError(0), options.getInt("interceptPrecision")) + ")"
	+ " [#chi^{2}/ndf=" + to_string_with_precision(linFitResult->Chi2 (), options.getInt("chi2precision")) + "/" + std::to_string(linFitResult->Ndf()) + "]";

	std::string 							linFitDef 									=	to_string_with_precision(linFitLine->GetParameter(1), 10) + "*x + " + to_string_with_precision(linFitLine->GetParameter(0), 10);
	TF1* 									linFitCopy 									=  new TF1((plotName+"_fitLinearCopy").c_str(), linFitDef.c_str(), histIsoVsPt->GetXaxis()->GetXmin(), histIsoVsPt->GetXaxis()->GetXmax());
	linFitCopy->SetLineColor(options.getTColFromHex("linearCol"));
	linFitCopy->SetLineWidth(options.getFloat("fitLnWidth"));
	linFitCopy->SetLineStyle(options.getInt("fitLnStl"));
	linFitCopy->Draw("A SAME L");
	legend.AddEntry(linFitCopy, linFitText.c_str() , "LPE");


	TF1* 									quadFitLine 								= 	new TF1((plotName+"_fitLinear").c_str(), "[0]+[1]*x+[2]*x*x", fitRange[0], fitRange[1]);
	TFitResultPtr 							quadFitResult 								=	isoPcVSpT->Fit(quadFitLine, options.getCSTR("fitOpt"));
	quadFitResult->Print("V");
	
	std::string 							quadFitText 								= 	"("
	+ to_string_with_precision(quadFitLine->GetParameter(2), options.getInt("c2Precision")) + "#pm" 
	+ to_string_with_precision(quadFitLine->GetParError(2), options.getInt("c2Precision")) + ")p_{T}^{2}+("
	+ to_string_with_precision(quadFitLine->GetParameter(1), options.getInt("c1Precision")) + "#pm"
	+ to_string_with_precision(quadFitLine->GetParError(1), options.getInt("c1Precision")) + ")p_{T}+("
	+ to_string_with_precision(quadFitLine->GetParameter(0), options.getInt("c0Precision")) + "#pm"
	+ to_string_with_precision(quadFitLine->GetParError(0), options.getInt("c0Precision")) + ")"
	+ " [#chi^{2}/ndf=" + to_string_with_precision(quadFitResult->Chi2 (), options.getInt("chi2precision")) + "/" + std::to_string(quadFitResult->Ndf()) + "]";
	
	std::string 							quadFitDef 									=	to_string_with_precision(quadFitLine->GetParameter(2), 10) + "*x*x +" + to_string_with_precision(quadFitLine->GetParameter(1), 10) + "*x + " + to_string_with_precision(quadFitLine->GetParameter(0), 10);
	TF1* 									quadFitCopy 									=  new TF1((plotName+"_fitQuadCopy").c_str(), quadFitDef.c_str(), histIsoVsPt->GetXaxis()->GetXmin(), histIsoVsPt->GetXaxis()->GetXmax());
	quadFitCopy->SetLineColor(options.getTColFromHex("quadCol"));
	quadFitCopy->SetLineWidth(options.getFloat("fitLnWidth"));
	quadFitCopy->SetLineStyle(options.getInt("fitLnStl"));
	quadFitCopy->Draw("A SAME L");
	legend.AddEntry(quadFitCopy, quadFitText.c_str() , "LPE");


	std::string outFileString = removeTrailingZeros(_etaMin) + ", " + removeTrailingZeros(_etaMax) + ", " +
	to_string_with_precision(linFitLine->GetParameter(0), 10) + ", " + to_string_with_precision(linFitLine->GetParError(0), 10) + ", " +
	to_string_with_precision(linFitLine->GetParameter(1), 10) + ", " + to_string_with_precision(linFitLine->GetParError(1), 10) + ", " +
	to_string_with_precision(quadFitLine->GetParameter(0), 10) + ", " + to_string_with_precision(quadFitLine->GetParError(0), 10) + ", " +
	to_string_with_precision(quadFitLine->GetParameter(1), 10) + ", " + to_string_with_precision(quadFitLine->GetParError(1), 10) + ", " +
	to_string_with_precision(quadFitLine->GetParameter(2), 10) + ", " + to_string_with_precision(quadFitLine->GetParError(2), 10);

	(*_pTscalingOutFile) << outFileString << endl;

	Float_t yMax = std::max(linFitLine->Eval(histIsoVsPt->GetXaxis()->GetBinUpEdge(histIsoVsPt->GetNbinsX())), quadFitLine->Eval(histIsoVsPt->GetXaxis()->GetBinUpEdge(histIsoVsPt->GetNbinsX())));

	histIsoVsPt->GetYaxis()->SetRangeUser(0., options.getFloat("yMaxRatio")*yMax);
	isoProj->GetXaxis()->SetRangeUser(0., options.getFloat("yMaxRatio")*yMax);
	pad0.cd();

	legend.Draw();
	
	padRwing.RedrawAxis();
	padRwing.RedrawAxis("G");
	
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
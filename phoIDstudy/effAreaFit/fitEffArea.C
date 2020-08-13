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
std::string 													optionsFile 			=	"fitEffAreaOptions.txt";
parseOptions 													options;
std::vector<std::string>										samples;
std::vector<std::string>										percentiles;
TH2F*															get2DHistogram(std::string _isoVar, Float_t _etaMin, Float_t _etaMax);
void 															initialize();
void 															getEffectiveArea(std::string _isoVar, Float_t _etaMin, Float_t _etaMax, ofstream* effAreaOutFile);

TH2F*															get2DHistogram(std::string _isoVar, Float_t _etaMin, Float_t _etaMax){

	std::string 								histName 									=	_isoVar + "_vs_rho_eta" + findAndReplaceAll(std::to_string(_etaMin) + "to" + std::to_string(_etaMax), ".", "p");
	std::string 								histTitle 									=	";"	+ varPlotInfo["rho"][0] + ";" + varPlotInfo[_isoVar][0];
	TH2F* 										histIsoVsRho								=	new TH2F(histName.c_str(), histTitle.c_str(), std::stoi(varPlotInfo["rho"][1]), std::stof(varPlotInfo["rho"][2]), std::stof(varPlotInfo["rho"][3]),
		std::stoi(varPlotInfo[_isoVar][1]), std::stof(varPlotInfo[_isoVar][2]), std::stof(varPlotInfo[_isoVar][3]));
	histIsoVsRho->SetContour(options.getInt("2dNcontours"));
	histIsoVsRho->Sumw2();

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
			histIsoVsRho->Fill(rho_, isolation_, weight);
		}
		closeTChain(tChain);
	}
	
	// histIsoVsRho->GetZaxis()->SetTitle(options.getCSTR("zTitle"));
	histIsoVsRho->GetXaxis()->CenterTitle();
	histIsoVsRho->GetYaxis()->CenterTitle();
	histIsoVsRho->GetZaxis()->CenterTitle();
	histIsoVsRho->GetZaxis()->RotateTitle();
	histIsoVsRho->GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histIsoVsRho->GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histIsoVsRho->GetZaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	histIsoVsRho->GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	histIsoVsRho->GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	histIsoVsRho->GetZaxis()->SetLabelSize(options.getDouble("pad0ZaxisLabelSize"));
	histIsoVsRho->GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	histIsoVsRho->GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	histIsoVsRho->GetZaxis()->SetTitleOffset(options.getDouble("pad0ZtitleOffset"));
	histIsoVsRho->GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	histIsoVsRho->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	return histIsoVsRho;
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
	percentiles																				=	options.getList("percentiles", ";");
	mkdir(options.get("writeDir"));
	
	std::vector<std::string> isoVars = options.getList("isoVar", ",");
	std::vector<std::string> etaBins = options.getList("etaBins", ";");
	
	for(std::string isoVar : isoVars){

		std::string effAreaFilePath = options.get("writeDir") + "/" +isoVar + ".txt";
		ofstream* effAreaOutFile = new ofstream(effAreaFilePath, std::ofstream::out);
		
		for(std::string iEtaBinStr : etaBins){
			
			std::vector<Float_t> iEtaMinMax = strToFloatList(iEtaBinStr, ",");
			
			getEffectiveArea(isoVar, iEtaMinMax[0], iEtaMinMax[1], effAreaOutFile);
			
		}
		
		effAreaOutFile->close();
		delete effAreaOutFile;
	}
	
};


void getEffectiveArea(std::string _isoVar, Float_t _etaMin, Float_t _etaMax, ofstream* effAreaOutFile){

	TH2F* 										histIsoVsRho								=	get2DHistogram(_isoVar, _etaMin, _etaMax);

	std::string 								plotName 									=	_isoVar + "_vs_rho_eta" + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p");
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

	TH1D* 										rhoProj 									= 	histIsoVsRho->ProjectionX((plotName+"_px").c_str());
	rhoProj->SetFillStyle(1001);
	rhoProj->SetFillColor(options.getTColFromHex("projXcol"));
	rhoProj->SetMinimum(0.);
	rhoProj->SetLineWidth(0.);
	rhoProj->GetXaxis()->SetTitle("");
	rhoProj->GetYaxis()->SetTitle("");
	rhoProj->GetXaxis()->SetLabelSize(0);
	rhoProj->GetYaxis()->SetLabelSize(0);
	rhoProj->SetStats(0);
	rhoProj->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));
	rhoProj->GetYaxis()->SetNdivisions(options.getInt("padUwingyNdivs"));
	rhoProj->Draw("HIST");
	// rhoProj->GetYaxis()->SetTickLength(0.);

	std::vector<Double_t> 						fitRange 									=	getNpercentMinInterval(rhoProj, 100. * options.getDouble("fitIntervalContent"));

	TArrow	rangeArrow(fitRange[0], options.getFloat("fitRangeArrowElevation")*rhoProj->GetMaximum(), fitRange[1], options.getFloat("fitRangeArrowElevation")*rhoProj->GetMaximum(),	options.getFloat("fitRangeArrowSize"), "<>");
	rangeArrow.SetLineStyle(options.getInt("fitRangeLineStyle"));
	rangeArrow.SetLineColor(options.getTColFromHex("fitRangeLineColor"));
	rangeArrow.SetLineWidth(options.getFloat("fitRangeLineWidth"));
	rangeArrow.SetFillColor(options.getTColFromHex("fitRangeLineColor"));
	rangeArrow.Draw();

	TLatex texWriter;
	texWriter.SetNDC(0);
	texWriter.SetTextAlign(21);
	std::string fitRangeInfo 									= 	"#bf{" + removeTrailingZeros(fitRange[0]) + " #leq #frac{#rho}{GeV} #leq " + removeTrailingZeros(fitRange[1]) + "} " + options.get("fitIntervalInfo");
	texWriter.SetTextSize(options.getFloat("fitRangeInfoTextSize"));
	texWriter.SetTextColor(options.getTColFromHex("fitRangeInfoTextColor"));
	Float_t fitRangeXpos = fitRange[0] + options.getFloat("fitRangeInfoDisplacement")*(fitRange[1] - fitRange[0]);
	texWriter.DrawLatex(fitRangeXpos, options.getFloat("fitRangeInfoElevation")*rhoProj->GetMaximum(), fitRangeInfo.c_str());

	padUwing.RedrawAxis();
	padUwing.RedrawAxis("G");
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");

	canvas.cd();
	texWriter.SetNDC(0);
	texWriter.SetTextAlign(11);
	texWriter.SetTextSize(options.getFloat("categoryInfoTextSize"));
	texWriter.SetTextColor(options.getTColFromHex("categoryInfoTextColor"));
	std::string 								etaTitle 									=	removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax);
	texWriter.DrawLatex(options.getFloat("categoryInfoX"), options.getFloat("categoryInfoY"), etaTitle.c_str());

	//// Y projection
	TPad padRwing("padRwing", "", options.getFloat("padRwingx1"), options.getFloat("padRwingy1"), options.getFloat("padRwingx2"), options.getFloat("padRwingy2"));
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

	TH1D* 										isoProj 									=	histIsoVsRho->ProjectionY((plotName+"_py").c_str());
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
	histIsoVsRho->Draw("colz");

	TLegend 									legend(options.getDouble("legx1"), options.getDouble("legy1"), options.getDouble("legx2"), options.getDouble("legy2"));
	legend.SetTextSize(options.getDouble("legTextSize"));
	legend.SetNColumns(options.getInt("legNcols"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetLineColor(options.getTColFromHex("legFillColor"));
	legend.SetBorderSize(options.getInt("legBorderWidth"));

	if(options.getInt("autoZrange")) histIsoVsRho->GetZaxis()->SetRangeUser(histIsoVsRho->GetMinimum(std::numeric_limits<Double_t>::min()) / options.getFloat("zMinRatio"), histIsoVsRho->GetMaximum() * options.getFloat("zMaxRatio"));

	Float_t yMax = -999.;
	for(std::string iPcStr : percentiles){
		
		Float_t 								iPercentile 								= 	std::stof(split_string(iPcStr, ",")[0])/100.;
		
		TH1D* 									iPcVSrho 									= 	histIsoVsRho->QuantilesX(iPercentile, (plotName+"_"+split_string(iPcStr, ",")[0]).c_str());
		iPcVSrho->SetLineColor(hex2rootColor(split_string(iPcStr, ",")[1]));
		iPcVSrho->SetLineWidth(options.getFloat("pcWidth"));
		iPcVSrho->SetMarkerColor(hex2rootColor(split_string(iPcStr, ",")[1]));
		iPcVSrho->SetMarkerSize(options.getFloat("pcMkrSz"));
		iPcVSrho->SetMarkerStyle(options.getInt("pcMkrStl"));

		TF1* 									iFitLine 									= 	new TF1((plotName+"_"+split_string(iPcStr, ",")[0]+"_fit").c_str(), "[0]*x+[1]", fitRange[0], fitRange[1]);
		
		TFitResultPtr 							iFitResult 									=	iPcVSrho->Fit(iFitLine, options.getCSTR("fitOpt"));
		iFitResult->Print("V");
		
		std::string 							iFitText 									= 	split_string(iPcStr, ",")[0] + "%: (" +
		to_string_with_precision(iFitLine->GetParameter(0), options.getInt("slopePrecision")) + "#pm" 
		+ to_string_with_precision(iFitLine->GetParError(0), options.getInt("slopePrecision")) + ")#rho + ("
		+ to_string_with_precision(iFitLine->GetParameter(1), options.getInt("interceptPrecision")) + "#pm"
		+ to_string_with_precision(iFitLine->GetParError(1), options.getInt("interceptPrecision")) + ")";
		+ " [#chi^{2}/ndf=" + to_string_with_precision(iFitResult->Chi2 (), options.getInt("chi2precision")) + "/" + std::to_string(iFitResult->Ndf()) + "]";
		
		std::string 							iLineDef 									=	to_string_with_precision(iFitLine->GetParameter(0), 10) + "*x + " + to_string_with_precision(iFitLine->GetParameter(1), 10);
		TF1* 									iLineCopy 									=  new TF1((plotName+"_"+split_string(iPcStr, ",")[0]+"_fit").c_str(), iLineDef.c_str(), histIsoVsRho->GetXaxis()->GetXmin(), histIsoVsRho->GetXaxis()->GetXmax());
		iLineCopy->SetLineColor(hex2rootColor(split_string(iPcStr, ",")[1]));
		iLineCopy->SetLineWidth(options.getFloat("fitLnWidth"));
		iLineCopy->SetLineStyle(options.getInt("fitLnStl"));
		
		legend.AddEntry(iPcVSrho, iFitText.c_str() , "LPE");
		
		iPcVSrho->Draw("A SAME PE");
		iLineCopy->Draw("A SAME L");

		std::string outFileString = split_string(iPcStr, ",")[0] + ", " + removeTrailingZeros(_etaMin) + ", " + removeTrailingZeros(_etaMax) + ", " + removeTrailingZeros(iFitLine->GetParameter(0));

		(*effAreaOutFile) << outFileString << endl;
		
		Float_t iFitLineYrange = iLineCopy->Eval(histIsoVsRho->GetXaxis()->GetBinUpEdge(histIsoVsRho->GetNbinsX()));
		
		if(iFitLineYrange > yMax){
			yMax = iFitLineYrange;
		}
		
		delete iFitLine;
		
	}
	
	histIsoVsRho->GetYaxis()->SetRangeUser(0., options.getFloat("yMaxRatio")*yMax);
	isoProj->GetXaxis()->SetRangeUser(0., options.getFloat("yMaxRatio")*yMax);
	pad0.cd();

	TLine fRa1ngeX1(fitRange[0], 0. ,fitRange[0], histIsoVsRho->GetYaxis()->GetXmax());
	fRa1ngeX1.SetLineStyle(options.getInt("fitRangeLineStyle0"));
	fRa1ngeX1.SetLineColor(options.getTColFromHex("fitRangeLineColor"));
	fRa1ngeX1.SetLineWidth(options.getFloat("fitRangeLineWidth0"));
	fRa1ngeX1.Draw();

	TLine fRa1ngeX2(fitRange[1], 0. ,fitRange[1], histIsoVsRho->GetYaxis()->GetXmax());
	fRa1ngeX2.SetLineStyle(options.getInt("fitRangeLineStyle0"));
	fRa1ngeX2.SetLineColor(options.getTColFromHex("fitRangeLineColor"));
	fRa1ngeX2.SetLineWidth(options.getFloat("fitRangeLineWidth0"));
	fRa1ngeX2.Draw();

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
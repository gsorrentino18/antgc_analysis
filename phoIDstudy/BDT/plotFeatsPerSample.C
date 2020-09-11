#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"

std::string 										optionsFile="plotFeatsPerSampleOptions.txt";
parseOptions 										options;
std::map<std::string,std::vector<std::string>>		varPlotInfo;
void 												makeFeatPlots();
void 												plotVar(std::string _plotVar);
Bool_t 												plotTest;
Bool_t 												plotTrain;
Bool_t 												plotValidation;

void plotVar(std::string _plotVar, Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax){

	if ( varPlotInfo.find(_plotVar) == varPlotInfo.end() ) {
		std::cout<<"Error! Variable "<<_plotVar<<" not found in file "<<options.get("inFile")<<std::endl;
		return;
	}

	UChar_t 									sSample 								=		options.getInt("sSample");
	UChar_t 									bSample 								=		options.getInt("bSample");
	std::string 								plotName 								=		options.get("sSample")+ "_" + options.get("bSample" ) + "_"+ _plotVar + "_eta_" + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p") + "_pT_" + findAndReplaceAll(removeTrailingZeros(_pTmin) + "to" + removeTrailingZeros(_pTmax), ".", "p");
	TChain*                                     featsTree      							=		openTChain(std::vector<std::string>({"/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root"}), "fullEB_Tree");
	TChain*                                     isoTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/isoShuffledTree.root"}), "fullEB_isoTree");
	TChain*                                     bdtTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/optimizedV1/BDTresults_0.root"}), "fullEB_BDT_Tree");
	featsTree->AddFriend(isoTree);
	featsTree->AddFriend(bdtTree);
	TTreeReader                             	inputTTreeReader(featsTree);
	TTreeReaderAnyValue<Bool_t>					isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Bool_t>					isTrain_										(inputTTreeReader, "isTrain");
	TTreeReaderAnyValue<Bool_t>					isValidation_									(inputTTreeReader, "isValidation");
	TTreeReaderAnyValue<Float_t>				phoPt_											(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Double_t>				xSecW_										(inputTTreeReader, "xSecW");
	TTreeReaderAnyValue<Float_t>				feat_											(inputTTreeReader, _plotVar);
	TTreeReaderAnyValue<UChar_t>				sampleIndex_									(inputTTreeReader, "sampleIndex");



	TH1D 										signalHist((plotName+"_signal").c_str(), "", std::stoi(varPlotInfo[_plotVar][1]), std::stof(varPlotInfo[_plotVar][2]), std::stof(varPlotInfo[_plotVar][3]));
	TH1D 										backgroundHist((plotName+"_background").c_str(), "", std::stoi(varPlotInfo[_plotVar][1]), std::stof(varPlotInfo[_plotVar][2]), std::stof(varPlotInfo[_plotVar][3]));

	signalHist.SetLineStyle(options.getInt("signalLineStyle"));
	signalHist.SetLineWidth(options.getInt("signalLineWidth"));
	signalHist.SetLineColor(options.getTColFromHex("signalColor"));
	backgroundHist.SetLineStyle(options.getInt("backgroundLineStyle"));
	backgroundHist.SetLineWidth(options.getInt("backgroundLineWidth"));
	backgroundHist.SetLineColor(options.getTColFromHex("backgroundColor"));	
	
	ULong64_t iEntry =0;
	while(inputTTreeReader.Next()){

		if(iEntry % 1000000 == 0)std::cout<<"Entry "<<iEntry<<std::endl;
		iEntry++;

		

		if(phoPt_ < _pTmin) continue;
		if(phoPt_ >= _pTmax) continue;
		Float_t 								absSCEta									=	std::abs(phoSCeta_);
		if(absSCEta <= _etaMin) continue;
		if(absSCEta > _etaMax) continue;

		Bool_t includePoint = (plotTest && !isTrain_) || (plotTrain && !isValidation_ && isTrain_) || (plotValidation && isValidation_);
		
		if(isSignal_){
			if(sSample != sampleIndex_) continue;
			signalHist.Fill(feat_, xSecW_);
		} else{
			if(bSample != sampleIndex_) continue;
			backgroundHist.Fill(feat_, xSecW_);
		}
	}
	
	normalizeHist(signalHist, 100.);
	normalizeHist(backgroundHist, 100.);

	TCanvas canvas((plotName + "_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);
	
	TPad pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad0.SetFillStyle(4000);
	pad0.SetFillColor(0);
	pad0.SetFrameFillStyle(4000);
	pad0.SetGrid(1,1);
	
	TLegend legend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	legend.SetTextSize(options.getDouble("legendTextSize"));
	legend.SetNColumns(options.getInt("legendNcols"));
	legend.SetBorderSize(0);
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	
	THStack 								hStack((plotName + "_stack").c_str(), "");
	hStack.Add(&signalHist, "HIST");
	hStack.Add(&backgroundHist, "HIST");

	canvas.Draw();
	canvas.cd();
	pad0.Draw();
	pad0.cd();
	
	hStack.Draw("NOSTACK HIST");
	
	hStack.GetXaxis()->SetTitle(varPlotInfo[_plotVar][0].c_str());
	hStack.GetXaxis()->CenterTitle();
	hStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	hStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	hStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));

	hStack.GetXaxis()->SetNoExponent(1);
	
	hStack.GetYaxis()->SetTitle(options.get("yTitle").c_str());
	hStack.GetYaxis()->CenterTitle();
	hStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	hStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	hStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));
	
	hStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	hStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	hStack.SetMaximum(options.getFloat("yScaleup")* hStack.GetMaximum("nostack"));
	
	if(options.getInt("logX") == 1) pad0.SetLogx();
	
	std::string 		legEn = "Prompt " + findAndReplaceAll(getFileName(vLookup(options.get("sSample"), options.get("sampleMap"), 0, 1)), ".root", "");
	legend.AddEntry(&signalHist, legEn.c_str(), "L");

	legEn = "Fake " + findAndReplaceAll(getFileName(vLookup(options.get("bSample"), options.get("sampleMap"), 0, 1)), ".root", "");
	legend.AddEntry(&backgroundHist, legEn.c_str(), "L");

	

	std::string 								legTitle 									=	removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax) + "  " + removeTrailingZeros(_pTmin) + "<p_{T}#leq" + removeTrailingZeros(_pTmax);
	legend.SetHeader(legTitle.c_str(), "C");
	legend.SetTextAlign(12);
	legend.Draw();
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	
	canvas.SaveAs((options.get("writeDir")+plotName+".png").c_str());
	canvas.SaveAs((options.get("writeDir")+plotName+".pdf").c_str());

	pad0.SetLogy();
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "/logY/" +plotName+".png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "/logY/"+plotName+".pdf").c_str());
	
	clearHeap();
};





void makeFeatPlots(){
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	options.parseIt(optionsFile, "==");

	CSVReader 									varFile(options.get("prettyNamesFile"), "==");
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
	
	std::vector<std::string> 					vars 										=	options.getList("vars", ",");
	std::vector<std::string> 					etaBins 									=	options.getList("etaBins", ";");
	std::vector<std::string> 					pTbins										= 	options.getList("pTbins", ";");
	
	plotTest 																				=	options.getInt("plotTest");
	plotTrain 																				=	options.getInt("plotTrain");
	plotValidation 																			=	options.getInt("plotValidation");

	mkdir(options.get("writeDir")+ "/logY/");

	for(std::string iVar : vars){

		std::cout<<"\n\nFeature: "<<iVar<<std::endl;

		for(std::string iPtBinStr : pTbins){
			std::vector<Float_t> iPtMinMax = strToFloatList(iPtBinStr, ",");
			std::cout<<"\tPt: "<<iPtMinMax[0]<<"-"<<iPtMinMax[1]<<std::endl;
			for(std::string iEtaBinStr : etaBins){
				std::vector<Float_t> iEtaMinMax = strToFloatList(iEtaBinStr, ",");
				std::cout<<"\t\tEta: "<<iEtaMinMax[0]<<"-"<<iEtaMinMax[1]<<std::endl;
				plotVar(iVar, iEtaMinMax[0], iEtaMinMax[1], iPtMinMax[0], iPtMinMax[1]);
			}
		}
	}
};
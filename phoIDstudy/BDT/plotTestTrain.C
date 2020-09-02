#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"

std::string 										optionsFile="plotTestTrainOptions.txt";
parseOptions 										options;
std::map<std::string,std::vector<std::string>>		varPlotInfo;
void 												initialize();
void 												plotBDT(Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax);
Bool_t 												plotTest;
Bool_t 												plotTrain;
Bool_t 												plotValidation;

void plotBDT(Float_t _etaMin, Float_t _etaMax, Float_t _pTmin, Float_t _pTmax){

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
	TTreeReaderAnyValue<Double_t>				bdtWeight_										(inputTTreeReader, "bdtWeightF");
	TTreeReaderAnyValue<Double_t>				bdtScore_										(inputTTreeReader, "bdtScore");

	TH1D 										testSignalHist("test_signal", "", std::stoi(varPlotInfo["bdtScore"][1]), std::stof(varPlotInfo["bdtScore"][2]), std::stof(varPlotInfo["bdtScore"][3]));
	TH1D 										trainSignalHist("test_signal", "", std::stoi(varPlotInfo["bdtScore"][1]), std::stof(varPlotInfo["bdtScore"][2]), std::stof(varPlotInfo["bdtScore"][3]));

	TH1D 										testBackgroundHist("test_background", "", std::stoi(varPlotInfo["bdtScore"][1]), std::stof(varPlotInfo["bdtScore"][2]), std::stof(varPlotInfo["bdtScore"][3]));
	TH1D 										trainBackgroundHist("train_background", "", std::stoi(varPlotInfo["bdtScore"][1]), std::stof(varPlotInfo["bdtScore"][2]), std::stof(varPlotInfo["bdtScore"][3]));

	testSignalHist.SetLineStyle(options.getInt("testLineStyle"));
	trainSignalHist.SetLineStyle(options.getInt("trainLineStyle"));
	testBackgroundHist.SetLineStyle(options.getInt("testLineStyle"));
	trainBackgroundHist.SetLineStyle(options.getInt("trainLineStyle"));

	testSignalHist.SetLineWidth(options.getInt("testLineWidth"));
	trainSignalHist.SetLineWidth(options.getInt("trainLineWidth"));
	testBackgroundHist.SetLineWidth(options.getInt("testLineWidth"));
	trainBackgroundHist.SetLineWidth(options.getInt("trainLineWidth"));
	
	testSignalHist.SetLineColor(options.getTColFromHex("testSignalColor"));
	trainSignalHist.SetLineColor(options.getTColFromHex("trainSignalColor"));
	testBackgroundHist.SetLineColor(options.getTColFromHex("testBackgroundColor"));
	trainBackgroundHist.SetLineColor(options.getTColFromHex("trainBackgroundColor"));

	testSignalHist.SetMarkerColor(options.getTColFromHex("testSignalColor"));
	trainSignalHist.SetMarkerColor(options.getTColFromHex("trainSignalColor"));
	testBackgroundHist.SetMarkerColor(options.getTColFromHex("testBackgroundColor"));
	trainBackgroundHist.SetMarkerColor(options.getTColFromHex("trainBackgroundColor"));

	testSignalHist.SetMarkerStyle(options.getInt("testSignalMkrStyle"));
	testBackgroundHist.SetMarkerStyle(options.getInt("testBackgroundMkrStyle"));

	testSignalHist.SetMarkerSize(options.getFloat("mkrSize"));
	testSignalHist.SetMarkerSize(options.getFloat("mkrSize"));

	ULong64_t iEntry =0;
	while(inputTTreeReader.Next()){

		if(iEntry % 1000000 == 0)std::cout<<"Entry "<<iEntry<<std::endl;
		iEntry++;

		if(phoPt_ < _pTmin) continue;
		if(phoPt_ >= _pTmax) continue;
		Float_t 								absSCEta									=	std::abs(phoSCeta_);
		if(absSCEta <= _etaMin) continue;
		if(absSCEta > _etaMax) continue;

		if(isTrain_ && !isValidation_){
			if(isSignal_) trainSignalHist.Fill(bdtScore_, bdtWeight_);
			else trainBackgroundHist.Fill(bdtScore_, bdtWeight_);
		} else if(!isTrain_){
			if(isSignal_) testSignalHist.Fill(bdtScore_, bdtWeight_);
			else testBackgroundHist.Fill(bdtScore_, bdtWeight_);
		}
	}
	
	normalizeHist(trainSignalHist, 100.);
	normalizeHist(testSignalHist, 100.);
	normalizeHist(trainBackgroundHist, 100.);
	normalizeHist(testBackgroundHist, 100.);

	TCanvas canvas("BDT_canvas", "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);
	
	TPad pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad0.SetFillStyle(4000);
	pad0.SetFillColor(0);
	pad0.SetFrameFillStyle(4000);
	pad0.SetGrid(1,1);
	
	TLegend sLegend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	sLegend.SetTextSize(options.getDouble("legendTextSize"));
	sLegend.SetNColumns(options.getInt("legendNcols"));
	sLegend.SetBorderSize(0);
	sLegend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	sLegend.SetFillStyle(options.getInt("legFillStyle"));
	
	THStack 								sStack("s_bdt_stack", "");
	sStack.Add(&trainSignalHist, "HIST");
	sStack.Add(&testSignalHist, "PE");

	canvas.Draw();
	canvas.cd();
	pad0.Draw();
	pad0.cd();
	
	sStack.Draw("NOSTACK HIST");
	
	sStack.GetXaxis()->SetTitle(varPlotInfo["bdtScore"][0].c_str());
	sStack.GetXaxis()->CenterTitle();
	sStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	sStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	sStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	sStack.GetXaxis()->SetNoExponent(1);
	sStack.GetYaxis()->SetTitle(options.get("yTitle").c_str());
	sStack.GetYaxis()->CenterTitle();
	sStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	sStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	sStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));
	sStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	sStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	sStack.SetMaximum(options.getFloat("yScaleup")* sStack.GetMaximum("nostack"));
	
	if(options.getInt("logX") == 1) pad0.SetLogx();
	
	std::string legTxt 		= "Prompt Train (" + to_string_with_precision(trainSignalHist.GetEntries(),0)+")";
	sLegend.AddEntry(&trainSignalHist, legTxt.c_str(), "L");
	legTxt 		= "Prompt Test (" + to_string_with_precision(testSignalHist.GetEntries(),0)+")";
	sLegend.AddEntry(&testSignalHist, legTxt.c_str(), "LPE");

	std::string 								legTitle 									=	removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax) + "  " + removeTrailingZeros(_pTmin) + "<p_{T}#leq" + removeTrailingZeros(_pTmax);
	sLegend.SetHeader(legTitle.c_str(), "C");
	sLegend.SetTextAlign(12);
	sLegend.Draw();
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();
	
	canvas.SaveAs((options.get("writeDir")+"signal_BDTscore.png").c_str());
	canvas.SaveAs((options.get("writeDir")+"signal_BDTscore.pdf").c_str());

	pad0.SetLogy();
	sStack.SetMaximum(options.getFloat("yScaleup")* sStack.GetMaximum("yScaleupLog"));
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "logY/signal_BDTscore.png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "logY/signal_BDTscore.pdf").c_str());


	pad0.Clear();
	pad0.SetLogy(0);

	THStack 								bStack("b_bdt_stack", "");
	bStack.Add(&trainBackgroundHist, "HIST");
	bStack.Add(&testBackgroundHist, "PE");

	TLegend bLegend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	bLegend.SetTextSize(options.getDouble("legendTextSize"));
	bLegend.SetNColumns(options.getInt("legendNcols"));
	bLegend.SetBorderSize(0);
	bLegend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	bLegend.SetFillStyle(options.getInt("legFillStyle"));

	legTxt 		= "Fake Train (" + to_string_with_precision(trainBackgroundHist.GetEntries(),0)+")";
	bLegend.AddEntry(&trainBackgroundHist, legTxt.c_str(), "L");
	legTxt 		= "Fake Test (" + to_string_with_precision(testBackgroundHist.GetEntries(),0)+")";
	bLegend.AddEntry(&testBackgroundHist, legTxt.c_str(), "LPE");
	bLegend.SetHeader(legTitle.c_str(), "C");
	bLegend.SetTextAlign(12);

	pad0.cd();
	bStack.Draw("NOSTACK HIST");
	bLegend.Draw();

	bStack.GetXaxis()->SetTitle(varPlotInfo["bdtScore"][0].c_str());
	bStack.GetXaxis()->CenterTitle();
	bStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	bStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	bStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	bStack.GetXaxis()->SetNoExponent(1);
	bStack.GetYaxis()->SetTitle(options.getCSTR("yTitle"));
	bStack.GetYaxis()->CenterTitle();
	bStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	bStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	bStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));
	bStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	bStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	bStack.SetMaximum(options.getFloat("yScaleup")* bStack.GetMaximum("nostack"));
	
	if(options.getInt("logX") == 1) pad0.SetLogx();

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();
	
	canvas.SaveAs((options.get("writeDir")+"background_BDTscore.png").c_str());
	canvas.SaveAs((options.get("writeDir")+"background_BDTscore.pdf").c_str());

	pad0.SetLogy();
	bStack.SetMaximum(options.getFloat("yScaleup")* sStack.GetMaximum("yScaleupLog"));
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "logY/background_BDTscore.png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "logY/background_BDTscore.pdf").c_str());

	clearHeap();
};





void initialize(){
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
	
	std::vector<std::string> 					etaBins 									=	options.getList("etaBins", ";");
	std::vector<std::string> 					pTbins										= 	options.getList("pTbins", ";");
	
	mkdir(options.get("writeDir")+ "/logY/");

	for(std::string iPtBinStr : pTbins){
		std::vector<Float_t> iPtMinMax = strToFloatList(iPtBinStr, ",");
		std::cout<<"\tPt: "<<iPtMinMax[0]<<"-"<<iPtMinMax[1]<<std::endl;
		for(std::string iEtaBinStr : etaBins){
			std::vector<Float_t> iEtaMinMax = strToFloatList(iEtaBinStr, ",");
			std::cout<<"\t\tEta: "<<iEtaMinMax[0]<<"-"<<iEtaMinMax[1]<<std::endl;
			plotBDT(iEtaMinMax[0], iEtaMinMax[1], iPtMinMax[0], iPtMinMax[1]);
		}
	}
};
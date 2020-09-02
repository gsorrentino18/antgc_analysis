#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"

std::string 										optionsFile="plotPtEtaOptions.txt";
parseOptions 										options;
std::map<std::string,std::vector<std::string>>		varPlotInfo;
void 												initialize();
void 												plotPtEta();


void plotPtEta(){

	std::vector<Double_t> 						pTbinning										= 	options.getDoubleList("pTbinning");
	std::vector<Double_t> 						etaBinning										= 	options.getDoubleList("etaBinning");

	Float_t 									pTmin = pTbinning[0];
	Float_t 									pTmax = pTbinning.back();

	Float_t 									etaMin = etaBinning[0];
	Float_t 									etaMax = etaBinning.back();

	TChain*                                     featsTree      							=		openTChain(std::vector<std::string>({"/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root"}), "fullEB_Tree");
	TChain*                                     isoTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/isoShuffledTree.root"}), "fullEB_isoTree");
	TChain*                                     bdtTree         						=		openTChain(std::vector<std::string>({"/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/optimizedV1/BDTresults_0.root"}), "fullEB_BDT_Tree");
	featsTree->AddFriend(isoTree);
	featsTree->AddFriend(bdtTree);
	TTreeReader                             	inputTTreeReader(featsTree);
	TTreeReaderAnyValue<Bool_t>					isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Float_t>				phoPt_											(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Float_t>				phoEta_											(inputTTreeReader, "phoEta");
	TTreeReaderAnyValue<Double_t>				bdtWeight_										(inputTTreeReader, "bdtWeightF");
	TTreeReaderAnyValue<Float_t>				xSecW_											(inputTTreeReader, "xSecW");


	TH1D 										signalPtHist("signalPtHist", "", pTbinning.size()-1, pTbinning.data());
	TH1D 										backgroundPtHist("backgroundPtHist", "", pTbinning.size()-1, pTbinning.data());
	TH1D 										backgroundPtRwHist("backgroundPtRwHist", "", pTbinning.size()-1, pTbinning.data());

	TH1D 										signalEtaHist("signalEtaHist", "", etaBinning.size()-1, etaBinning.data());
	TH1D 										backgroundEtaHist("backgroundEtaHist", "", etaBinning.size()-1, etaBinning.data());
	TH1D 										backgroundEtaRwHist("backgroundEtaRwHist", "", etaBinning.size()-1, etaBinning.data());

	signalPtHist.SetLineStyle(options.getInt("lineStyle"));
	backgroundPtHist.SetLineStyle(options.getInt("lineStyle"));
	backgroundPtRwHist.SetLineStyle(options.getInt("lineStyle"));
	signalEtaHist.SetLineStyle(options.getInt("lineStyle"));
	backgroundEtaHist.SetLineStyle(options.getInt("lineStyle"));
	backgroundEtaRwHist.SetLineStyle(options.getInt("lineStyle"));
	
	signalPtHist.SetLineColor(options.getTColFromHex("signalColor"));
	backgroundPtHist.SetLineColor(options.getTColFromHex("backgroundColor"));
	backgroundPtRwHist.SetLineColor(options.getTColFromHex("backgroundRwColor"));
	signalEtaHist.SetLineColor(options.getTColFromHex("signalColor"));
	backgroundEtaHist.SetLineColor(options.getTColFromHex("backgroundColor"));
	backgroundEtaRwHist.SetLineColor(options.getTColFromHex("backgroundRwColor"));

	signalPtHist.SetLineWidth(options.getFloat("lineWidth"));
	backgroundPtHist.SetLineWidth(options.getFloat("lineWidth"));
	backgroundPtRwHist.SetLineWidth(options.getFloat("lineWidth"));
	signalEtaHist.SetLineWidth(options.getFloat("lineWidth"));
	backgroundEtaHist.SetLineWidth(options.getFloat("lineWidth"));
	backgroundEtaRwHist.SetLineWidth(options.getFloat("lineWidth"));
	
	signalPtHist.SetMarkerColor(options.getTColFromHex("signalColor"));
	backgroundPtHist.SetMarkerColor(options.getTColFromHex("backgroundColor"));
	backgroundPtRwHist.SetMarkerColor(options.getTColFromHex("backgroundRwColor"));
	signalEtaHist.SetMarkerColor(options.getTColFromHex("signalColor"));
	backgroundEtaHist.SetMarkerColor(options.getTColFromHex("backgroundColor"));
	backgroundEtaRwHist.SetMarkerColor(options.getTColFromHex("backgroundRwColor"));

	backgroundPtRwHist.SetMarkerStyle(options.getInt("backgroundRwMkrStyle"));
	backgroundEtaRwHist.SetMarkerStyle(options.getInt("backgroundRwMkrStyle"));

	backgroundPtRwHist.SetMarkerSize(options.getFloat("mkrSize"));
	backgroundEtaRwHist.SetMarkerSize(options.getFloat("mkrSize"));

	ULong64_t iEntry =0;
	while(inputTTreeReader.Next()){

		if(iEntry % 1000000 == 0)std::cout<<"Entry "<<iEntry<<std::endl;
		iEntry++;

		if(phoPt_ < pTmin) continue;
		if(phoPt_ >= pTmax) continue;
		if(phoEta_ < etaMin) continue;
		if(phoEta_ > etaMax) continue;

		if(isSignal_){
			signalPtHist.Fill(phoPt_, xSecW_);
			signalEtaHist.Fill(phoEta_, xSecW_);
		} else{

			backgroundPtHist.Fill(phoPt_, xSecW_);
			backgroundEtaHist.Fill(phoEta_, xSecW_);

			backgroundPtRwHist.Fill(phoPt_, bdtWeight_);
			backgroundEtaRwHist.Fill(phoEta_, bdtWeight_);
		}
	}
	
	normalizeHist(signalPtHist, 100.);
	normalizeHist(signalEtaHist, 100.);
	normalizeHist(backgroundPtHist, 100.);
	normalizeHist(backgroundEtaHist, 100.);
	normalizeHist(backgroundPtRwHist, 100.);
	normalizeHist(backgroundEtaRwHist, 100.);

	TCanvas canvas("BDT_canvas", "", options.getDouble("canvasX"), options.getDouble("canvasY"));
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
	
	THStack 								ptStack("s_bdt_stack", "");
	ptStack.Add(&signalPtHist, "HIST");
	ptStack.Add(&backgroundPtHist, "HIST");
	ptStack.Add(&backgroundPtRwHist, "P");

	canvas.Draw();
	canvas.cd();
	pad0.Draw();
	pad0.cd();
	
	ptStack.Draw("NOSTACK HIST");
	
	ptStack.GetXaxis()->SetTitle("p_{T} (GeV)");
	ptStack.GetXaxis()->CenterTitle();
	ptStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	ptStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	ptStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	ptStack.GetXaxis()->SetNoExponent(1);
	ptStack.GetYaxis()->SetTitle(options.get("yTitle").c_str());
	ptStack.GetYaxis()->CenterTitle();
	ptStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	ptStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	ptStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));
	ptStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	ptStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	ptStack.SetMaximum(options.getFloat("yScaleup")* ptStack.GetMaximum("nostack"));
	
	if(options.getInt("logX") == 1) pad0.SetLogx();
	
	std::string legTxt 		= "Prompt (" + to_string_with_precision(signalPtHist.GetEntries(),0)+")";
	legend.AddEntry(&signalPtHist, legTxt.c_str(), "L");
	legTxt 		= "Fake (" + to_string_with_precision(backgroundPtHist.GetEntries(),0)+")";
	legend.AddEntry(&backgroundPtHist, legTxt.c_str(), "L");
	legTxt 		= "Fake Reweighted (" + to_string_with_precision(backgroundPtRwHist.GetEntries(),0)+")";
	legend.AddEntry(&backgroundPtRwHist, legTxt.c_str(), "LPE");
	legend.Draw();
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();
	
	canvas.SaveAs((options.get("writeDir")+"ptRweighting.png").c_str());
	canvas.SaveAs((options.get("writeDir")+"ptRweighting.pdf").c_str());

	pad0.SetLogy();
	ptStack.SetMaximum(options.getFloat("yScaleup")* ptStack.GetMaximum("yScaleupLog"));
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "logY/ptRweighting.png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "logY/ptRweighting.pdf").c_str());


	pad0.Clear();
	pad0.SetLogy(0);

	THStack 								bStack("b_bdt_stack", "");
	bStack.Add(&signalEtaHist, "HIST");
	bStack.Add(&backgroundEtaHist, "HIST");
	bStack.Add(&backgroundEtaRwHist, "P");


	pad0.cd();
	bStack.Draw("NOSTACK HIST");
	legend.Draw();

	bStack.GetXaxis()->SetTitle("#eta");
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
	
	canvas.SaveAs((options.get("writeDir")+"etaReweighting.png").c_str());
	canvas.SaveAs((options.get("writeDir")+"etaReweighting.pdf").c_str());

	pad0.SetLogy();
	bStack.SetMaximum(options.getFloat("yScaleup")* ptStack.GetMaximum("yScaleupLog"));
	
	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "logY/etaReweighting.png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "logY/etaReweighting.pdf").c_str());

	clearHeap();
};





void initialize(){
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	options.parseIt(optionsFile, "==");
	mkdir(options.get("writeDir")+ "/logY/");
	plotPtEta();
};
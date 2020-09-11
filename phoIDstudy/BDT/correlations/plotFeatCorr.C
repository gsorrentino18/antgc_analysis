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

std::string 													optionsFile 			=	"plotFeatCorrOptions.txt";
parseOptions 													options;

std::vector<std::string>					feats;

TH2F* getHist(){
	Bool_t 										plotSignal=options.getInt("plotSignal");
	Float_t 									pTmin 	= options.getFloatList("pTrange")[0];
	Float_t 									pTmax 	= options.getFloatList("pTrange")[1];
	Float_t 									etaMin 	= options.getFloatList("etaRange")[0];
	Float_t 									etaMax 	= options.getFloatList("etaRange")[1];

	TChain*										fTree 										=	openTChain((std::vector<std::string>){options.get("featsFile")}, options.get("featsTree"));
	TChain*										iTree 										=	openTChain((std::vector<std::string>){options.get("indexFile")}, options.get("indexTree"));
	fTree->AddFriend(iTree);
	TTreeReader                             	inputTTreeReader(fTree);
	TTreeReaderAnyValue<Bool_t> 				isSignal 										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Bool_t> 				isTrain 										(inputTTreeReader, "isTrain");
	TTreeReaderAnyValue<Double_t>				flatPtEtaRwNoXsec								(inputTTreeReader, "flatPtEtaRwNoXsec");
	TTreeReaderAnyValue<Float_t>				phoSCeta										(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<Float_t>				phoPt											(inputTTreeReader, "phoPt");
	
	std::vector<TTreeReaderAnyValue<Float_t>*> 	featBranches;
	for(const std::string & iFeat : feats){
		featBranches.push_back(new TTreeReaderAnyValue<Float_t>(inputTTreeReader, iFeat));
	}

	correlationMatix 							fCorr(featBranches.size());
	
	while(inputTTreeReader.Next()){
		if(!isTrain) continue;
		if(plotSignal != isSignal) continue;
		if(phoPt < pTmin) continue;
		if(phoPt >= pTmax) continue;
		Float_t 								absSCEta									=	std::abs(phoSCeta);
		if (absSCEta <= etaMin) continue;
		if (absSCEta > etaMax) continue;
		std::vector<Double_t> featVals;
		for(TTreeReaderAnyValue<Float_t> * iFeatBr : featBranches){
			featVals.push_back(iFeatBr->get());
		}
		
		fCorr.addPoint(featVals, flatPtEtaRwNoXsec);
	}
	
	TH2F* 										corrMatHist 								=	fCorr.getCorrelationHist();
	
	for(TTreeReaderAnyValue<Float_t>* iFeat : featBranches){
		delete iFeat;
	}

	closeTChain(fTree);
	closeTChain(iTree);

	return corrMatHist;
};

void plotCorrelations() {
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	options.parseIt(optionsFile, "==");
	gStyle->SetPalette(options.getInt("colorPalette"));

	gStyle->SetPaintTextFormat(options.getCSTR("decPlaces"));

	std::map<std::string,std::vector<std::string>>					varPlotInfo;
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

	if(options.getInt("plotSignal")){
		feats = options.getList("sigFeatList");

	} else {
		feats = options.getList("bgFeatList");
	}
	
	TH2F* 					corrMatHist;

	std::string 								writePath 									= 	options.get("writeDir") + "/" + options.get("plotSignal") +"_featCorrelations";
	std::string 								histName 									=	"correlations_" + options.get("plotSignal");

	if(options.getInt("useSaved") && file_exists(writePath+".root") ){
		corrMatHist 																		=	(TH2F*)	getHistFromFile(histName, writePath+".root");
	} else{
		corrMatHist																		= 	getHist();
		corrMatHist->SetName(histName.c_str());
		writeToFile(corrMatHist, writePath+".root", "RECREATE", 1);
	}


	ofstream 								featCorrDump(writePath+".txt", std::ofstream::out);
	for(Int_t iBinX = 0; iBinX <= corrMatHist->GetNbinsX(); iBinX++){
		if(iBinX == 0) featCorrDump<<",";
		else  featCorrDump<<feats[iBinX-1]<<",";
		for(Int_t iBinY = 1; iBinY <= corrMatHist->GetNbinsY(); iBinY++){
			if(iBinX == 0){
				featCorrDump<<feats[iBinY-1]<<(iBinY < corrMatHist->GetNbinsY() ? "," : "");
			} else{
				featCorrDump<<to_string_with_precision(corrMatHist->GetBinContent(iBinX, iBinY), 10)<<(iBinY < corrMatHist->GetNbinsY() ? "," : "");
			}
		}		
		featCorrDump<<std::endl;
	}
	featCorrDump.close();

	TCanvas 									canvas("canvas", "", options.getDouble("canvasX"), options.getDouble("canvasY"));
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
	if(options.getInt("logZ")) pad0.SetLogz();

	TH2F* corrMatHistAbs = (TH2F*) getAbsHist(corrMatHist);	
	corrMatHistAbs->SetContour(options.getInt("2dNcontours"));
	corrMatHistAbs->Draw(options.getCSTR("drawOpt1"));


	corrMatHistAbs->GetXaxis()->CenterTitle();
	corrMatHistAbs->GetYaxis()->CenterTitle();
	corrMatHistAbs->GetZaxis()->CenterTitle();
	corrMatHistAbs->GetZaxis()->RotateTitle();
	corrMatHistAbs->GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	corrMatHistAbs->GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	corrMatHistAbs->GetZaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	corrMatHistAbs->GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	corrMatHistAbs->GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	corrMatHistAbs->GetZaxis()->SetLabelSize(options.getDouble("pad0ZaxisLabelSize"));
	corrMatHistAbs->GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	corrMatHistAbs->GetXaxis()->SetLabelOffset(options.getDouble("pad0XlabelOffset"));
	corrMatHistAbs->GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	corrMatHistAbs->GetZaxis()->SetTitleOffset(options.getDouble("pad0ZtitleOffset"));
	corrMatHistAbs->GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	corrMatHistAbs->GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	corrMatHistAbs->GetXaxis()->SetAlphanumeric();
	corrMatHistAbs->GetYaxis()->SetAlphanumeric();

	for(UInt_t iBin = 0; iBin < feats.size(); iBin++){
		std::string binLabel = varPlotInfo[feats[iBin]][0];
		binLabel = findAndReplaceAll(binLabel, "(GeV)", "");
		trim(binLabel);
		corrMatHistAbs->GetYaxis()->SetBinLabel(iBin+1, binLabel.c_str());
		// corrMatHistAbs->GetXaxis()->SetBinLabel(iBin+1, binLabel.c_str());
		corrMatHistAbs->GetXaxis()->ChangeLabel(iBin+1, options.getFloat("labAngle"), options.getFloat("xlabSize"),options.getInt("labAlign"),-1,-1,binLabel);
	}

	corrMatHistAbs->GetXaxis()->ChangeLabel(feats.size()+1, options.getFloat("labAngle"), 0.,options.getInt("labAlign"),-1,-1,"");

	corrMatHistAbs->GetXaxis()->CenterLabels();


	corrMatHist->SetMarkerColor(options.getTColFromHex("tCol"));
	corrMatHist->Draw(options.getCSTR("drawOpt2"));

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

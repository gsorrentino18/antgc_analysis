#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"

std::string 										optionsFile									=	"getAltSampleEfficiencyV2Options.txt";
parseOptions 										options;
UChar_t 											idBit 										=	3;
Bool_t 												useSaved									=	1;
void 												initialize();
TGraphAsymmErrors* getPtEfficiency(std::string _plotName, std::string _inDir, std::string _inSamples, std::vector<Double_t> _ptBinning, Float_t _etaMin, Float_t _etaMax);
void	 	drawPtEfficiency(Bool_t signal, std::vector<Double_t> _ptBinning, Float_t _etaMin, Float_t _etaMax);

TGraphAsymmErrors* getEtaEfficiency(std::string _plotName, std::string _inDir, std::string _inSamples, std::vector<Double_t> _etaBinning, Float_t _ptMin, Float_t _ptMax);
void	 	drawEtaEfficiency(Bool_t signal, std::vector<Double_t> _etaBinning, Float_t _ptMin, Float_t _ptMax);

void initialize(){
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	options.parseIt(optionsFile, "==");

	mkdir(options.get("writeDir"));

	idBit 																						=	options.getInt("idBit");
	useSaved																					=	options.getInt("useSaved");
	std::vector<Double_t> 						pTbinning										= 	options.getDoubleList("pTbinning");
	std::vector<Double_t> 						etaRange										= 	options.getDoubleList("etaRange");
	std::vector<Double_t> 						ptRange											= 	options.getDoubleList("ptRange");
	std::vector<Double_t> 						etaBinning										= 	options.getDoubleList("etaBinning");

	drawPtEfficiency(1, pTbinning, etaRange[0], etaRange[1]);
	// drawPtEfficiency(0, pTbinning, etaRange[0], etaRange[1]);

	// drawEtaEfficiency(1, etaBinning, ptRange[0], ptRange[1]);
	// drawEtaEfficiency(0, etaBinning, ptRange[0], ptRange[1]);
	
	std::cout<<"Done"<<std::endl;
};


TGraphAsymmErrors* getEtaEfficiency(std::string _plotName, std::string _inDir, std::string _inSamples, std::vector<Double_t> _etaBinning, Float_t _ptMin, Float_t _ptMax){

	Float_t 									etaMin 									=		_etaBinning[0];
	Float_t 									etaMax 									=		_etaBinning.back();

	TEfficiency									eff(_plotName.c_str(), "", _etaBinning.size()-1, _etaBinning.data());
	eff.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));

	std::vector<std::string> 					inSamplesV								= 		split_string(_inSamples, ",");

	for(std::string iSample : inSamplesV){

		Double_t 								xSection 			=	std::stod(vLookup(iSample, options.get("xSectionMap"), 0, 2));
		Double_t 								sumGenWeight		=	std::stod(vLookup(iSample, options.get("xSectionMap"), 0, 7));

		TChain*                                 inTree      							=		openTChain(std::vector<std::string>({_inDir + "/" + iSample + ".root"}), options.get("inTreeName"));
		TTreeReader                             inputTTreeReader(inTree);
		TTreeReaderAnyValue<UChar_t>   			phoPFClusIDbits(inputTTreeReader, "phoPFClusIDbits");
		TTreeReaderAnyValue<Float_t>            puWeight(inputTTreeReader, "puWeight");
		TTreeReaderAnyValue<Float_t>            genWeight(inputTTreeReader, "genWeight");
		TTreeReaderAnyValue<Float_t>			phoPt(inputTTreeReader, "phoPt");
		TTreeReaderAnyValue<Float_t>			phoSCeta(inputTTreeReader, "phoSCeta");

		while(inputTTreeReader.Next()){


			UChar_t 							idBits 										=	phoPFClusIDbits;
			Double_t							evWeight 									=	xSection * puWeight * genWeight/sumGenWeight;

			// std::cout<<getBit(idBits, idBit)<<" "<<evWeight<<" "<<phoSCeta<<" "<<phoPt<<std::endl;

			if (phoSCeta <= etaMin) continue;
			if (phoSCeta >= etaMax) continue;

			if(phoPt < _ptMin) continue;
			if(phoPt > _ptMax) continue;


			eff.FillWeighted(getBit(idBits, idBit), evWeight, phoSCeta);
		}

		closeTChain(inTree);
	}

	TGraphAsymmErrors* 							effGraph 									=	eff.CreateGraph(options.getCSTR("gOpt"));

	return effGraph;
}

void	 	drawEtaEfficiency(Bool_t signal, std::vector<Double_t> _etaBinning, Float_t _ptMin, Float_t _ptMax){

	TMultiGraph 								stack;

	TLegend legend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	legend.SetTextSize(options.getDouble("legendTextSize"));
	legend.SetNColumns(options.getInt("legendNcols"));
	legend.SetBorderSize(0);
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	std::string 								legTitle 								=	(signal ? "Prompt  " : "Fake  ") +  removeTrailingZeros(_ptMin) + " GeV #leqp_{T}<" +removeTrailingZeros(_ptMax) + " GeV";
	legend.SetHeader(legTitle.c_str(), "C");
	legend.SetTextAlign(12);
	
	std::string 								inDir 	= signal ? options.get("signalInDir") : options.get("backgroundInDir");
	std::vector<std::string>					samples = signal ? options.getList("signalSamples", ";") :options.getList("backgroundSamples", ";");
	for(std::string iSample : samples){

		std::vector<string> 					iSampleParts 									= 	split_string(iSample, ":");
		TGraphAsymmErrors*						iSampleEffGraph;

		std::string 							writePath 										=	options.get("writeDir") + "/" + std::to_string(signal) + "_eta_" + iSampleParts[2] + ".root";

		if(useSaved && file_exists(writePath)){
			iSampleEffGraph																		=	(TGraphAsymmErrors*) getObjectFromFile(iSampleParts[2], writePath);
		} else{
			iSampleEffGraph																		= 	getEtaEfficiency(iSampleParts[2], inDir, iSampleParts[3], _etaBinning, _ptMin, _ptMax);
			iSampleEffGraph->SetName(iSampleParts[2].c_str());
			writeToFile(iSampleEffGraph, writePath, "RECREATE", 1);
		}
		
		iSampleEffGraph->SetLineStyle(options.getInt("lineStyle"));
		iSampleEffGraph->SetLineWidth(options.getInt("lineWidth"));
		iSampleEffGraph->SetMarkerSize(options.getFloat("mkrSize"));
		iSampleEffGraph->SetMarkerStyle(options.getInt("mkrStyle"));
		iSampleEffGraph->SetLineColor(hex2rootColor(iSampleParts[0]));
		iSampleEffGraph->SetMarkerColor(hex2rootColor(iSampleParts[0]));

		stack.Add(iSampleEffGraph, options.getCSTR("dOpt"));
		legend.AddEntry(iSampleEffGraph, iSampleParts[1].c_str(), "LPE");
	}

	TCanvas canvas((std::to_string(signal) + "_pTeff_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);
	
	TPad pad("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad.SetFillStyle(4000);
	pad.SetFillColor(0);
	pad.SetFrameFillStyle(4000);
	pad.SetGrid(1,1);

	canvas.Draw();
	canvas.cd();
	pad.Draw();
	pad.cd();

	stack.Draw(options.getCSTR("dOpt"));
	legend.Draw();

	stack.GetXaxis()->SetTitle("#eta^{SC}");
	stack.GetYaxis()->SetTitle("Efficiency");

	stack.GetXaxis()->CenterTitle();
	stack.GetYaxis()->CenterTitle();
	stack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	stack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	stack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	stack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	stack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	stack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	stack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	stack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	stack.SetMaximum(options.getFloat("yScaleupSig")*stack.GetYaxis()->GetXmax());

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ std::to_string(signal) + "_eta_efficiency.png").c_str());
	canvas.SaveAs((options.get("writeDir")+ std::to_string(signal) + "_eta_efficiency.pdf").c_str());

	clearHeap();
};




void	 drawPtEfficiency(Bool_t signal, std::vector<Double_t> _ptBinning, Float_t _etaMin, Float_t _etaMax){

	TMultiGraph 								stack;

	TLegend legend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	legend.SetTextSize(options.getDouble("legendTextSize"));
	legend.SetNColumns(options.getInt("legendNcols"));
	legend.SetBorderSize(0);
	legend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	legend.SetFillStyle(options.getInt("legFillStyle"));
	std::string 								legTitle 								=	(signal ? "Prompt  " : "Fake  ") + removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax);
	legend.SetHeader(legTitle.c_str(), "C");
	legend.SetTextAlign(12);
	
	std::string 								inDir 	= signal ? options.get("signalInDir") : options.get("backgroundInDir");
	std::vector<std::string>					samples = signal ? options.getList("signalSamples", ";") :options.getList("backgroundSamples", ";");
	for(std::string iSample : samples){


		std::vector<string> 					iSampleParts 									= 	split_string(iSample, ":");
		TGraphAsymmErrors*						iSampleEffGraph;

		std::string 							writePath 										=	options.get("writeDir") + "/" + std::to_string(signal) + "_pT_" + iSampleParts[2] + ".root";

		if(useSaved && file_exists(writePath)){
			iSampleEffGraph																		=	(TGraphAsymmErrors*) getObjectFromFile(iSampleParts[2], writePath);
		} else{
			iSampleEffGraph																		= 	getPtEfficiency(iSampleParts[2], inDir, iSampleParts[3], _ptBinning, _etaMin, _etaMax);
			iSampleEffGraph->SetName(iSampleParts[2].c_str());
			writeToFile(iSampleEffGraph, writePath, "RECREATE", 1);
		}
		
		iSampleEffGraph->SetLineStyle(options.getInt("lineStyle"));
		iSampleEffGraph->SetLineWidth(options.getInt("lineWidth"));
		iSampleEffGraph->SetMarkerSize(options.getFloat("mkrSize"));
		iSampleEffGraph->SetMarkerStyle(options.getInt("mkrStyle"));
		iSampleEffGraph->SetLineColor(hex2rootColor(iSampleParts[0]));
		iSampleEffGraph->SetMarkerColor(hex2rootColor(iSampleParts[0]));

		stack.Add(iSampleEffGraph, options.getCSTR("dOpt"));
		legend.AddEntry(iSampleEffGraph, iSampleParts[1].c_str(), "LPE");
	}

	TCanvas canvas((std::to_string(signal) + "_pTeff_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);
	
	TPad pad("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad.SetFillStyle(4000);
	pad.SetFillColor(0);
	pad.SetFrameFillStyle(4000);
	pad.SetGrid(1,1);

	canvas.Draw();
	canvas.cd();
	pad.Draw();
	pad.cd();

	stack.Draw(options.getCSTR("dOpt"));
	legend.Draw();

	stack.GetXaxis()->SetTitle("p_{T} (GeV)");
	stack.GetYaxis()->SetTitle("Efficiency");

	stack.GetXaxis()->CenterTitle();
	stack.GetYaxis()->CenterTitle();
	stack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	stack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	stack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	stack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	stack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	stack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	stack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	stack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	stack.SetMaximum(options.getFloat("yScaleupSig")*stack.GetYaxis()->GetXmax());

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ std::to_string(signal) + "_pT_efficiency.png").c_str());
	canvas.SaveAs((options.get("writeDir")+ std::to_string(signal) + "_pT_efficiency.pdf").c_str());

	clearHeap();
};

TGraphAsymmErrors* getPtEfficiency(std::string _plotName, std::string _inDir, std::string _inSamples, std::vector<Double_t> _ptBinning, Float_t _etaMin, Float_t _etaMax){

	Float_t 									ptMin 									=		_ptBinning[0];
	Float_t 									ptMax 									=		_ptBinning.back();

	TEfficiency									eff(_plotName.c_str(), "", _ptBinning.size()-1, _ptBinning.data());
	eff.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));


	std::vector<std::string> 					inSamplesV								= 		split_string(_inSamples, ",");

	for(std::string iSample : inSamplesV){

		Double_t 								xSection 			=	std::stod(vLookup(iSample, options.get("xSectionMap"), 0, 2));
		Double_t 								sumGenWeight		=	std::stod(vLookup(iSample, options.get("xSectionMap"), 0, 7));

		TChain*                                 inTree      							=		openTChain(std::vector<std::string>({_inDir + "/" + iSample + ".root"}), options.get("inTreeName"));
		TChain*                                 inTree2      							=		openTChain(std::vector<std::string>({"data/" + iSample + ".root"}), "fullEB_BDT_Tree");
		inTree->AddFriend(inTree2);
		TTreeReader                             inputTTreeReader(inTree);
		TTreeReaderAnyValue<UChar_t>   			phoPFClusIDbits(inputTTreeReader, "phoPFClusIDbits");
		TTreeReaderAnyValue<Float_t>            puWeight(inputTTreeReader, "puWeight");
		TTreeReaderAnyValue<Float_t>            genWeight(inputTTreeReader, "genWeight");
		TTreeReaderAnyValue<Float_t>			phoPt(inputTTreeReader, "phoPt");
		TTreeReaderAnyValue<Float_t>			phoSCeta(inputTTreeReader, "phoSCeta");

		TTreeReaderAnyValue<Double_t>				bdtScore_										(inputTTreeReader, "bdtScore");
		TTreeReaderAnyValue<Float_t>				phoHoverE_										(inputTTreeReader, "phoHoverE");
		TTreeReaderAnyValue<Float_t>				phoPFECALClusIsoCorr_							(inputTTreeReader, "phoPFECALClusIsoCorr");
		TTreeReaderAnyValue<Float_t>				phoPFHCALClusIsoCorr_							(inputTTreeReader, "phoPFHCALClusIsoCorr");
		TTreeReaderAnyValue<Float_t>				phoTkrIsoCorr_									(inputTTreeReader, "phoTkrIsoCorr");

		while(inputTTreeReader.Next()){

			Float_t 							absSCEta									=	std::abs(phoSCeta);
			if (absSCEta < _etaMin) continue;
			if (absSCEta >= _etaMax) continue;

			if(phoPt < ptMin) continue;
			if(phoPt > ptMax) continue;

			
			Double_t							evWeight 									=	xSection * puWeight * genWeight/sumGenWeight;

			// UChar_t 							idBits 										=	phoPFClusIDbits;
			// eff.FillWeighted(getBit(idBits, idBit), evWeight, phoPt);

			Bool_t pass95 = (bdtScore_ >= 8.4858950660326768e-02) && (phoHoverE_ < 4.5018283831009365e-02) && (phoPFECALClusIsoCorr_ < 3.4103181645009730e+00) && (phoPFHCALClusIsoCorr_ < 9.2786995496122131e+00) && (phoTkrIsoCorr_ < 4.2117599454665431e+00);
			eff.FillWeighted(pass95, evWeight, phoPt);
		}

		closeTChain(inTree);
	}

	TGraphAsymmErrors* 							effGraph 									=	eff.CreateGraph(options.getCSTR("gOpt"));

	return effGraph;
};
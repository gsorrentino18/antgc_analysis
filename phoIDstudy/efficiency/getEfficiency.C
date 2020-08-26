#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"

std::string 										optionsFile="getEfficiencyOptions.txt";
parseOptions 										options;
void getPtEfficiency(std::vector<Float_t> _ptBinning, Float_t _etaMin, Float_t _etaMax);
void getEtaEfficiency(std::vector<Float_t> _etaBinning, Float_t _ptMin, Float_t _ptMax);
void getNvtxEfficiency(Float_t _ptMin, Float_t _ptMax, Float_t _etaMin, Float_t _etaMax);
void initialize();

void getPtEfficiency(const std::vector<Double_t> _ptBinning, Float_t _etaMin, Float_t _etaMax){

	Float_t 									ptMin 									=		_ptBinning[0];
	Float_t 									ptMax 									=		_ptBinning.back();

	TChain*                                     featsTree      							=		openTChain(std::vector<std::string>({options.get("featsFile")}), options.get("featsTree"));
	TChain*                                     bdtTree         						=		openTChain(std::vector<std::string>({options.get("bdtFile")}), options.get("bdtTree"));
	TChain*                                     idTree      							=		openTChain(std::vector<std::string>({options.get("idFile")}), options.get("idTree"));
	featsTree->AddFriend(bdtTree);
	featsTree->AddFriend(idTree);
	TTreeReader                             	inputTTreeReader(featsTree);
	TTreeReaderAnyValue<Double_t>				bdtWeightF_										(inputTTreeReader, "bdtWeightF");
	TTreeReaderAnyValue<Bool_t>  				pass95_											(inputTTreeReader, "pass95");
	TTreeReaderAnyValue<Bool_t>  				pass90_											(inputTTreeReader, "pass90");
	TTreeReaderAnyValue<Bool_t>  				pass80_											(inputTTreeReader, "pass80");
	TTreeReaderAnyValue<Bool_t>  				pass70_											(inputTTreeReader, "pass70");
	TTreeReaderAnyValue<Bool_t>  				isTrain_										(inputTTreeReader, "isTrain");
	TTreeReaderAnyValue<Bool_t>  				isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Float_t>				phoPt_											(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");

	std::string 								plotame 								=	"Efficiency_vs_pT_eta" + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p");
	TEfficiency									sigEff95(("signal95_"+plotame).c_str(), "", _ptBinning.size()-1, _ptBinning.data());
	TEfficiency									sigEff90(("signal90_"+plotame).c_str(), "", _ptBinning.size()-1, _ptBinning.data());
	TEfficiency									sigEff80(("signal80_"+plotame).c_str(), "", _ptBinning.size()-1, _ptBinning.data());
	TEfficiency									sigEff70(("signal70_"+plotame).c_str(), "", _ptBinning.size()-1, _ptBinning.data());
	
	TEfficiency									bgEff95(("background95_"+plotame).c_str(), "", _ptBinning.size()-1, _ptBinning.data());
	TEfficiency									bgEff90(("background90_"+plotame).c_str(), "", _ptBinning.size()-1, _ptBinning.data());
	TEfficiency									bgEff80(("background90_"+plotame).c_str(), "", _ptBinning.size()-1, _ptBinning.data());
	TEfficiency									bgEff70(("background90_"+plotame).c_str(), "", _ptBinning.size()-1, _ptBinning.data());

	sigEff95.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff90.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff80.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff70.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));

	bgEff95.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff90.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff80.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff70.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));	

	while(inputTTreeReader.Next()){

		if(isTrain_) continue;
		
		Float_t 							absSCEta									=	std::abs(phoSCeta_);
		if (absSCEta <= _etaMin) continue;
		if (absSCEta > _etaMax) continue;

		if(phoPt_ < ptMin) continue;
		if(phoPt_ > ptMax) continue;

		if(isSignal_){
			sigEff95.FillWeighted(pass95_, bdtWeightF_, phoPt_);
			sigEff90.FillWeighted(pass90_, bdtWeightF_, phoPt_);
			sigEff80.FillWeighted(pass80_, bdtWeightF_, phoPt_);
			sigEff70.FillWeighted(pass70_, bdtWeightF_, phoPt_);
		} else{
			bgEff95.FillWeighted(pass95_, bdtWeightF_, phoPt_);
			bgEff90.FillWeighted(pass90_, bdtWeightF_, phoPt_);
			bgEff80.FillWeighted(pass80_, bdtWeightF_, phoPt_);
			bgEff70.FillWeighted(pass70_, bdtWeightF_, phoPt_);
		}
	}

	closeTChain(featsTree);
	closeTChain(bdtTree);
	closeTChain(idTree);

	TGraphAsymmErrors* 							sigEff95Graph 							=	sigEff95.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff90Graph 							=	sigEff90.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff80Graph 							=	sigEff80.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff70Graph 							=	sigEff70.CreateGraph(options.getCSTR("gOpt"));

	sigEff95Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff95Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff95Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff95Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff90Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff90Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff90Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff90Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff80Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff80Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff80Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff80Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff70Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff70Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff70Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff70Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff95Graph->SetLineColor(options.getTColFromHex("eff95Col"));
	sigEff90Graph->SetLineColor(options.getTColFromHex("eff90Col"));
	sigEff80Graph->SetLineColor(options.getTColFromHex("eff80Col"));
	sigEff70Graph->SetLineColor(options.getTColFromHex("eff70Col"));

	sigEff95Graph->SetMarkerColor(options.getTColFromHex("eff95Col"));
	sigEff90Graph->SetMarkerColor(options.getTColFromHex("eff90Col"));
	sigEff80Graph->SetMarkerColor(options.getTColFromHex("eff80Col"));
	sigEff70Graph->SetMarkerColor(options.getTColFromHex("eff70Col"));

	TMultiGraph 								sStack;
	sStack.Add(sigEff95Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff90Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff80Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff70Graph, options.getCSTR("dOpt"));
	
	TLegend sLegend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	sLegend.SetTextSize(options.getDouble("legendTextSize"));
	sLegend.SetNColumns(options.getInt("legendNcols"));
	sLegend.SetBorderSize(0);
	sLegend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	sLegend.SetFillStyle(options.getInt("legFillStyle"));
	std::string 								sLegTitle 								=	"Prompt  " + removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax);
	sLegend.SetHeader(sLegTitle.c_str(), "C");
	sLegend.SetTextAlign(12);
	sLegend.AddEntry(sigEff95Graph, "Very Loose", "LPE");
	sLegend.AddEntry(sigEff90Graph, "Loose", "LPE");
	sLegend.AddEntry(sigEff80Graph, "Medium", "LPE");
	sLegend.AddEntry(sigEff70Graph, "Tight", "LPE");

	TCanvas canvas((plotame + "_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
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

	sStack.Draw(options.getCSTR("dOpt"));
	sLegend.Draw();

	sStack.GetXaxis()->SetTitle("p_{T} (GeV)");
	sStack.GetYaxis()->SetTitle("Efficiency");

	sStack.GetXaxis()->CenterTitle();
	sStack.GetYaxis()->CenterTitle();
	sStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	sStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	sStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	sStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	sStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	sStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	sStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	sStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	sStack.SetMaximum(options.getFloat("yScaleupSig")*sStack.GetYaxis()->GetXmax());

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "signal_" + plotame+".png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "signal_" + plotame+".pdf").c_str());

	TGraphAsymmErrors* 							bgEff95Graph 							=	bgEff95.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff90Graph 							=	bgEff90.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff80Graph 							=	bgEff80.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff70Graph 							=	bgEff70.CreateGraph(options.getCSTR("gOpt"));

	bgEff95Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff95Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff95Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff95Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff90Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff90Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff90Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff90Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff80Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff80Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff80Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff80Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff70Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff70Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff70Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff70Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff95Graph->SetLineColor(options.getTColFromHex("eff95Col"));
	bgEff90Graph->SetLineColor(options.getTColFromHex("eff90Col"));
	bgEff80Graph->SetLineColor(options.getTColFromHex("eff80Col"));
	bgEff70Graph->SetLineColor(options.getTColFromHex("eff70Col"));

	bgEff95Graph->SetMarkerColor(options.getTColFromHex("eff95Col"));
	bgEff90Graph->SetMarkerColor(options.getTColFromHex("eff90Col"));
	bgEff80Graph->SetMarkerColor(options.getTColFromHex("eff80Col"));
	bgEff70Graph->SetMarkerColor(options.getTColFromHex("eff70Col"));

	TMultiGraph 								bStack;
	bStack.Add(bgEff95Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff90Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff80Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff70Graph, options.getCSTR("dOpt"));

	TLegend bLegend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	bLegend.SetTextSize(options.getDouble("legendTextSize"));
	bLegend.SetNColumns(options.getInt("legendNcols"));
	bLegend.SetBorderSize(0);
	bLegend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	bLegend.SetFillStyle(options.getInt("legFillStyle"));
	std::string 								bLegTitle 								=	"Fake  " + removeTrailingZeros(_etaMin) + "#leq|#eta|<" +removeTrailingZeros(_etaMax);
	bLegend.SetHeader(bLegTitle.c_str(), "C");
	bLegend.SetTextAlign(12);
	bLegend.AddEntry(bgEff95Graph, "Very Loose", "LPE");
	bLegend.AddEntry(bgEff90Graph, "Loose", "LPE");
	bLegend.AddEntry(bgEff80Graph, "Medium", "LPE");
	bLegend.AddEntry(bgEff70Graph, "Tight", "LPE");

	pad.Clear();
	pad.cd();


	bStack.Draw((options.get("dOpt")).c_str());
	bStack.SetMaximum(options.getFloat("yScaleupBG")*bStack.GetYaxis()->GetXmax());

	bStack.Draw(options.getCSTR("dOpt"));
	bLegend.Draw();

	bStack.GetXaxis()->SetTitle("p_{T} (GeV)");
	bStack.GetYaxis()->SetTitle("Efficiency");

	bStack.GetXaxis()->CenterTitle();
	bStack.GetYaxis()->CenterTitle();
	bStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	bStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	bStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	bStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	bStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	bStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	bStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	bStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "background_" + plotame+".png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "background_" + plotame+".pdf").c_str());
	
	clearHeap();
	
};


void getEtaEfficiency(const std::vector<Double_t> _etaBinning, Float_t _ptMin, Float_t _ptMax){

	Float_t 									etaMin 									=		_etaBinning[0];
	Float_t 									etaMax 									=		_etaBinning.back();

	TChain*                                     featsTree      							=		openTChain(std::vector<std::string>({options.get("featsFile")}), options.get("featsTree"));
	TChain*                                     bdtTree         						=		openTChain(std::vector<std::string>({options.get("bdtFile")}), options.get("bdtTree"));
	TChain*                                     idTree      							=		openTChain(std::vector<std::string>({options.get("idFile")}), options.get("idTree"));
	featsTree->AddFriend(bdtTree);
	featsTree->AddFriend(idTree);
	TTreeReader                             	inputTTreeReader(featsTree);
	TTreeReaderAnyValue<Double_t>				bdtWeightF_										(inputTTreeReader, "bdtWeightF");
	TTreeReaderAnyValue<Bool_t>  				pass95_											(inputTTreeReader, "pass95");
	TTreeReaderAnyValue<Bool_t>  				pass90_											(inputTTreeReader, "pass90");
	TTreeReaderAnyValue<Bool_t>  				pass80_											(inputTTreeReader, "pass80");
	TTreeReaderAnyValue<Bool_t>  				pass70_											(inputTTreeReader, "pass70");
	TTreeReaderAnyValue<Bool_t>  				isTrain_										(inputTTreeReader, "isTrain");
	TTreeReaderAnyValue<Bool_t>  				isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Float_t>				phoPt_											(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");

	std::string 								plotame 								=	"Efficiency_vs_Eta_pt" + findAndReplaceAll(removeTrailingZeros(_ptMin) + "to" + removeTrailingZeros(_ptMax), ".", "p");
	TEfficiency									sigEff95(("signal95_"+plotame).c_str(), "", _etaBinning.size()-1, _etaBinning.data());
	TEfficiency									sigEff90(("signal90_"+plotame).c_str(), "", _etaBinning.size()-1, _etaBinning.data());
	TEfficiency									sigEff80(("signal80_"+plotame).c_str(), "", _etaBinning.size()-1, _etaBinning.data());
	TEfficiency									sigEff70(("signal70_"+plotame).c_str(), "", _etaBinning.size()-1, _etaBinning.data());
	
	TEfficiency									bgEff95(("background95_"+plotame).c_str(), "", _etaBinning.size()-1, _etaBinning.data());
	TEfficiency									bgEff90(("background90_"+plotame).c_str(), "", _etaBinning.size()-1, _etaBinning.data());
	TEfficiency									bgEff80(("background90_"+plotame).c_str(), "", _etaBinning.size()-1, _etaBinning.data());
	TEfficiency									bgEff70(("background90_"+plotame).c_str(), "", _etaBinning.size()-1, _etaBinning.data());

	sigEff95.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff90.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff80.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff70.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));

	bgEff95.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff90.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff80.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff70.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));	

	while(inputTTreeReader.Next()){

		if(isTrain_) continue;
		
		if (phoSCeta_ < etaMin) continue;
		if (phoSCeta_ > etaMax) continue;

		if(phoPt_ < _ptMin) continue;
		if(phoPt_ > _ptMax) continue;

		if(isSignal_){
			sigEff95.FillWeighted(pass95_, bdtWeightF_, phoSCeta_);
			sigEff90.FillWeighted(pass90_, bdtWeightF_, phoSCeta_);
			sigEff80.FillWeighted(pass80_, bdtWeightF_, phoSCeta_);
			sigEff70.FillWeighted(pass70_, bdtWeightF_, phoSCeta_);
		} else{
			bgEff95.FillWeighted(pass95_, bdtWeightF_, phoSCeta_);
			bgEff90.FillWeighted(pass90_, bdtWeightF_, phoSCeta_);
			bgEff80.FillWeighted(pass80_, bdtWeightF_, phoSCeta_);
			bgEff70.FillWeighted(pass70_, bdtWeightF_, phoSCeta_);
		}
	}

	closeTChain(featsTree);
	closeTChain(bdtTree);
	closeTChain(idTree);

	TGraphAsymmErrors* 							sigEff95Graph 							=	sigEff95.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff90Graph 							=	sigEff90.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff80Graph 							=	sigEff80.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff70Graph 							=	sigEff70.CreateGraph(options.getCSTR("gOpt"));

	sigEff95Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff95Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff95Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff95Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff90Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff90Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff90Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff90Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff80Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff80Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff80Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff80Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff70Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff70Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff70Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff70Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff95Graph->SetLineColor(options.getTColFromHex("eff95Col"));
	sigEff90Graph->SetLineColor(options.getTColFromHex("eff90Col"));
	sigEff80Graph->SetLineColor(options.getTColFromHex("eff80Col"));
	sigEff70Graph->SetLineColor(options.getTColFromHex("eff70Col"));

	sigEff95Graph->SetMarkerColor(options.getTColFromHex("eff95Col"));
	sigEff90Graph->SetMarkerColor(options.getTColFromHex("eff90Col"));
	sigEff80Graph->SetMarkerColor(options.getTColFromHex("eff80Col"));
	sigEff70Graph->SetMarkerColor(options.getTColFromHex("eff70Col"));

	TMultiGraph 								sStack;
	sStack.Add(sigEff95Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff90Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff80Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff70Graph, options.getCSTR("dOpt"));
	
	TLegend sLegend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	sLegend.SetTextSize(options.getDouble("legendTextSize"));
	sLegend.SetNColumns(options.getInt("legendNcols"));
	sLegend.SetBorderSize(0);
	sLegend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	sLegend.SetFillStyle(options.getInt("legFillStyle"));
	std::string 								sLegTitle 								=	"Prompt  " + removeTrailingZeros(_ptMin) + " GeV #leqp_{T}<" +removeTrailingZeros(_ptMax) + " GeV";
	sLegend.SetHeader(sLegTitle.c_str(), "C");
	sLegend.SetTextAlign(12);
	sLegend.AddEntry(sigEff95Graph, "Very Loose", "LPE");
	sLegend.AddEntry(sigEff90Graph, "Loose", "LPE");
	sLegend.AddEntry(sigEff80Graph, "Medium", "LPE");
	sLegend.AddEntry(sigEff70Graph, "Tight", "LPE");

	TCanvas canvas((plotame + "_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
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

	sStack.Draw(options.getCSTR("dOpt"));
	sLegend.Draw();

	sStack.GetXaxis()->SetTitle("#eta^{SC}");
	sStack.GetYaxis()->SetTitle("Efficiency");

	sStack.GetXaxis()->CenterTitle();
	sStack.GetYaxis()->CenterTitle();
	sStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	sStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	sStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	sStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	sStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	sStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	sStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	sStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	sStack.SetMaximum(options.getFloat("yScaleupSig")*sStack.GetYaxis()->GetXmax());

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "signal_" + plotame+".png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "signal_" + plotame+".pdf").c_str());

	TGraphAsymmErrors* 							bgEff95Graph 							=	bgEff95.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff90Graph 							=	bgEff90.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff80Graph 							=	bgEff80.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff70Graph 							=	bgEff70.CreateGraph(options.getCSTR("gOpt"));

	bgEff95Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff95Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff95Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff95Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff90Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff90Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff90Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff90Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff80Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff80Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff80Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff80Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff70Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff70Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff70Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff70Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff95Graph->SetLineColor(options.getTColFromHex("eff95Col"));
	bgEff90Graph->SetLineColor(options.getTColFromHex("eff90Col"));
	bgEff80Graph->SetLineColor(options.getTColFromHex("eff80Col"));
	bgEff70Graph->SetLineColor(options.getTColFromHex("eff70Col"));

	bgEff95Graph->SetMarkerColor(options.getTColFromHex("eff95Col"));
	bgEff90Graph->SetMarkerColor(options.getTColFromHex("eff90Col"));
	bgEff80Graph->SetMarkerColor(options.getTColFromHex("eff80Col"));
	bgEff70Graph->SetMarkerColor(options.getTColFromHex("eff70Col"));

	TMultiGraph 								bStack;
	bStack.Add(bgEff95Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff90Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff80Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff70Graph, options.getCSTR("dOpt"));

	TLegend bLegend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	bLegend.SetTextSize(options.getDouble("legendTextSize"));
	bLegend.SetNColumns(options.getInt("legendNcols"));
	bLegend.SetBorderSize(0);
	bLegend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	bLegend.SetFillStyle(options.getInt("legFillStyle"));
	std::string 								bLegTitle 								=	"Fake  " + removeTrailingZeros(_ptMin) + " GeV #leqp_{T}<" +removeTrailingZeros(_ptMax) + " GeV";
	bLegend.SetHeader(bLegTitle.c_str(), "C");
	bLegend.SetTextAlign(12);
	bLegend.AddEntry(bgEff95Graph, "Very Loose", "LPE");
	bLegend.AddEntry(bgEff90Graph, "Loose", "LPE");
	bLegend.AddEntry(bgEff80Graph, "Medium", "LPE");
	bLegend.AddEntry(bgEff70Graph, "Tight", "LPE");	

	pad.Clear();
	pad.cd();

	bStack.Draw((options.get("dOpt")).c_str());
	bStack.SetMaximum(options.getFloat("yScaleupBG")*bStack.GetYaxis()->GetXmax());

	bStack.Draw(options.getCSTR("dOpt"));
	bLegend.Draw();

	bStack.GetXaxis()->SetTitle("#eta^{SC}");
	bStack.GetYaxis()->SetTitle("Efficiency");

	bStack.GetXaxis()->CenterTitle();
	bStack.GetYaxis()->CenterTitle();
	bStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	bStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	bStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	bStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	bStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	bStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	bStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	bStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "background_" + plotame+".png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "background_" + plotame+".pdf").c_str());
	
	clearHeap();
	
};


void getNvtxEfficiency(Float_t _ptMin, Float_t _ptMax, Float_t _etaMin, Float_t _etaMax){

	TChain*                                     featsTree      							=		openTChain(std::vector<std::string>({options.get("featsFile")}), options.get("featsTree"));
	TChain*                                     bdtTree         						=		openTChain(std::vector<std::string>({options.get("bdtFile")}), options.get("bdtTree"));
	TChain*                                     idTree      							=		openTChain(std::vector<std::string>({options.get("idFile")}), options.get("idTree"));
	featsTree->AddFriend(bdtTree);
	featsTree->AddFriend(idTree);
	TTreeReader                             	inputTTreeReader(featsTree);
	TTreeReaderAnyValue<Double_t>				bdtWeightF_										(inputTTreeReader, "bdtWeightF");
	TTreeReaderAnyValue<Bool_t>  				pass95_											(inputTTreeReader, "pass95");
	TTreeReaderAnyValue<Bool_t>  				pass90_											(inputTTreeReader, "pass90");
	TTreeReaderAnyValue<Bool_t>  				pass80_											(inputTTreeReader, "pass80");
	TTreeReaderAnyValue<Bool_t>  				pass70_											(inputTTreeReader, "pass70");
	TTreeReaderAnyValue<Bool_t>  				isTrain_										(inputTTreeReader, "isTrain");
	TTreeReaderAnyValue<Bool_t>  				isSignal_										(inputTTreeReader, "isSignal");
	TTreeReaderAnyValue<Float_t>				phoPt_											(inputTTreeReader, "phoPt");
	TTreeReaderAnyValue<Float_t>				phoSCeta_										(inputTTreeReader, "phoSCeta");
	TTreeReaderAnyValue<UChar_t>				nVtx_											(inputTTreeReader, "nVtx");

	std::string 								plotame 								=	"Efficiency_vs_Nvtx_eta" + findAndReplaceAll(removeTrailingZeros(_etaMin) + "to" + removeTrailingZeros(_etaMax), ".", "p")+ "_pt" + findAndReplaceAll(removeTrailingZeros(_ptMin) + "to" + removeTrailingZeros(_ptMax), ".", "p");
	TEfficiency									sigEff95(("signal95_"+plotame).c_str(), "", std::stoi(options.getList("nVtxBinning")[0]), options.getFloatList("nVtxBinning")[1], options.getFloatList("nVtxBinning")[2]);
	TEfficiency									sigEff90(("signal90_"+plotame).c_str(), "", std::stoi(options.getList("nVtxBinning")[0]), options.getFloatList("nVtxBinning")[1], options.getFloatList("nVtxBinning")[2]);
	TEfficiency									sigEff80(("signal80_"+plotame).c_str(), "", std::stoi(options.getList("nVtxBinning")[0]), options.getFloatList("nVtxBinning")[1], options.getFloatList("nVtxBinning")[2]);
	TEfficiency									sigEff70(("signal70_"+plotame).c_str(), "", std::stoi(options.getList("nVtxBinning")[0]), options.getFloatList("nVtxBinning")[1], options.getFloatList("nVtxBinning")[2]);
	
	TEfficiency									bgEff95(("background95_"+plotame).c_str(), "", std::stoi(options.getList("nVtxBinning")[0]), options.getFloatList("nVtxBinning")[1], options.getFloatList("nVtxBinning")[2]);
	TEfficiency									bgEff90(("background90_"+plotame).c_str(), "", std::stoi(options.getList("nVtxBinning")[0]), options.getFloatList("nVtxBinning")[1], options.getFloatList("nVtxBinning")[2]);
	TEfficiency									bgEff80(("background90_"+plotame).c_str(), "", std::stoi(options.getList("nVtxBinning")[0]), options.getFloatList("nVtxBinning")[1], options.getFloatList("nVtxBinning")[2]);
	TEfficiency									bgEff70(("background90_"+plotame).c_str(), "", std::stoi(options.getList("nVtxBinning")[0]), options.getFloatList("nVtxBinning")[1], options.getFloatList("nVtxBinning")[2]);

	sigEff95.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff90.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff80.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	sigEff70.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));

	bgEff95.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff90.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff80.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));
	bgEff70.SetStatisticOption((TEfficiency::EStatOption)options.getInt("statOption"));	

	while(inputTTreeReader.Next()){

		if(isTrain_) continue;
		
		Float_t 							absSCEta									=	std::abs(phoSCeta_);
		if (absSCEta <= _etaMin) continue;
		if (absSCEta > _etaMax) continue;

		if(phoPt_ < _ptMin) continue;
		if(phoPt_ > _ptMax) continue;

		if(isSignal_){
			sigEff95.FillWeighted(pass95_, bdtWeightF_, 0.2+nVtx_);
			sigEff90.FillWeighted(pass90_, bdtWeightF_, 0.2+nVtx_);
			sigEff80.FillWeighted(pass80_, bdtWeightF_, 0.2+nVtx_);
			sigEff70.FillWeighted(pass70_, bdtWeightF_, 0.2+nVtx_);
		} else{
			bgEff95.FillWeighted(pass95_, bdtWeightF_, 0.2+nVtx_);
			bgEff90.FillWeighted(pass90_, bdtWeightF_, 0.2+nVtx_);
			bgEff80.FillWeighted(pass80_, bdtWeightF_, 0.2+nVtx_);
			bgEff70.FillWeighted(pass70_, bdtWeightF_, 0.2+nVtx_);
		}
	}

	closeTChain(featsTree);
	closeTChain(bdtTree);
	closeTChain(idTree);

	TGraphAsymmErrors* 							sigEff95Graph 							=	sigEff95.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff90Graph 							=	sigEff90.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff80Graph 							=	sigEff80.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							sigEff70Graph 							=	sigEff70.CreateGraph(options.getCSTR("gOpt"));

	sigEff95Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff95Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff95Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff95Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff90Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff90Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff90Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff90Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff80Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff80Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff80Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff80Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff70Graph->SetLineStyle(options.getInt("signalLineStyle"));
	sigEff70Graph->SetLineWidth(options.getInt("signalLineWidth"));
	sigEff70Graph->SetMarkerSize(options.getFloat("signalMkrSize"));
	sigEff70Graph->SetMarkerStyle(options.getInt("signalMkrStyle"));

	sigEff95Graph->SetLineColor(options.getTColFromHex("eff95Col"));
	sigEff90Graph->SetLineColor(options.getTColFromHex("eff90Col"));
	sigEff80Graph->SetLineColor(options.getTColFromHex("eff80Col"));
	sigEff70Graph->SetLineColor(options.getTColFromHex("eff70Col"));

	sigEff95Graph->SetMarkerColor(options.getTColFromHex("eff95Col"));
	sigEff90Graph->SetMarkerColor(options.getTColFromHex("eff90Col"));
	sigEff80Graph->SetMarkerColor(options.getTColFromHex("eff80Col"));
	sigEff70Graph->SetMarkerColor(options.getTColFromHex("eff70Col"));

	TMultiGraph 								sStack;
	sStack.Add(sigEff95Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff90Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff80Graph, options.getCSTR("dOpt"));
	sStack.Add(sigEff70Graph, options.getCSTR("dOpt"));
	
	TLegend sLegend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	sLegend.SetTextSize(options.getDouble("legendTextSize"));
	sLegend.SetNColumns(options.getInt("legendNcols"));
	sLegend.SetBorderSize(0);
	sLegend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	sLegend.SetFillStyle(options.getInt("legFillStyle"));
	std::string 								sLegTitle 								=	"Prompt";
	sLegend.SetHeader(sLegTitle.c_str(), "C");
	sLegend.SetTextAlign(12);
	sLegend.AddEntry(sigEff95Graph, "Very Loose", "LPE");
	sLegend.AddEntry(sigEff90Graph, "Loose", "LPE");
	sLegend.AddEntry(sigEff80Graph, "Medium", "LPE");
	sLegend.AddEntry(sigEff70Graph, "Tight", "LPE");

	TCanvas canvas((plotame + "_canvas").c_str(), "", options.getDouble("canvasX"), options.getDouble("canvasY"));
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

	sStack.Draw(options.getCSTR("dOpt"));
	sLegend.Draw();

	sStack.GetXaxis()->SetTitle("N_{vtx}");
	sStack.GetYaxis()->SetTitle("Efficiency");

	sStack.GetXaxis()->CenterTitle();
	sStack.GetYaxis()->CenterTitle();
	sStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	sStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	sStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	sStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	sStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	sStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	sStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	sStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	sStack.SetMaximum(options.getFloat("yScaleupSig")*sStack.GetYaxis()->GetXmax());

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "signal_" + plotame+".png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "signal_" + plotame+".pdf").c_str());

	TGraphAsymmErrors* 							bgEff95Graph 							=	bgEff95.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff90Graph 							=	bgEff90.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff80Graph 							=	bgEff80.CreateGraph(options.getCSTR("gOpt"));
	TGraphAsymmErrors* 							bgEff70Graph 							=	bgEff70.CreateGraph(options.getCSTR("gOpt"));

	bgEff95Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff95Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff95Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff95Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff90Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff90Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff90Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff90Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff80Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff80Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff80Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff80Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff70Graph->SetLineStyle(options.getInt("backgroundLineStyle"));
	bgEff70Graph->SetLineWidth(options.getInt("backgroundLineWidth"));
	bgEff70Graph->SetMarkerSize(options.getFloat("backgroundMkrSize"));
	bgEff70Graph->SetMarkerStyle(options.getInt("backgroundMkrStyle"));

	bgEff95Graph->SetLineColor(options.getTColFromHex("eff95Col"));
	bgEff90Graph->SetLineColor(options.getTColFromHex("eff90Col"));
	bgEff80Graph->SetLineColor(options.getTColFromHex("eff80Col"));
	bgEff70Graph->SetLineColor(options.getTColFromHex("eff70Col"));

	bgEff95Graph->SetMarkerColor(options.getTColFromHex("eff95Col"));
	bgEff90Graph->SetMarkerColor(options.getTColFromHex("eff90Col"));
	bgEff80Graph->SetMarkerColor(options.getTColFromHex("eff80Col"));
	bgEff70Graph->SetMarkerColor(options.getTColFromHex("eff70Col"));

	TMultiGraph 								bStack;
	bStack.Add(bgEff95Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff90Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff80Graph, options.getCSTR("dOpt"));
	bStack.Add(bgEff70Graph, options.getCSTR("dOpt"));

	TLegend bLegend(options.getDouble("legendx1"), options.getDouble("legendy1"), options.getDouble("legendx2"), options.getDouble("legendy2"));
	bLegend.SetTextSize(options.getDouble("legendTextSize"));
	bLegend.SetNColumns(options.getInt("legendNcols"));
	bLegend.SetBorderSize(0);
	bLegend.SetFillColorAlpha(options.getTColFromHex("legFillColor"), options.getFloat("legFillColorAlpha"));
	bLegend.SetFillStyle(options.getInt("legFillStyle"));
	std::string 								bLegTitle 								=	"Fake";
	bLegend.SetHeader(bLegTitle.c_str(), "C");
	bLegend.SetTextAlign(12);
	bLegend.AddEntry(bgEff95Graph, "Very Loose", "LPE");
	bLegend.AddEntry(bgEff90Graph, "Loose", "LPE");
	bLegend.AddEntry(bgEff80Graph, "Medium", "LPE");
	bLegend.AddEntry(bgEff70Graph, "Tight", "LPE");

	pad.Clear();
	pad.cd();


	bStack.Draw((options.get("dOpt")).c_str());
	bStack.SetMaximum(options.getFloat("yScaleupBG")*bStack.GetYaxis()->GetXmax());

	bStack.Draw(options.getCSTR("dOpt"));
	bLegend.Draw();

	bStack.GetXaxis()->SetTitle("N_{vtx}");
	bStack.GetYaxis()->SetTitle("Efficiency");

	bStack.GetXaxis()->CenterTitle();
	bStack.GetYaxis()->CenterTitle();
	bStack.GetXaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	bStack.GetYaxis()->SetTitleSize(options.getDouble("pad0axisTitleSize"));
	bStack.GetXaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	bStack.GetYaxis()->SetLabelSize(options.getDouble("pad0axisLabelSize"));
	bStack.GetXaxis()->SetTitleOffset(options.getDouble("pad0XtitleOffset"));
	bStack.GetYaxis()->SetTitleOffset(options.getDouble("pad0YtitleOffset"));	
	bStack.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	bStack.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "background_" + plotame+".png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "background_" + plotame+".pdf").c_str());
	
	clearHeap();

};

void initialize(){
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	options.parseIt(optionsFile, "==");

	mkdir(options.get("writeDir"));

	std::vector<Double_t> 						pTbinning										= 	options.getDoubleList("pTbinning");
	std::vector<Float_t> 						etaRange										= 	options.getFloatList("etaRange");
	std::vector<Float_t> 						ptRange											= 	options.getFloatList("ptRange");

	
	std::vector<Double_t> 						etaBinning										= 	options.getDoubleList("etaBinning");


	// getPtEfficiency(pTbinning, etaRange[0], etaRange[1]);
	// getEtaEfficiency(etaBinning, ptRange[0], ptRange[1]);

	getNvtxEfficiency(ptRange[0], ptRange[1], etaRange[0], etaRange[1]);

	
	std::cout<<"Done"<<std::endl;
};
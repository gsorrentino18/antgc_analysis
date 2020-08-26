#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"
std::string 										optionsFile="/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/training/plotFeatImpOptions.txt";
parseOptions 										options;
std::map<std::string,std::vector<std::string>>		varPlotInfo;

void plotFeatImps(){
	gROOT->SetBatch();
	options.parseIt(optionsFile, "==");

	CSVReader 											varFile(options.get("brSpecFile"), "==");
	std::vector<std::vector<std::string>>				varDat 		= varFile.getData();
	for(const std::vector<std::string> & iRow : varDat){
		if(iRow.size()==3) {
			std::vector<std::string> iVarPlotInfo = {iRow[1]};
			std::vector<std::string> binInfo = split_string(iRow[2]);
			iVarPlotInfo.insert(iVarPlotInfo.end(),binInfo.begin(), binInfo.end());
			varPlotInfo[iRow[0]]		= iVarPlotInfo;
		}
	};

	CSVReader 											featsFile(options.get("featsPath"), ",");
	std::vector<std::vector<std::string>>				featsDat 		= featsFile.getData();

	Int_t nBins = featsDat.size();

	TH1F featBar("featBar", "", nBins, 0., nBins);
	for(Int_t iRow = 0; iRow < nBins; iRow++){
		Int_t iBin = nBins - iRow;
		featBar.Fill(-0.1+iBin, std::stof(featsDat[iRow][1]));

		featBar.GetXaxis()->SetBinLabel(iBin, varPlotInfo[featsDat[iRow][0]][0].c_str());
	}

	featBar.GetXaxis()->SetLabelSize(options.getDouble("pad0YaxisLabelSize"));
	featBar.GetYaxis()->SetLabelSize(options.getDouble("pad0XaxisLabelSize"));
	featBar.GetYaxis()->SetNdivisions(options.getInt("pad0yNdivs"));
	featBar.GetXaxis()->SetNdivisions(options.getInt("pad0xNdivs"));
	featBar.GetYaxis()->SetTitle("Feature Importance");
	featBar.GetYaxis()->SetTitleOffset(options.getFloat("pad0XtitleOffset"));
	featBar.GetYaxis()->SetTitleSize(options.getFloat("pad0axisTitleSize"));

	featBar.SetMaximum(options.getFloat("maxFactor")*featBar.GetMaximum());
	

	featBar.SetFillColor(options.getTColFromHex("col"));
	featBar.SetBarWidth(options.getFloat("barW"));
	featBar.SetBarOffset(options.getFloat("barSep"));
	featBar.SetStats(0);

	TCanvas canvas("canvas", "", options.getDouble("canvasX"), options.getDouble("canvasY"));
	canvas.SetFillStyle(4000);
	
	TPad pad0("pad0", "", options.getDouble("pad0x1"), options.getDouble("pad0y1"), options.getDouble("pad0x2"), options.getDouble("pad0y2"));
	pad0.SetMargin(options.getDouble("pad0marginL"), options.getDouble("pad0marginR"), options.getDouble("pad0marginB"), options.getDouble("pad0marginT"));
	pad0.SetFillStyle(4000);
	pad0.SetFillColor(0);
	pad0.SetFrameFillStyle(4000);
	pad0.SetGrid(1,1);

	canvas.Draw();
	canvas.cd();
	pad0.Draw();
	pad0.cd();

	featBar.Draw(options.getCSTR("drawOption"));

	mkdir(options.get("writeDir"));

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	canvas.RedrawAxis();
	canvas.Update();
	canvas.Modified();

	canvas.SaveAs((options.get("writeDir")+ "FeatImportances.png").c_str());
	canvas.SaveAs((options.get("writeDir")+ "FeatImportances.pdf").c_str());
	
	clearHeap();
};

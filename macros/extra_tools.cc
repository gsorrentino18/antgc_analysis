//////////////////////////////////////////////////////////////////////
//	Author: Mohammad Abrar Wadud, Univeristy of Minnesota			//
//	Date: Jan/03/2020												//
//////////////////////////////////////////////////////////////////////

#ifndef EXTRATOOLS_H
#define EXTRATOOLS_H

// #include "KDEProducer1D.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include "TKey.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TAxis.h"
#include "TObjArray.h"
#include "TColor.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TCollection.h"
#include "TKey.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TColor.h"
#include "RtypesCore.h"
#include <boost/algorithm/string.hpp>
#include "boost/algorithm/string/replace.hpp"
#include <boost/algorithm/string_regex.hpp>
#include <boost/type_index.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/lexical_cast/bad_lexical_cast.hpp"
#include <sys/stat.h>
#include "iostream"
#include "sstream"
#include "fstream"
#include "vector"
#include "math.h"
#include <algorithm>
#include <cctype>
#include <locale>
#include <numeric>
#include <regex>
#include <dirent.h>
#include <errno.h>
#include <dirent.h>
#include <stdio.h>
#include <iterator>
#include <string>
#include <chrono>
#include <ctime>
#include "limits"
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>



/*************************************************************Declarations*************************************************************/
Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t deltaPhi(Double_t phi1, Double_t phi2);
template <typename anytype1, typename anytype2>
void setBit(anytype1 & _container, anytype2 _pos, Bool_t _bitVal);
template <typename anytype1, typename anytype2>
Bool_t getBit(anytype1 _container, anytype2 _pos);
std::string removeNonAlpha(std::string word);
template <class any_number>
std::string removeTrailingZeros(any_number number);
Bool_t file_exists(std::string fileName);
Int_t mkdir(std::string dir_path);
std::map<std::string, Double_t> load_xsecs(std::string filepath);
std::vector<std::string> getObjectList(std::string filepath, std::string objtype, std::vector<std::string> exclusion_list={});
Bool_t match(std::string _needle, std::string _haystack);
std::string ReadNthLine(std::string filename, int N);
UInt_t countLines(std::string filename);
std::vector<std::string> split_string(std::string _string, std::string _delimiter=",", Bool_t _trim=1);
std::string get_cell(std::string filename, UInt_t row, UInt_t column, std::string _delimiter=",");
std::string getFileName(std::string _filepath);
std::string getDirPath(std::string _somePath);
std::vector<std::string> getNonemptyLines(std::string filepath);
std::vector<std::string> getNonemptyLinesWithFilenameKeyword(std::string filepath, std::string keyword, std::string exclude="");
std::vector<std::string> getLinesRegex(std::string _filepath, std::string _regexStr);
TH1* getHistFromFile(std::string _histName, std::string _filename, Bool_t _verbose=1, TDirectory *_dir = nullptr);
TObject *getObjectFromFile(std::string _objectName, std::string _filename);
TH1* rebinNHist(TH1* _hist, Int_t _N);
TH1* rebinHist(TH1* _hist, std::vector<Double_t> _newBins);
std::vector<Double_t> getGoodBins(TH1* _hist, Double_t _statUnc, Double_t _reScale = -999., Int_t _nbinPar=10);
Double_t sumNextNbins(TH1* _hist, Int_t _n, Int_t _curr);
void copyHistAtts(TH1* _source, TH1* _mock);
Double_t getSumW(std::string _cutflowfile);
void ltrim(std::string &s);
void rtrim(std::string &s);
void trim(std::string &s);
std::string ltrim_copy(std::string s);
std::string rtrim_copy(std::string s);
std::string trim_copy(std::string s);
Double_t ams(Double_t _s, Double_t _b);
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 0);
std::string getUnit(TH1* _hist);
std::string getUnit(std::string ystring);
std::string eraseUnit(std::string ystring);
std::string first_numberstring(std::string const & str);
std::vector<Double_t> getXbins(TH1 *hist);
void closeTChain(TChain * & _chain);
void setFrameColor(TAxis* _axis, std::string _color);
void setFrameColor(TH1* _hist, std::string _color);
void setFrameColor(THStack* _stack, std::string _color);
TTree *loadTreeFromFile(std::string _treeName, std::string _filename);
Char_t isDirectory(std::string filePath, Bool_t _verbose=0);
void addPointToGraph(TGraph & _graph, Double_t _x, Double_t _y);
void graphStats(TGraphAsymmErrors* graph, Double_t &mean, Double_t &stdev);
TH1D* graph2hist(TGraphAsymmErrors* graph, UInt_t ndivs, Double_t ylow, Double_t yhigh);
Int_t writeToFile(TObject *_object, std::string _filePath, std::string _mode = "RECREATE", Bool_t _verbose=1);
std::string findAndReplaceAll(std::string data, std::string toSearch, std::string replaceStr);
TChain *openTChain(std::string _chainListFile, std::string _treeName="", Bool_t _verbose = 1);
TChain *openTChain(std::vector<std::string> _chainList, std::string _treeName="", Bool_t _verbose = 1);
TChain *openTChainWithFilesInDir(std::string _dirPath, std::string _treeName="");
std::vector<std::pair<std::string, std::string>> getBranchList(std::string _treePath, std::string _treeName="", Bool_t _verbose=1);
std::vector<std::string> listFilesInDir(std::string _dirPath, std::string _regexStr="", Bool_t _verb=0);
Bool_t matchRegex(std::string _str,std::string _regexStr);
std::string getTreeNameInFile(std::string _filePath);
TH1F *mergeBins(std::string _fileList, std::string _histName, std::string _sumWeightsHistname, std::string _xsecMap, Int_t _nameCol=0, Int_t _xSecCol=2, std::string _path="");
TH1F *mergeBins(std::vector<std::string> _fileList, std::string _histName, std::string _sumWeightsHistname, std::string _xsecMap, Int_t _nameCol=0, Int_t _xSecCol=2, std::string _path="");
std::string vLookup(std::string _lookupKey, std::string _inFile, Int_t _lookupCol, Int_t _valCol, Bool_t _regex=0);
Bool_t isROOTfile(std::string _filepath);
std::vector<Float_t> getXlimits(std::vector<TH1*> _hists, Float_t _binThreshold=0.);
void clearHeap();
Double_t weightedYmean(TH1 *_hist);
Bool_t branchExists(std::string _branchName, TTree *_tree);
Float_t getMean(std::vector<Float_t> _set);
template <class ObjType>
ObjType copyObjectDeleteSrc(ObjType *_original);
template<typename T1, typename T2>
Int_t findSecondaryIndex(T1 searchIndex, std::vector<T2> container);
Short_t findSecondaryIndex(Short_t searchIndex, std::vector<Short_t> container);
std::string getCurrentTime();
template <typename anytype>
void eraseElement(std::vector<anytype> & _vec, UInt_t _Index2Erase);
Double_t getCategoryBoundary(TH1*_signal, TH1*_background);
std::vector<Double_t> vecString2vecDouble(std::vector<std::string> _numStrings);
Bool_t stringIsNumber(std::string _isThisANumber);
TDirectory *mkTFileDir(TFile *_file, std::string _dir);
Int_t caselessSubstringPos( std::string str1, std::string str2);
Bool_t findStringIC(const std::string & strHaystack, const std::string & strNeedle);
Bool_t stringIsEmpty(std::string str);
template<typename T>
UChar_t setLineAtts(T* graph, std::string _atts, std::string _delimiter=",");
template<typename T>
UChar_t setMarkerAtts(T* graph, std::string _atts, std::string _delimiter=",");
template<typename T>
UChar_t setFillAtts(T* graph, std::string _atts, std::string _delimiter=",");
UChar_t setAxisAtts(TAxis* _axis, std::string _atts, std::string _delimiter=",");
template<typename T>
UChar_t setPadAtts(T* _pad, std::string _atts, std::string _delimiter=",");
template<typename T>
void setFrameAtts(T* _hist, std::string _style, std::string _del1=";", std::string _del2=",");
template<typename T>
void setHistStyle(T* _hist, std::string _style, std::string _del1=";", std::string _del2=",");
std::string stringToUpper(std::string _str);
std::string stringToLower(std::string _str);
template<typename T1, typename T2>
Int_t findIndex(const std::vector<T1> & _haystack, T2 _needle);
Bool_t objectExists(TFile *_File, std::string _objectName);
Bool_t objectExists(std::string _FilePath, std::string _objectName);
std::vector<Float_t> strToFloatList(std::string _listString, std::string _delimiter=",");
std::vector<Double_t> strToDoubleList(std::string _listString, std::string _delimiter=",");
template<typename T>
void setAxes(T* _graph, std::string _axesOptions){};
Int_t hex2rootColor(std::string _hexCode);
std::string sysexec(std::string cmd);
std::vector<Double_t> getNpercentMinInterval(TH1 *_hist, Float_t _N);
std::vector<Double_t> getNpercentLinearInterval(TH1 *_hist, Float_t _N);
std::vector<std::string> prefixVecString(std::vector<std::string> _vecStr, std::string _prefix, Int_t _insPos = 0);
void normalizeHist(TH1* _hist, Double_t _norm=100.);
void normalizeHist(TH1& _hist, Double_t _norm=100.);
std::vector<Int_t> strToIntList(std::string _listString, std::string _delimiter=",");
Double_t spearmanR(TH2 *_hist, Int_t _Nbins = 1000);
Float_t relInvMass(Float_t _pT1, Float_t _eta1, Float_t _phi1, Float_t _pT2, Float_t _eta2, Float_t _phi2);
template<typename T>
void remove_intersection(std::vector<T>& a, std::vector<T>& b);
TH1* getSqrtHist(TH1* _hist);
void removeNegativeBins(TH1* _hist);
Float_t foldedPhi(Float_t _phi);
bool isInteger(const std::string & s);
Char_t RGB2ColPlt(std::string _CSVFile, Int_t _firstCol);
TH2* getAbsHist(TH2* _hist);

// void setTFileDir(std::vector<TH1*> _hists);

// void setTFileDir(std::vector<TH1&> _hists, TFile & _file, std::string _dir = ""){
// 	_file.mkdir(_dir.c_str());
// 	for(TH1 & iHist : _hists){
// 		iHist.SetDirectory(_file.GetDirectory(_dir.c_str()));
// 	}
// };



template <typename anytype>
struct TTreeReaderAnyValue;
template <typename anytype>
struct TTreeReaderVectorValue;
template <typename anytype>
struct TTreeReaderArrayValue;
struct TTreePathok;
struct plot_variable;
struct histogram_template;
struct twoDhistogram_template;
struct bit_histogram_template;
struct bit_twoDhistogram_template;
template <typename anytype>
struct vector_association;
struct BinCollection;
struct signal_atts;
struct sample;
struct Profile2D;
struct parseOptions;
struct bitBar2D;
struct effectiveAreaMap;
struct isoPtScalingMap;
struct coronaCorrections;
struct BGset;
template <typename anytype>
struct smarter_ptr;
template <typename anytype>
struct indexer;
struct DinkyHister;
struct logStream;
struct correlationMatix;


class CSVReader;
class PileupReWeighting;


const std::map<std::string, Double_t> xsec_unit_map = {
	{"fb", 1.0e-3},
	{"pb", 1.0},
	{"nb", 1.0e3}
};

// struct TTreePathok{
// 	TTreePathok(){};
// 	TTreePathok(TChain* _tree){
// 		tree =  _tree;
// 		pathok.SetTree(tree);
// 	};

// 	TTreePathok(std::string _fileListPath, std::string _treeName){
// 		tree =  openTChain(_fileListPath, _treeName);
// 		pathok.SetTree(tree);
// 	};

// 	void setEntry(Long64_t _ientry){
// 		currentEntry = _ientry;
// 	};

// 	TChain* tree = nullptr;
// 	TTreeReader pathok;
// 	Long64_t	currentEntry = -999;

// 	template <typename anytype>
// 	anytype getSingleton(std::string _branchName){
// 		TTreeReaderAnyValue<anytype> singleton(pathok, _branchName);
// 		pathok.SetEntry(currentEntry);
// 		return (anytype) singleton;
// 	};	
// };


/*************************************************************Definitions*************************************************************/
struct eventCategory {
	std::string name;
	Bool_t *Pass = nullptr;
	eventCategory(Bool_t &_Pass, std::string _name): Pass(&_Pass), name(_name){
	};
	eventCategory(){};
	~eventCategory(){};
};

template <typename anytype1, typename anytype2>
void setBit(anytype1 & _container, anytype2 _pos, Bool_t _bitVal){
	_container ^= (-_bitVal ^ _container) & (1UL << _pos);
};


template <typename anytype1, typename anytype2>
Bool_t getBit(anytype1 _container, anytype2 _pos){
	return (_container>>_pos) & 1;
};


template <typename anytype>
struct TTreeReaderAnyValue{
	TTreeReaderValue<anytype> *val = nullptr;
	TTreeReaderAnyValue(TTreeReader & ttreereader, std::string branchname){
		set(ttreereader, branchname);
	};
	TTreeReaderAnyValue(){};
	~TTreeReaderAnyValue(){
		delete val;
		val = nullptr;
	};
	void set(TTreeReader & ttreereader, std::string branchname){
		// if(!branchExists(branchname, ttreereader.GetTree())){
		// 	std::cout<<"\t\tError! Branch "<<branchname<<" does not exist in TTree "<<ttreereader.GetTree()->GetName()<<std::endl;
		// 	exit(EXIT_FAILURE);
		// }
		delete val;
		val = new TTreeReaderValue<anytype>(ttreereader,branchname.c_str());
		std::cout<<"\t\tInitialized branch <"<< boost::typeindex::type_id<anytype>().pretty_name() <<"> "<<branchname<<std::endl;
	};
	operator anytype () const {
		return **val;
	}
	anytype get(){
		return **val;
	};
};


template <typename anytype>
struct TTreeReaderVectorValue: TTreeReaderAnyValue<std::vector<anytype>>{
	TTreeReaderVectorValue(TTreeReader & ttreereader, std::string branchname):TTreeReaderAnyValue<std::vector<anytype>>(ttreereader, branchname){
	};
	TTreeReaderVectorValue(){};
	~TTreeReaderVectorValue(){
	};
	anytype operator[](UInt_t index){
		return (*(this->val))->at(index);
	};
	anytype at(UInt_t index){
		return (*(this->val))->at(index);
	};
	UInt_t size(){
		return (*(this->val))->size();
	};

	//define begin & end function for iterator
};


template <typename anytype>
struct TTreeReaderArrayValue{
	TTreeReaderArray<anytype> *val = nullptr;
	TTreeReaderArrayValue(TTreeReader & ttreereader, std::string branchname){
		set(ttreereader, branchname);
	};
	TTreeReaderArrayValue(){};
	~TTreeReaderArrayValue(){
		delete val;
		val = nullptr;
	};
	void set(TTreeReader & ttreereader, std::string branchname){
		// if(!branchExists(branchname, ttreereader.GetTree())){
		// 	std::cout<<"\t\tError! Branch "<<branchname<<" does not exist in TTree "<<ttreereader.GetTree()->GetName()<<std::endl;
		// 	exit (EXIT_FAILURE);
		// }
		val = new TTreeReaderArray<anytype>(ttreereader,branchname.c_str());
		std::cout<<"\t\tInitialized array branch <"<< boost::typeindex::type_id<anytype>().pretty_name() <<"> "<<branchname<<std::endl;
	};
	anytype operator[](UInt_t index){
		return (*(this->val)).At(index);
	};
	anytype at(UInt_t index){
		return (*(this->val)).At(index);
	};
	UInt_t size(){
		return (*(this->val)).GetSize();
	};
	//returns pointer to 0th position
	operator anytype* () const{
		return &(*(*(this->val)).begin());
	};
};


template <typename anytype>
struct vector_association{
	std::vector<anytype> *vec = nullptr;
	anytype *var = nullptr;
	vector_association(anytype *_var, std::vector<anytype> *_vec): var(_var), vec(_vec){};
	vector_association(){};
	~vector_association(){};
	anytype operator[](UInt_t index){
		return vec->at(index);
	}
	void push_back(){
		vec->push_back(*var);
	};
	void clear(){
		vec->clear();
	};
};


struct plot_variable {
	Float_t *xptr = nullptr;
	Bool_t destroyPtr = 0;
	Float_t xmin;
	Float_t xmax;
	ULong_t nbins;
	Double_t *xBins = nullptr;
	std::string xtitle;
	std::string xunit;

	plot_variable(Float_t _xmin, Float_t _xmax, ULong_t _nbins,	std::string _xtitle = "", std::string _xunit = ""){
		set(_xmin, _xmax, _nbins, _xtitle, _xunit);
	};

	plot_variable(const Double_t *_xBins, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = ""){
		set(_xBins, _nbins, _xtitle, _xunit);
	};

	plot_variable(Float_t &_xptr, Float_t _xmin, Float_t _xmax, ULong_t _nbins,	std::string _xtitle = "", std::string _xunit = ""){
		set(_xptr, _xmin, _xmax, _nbins, _xtitle, _xunit);
	};

	plot_variable(Float_t &_xptr, const Double_t *_xBins, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = ""){
		set(_xptr, _xBins, _nbins, _xtitle, _xunit);
	};

	plot_variable(){};

	void set(Float_t _xmin, Float_t _xmax, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = ""){
		xptr = new Float_t;
		destroyPtr = 1;
		xmin = _xmin;
		xmax = _xmax;
		nbins = _nbins;
		xtitle = _xtitle;
		xunit = _xunit;
	};

	void set(const Double_t *_xBins, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = ""){
		xptr = new Float_t;
		destroyPtr = 1;
		nbins = _nbins;
		xtitle = _xtitle;
		xunit = _xunit;
		xBins = new Double_t[nbins+1];
		for(ULong_t i = 0; i <= _nbins; i++){
			xBins[i] = _xBins[i];
		}
	};

	void set(Float_t &_xptr, Float_t _xmin, Float_t _xmax, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = ""){
		xptr = &_xptr;
		xmin = _xmin;
		xmax = _xmax;
		nbins = _nbins;
		xtitle = _xtitle;
		xunit = _xunit;
	};

	void set(Float_t &_xptr, const Double_t *_xBins, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = ""){
		xptr = &_xptr;
		nbins = _nbins;
		xtitle = _xtitle;
		xunit = _xunit;
		xBins = new Double_t[nbins+1];
		for(ULong_t i = 0; i <= _nbins; i++){
			xBins[i] = _xBins[i];
		}
	};
	
	~plot_variable(){
		if(destroyPtr) delete xptr;
		delete [] xBins;
	};

	operator Float_t (){
		return *xptr;
	};

	void operator = (const Float_t & assignVal){ 
		*xptr = assignVal;
	};
};


struct histogram_template {
	const plot_variable *var = nullptr;
	const Float_t * varPtr = nullptr;
	TH1F* hist = nullptr;
	std::string histTitle;
	std::string histName;
	histogram_template(const plot_variable &_var, std::string _histTitle="", std::string _histName="", TH1F* _hist=nullptr){
		set(_var, _histTitle, _histName, _hist);
	};

	histogram_template(){};

	void set(const plot_variable &_var, std::string _histTitle="", std::string _histName="", TH1F* _hist=nullptr){
		var = &_var;
		varPtr = _var.xptr;
		histTitle = _histTitle;
		hist = _hist;
		if(hist) isUser = 1;
	};

	void initializehist(std::string name_prefix="", std::string title_prefix="", Bool_t Yunit = 1, TFile *_file = nullptr, std::string _TFileDir=""){
		if(isUser) std::cout<<"\tWarning: Data member .hist will not point to user-set histogram!"<<std::endl;
		std::string binwidth_string = "";
		if(var->xBins == nullptr && Yunit){
			binwidth_string = removeTrailingZeros((var->xmax - var->xmin)/(Float_t)var->nbins);
			if(!binwidth_string.compare("1") && !var->xunit.empty()) binwidth_string.pop_back();
			else binwidth_string += " ";
			binwidth_string = "/"+binwidth_string;
			binwidth_string += var->xunit;
			rtrim(binwidth_string);
		}
		if(histName.empty()) histName = removeNonAlpha(name_prefix) + "_1D_" + removeNonAlpha(var->xtitle);
		trim(histName);
		trim(histTitle);
		std::string titles = title_prefix + histTitle + ";" + var->xtitle + (var->xunit.empty()?"":(" ["+var->xunit+"]")) + ";" + "# of events"+binwidth_string;
		if(var->xBins != nullptr) hist = new TH1F(histName.c_str(), titles.c_str(), var->nbins, var->xBins);
		else hist = new TH1F(histName.c_str(), titles.c_str(), var->nbins, var->xmin, var->xmax);
		hist->GetXaxis()->CenterTitle();
		hist->GetYaxis()->CenterTitle();
		hist->Sumw2();
		if(_file){
			mkTFileDir(_file, _TFileDir);
			hist->SetDirectory(_file->GetDirectory(_TFileDir.c_str()));
		} 
		isUser = 0;
		std::cout<<"\t\tInitialized TH1F "<<histName<<std::endl;
	};

	void fill(Double_t weight = 1.0){
		if(!hist){
			std::cout<<"Cannot fill! TH1F ("<<var->xtitle <<") is uninitialized!";
			return;
		}
		hist->Fill(*varPtr, weight);
	};

	~histogram_template(){
		// if(hist) delete hist;
	}
private:
	Bool_t isUser = 0;
};


struct twoDhistogram_template {
	const plot_variable *xvar = nullptr;
	const plot_variable *yvar = nullptr;

	const Float_t * xPtr = nullptr;
	const Float_t * yPtr = nullptr;

	TH2F* hist = nullptr;
	std::string histTitle;
	std::string histName;

	twoDhistogram_template(const plot_variable &_xvar, const plot_variable &_yvar, std::string _histTitle = "", std::string _histName="", TH2F* _hist = nullptr){
		set(_xvar, _yvar, _histTitle, _histName, _hist);
	};

	twoDhistogram_template(){};

	void set(const plot_variable &_xvar, const plot_variable &_yvar, std::string _histTitle = "", std::string _histName="", TH2F* _hist = nullptr){
		xvar = &_xvar;
		yvar = &_yvar;
		xPtr = _xvar.xptr;
		yPtr = _yvar.xptr;
		histTitle = _histTitle;
		hist = _hist;
		if(hist) isUser = 1;
	};

	void fill(Double_t weight = 1.0){
		if(!hist){
			std::cout<<"Cannot fill! TH2F ("<<xvar->xtitle << " VS "<< yvar->xtitle<<") is uninitialized!";
			return;
		}
		hist->Fill(*xPtr, *yPtr, weight);
	};

	void initializehist(std::string name_prefix = "", std::string title_prefix="", TFile *_file = nullptr, std::string _TFileDir=""){
		if(isUser) std::cout<<"\tWarning: Data member .hist will not point to user-set histogram!"<<std::endl;
		if(histName.empty())histName = removeNonAlpha(name_prefix) + "_2D_" + removeNonAlpha(xvar->xtitle +"\\ VS\\ " + yvar->xtitle);
		trim(histName);
		trim(histTitle);
		std::string titles = title_prefix + histTitle + ";" + xvar->xtitle + (xvar->xunit.empty()?"":(" ["+xvar->xunit+"]")) + ";" + yvar->xtitle + (yvar->xunit.empty()?"":(" ["+yvar->xunit+"]")) ;
		if(xvar->xBins != nullptr && yvar->xBins == nullptr){
			hist = new TH2F(histName.c_str(), titles.c_str(), xvar->nbins, xvar->xBins, yvar->nbins, yvar->xmin, yvar->xmax);
		}
		else if(xvar->xBins == nullptr && yvar->xBins != nullptr){
			hist = new TH2F(histName.c_str(), titles.c_str(), xvar->nbins, xvar->xmin, xvar->xmax, yvar->nbins, yvar->xBins);
		}
		else if(xvar->xBins != nullptr && yvar->xBins != nullptr){
			hist = new TH2F(histName.c_str(), titles.c_str(), xvar->nbins, xvar->xBins, yvar->nbins, yvar->xBins);
		}
		else{
			hist = new TH2F(histName.c_str(), titles.c_str(), xvar->nbins, xvar->xmin, xvar->xmax, yvar->nbins, yvar->xmin, yvar->xmax);
		}
		hist->GetXaxis()->CenterTitle();
		hist->GetYaxis()->CenterTitle();
		hist->SetTitle(histTitle.c_str());
		hist->Sumw2();

		if(_file){
			mkTFileDir(_file, _TFileDir);
			hist->SetDirectory(_file->GetDirectory(_TFileDir.c_str()));
		} 

		isUser = 0;
		std::cout<<"\t\tInitialized TH2F "<<histName<<std::endl;
	};

	~twoDhistogram_template(){
		// if(hist) delete hist;
	}
private:
	Bool_t isUser = 0;
};


struct bit_histogram_template: histogram_template {
	Bool_t *selector = nullptr;
	std::string prefix;
	bit_histogram_template(){};

	bit_histogram_template(const plot_variable &_var, Bool_t *_selector, std::string _prefix){
		bitHistSet(_var, _selector, _prefix);
	};

	void bitHistSet(const plot_variable &_var, Bool_t *_selector, std::string _prefix){
		selector = _selector;
		prefix = _prefix;
		set(_var);
	};

	void bitHistInit(){
		initializehist(prefix, prefix);
	};

	void fillBit(Double_t _weight=1.){
		if(*selector){
			fill(_weight);
		} else{
			hist->Fill(std::numeric_limits<Float_t>::max(), _weight);
		}
	};
};


struct bit_twoDhistogram_template: twoDhistogram_template{
	Bool_t *selector = nullptr;
	std::string prefix;
	bit_twoDhistogram_template(){};

	bit_twoDhistogram_template(const plot_variable &_xvar, const plot_variable &_yvar, Bool_t *_selector, std::string _prefix){
		bitHistSet(_xvar, _yvar, _selector, _prefix);
	};

	void bitHistSet(const plot_variable &_xvar, const plot_variable &_yvar, Bool_t *_selector, std::string _prefix){
		selector = _selector;
		prefix = _prefix;
		set(_xvar, _yvar);
	};

	void bitHistInit(){
		initializehist(prefix, prefix);
	};

	void fillBit(Double_t _weight=1.){
		if(*selector){
			fill(_weight);
		} else{
			hist->Fill(*(xvar->xptr), std::numeric_limits<Float_t>::max(), _weight);
		}
	};
};

struct signal_atts {
	signal_atts(std::string _couplingname, std::string _legend, std::string _color, Int_t _markerstyle) : couplingname(_couplingname), legend(_legend), color(_color), markerstyle(_markerstyle){
	}

	signal_atts(){
	}
	std::string couplingname;
	std::string legend;
	std::string color;
	Int_t markerstyle;
	// std::string operator[]{
	// 	return couplingname;
	// };
};


struct sample {

	sample(std::string _ntuple, std::string _legend = "", Int_t _marker = 20, std::string _color = "#252525", Bool_t _drawLine = 0, TFile *_file = NULL, Float_t _luminosity=-999){
		set(_ntuple, _legend, _marker, _color, _drawLine,_file, _luminosity);
	}

	sample(){
	}
	std::string ntuple;
	std::string legend;
	Int_t marker;
	Int_t lineStyle = 1;
	std::string color;
	Bool_t drawLine = 0;
	TFile *file;
	Float_t luminosity = 0.;

	void set(std::string _ntuple, std::string _legend = "", Int_t _marker = 20, std::string _color = "#252525", Bool_t _drawLine = 0, TFile *_file = NULL, Float_t _luminosity=-999){
		ntuple= _ntuple;
		legend = _legend;
		marker = _marker;
		color = _color;
		drawLine = _drawLine;
		file = _file;
		luminosity = _luminosity;

		std::cout<<"\tInitialized sample: "<<std::endl<<
		"\t\tntuple: "<<ntuple<<std::endl<<
		"\t\tlegend: "<<legend<<std::endl<<
		"\t\tmarker: "<<marker<<std::endl<<
		"\t\tcolor: "<<color<<std::endl<<
		"\t\tdrawLine: "<<drawLine<<std::endl<<
		"\t\tfile: "<<file<<std::endl<<
		"\t\tluminosity: "<<luminosity<<std::endl;
	};

	void assignAtt(TH1 *_hist, Float_t _markerSize=1.5, Float_t _lineWidth=2.){
		if(marker>0){
			_hist->SetMarkerStyle(marker);
			_hist->SetMarkerSize(_markerSize);
			_hist->SetMarkerColor(TColor::GetColor(color.c_str()));
			_hist->SetLineColor(TColor::GetColor(color.c_str()));
			_hist->SetLineWidth(_lineWidth);
		} else{
			_hist->SetFillStyle((-marker));
			_hist->SetFillColor(TColor::GetColor(color.c_str()));
			_hist->SetLineColor(TColor::GetColor(color.c_str()));
			_hist->SetLineWidth(_lineWidth);
		}
		_hist->SetLineStyle(lineStyle);
	};

	void setLineStyle(Int_t _lineStyle){
		lineStyle = _lineStyle;
	};
};


struct Profile2D{
	TH2D *hist = nullptr;
	std::vector<Double_t> _bin_entries;
	UInt_t nBinsTot;
	Profile2D(std::string _name, std::string _title, UInt_t _nBinsX, Double_t _xMin, Double_t _xMax, UInt_t _nBinsY, Double_t _yMin, Double_t _yMax){
		set(_name, _title, _nBinsX, _xMin, _xMax, _nBinsY, _yMin, _yMax);
	};
	Profile2D(){};
	~Profile2D(){};
	void set(std::string _name, std::string _title, UInt_t _nBinsX, Double_t _xMin, Double_t _xMax, UInt_t _nBinsY, Double_t _yMin, Double_t _yMax){
		hist = new TH2D(_name.c_str(), _title.c_str(), _nBinsX, _xMin, _xMax, _nBinsY, _yMin, _yMax);
		hist->GetXaxis()->CenterTitle();
		hist->GetYaxis()->CenterTitle();
		nBinsTot = (_nBinsX+2) * (_nBinsY+2);
		_bin_entries.reserve(nBinsTot);
		for(UInt_t i = 0; i < nBinsTot; i++){
			_bin_entries.push_back(0.);
		}
	};
	void fill(Double_t _x, Double_t _y, Double_t _Zvalue, Double_t _weight=1.0){
		hist->Fill(_x, _y, _Zvalue * _weight);
		UInt_t whichBin = hist->FindBin(_x, _y);
		_bin_entries[whichBin] += _weight;
	}
	TH2D * getProfile(){
		std::string profile_name = (std::string) hist->GetName() + "_profile";
		std::string profile_title = "Profile\\ " + (std::string) hist->GetTitle();
		TH2D *profile = (TH2D*) hist->Clone(profile_name.c_str());
		profile->Reset();
		profile->GetXaxis()->CenterTitle();
		profile->GetYaxis()->CenterTitle();
		profile->SetTitle(profile_title.c_str());
		Double_t sumEntries = std::accumulate(_bin_entries.begin(), _bin_entries.end(), 0);
		for(UInt_t i = 0; i < nBinsTot; i++){
			Double_t _mean = _bin_entries[i] > 0. ? hist->GetBinContent(i)/_bin_entries[i] : 0. ;
			if(_bin_entries[i] > 0.) profile->SetBinContent(i, _mean);
		}
		return profile;
	}
	operator TH2D* (){
		return hist;
	}
};

template <typename anytype>
struct indexer{
	indexer(){};
	indexer(std::vector<anytype> _list){
		init(_list);
	};
	~indexer(){};

	void init(std::vector<anytype> _list){
		for(UInt_t iIndex = 0; iIndex < _list.size(); iIndex++){
			theIndex[_list[iIndex]]	= iIndex;
		}
	};

	UInt_t size(){
		return theIndex.size();
	};

	std::map<anytype, Int_t> theIndex;

	Int_t operator[](anytype needle){
		if ( theIndex.find(needle) != theIndex.end()) return theIndex[needle];
		else return -999;	
	};
};

class CSVReader{
private:
	std::string fileName;
	std::string delimeter;

public:
	CSVReader(std::string filename, std::string delm = ",") :
	fileName(filename), delimeter(delm){};

	std::vector<std::vector<std::string> > getData(){
		std::ifstream file(fileName);
		std::vector<std::vector<std::string> > dataList;
		std::string line = "";
		while (getline(file, line)){
			std::vector<std::string> vec = split_string(line, delimeter);
			// boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
			// boost::algorithm::split_regex(vec, line, boost::regex(delimeter));
			if(vec.empty()) continue;
			dataList.push_back(vec);
		}
		// Close the File
		file.close();

		return dataList;
	};
};


struct effectiveAreaMap{
	std::map<Float_t, Float_t> effectiveAreas;
	std::string mapFile;

	effectiveAreaMap(){};

	effectiveAreaMap(std::string mapFile, Bool_t _verbose=1, std::string _delimiter="        ", Int_t _firstLine = 2){
		init(mapFile, _verbose, _delimiter, _firstLine);
	};

	void init(std::string _mapFile, Bool_t _verbose=1, std::string _delimiter="        ", Int_t _firstLine = 2){
		if(!file_exists(_mapFile)){
			std::cout<<"\tError! File does not exist! "<<mapFile<<std::endl;
			exit(EXIT_FAILURE);
		}

		mapFile = _mapFile;
		CSVReader readFile(mapFile, _delimiter);
		std::vector<std::vector<std::string>> effAreaData 		= 	readFile.getData();

		for(UInt_t iRow = _firstLine; iRow < effAreaData.size(); iRow++){
			Float_t uBound 										= 	std::stof(effAreaData[iRow][1]);
			Float_t effArea 									= 	std::stof(effAreaData[iRow][2]);
			effectiveAreas[uBound] 								= 	effArea;
		}

		if(_verbose){
			std::cout<<"\tEffective areas loaded from file: "<< mapFile<<std::endl;
			std::cout<<"\t\t\t|eta|\t\tEA:"<<std::endl;
			for(std::map<Float_t, Float_t>::iterator iEl 		= 	effectiveAreas.begin(); iEl != effectiveAreas.end(); iEl++){
				std::cout<<"\t\t\t"<<iEl->first<<"\t\t"<<iEl->second<<std::endl;
			}
		}
	};

	Float_t getEffectiveArea(Float_t absEta){
		if(effectiveAreas.empty()) return 0.;
		std::map<Float_t, Float_t>::iterator iBin = effectiveAreas.lower_bound(absEta);
		if(iBin == effectiveAreas.end()){
			std::cout<<"Error! Eta "<<absEta<<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}
		return iBin->second;
	};

	Float_t getEffectiveAreaAbs(Float_t eta){
		if(effectiveAreas.empty()) return 0.;
		std::map<Float_t, Float_t>::iterator iBin = effectiveAreas.lower_bound(std::abs(eta));
		if(iBin == effectiveAreas.end()){
			std::cout<<"Error! Eta "<<eta<<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}
		return iBin->second;
	};
};

struct isoPtScalingMap{

	std::map<Float_t, std::pair<Float_t,Float_t>> ptScalingCoeffs;
	std::map<Float_t, std::pair<Float_t,Float_t>>::iterator etaBin;
	Bool_t isQuadratic;
	std::string mapFile;

	isoPtScalingMap(){};

	isoPtScalingMap(std::string _mapFile, Bool_t _isQuadratic, Bool_t _verbose=1, std::string _delimiter="        ", Int_t _firstLine = 2){
		init(_mapFile, _isQuadratic, _verbose, _delimiter, _firstLine);
	};

	void init(std::string _mapFile, Bool_t _isQuadratic, Bool_t _verbose=1, std::string _delimiter="        ", Int_t _firstLine = 2){
		if(!file_exists(_mapFile)){
			std::cout<<"\tError! Isolation scaling file does not exist! "<<_mapFile<<std::endl;
			exit(EXIT_FAILURE);
		}
		mapFile = _mapFile;
		isQuadratic = _isQuadratic;

		CSVReader readFile(_mapFile, _delimiter);
		std::vector<std::vector<std::string>> ptScalingData 		= 	readFile.getData();

		for(UInt_t iRow = _firstLine; iRow < ptScalingData.size(); iRow++){
			Float_t uBound 										= 	std::stof(ptScalingData[iRow][1]);
			Float_t c1 											= 	std::stof(ptScalingData[iRow][3]);
			Float_t c2 											= 	_isQuadratic ? std::stof(ptScalingData[iRow][4]) : 0.;
			ptScalingCoeffs[uBound] 								= 	std::make_pair(c1,c2);
		}

		if(_verbose){
			std::cout<<"\tIsolation pT scalingsa loaded from file: "<< _mapFile<<std::endl;
			std::cout<<"\t\t\t|eta|\t\t\tc1:\t\t\tc2"<<std::endl;
			for(std::map<Float_t, std::pair<Float_t,Float_t>>::iterator iEl	= 	ptScalingCoeffs.begin(); iEl != ptScalingCoeffs.end(); iEl++){
				std::cout<<"\t\t\t"<<iEl->first<<"\t\t\t"<<iEl->second.first<<"\t\t"<<iEl->second.second<<std::endl;
			}
		}
	};

	Float_t getPtScaling(Float_t _absEta, Float_t _pT){
		etaBin = ptScalingCoeffs.lower_bound(_absEta);

		if(etaBin == ptScalingCoeffs.end()){
			std::cout<<"Error! Eta "<< _absEta <<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}

		Float_t scalingCorrection = (etaBin->second.first)*_pT + (isQuadratic ? (etaBin->second.second)*_pT*_pT : 0.);
		return scalingCorrection;
	};

	Float_t getPtScalingAbs(Float_t _eta, Float_t _pT){
		etaBin = ptScalingCoeffs.lower_bound(std::abs(_eta));

		if(etaBin == ptScalingCoeffs.end()){
			std::cout<<"Error! Eta "<< _eta <<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}

		Float_t scalingCorrection = (etaBin->second.first)*_pT + (isQuadratic ? (etaBin->second.second)*_pT*_pT : 0.);
		return scalingCorrection;
	};

};


struct coronaCorrections{
	std::vector<std::vector<Float_t>>							percentiles;
	std::vector<std::vector<std::vector<Float_t>>> 				percentileLines;
	std::vector<Float_t> 										etaLimits;
	

	coronaCorrections(std::string mapFile, std::string pDelim = ";", std::string sDelim=",", Bool_t verbose = 1){
		init(mapFile, pDelim, sDelim, verbose);
	};

	coronaCorrections(){};

	void init(std::string mapFile, std::string pDelim = ";", std::string sDelim=",", Bool_t verbose = 1){
		percentileLines.clear();
		etaLimits.clear();
		percentiles.clear();

		CSVReader readFile(mapFile, pDelim);
		std::vector<std::vector<std::string>> effAreaData 		= 	readFile.getData();

		if(verbose) std::cout<<"\tLoading Corona corrections from file "<<mapFile<<std::endl;

		for(UInt_t i = 0; i < effAreaData.size(); i++){

			etaLimits.push_back(strToFloatList(effAreaData[i][0], sDelim)[1]);

			std::vector<std::vector<Float_t>>	etaLines;
			std::vector<Float_t>				iPercentiles;

			for(UInt_t j = 1; j < effAreaData[i].size(); j++){
				std::vector<Float_t> lineDef =	strToFloatList(effAreaData[i][j], sDelim);
				
				if(lineDef.size() != 3) continue;

				iPercentiles.push_back(lineDef[0]);
				etaLines.push_back({lineDef[1], lineDef[2]});
			};

			percentileLines.push_back(etaLines);
			percentiles.push_back(iPercentiles);
		};

		if(verbose) printData();
	};

	void printData(){
		for(UInt_t iEta = 0; iEta< etaLimits.size(); iEta++){
			std::cout<<"eta < "<<etaLimits[iEta]<<std::endl;
			for(UInt_t iP = 0; iP < percentileLines[iEta].size(); iP++){
				std::cout<<"\t"<<percentiles[iEta][iP]*100.<<"%-tile"<<std::endl;
				for(UInt_t iData = 0; iData < percentileLines[iEta][iP].size(); iData++){
					std::cout<<"\t\t"<<percentileLines[iEta][iP][iData]<<"\t";
				}
				std::cout<<std::endl;
			}
		}
	};

	Float_t getCorrection(Float_t _absEta, Float_t _iso, Float_t _rho, Bool_t verbose = 0){
		Int_t etaIndex = -999;
		for(UInt_t iEta = 0; iEta < etaLimits.size(); iEta++){
			if(_absEta < etaLimits[iEta]){
				etaIndex = iEta;
				break;
			};
		};
		
		if(etaIndex < 0) return 0.;

		Bool_t 	foundCorrection 		= 		0;
		Int_t 	corrP 					= 		0;
		
		for(Int_t iPercentile = percentileLines[etaIndex].size() - 1; iPercentile > -1 ; iPercentile--){
			Float_t yLine 				=		percentileLines[etaIndex][iPercentile][0]*_rho + percentileLines[etaIndex][iPercentile][1];
			if(verbose) std::cout<<percentiles[etaIndex][iPercentile]*100<<"%-tile \t yLine = "<<yLine<<std::endl;
			if(_iso < yLine) {
				corrP 			= 		iPercentile;
				break;
			};
		}

		if(verbose) std::cout<<"eta = "<< etaLimits[etaIndex]<<"\t % = "<< percentiles[etaIndex][corrP]*100<<"\t m = "<<percentileLines[etaIndex][corrP][0]<<"\t c = "<<percentileLines[etaIndex][corrP][1]<<std::endl;

		Float_t corrrection 				= 		percentileLines[etaIndex][corrP][0]*_rho;
		
		return corrrection;
	};

};


struct parseOptions {

	std::string optFile;
	std::map<std::string, std::string> optMap;
	Bool_t isParsed = false;

	parseOptions(std::string _optFile, std::string _delimiter=",", Bool_t _verbose = 1): optFile(_optFile){
		parseIt(optFile, _delimiter, _verbose);

	};

	parseOptions(){};

	Bool_t keyExists(std::string _key){
		if(optMap.find(_key) == optMap.end()){
			return 0;
		} else{
			return 1;
		}
	}

	void parseIt(std::string _optFile, std::string _delimiter=",", Bool_t _verbose = 1){
		if(!file_exists(_optFile)){
			std::cout<<"Error! Options file "<<_optFile<<" not found!"<<std::endl;
		}
		optMap.clear();
		if(_verbose) std::cout<<"\t\tOptions parsed from file "<<_optFile<<"... "<<std::endl;
		CSVReader _csvFile(_optFile, _delimiter);
		std::vector<std::vector<std::string>> _data = _csvFile.getData();

		std::string lastOption="";
		for(UInt_t i = 0; i < _data.size(); i++){
			std::string _optName = _data[i][0];
			trim(_optName);
			if(_optName.empty()) continue;
			if(match("#*", _optName)) continue;
			if((_data[i].size()==1)){
				optMap[lastOption] += _data[i][0];
				continue;
			};
			std::string _optVal = _data[i][1];
			optMap[_optName] = _optVal;
			lastOption = _optName;
		}

		isParsed = true;

		if(_verbose){
			for(std::map<std::string, std::string>::iterator optKey = optMap.begin(); optKey != optMap.end(); optKey++){
				std::cout<<"\t\t\t"<<optKey->first<<"\t\t"<<optKey->second<<std::endl;
			}

			std::cout<<"\tParsed!"<<std::endl;
		};
	};

	Float_t getFloat(std::string _opt){
		return std::stof(get(_opt));
	};

	Double_t getDouble(std::string _opt){
		return std::stod(get(_opt));
	};

	Int_t getInt(std::string _opt){
		return std::stoi(get(_opt));
	};

	std::string get(std::string _opt){
		if(!keyExists(_opt)){
			cout<<"Error! Key "<<_opt<<" not found!"<<std::endl;
		}
		return optMap.at(_opt);
	};

	const char* getCSTR(std::string _opt){
		if(!keyExists(_opt)){
			cout<<"Error! Key "<<_opt<<" not found!"<<std::endl;
		}
		return optMap.at(_opt).c_str();
	};

	std::vector<std::string> getList(std::string _opt, std::string _delimiter = ","){
		return split_string(get(_opt), _delimiter, 1);
	};

	std::vector<Float_t> getFloatList(std::string _opt, std::string _delimiter = ","){
		return  strToFloatList(get(_opt), _delimiter);
	};

	std::vector<Double_t> getDoubleList(std::string _opt, std::string _delimiter = ","){
		return  strToDoubleList(get(_opt), _delimiter);
	};

	std::vector<Int_t> getIntList(std::string _opt, std::string _delimiter = ","){
		return  strToIntList(get(_opt), _delimiter);
	};

	std::vector<Int_t> getTColListFromHexList(std::string _opt, std::string _delimiter = ","){

		std::vector<std::string> 	hexList 		=	split_string(get(_opt), _delimiter);
		std::vector<Int_t> 			tColList;

		for(std::string & hexCol : hexList){
			tColList.push_back(hex2rootColor(hexCol));
		};

		return  tColList;
	};

	Int_t getTColFromHex(std::string _opt){
		return hex2rootColor(get(_opt));
	};
};


struct bitBar2D {
	TH2F hist;
	std::vector<const Bool_t*> xVars;
	const Float_t *yVar = nullptr;

	bitBar2D(){};
	bitBar2D(std::vector<std::pair<Bool_t*, std::string>> _xVars, const plot_variable & _yVar){
		init(_xVars, _yVar);
	};

	void init(std::vector<std::pair<Bool_t*, std::string>> _xVars, const plot_variable & _yVar){
		yVar = _yVar.xptr;
		Int_t nxBins = _xVars.size();

		std::string histName = removeNonAlpha(_yVar.xtitle) + "__vs__";
		std::string histTitle = _yVar.xtitle + " vs (";
		for(UInt_t i = 0; i < _xVars.size(); i++){
			histName += removeNonAlpha(_xVars[i].second)+"_";
			histTitle += _xVars[i].second + ", ";
			xVars.push_back(_xVars[i].first);
		}
		histName.pop_back();
		histTitle.pop_back();
		histTitle.pop_back();
		histTitle += ")";

		histTitle += ";;" + _yVar.xtitle;
		histTitle += _yVar.xunit.empty() ? "" : ("[" + _yVar.xunit + "]");
		if(_yVar.xBins == nullptr) hist = TH2F(histName.c_str(), histTitle.c_str(), nxBins, 0., (Float_t)nxBins, _yVar.nbins, _yVar.xmin, _yVar.xmax);
		else hist = TH2F(histName.c_str(), histTitle.c_str(), nxBins, 0., (Float_t)nxBins, _yVar.nbins, _yVar.xBins);

		for(UInt_t i = 0; i < _xVars.size(); i++){
			hist.GetXaxis()->SetBinLabel(i+1, _xVars[i].second.c_str());
		}

		hist.GetYaxis()->CenterTitle();

		std::cout<<"\t\t\tInitialized bitBar2D "<<histTitle<<std::endl;
	}

	void fill(Float_t weight = 1.){
		for(UInt_t i = 0; i < xVars.size(); i++){
			hist.Fill((Float_t)*xVars[i]+0.01, *yVar, (Float_t)*xVars[i] * weight);
		}
	}
};

struct BGset{
	BGset(){};

	BGset(std::string _filenames, std::string _legend, std::string _hexColor, std::string _separator="", std::string _xsecFacator="1"){
		init(_filenames, _legend, _hexColor, _separator, _xsecFacator);
	};

	std::vector<std::string> fileNames;
	std::string legend;
	std::string hexColor;
	Double_t xsecFacator;

	void init(std::string _filenames, std::string _legend, std::string _hexColor, std::string _separator=",", std::string _xsecFacator="1"){
		fileNames = split_string(_filenames, _separator);
		legend = _legend;
		hexColor = _hexColor;
		trim(_xsecFacator);
		xsecFacator=std::stod(_xsecFacator);

		std::cout<<"\tInitialized background set:"<<std::endl<<"\t\tFiles:"<<std::endl;;

		for(std::string fileName : fileNames){
			std::cout<<"\t\t\t"<<fileName<<std::endl;
		}

		std::cout<<"\t\tLegend: "<<legend<<std::endl;
		std::cout<<"\t\tColor: "<<hexColor<<std::endl;
		std::cout<<"\t\txSec k-factor: "<<xsecFacator<<std::endl;
	};

	void assignAttFill(TH1 *_hist, TLegend *_legend = nullptr){
		_hist->SetFillColor(TColor::GetColor(hexColor.c_str()));
		if(_legend){
			TLegendEntry *legEntry = _legend->AddEntry(_hist, legend.c_str(), "F");
			legEntry->SetTextColor(TColor::GetColor(hexColor.c_str()));
		}
	};

	void assignAttLine(TH1 *_hist, TLegend *_legend = nullptr){
		_hist->SetLineColor(TColor::GetColor(hexColor.c_str()));
		if(_legend){
			TLegendEntry *legEntry = _legend->AddEntry(_hist, legend.c_str(), "LFPE");
			legEntry->SetTextColor(TColor::GetColor(hexColor.c_str()));
		}
	};
};


template <typename anytype>
struct smarter_ptr{
	anytype * thePtr = nullptr;
	
	smarter_ptr(){};
	
	smarter_ptr(anytype *_thePtr):thePtr(_thePtr){};

	~smarter_ptr(){
		delete thePtr;
		thePtr = nullptr;
	};

	anytype * get(){
		return thePtr;
	};

	anytype* operator->(){
		return thePtr;
	}

	smarter_ptr & operator = (const anytype * & _assign){
		thePtr = _assign;
	};

	operator anytype*(){
		return thePtr;
	};

	// anytype operator -> anytype(){
	// 	return *thePtr;
	// };
};

class PileupReWeighting {
public:

	PileupReWeighting( ){ };

	PileupReWeighting( std::string mcFile, std::string dataFile, std::string mcHistName, std::string dataHistName){
		init(mcFile, dataFile, mcHistName, dataHistName);
	};

	void init( std::string mcFile, std::string dataFile, std::string mcHistName, std::string dataHistName){

		std::cout<<"\tCreating Pileup Reweighter with "<<std::endl<<
		"\t\t\tdata histogram "<< dataHistName<<" from file "<<dataFile<<std::endl<<
		"\t\t\tmc histogram "<< mcHistName<<" from file "<<mcFile<<std::endl;

		weights_ = (TH1*) getHistFromFile(dataHistName, dataFile);
		MC_distr_ = (TH1*) getHistFromFile(mcHistName, mcFile);

		int NBins = weights_->GetNbinsX();

		weights_->Scale( 1000000.0/ weights_->Integral(0, NBins+1));
		MC_distr_->Scale( 1000000.0/ MC_distr_->Integral(0, NBins+1));

		weights_->SetName("pileupWeights");
		weights_->Divide(MC_distr_);
		MC_distr_->Delete();

		std::cout << "\t\tPileup weights: " << std::endl;

		for(int ibin = 1; ibin<NBins+1; ++ibin){
			std::cout << "\t\t\t" << ibin-1 << "\t" << weights_->GetBinContent(ibin) << std::endl;
		}
		std::cout<<"\tPileup weights calculated!"<<std::endl;
	};

	Double_t weight( Float_t n_int ){
		Int_t bin = weights_->GetXaxis()->FindBin(n_int);
		return weights_->GetBinContent(bin);
	}

	~PileupReWeighting(){
		delete weights_;
	};

protected:
	TH1* weights_	= nullptr;
	TH1* MC_distr_	= nullptr;
};


Double_t deltaR(Double_t _eta1, Double_t _phi1, Double_t _eta2, Double_t _phi2){
	Double_t _deltaEta = _eta1 - _eta2;
	Double_t _deltaPhi = deltaPhi(_phi1, _phi2);
	Double_t _deltaR = std::sqrt(_deltaEta*_deltaEta + _deltaPhi*_deltaPhi);
	return _deltaR;
};


Double_t deltaPhi(Double_t phi1, Double_t phi2){
	Double_t tmpDeltaPhi = std::abs(phi2-phi1);
	Double_t minDeltaPhi = (tmpDeltaPhi > TMath::Pi()) ? (2*TMath::Pi() - tmpDeltaPhi) : tmpDeltaPhi;
	return minDeltaPhi;
};


std::string removeNonAlpha(std::string word){
	word.erase(std::remove_if(word.begin(), word.end(),
		[](char ch){
			return !::iswalnum(ch);
		}), word.end());
	return word;
};

std::string removeNonAlphaSmart(std::string word, std::string texts2eraseCSV){

	word.erase(std::remove_if(word.begin(), word.end(),
		[](char ch){
			Bool_t doRemove =(!std::iswalnum(ch)) && (ch != '.') && (ch != '_');
			return doRemove;
		}), word.end());

	word = findAndReplaceAll(word, ".", "p");

	std::vector<std::string> strings2erase 	=	split_string(texts2eraseCSV);

	for(std::string & iErase : strings2erase){
		word = findAndReplaceAll(word, iErase, "");
	}

	return word;
};


template <class any_number>
std::string removeTrailingZeros(any_number number){
	std::string str = std::to_string (number);
	str.erase( str.find_last_not_of('0') + 1, std::string::npos);
	str.erase(str.find_last_not_of('.') + 1, std::string::npos);
	if(str.length()>0 && !str.substr(str.length()-1,1).compare(".")) str.pop_back();
	return str;
};


Bool_t file_exists(std::string fileName){
	// std::ifstream infile(fileName);
	// return infile.good();
	if(isDirectory(fileName)==1) return 0;
	struct stat buffer;
	return (stat (fileName.c_str(), &buffer) == 0);
};


Int_t mkdir(std::string dir_path){
	std::string command = "mkdir -p " + dir_path;
	const int dir_err = system(command.c_str());

	// if(checkIfDirectory(dir_path)){
	// 	std::cout<<"Directory "<<dir_path<<" already exists"<<std::endl;
	// 	return 1;
	// }

	if (-1 == dir_err){
		printf("Error creating directory!");
	}
	else std::cout <<"Created directory: " << dir_path << std::endl;
	return dir_err;
};


std::map<std::string, Double_t> load_xsecs(std::string filepath){
	std::map<std::string, Double_t> value_map;
	if(!file_exists(filepath)){
		std::cout<<"Cannot load file "<<filepath<<std::endl;
		return value_map;
	}
	CSVReader reader(filepath);
	std::vector<std::vector<std::string>> data_matrix = reader.getData();
	for(auto & row : data_matrix){
		std::string signal_name = row[0];
		Double_t xsec_val = std::stod(row[1]);
		std::string unit_name = row[3];
		value_map[signal_name] = xsec_val * xsec_unit_map.at(unit_name);
	}

	std::cout<<"Xsections [in pb] read from file "<<filepath<<std::endl;
	for(auto & it : value_map){
		std::cout<<"\t "<<it.first<<" \t "<<it.second<<std::endl;
	}
	return value_map;
};


std::vector<std::string> getObjectList(std::string filepath, std::string objtype, std::vector<std::string> exclusion_list){
	std::cout<<"Making object list..."<<std::endl;
	std::vector<std::string> obj_names;
	TFile rootfile(filepath.c_str(), "READ");
	TIter next(rootfile.GetListOfKeys());
	TKey *key;
	while ((key = (TKey*) next())){
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom(objtype.c_str())) continue;
		TObject * g = key->ReadObj();
		std::string name = g->GetName();
		Bool_t to_exclude = (std::find(exclusion_list.begin(), exclusion_list.end(), name) != exclusion_list.end());
		if(to_exclude){
			std::cout << "\t " << "Excluding "<< objtype<<" : "<< name << std::endl;
			continue;
		}
		std::cout << "\t " << "Added "<< objtype<<" : "<< name << std::endl;
		obj_names.push_back(name);
		g->Delete();
	}
	rootfile.Close();
	std::sort( obj_names.begin(), obj_names.end() );
	obj_names.erase( std::unique( obj_names.begin(), obj_names.end() ), obj_names.end());
	return obj_names;
};


Bool_t match(std::string _needle, std::string _haystack){
	char const *needle = _needle.c_str();
	char const *haystack = _haystack.c_str();
	for (; *needle != '\0'; ++needle){
		switch (*needle){
			case '?':{
				if (*haystack == '\0')	return false;
				++haystack;
				break;
			}
			case '*':{
				if (needle[1] == '\0')	return true;
				size_t max = strlen(haystack);
				for (size_t i = 0; i < max; i++)
					if (match(needle + 1, haystack + i)) return true;
				return false;
			}
			default:
			if (*haystack != *needle) return false;
			++haystack;
		}
	}
	return *haystack == '\0';
};


std::string ReadNthLine(std::string filename, int N){
	std::ifstream in(filename.c_str());
	std::string s;
   //for performance
	s.reserve(200);
   //skip N lines
	for(int i = 0; i < N; ++i){
		std::getline(in, s);
	}
	std::getline(in,s);
	return s;
};


UInt_t countLines(std::string filename){
	std::ifstream myfile(filename);
	// new lines will be skipped unless we stop it from happening:
	myfile.unsetf(std::ios_base::skipws);
	// count the newlines with an algorithm specialized for counting:
	UInt_t line_count = std::count(std::istream_iterator<char>(myfile),	std::istream_iterator<char>(), '\n');
	return line_count;
};


// std::vector<std::string> split_string(std::string _string, std::string _delimiter){
// 	std::vector<string> _split_string;
// 	boost::split(_split_string,_string,boost::is_any_of(_delimiter));
// 	return _split_string;
// };


std::vector<std::string> split_string(std::string _string, std::string _delimiter, Bool_t _trim){
	size_t pos = 0;
	std::string token;
	std::vector<std::string> res;
	while ((pos = _string.find(_delimiter)) != std::string::npos){
		token = _string.substr(0, pos);
		if(_trim) trim(token);
		_string.erase(0, pos + _delimiter.length());
		res.push_back(token);
	}
	if(_trim) trim(_string);
	res.push_back(_string);
	return res;
};


std::string get_cell(std::string filename, UInt_t row, UInt_t column, std::string _delimiter){
	return split_string(ReadNthLine(filename, row), _delimiter)[column];
};


std::vector<std::string> splitpath(const std::string& str, const std::set<char> delimiters){
	std::vector<std::string> result;

	char const* pch = str.c_str();
	char const* start = pch;
	for(; *pch; ++pch)
	{
		if (delimiters.find(*pch) != delimiters.end())
		{
			if (start != pch)
			{
				std::string str(start, pch);
				result.push_back(str);
			}
			else
			{
				result.push_back("");
			}
			start = pch + 1;
		}
	}
	result.push_back(start);

	return result;
}


std::string getFileName(std::string _filepath){
	constexpr char sep = '/';
	size_t i = _filepath.rfind(sep, _filepath.length());
	if (i != string::npos){
		return (_filepath.substr(i+1, _filepath.length() - i));
	}
	return(_filepath);
};


std::string getDirPath(std::string _somePath){
	std::string directory;
	const size_t last_slash_idx = _somePath.rfind('/');
	if (std::string::npos != last_slash_idx){
		directory = _somePath.substr(0, last_slash_idx);
	}
	return directory;
};


std::vector<std::string> getNonemptyLines(std::string filepath){
	std::ifstream infile(filepath);
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(infile, line))
	{
		trim(line);
		if(!line.empty()) lines.push_back(line);
	}
	return lines;
};


std::vector<std::string> getNonemptyLinesWithFilenameKeyword(std::string filepath, std::string keyword, std::string exclude){
	std::ifstream infile(filepath);
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(infile, line))
	{
		std::string filename = getFileName(line);
		if(filename.empty()) continue;
		if(!match(keyword.c_str(), filename.c_str())) continue;
		if(!exclude.empty() && match(exclude.c_str(), filename.c_str())) continue;
		lines.push_back(line);
	}
	return lines;
};


std::vector<std::string> getLinesRegex(std::string _filepath, std::string _regexStr){
	regex _regex(_regexStr);
	std::ifstream infile(_filepath);
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(infile, line)){
		if(!std::regex_match(line, _regex)) continue;
		lines.push_back(line);
	}
	return lines;
};


TH1* getHistFromFile(std::string _histName, std::string _filename, Bool_t _verbose, TDirectory *_dir){

	TFile _file(_filename.c_str(), "READ");
	if(!objectExists(&_file, _histName)){
		_file.Close();
		std::cout<<"Error! Histogram "<<_histName<<" not found in file "<<_filename<<std::endl;
		return nullptr;
	}

	TH1* _hist = (TH1*) _file.Get(_histName.c_str());
	if(!_dir)_hist->SetDirectory(0);
	else _hist->SetDirectory(_dir);
	_file.Close();
	if(_verbose)std::cout<<"\t\tLoaded histogram "<<_histName<<" (N="<<_hist->GetEntries() <<") from file "<<_filename<<std::endl;
	if(!_hist) std::cout<<"Error! Histogram "<< _histName<<" not found in file "<<_filename<<std::endl;
	return _hist;
};


TObject *getObjectFromFile(std::string _objectName, std::string _filename){
	TFile _file(_filename.c_str(), "READ");
	TObject* _object = (TObject*) _file.Get(_objectName.c_str());
	// _object->SetDirectory(0);
	_file.Close();
	if(!_object){
		std::cout<<"Error reading "<<_objectName<< " from file "<<_filename<<std::endl;
		return nullptr;
	}
	std::cout<<"\t\tLoaded TObject "<<_objectName<<" from file "<<_filename<<std::endl;
	return _object;
};


TH1* rebinHist(TH1* _hist, Double_t _statUnc){
	std::vector<Double_t> _goodBins = getGoodBins(_hist, _statUnc);
	if(_goodBins.size()<2) return _hist;
	std::string _newname = "rebinned_" + (std::string)_hist->GetName();
	TH1* _rebinnedHist = (TH1*) _hist->Rebin(_goodBins.size()-1, _newname.c_str(), _goodBins.data());
	std::cout<<"\t\t\tRebinned "<<_hist->GetName()<<". New bins:";
	for(Double_t iBin : _goodBins){
		std::cout<<"\t"<<iBin;
	}
	std::cout<<std::endl;	
	// _rebinnedHist->Scale(1.,"width");
	return _rebinnedHist;
};


TH1* rebinHist(TH1* _hist, std::vector<Double_t> _newBins){
	std::string _newname = "rebinned_" + (std::string)_hist->GetName();
	TH1 *_hist_rebinned = (TH1*) _hist->Rebin(_newBins.size()-1, _newname.c_str(), _newBins.data());
	std::cout<<"\t\t\tRebinned "<<_hist->GetName()<<". New bins:";
	
	for(Double_t iBin : _newBins){
		std::cout<<"\t"<<iBin;
	}

	_hist->Delete();
	return _hist_rebinned;
};


TH1* rebinNHist(TH1* _hist, Int_t _N){
	std::string _newname = "rebinned_" + std::to_string(_N) + "_" + (std::string)_hist->GetName();
	TH1 *_hist_rebinned = (TH1*) _hist->Rebin(_N, _newname.c_str());
	_hist->Delete();
	return _hist_rebinned;
};

Double_t sumNextNbins(TH1* _hist, Int_t _n, Int_t _curr){
	Double_t _NbinSum = 0.;
	Int_t _countEnd = (_hist->GetNbinsX() <= (_curr + _n)) ? _hist->GetNbinsX() : _curr + _n;
	for(Int_t i = _curr; i < _countEnd+1; i++){
		_NbinSum += _hist->GetBinContent(i);
	}
	return _NbinSum;
};


std::vector<Double_t> getGoodBins(TH1* _hist, Double_t _statUnc, Double_t _reScale, Int_t _nbinPar){

	TH1 *reScaled = (TH1*)_hist->Clone("reScaled");
	if(_reScale > 0.){
		reScaled->Scale(_reScale / reScaled->Integral(0, reScaled->GetNbinsX()));
		_hist = reScaled;
	}
	UInt_t _nBins = _hist->GetXaxis()->GetNbins();
	if (_nBins == 0){
		std::cout<<"\t Hist "<<_hist->GetName()<<" has no bins!"<<std::endl;
		return {};
	}
	_hist->Sumw2();
	std::vector<Double_t> good_bins;
	good_bins.push_back(_hist->GetXaxis()->GetBinLowEdge(1));
	Double_t _runningBinSum = 0.;
	UInt_t _first50pcBins = std::ceil(0.5 *_nBins);
	for(UInt_t i = 1; i < _nBins+1; i++){
		Double_t _bincontent = _hist->GetBinContent(i);
		Double_t _next10bincontent = sumNextNbins(_hist,_nbinPar,i);
		if(_bincontent != _bincontent) return {};
		Double_t _binUpedge = _hist->GetXaxis()->GetBinUpEdge(i);
		Double_t _binLowedge = _hist->GetXaxis()->GetBinLowEdge(i);
		//(i <=_first50pcBins && _bincontent == 0.)
		// if((_next10bincontent == 0. && i < _nBins) || (_bincontent == 0.)){
		if((_next10bincontent == 0. && i < _nBins)){
			if(good_bins.size()>0){
				Double_t _lastElement = good_bins.back();
				if(_binLowedge -_lastElement > 0.000000001) good_bins.push_back(_binLowedge);
			}
			good_bins.push_back(_binUpedge);
			_runningBinSum = 0.;
			continue;
		}
		_runningBinSum += _bincontent;
		if(1/std::sqrt(_runningBinSum) < _statUnc){
			good_bins.push_back(_binUpedge);
			_runningBinSum =0;
		}
		if(i < _nBins) continue;
		Double_t _lastElement = good_bins.back();
		if(_binUpedge -_lastElement < 0.00000001) continue;
		good_bins.push_back(_binUpedge);
		_runningBinSum =0;
		if(1/std::sqrt(_runningBinSum) > _statUnc) continue;
		if(good_bins.size()==1) continue;
		//erase second last element
		good_bins.erase(good_bins.begin() + good_bins.size()-2);
	};
	reScaled->Delete();
	return good_bins;
};



Double_t getSumW(std::string _cutflowfile){
	std::string sumWline = ReadNthLine(_cutflowfile, 3);
	std::string sumWstring = (split_string(sumWline, ",")).at(0);
	Double_t _sumW = std::stod(sumWstring);
	std::cout<<"\tFrom "<<_cutflowfile<<" sumW="<<_sumW<<std::endl;
	return _sumW;
};


// trim from start (in place)
void ltrim(std::string &s){
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch){
		return !std::isspace(ch);
	}));
}


// trim from end (in place)
void rtrim(std::string &s){
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch){
		return !std::isspace(ch);
	}).base(), s.end());
}


// trim from both ends (in place)
void trim(std::string &s){
	ltrim(s);
	rtrim(s);
}


// trim from start (copying)
std::string ltrim_copy(std::string s){
	ltrim(s);
	return s;
}


// trim from end (copying)
std::string rtrim_copy(std::string s){
	rtrim(s);
	return s;
}


// trim from both ends (copying)
std::string trim_copy(std::string s){
	trim(s);
	return s;
};


Double_t ams(Double_t _s, Double_t _b){
	return std::sqrt(2*(1+_s+_b)*std::log(1+_s/_b)-2*_s);
};


template <typename T>
std::string to_string_with_precision(const T a_value, const int n){
	std::ostringstream out;
	out << std::fixed << std::setprecision(n) << a_value;
	return out.str();
};


std::string getUnit(TH1* _hist){
	std::string ystring = _hist->GetYaxis()->GetTitle();
	return getUnit(ystring);
};


std::string getUnit(std::string ystring){
	size_t lastSlash = ystring.find_last_of("/");
	if(lastSlash == std::string::npos) return "";
	std::string _unit = ystring.substr(lastSlash+1);
	trim(_unit);
	std::string number = first_numberstring(_unit);
	std::string unit_text = _unit;
	boost::replace_all(unit_text, number, "");
	trim(unit_text);
	if(unit_text.empty()) return std::to_string(1);
	else return unit_text;
};


std::string eraseUnit(std::string ystring){
	size_t lastSlash = ystring.find_last_of("/");
	if(lastSlash == std::string::npos) return ystring;
	std::string _unit = ystring.substr(lastSlash);
	boost::replace_all(ystring, _unit, "");
	return ystring;
}


std::string first_numberstring(std::string const & str){
	std::size_t const n = str.find_first_of("0123456789.");
	if (n != std::string::npos)
	{
		std::size_t const m = str.find_first_not_of("0123456789.", n);
		return str.substr(n, m != std::string::npos ? m-n : m);
	}
	return std::string();
};


void closeTChain(TChain * & _chain){
	if(!_chain){
		std::cout<<"Error! TChain is null!"<<std::endl;
	}
	TFile *file = _chain->GetCurrentFile();
	_chain->SetDirectory(0);
	if(file) delete file;
	_chain = nullptr;
};

void closeTTree(TTree * & _chain){
	if(!_chain){
		std::cout<<"Error! TChain is null!"<<std::endl;
	}
	TFile *file = _chain->GetCurrentFile();
	_chain->SetDirectory(0);
	if(file) file->Close();
	_chain = nullptr;
};


std::vector<Double_t> getXbins(TH1* hist){
	UInt_t nbins = hist->GetNbinsX();
	if (nbins == 0){
		std::cout << "Error : No Bins!" << std::endl;
		return {};
	};
	std::vector<Double_t> bins;
	TAxis* axis = hist->GetXaxis();
	bins.push_back(axis->GetBinLowEdge(1));
	for (UInt_t i = 1; i < nbins + 1; i++){
		bins.push_back(axis->GetBinUpEdge(i));
	}
	return bins;
};


void setFrameColor(TAxis* _axis, std::string _color){
	Int_t color_val = TColor::GetColor(_color.c_str());
	_axis->SetAxisColor(color_val);
	_axis->SetLabelColor(color_val);
	_axis->SetTitleColor(color_val);
};


void setFrameColor(TH1* _hist, std::string _color){
	setFrameColor(_hist->GetXaxis(), _color);
	setFrameColor(_hist->GetYaxis(), _color);
};


void setFrameColor(THStack* _stack, std::string _color){
	setFrameColor(_stack->GetXaxis(), _color);
	setFrameColor(_stack->GetYaxis(), _color);
};


TTree *loadTreeFromFile(std::string _treeName, std::string _filename){
	TFile *_file = new TFile(_filename.c_str(), "READ");
	if(!_file){
		std::cout<<"Error! Cannot read file "<<_filename<<std::endl;
		return nullptr;
	}

	TTree *_tree = (TTree*) _file->Get(_treeName.c_str());
	if(!_tree){
		std::cout<<"Error! Cannot read TTree "<<_file<<" from successfully loaded TFile "<<_filename <<std::endl;
		return nullptr;
	}

	std::cout<<"Successfully loaded TTree "<<_treeName<<" from TFile "<<_filename <<std::endl;
	std::cout<<"Remember to call TTree*->GetDirectory()->Close() once the tree is not needed any more!" <<std::endl;

	return _tree;
};


Char_t isDirectory(std::string filePath, Bool_t _verbose){
	DIR* dir = opendir(filePath.c_str());
	if (dir){
		closedir(dir);
		return 1;
	} else if (ENOENT == errno){
		return 0;
	} else{
		if(_verbose) std::cout<<"Error checking path : "<<filePath<<std::endl;
		return -1;
	}
}


void addPointToGraph(TGraph & _graph, Double_t _x, Double_t _y){
	_graph.SetPoint(_graph.GetN(), _x, _y);
};


void graphStats(TGraphAsymmErrors* graph, Double_t &mean, Double_t &stdev){
	UInt_t Npoints = graph->GetN();
	std::vector<Double_t> xVals(graph->GetY(), graph->GetY() + Npoints), xErrsL(graph->GetEYlow(), graph->GetEYlow() + Npoints), xErrsH(graph->GetEYhigh(), graph->GetEYhigh() + Npoints);
	Double_t W, sumW = 0, sumW2 = 0, sumXW = 0, sumDev2W = 0;
	for (UInt_t i = 0; i < Npoints; i++){
		if (xErrsL[i] * xErrsH[i] > 0){
			W = 1 / (xErrsL[i] * xErrsH[i]);
			sumW += W;
			sumXW += xVals[i] * W;
		}
	}
	if (sumW > 0){
		mean = sumXW / sumW;
		for (UInt_t i = 0; i < Npoints; i++){
			if (xErrsL[i] * xErrsH[i] > 0){
				W = 1 / (xErrsL[i] * xErrsH[i]);
				sumW2 += W * W;
				sumDev2W += (xVals[i] - mean) * (xVals[i] - mean) * W;
			}
		}
		stdev = TMath::Sqrt(sumDev2W / (sumW - sumW2 / sumW));
	} else {
		mean = -999;
		stdev = -999;
	}
};


TH1D* graph2hist(TGraphAsymmErrors* graph, UInt_t ndivs, Double_t ylow, Double_t yhigh){
	std::string name = graph->GetName();
	name.append("_hist_" + to_string_with_precision(rand() % 1000, 0));
	UInt_t Npoints = graph->GetN();

	std::vector<Double_t> xVals(graph->GetY(), graph->GetY() + Npoints), xErrsL(graph->GetEYlow(), graph->GetEYlow() + Npoints), xErrsH(graph->GetEYhigh(), graph->GetEYhigh() + Npoints);

	TH1D* hist = new TH1D("", "", ndivs, ylow, yhigh);
	hist->SetFillColor(graph->GetLineColor());
	hist->SetLineColor(graph->GetLineColor());
	hist->SetMarkerColor(graph->GetLineColor());

	Double_t w = 1;
	for (UInt_t i = 0; i < Npoints; i++){
		if (!(xErrsL[i] <= 0  || xErrsH[i] <= 0))w = 1 / (xErrsL[i] * xErrsH[i]);
		else continue;
		hist->Fill(xVals[i], w);
	}
	Npoints = hist->GetEffectiveEntries();
	hist->Scale(Npoints / hist->Integral());
	return hist;
};


Int_t writeToFile(TObject *_object, std::string _filePath, std::string _mode, Bool_t _verbose){

	mkdir(getDirPath(_filePath));

	TFile _outfile(_filePath.c_str(), _mode.c_str());

	if(_outfile.IsZombie()){
		std::cout<<"\t\tFailed to open TFile "<<_filePath<<" in the mode "<<_mode<<std::endl;
		return -1;
	}

	_outfile.cd();
	_object->Write(_object->GetName());
	_outfile.Close();
	if(_verbose) std::cout<<"\t\tWritten "<< _object->GetName()<<" to file "<<_filePath<<std::endl;
	return 0;
};


std::string findAndReplaceAll(std::string data, std::string toSearch, std::string replaceStr){
	// Get the first occurrence
	size_t pos = data.find(toSearch);

	// Repeat till end is reached
	while( pos != std::string::npos)
	{
		// Replace this occurrence of Sub String
		data.replace(pos, toSearch.size(), replaceStr);
		// Get the next occurrence from the current position
		pos =data.find(toSearch, pos + replaceStr.size());
	}

	return data;
};


TChain *openTChain(std::string _chainListFile, std::string _treeName, Bool_t _verbose){
	if(!file_exists(_chainListFile)){
		std::cout<<"Error! Chain list file does not exist:"<<_chainListFile<<std::endl;
		exit(EXIT_FAILURE);
	}
	if(_verbose) std::cout<<"\tMaking TChain with root files listed in "<<_chainListFile<<std::endl;
	std::vector<std::string> _ntuples = getNonemptyLines(_chainListFile);
	if(_ntuples.size() == 0){
		std::cout<<"Error! File list "<<_chainListFile<<" is invalid!"<<std::endl;
		return nullptr;
	}

	return openTChain(_ntuples, _treeName, _verbose);
};


TChain *openTChain(std::vector<std::string> _chainList, std::string _treeName, Bool_t _verbose){
	if(_treeName.empty()){
		TFile _testF(_chainList[0].c_str(), "READ");
		if (!_testF.GetListOfKeys()){
			std::cout<<"Error! Cannot open TChain. No key found in file "<<_chainList[0]<<std::endl;
			exit(EXIT_FAILURE);
		}
		TIter lNext(_testF.GetListOfKeys()) ;
		TObject* lObj ;
		Char_t treeCounter = 0;
		while((lObj = (TObject*)lNext())){
			if(treeCounter > 1){
				std::cout<<"Error! No TTree name provided & there are >1 trees in the file "<<_chainList[0]<<std::endl;
				return nullptr;
			}
			if(lObj->InheritsFrom(TTree::Class())){
				_treeName = lObj->GetName();
				treeCounter++;
			}
		}
		_testF.Close();
	}

	if(_treeName.empty()){
		std::cout<<"Error! No TTree name found!"<<std::endl;
		return nullptr;
	}

	TChain *_bChain = new TChain(_treeName.c_str());
	if(_verbose) std::cout<<"\tMaking TChain from trees "<<_treeName<<std::endl;
	for(auto & _ntuple : _chainList){
		if(_ntuple.find("root:") != std::string::npos){
		} else if(!file_exists(_ntuple)){
			std::cout<<"Error! File does not exist "<<_ntuple<<std::endl;
			closeTChain(_bChain);
			return nullptr;
		};
		_bChain->Add(_ntuple.c_str());
		if(_verbose) std::cout << "\tAdded file "<< _ntuple <<std::endl;
	}
	if(_verbose) std::cout<<"\tBuilt TChain!"<<std::endl;

	return _bChain;
};


TChain *openTChainWithFilesInDir(std::string _dirPath, std::string _treeName){
	if(isDirectory(_dirPath)!=1){
		std::cout<<"Error! Directory "<<_dirPath<<" doesn't exist"<<std::endl;
		return nullptr;
	}

	std::cout<<"Making TChain from trees "<<_treeName<<" in directory "<<_dirPath<<std::endl;
	std::vector<std::string> _fileList = listFilesInDir(_dirPath, ".*\.root", 1);
	// TChain *_bChain = openTChain(_fileList, _treeName);
	return openTChain(_fileList, _treeName);
};


std::vector<std::pair<std::string, std::string>> getBranchList(std::string _treePath, std::string _treeName, Bool_t _verbose){
	TChain *bChain = openTChain(std::vector<std::string> {_treePath}, _treeName);
	if(!bChain){
		std::cout<<"Error! Could no create TChain!"<<std::endl;
		return {};
	}

	TObjArray *bArray = (TObjArray*)bChain->GetListOfBranches()->Clone();
	bArray->SetOwner(kFALSE);
	// bArray->Sort();

	std::vector<std::pair<std::string, std::string>> bList;

	std::cout<<"List of branches in TTree "<<_treeName<<" in file "<<_treePath<<":"<<std::endl;
	TIter iBr(bArray);
	TObject* bObj = nullptr;
	while ((bObj = (TObject*)iBr())){
		std::string bName = bObj->GetName();
		std::string bType = bChain->GetLeaf(bName.c_str())->GetTypeName();
		bList.emplace_back(bType, bName);
		if(_verbose) std::cout<<"\t"<<bType<<",\t"<<bName<<std::endl;
	}

	closeTChain(bChain);

	return bList;
};


std::vector<std::string> listFilesInDir(std::string _dirPath, std::string _regexStr, Bool_t _verb){

	if(!isDirectory(_dirPath)){
		std::cout<<"Error! Directory ("<<_dirPath<<") does not exist!"<<std::endl;
		return {};
	}

	std::vector<std::string> _fileList;

	DIR           *dirp;
	struct dirent *directory;

	regex regex_;

	if(!_regexStr.empty()){
		regex regextmp_(_regexStr);
		regex_ = regextmp_;
	}

	if(_verb) std::cout <<"Files in directory "<< _dirPath<<" matching  \""<<_regexStr<< "\" :" << std::endl;
	dirp = opendir(_dirPath.c_str());
	if (dirp){
		while ((directory = readdir(dirp)) != NULL){
			if(!_regexStr.empty()){
				if(!std::regex_match(directory->d_name, regex_)) continue;
			}
			if(_verb) std::cout <<"\t\t"<< getFileName(directory->d_name) << std::endl;
			std::string _filePath = _dirPath + "/" + directory->d_name;
			_fileList.push_back(_filePath);
		}
		closedir(dirp);
	}

	return _fileList;
};


Bool_t matchRegex(std::string _str,std::string _regexStr){
	regex regex_(_regexStr);
	return std::regex_match(_str, regex_);
};


std::string getTreeNameInFile(std::string _filePath){
	std::string _treeName="";
	TFile _testF(_filePath.c_str(), "READ");
	if (!_testF.GetListOfKeys()){
		std::cout<<"Error! Cannot open TChain. No key found in file "<<_filePath<<std::endl;
	}
	TIter lNext(_testF.GetListOfKeys()) ;
	TObject* lObj ;
	Char_t treeCounter = 0;
	while((lObj = (TObject*)lNext())){
		if(treeCounter > 1){
			std::cout<<"Error! There are >1 trees in the file "<<_filePath<<std::endl;
			return "";
		}
		if(lObj->InheritsFrom(TTree::Class())){
			_treeName = lObj->GetName();
			treeCounter++;
		}
	}
	_testF.Close();
	return _treeName;
};


// 	mergeBins
// 	Adds histograms for MC samples of the same process produced in different phase space bins
// 	_inFiles : list of root files containing histograms
// 	_histName : the histogram to be merged
// 	_xsecMap : CSV file listing samples and cross sections
// 	root file name is used to look up cross section
TH1F *mergeBins(std::vector<std::string> _inFiles, std::string _histName, std::string _sumWeightsHistname, std::string _xsecMap, Int_t _nameCol, Int_t _xSecCol, std::string _path){

	TH1F *_mergedHist = (TH1F*) getHistFromFile(_histName, _path + _inFiles[0]);
	_mergedHist->Reset("ICESM");
	_mergedHist->SetDirectory(0);
	_mergedHist->Sumw2();
	std::string _newName = _mergedHist->GetName();
	_newName += "_mergedBins";
	_mergedHist->SetName(_newName.c_str());

	for(const auto & _file : _inFiles){
		std::string _sampleName = getFileName(_file);
		_sampleName = findAndReplaceAll(_sampleName, ".root", "");

		TH1F *_sumWHist = (TH1F*) getHistFromFile(_sumWeightsHistname, _path + _file);
		Double_t _sumW = _sumWHist->GetBinContent(1);
		_sumWHist->Delete();

		Double_t _xSection = std::stod(vLookup(_sampleName, _xsecMap, _nameCol, _xSecCol));

		TH1F *_binHist = (TH1F*) getHistFromFile(_histName, _path + _file);
		_binHist->Scale(_xSection/_sumW);

		Double_t _integral = _binHist->Integral();

		_mergedHist->Add(_binHist);
		_binHist->Delete();

		std::cout<<"\t\t"<<_file<<":\t xSection = "<<_xSection<<"\tSumW = "<<_sumW<<"\tIntegral = "<< _integral<<std::endl;

	}

	std::cout<<"\t\t\tMerged all bins! Integral = "<<_mergedHist->Integral()<<std::endl;

	return _mergedHist;
};


// 	mergeBins
// 	Adds histograms for MC samples of the same process produced in different phase space bins
// 	_fileList : list of root files containing histograms
// 	_histName : the histogram to be merged
// 	_xsecMap : CSV file listing samples and cross sections
// 	root file name is used to look up cross section
TH1F *mergeBins(std::string _fileList, std::string _histName, std::string _sumWeightsHistname, std::string _xsecMap, Int_t _nameCol, Int_t _xSecCol, std::string _path){
	std::cout<<"Merging TH1F "<<_histName<<" using file list "<< _fileList<<std::endl;
	std::vector<std::string> _inFiles;
	if(!isROOTfile(_fileList))_inFiles = getNonemptyLines(_fileList);
	else _inFiles = {_fileList};

	std::sort(_inFiles.begin(), _inFiles.end());

	return mergeBins(_inFiles, _histName, _sumWeightsHistname, _xsecMap, _nameCol, _xSecCol, _path);
};


std::string vLookup(std::string _lookupKey, std::string _inFile, Int_t _lookupCol, Int_t _valCol, Bool_t _regex){
	CSVReader reader(_inFile);
	std::vector<std::vector<std::string>> data_matrix = reader.getData();

	Int_t _matchedRow = -999;
	for(UInt_t i = 0; i < data_matrix.size(); i++){
		std::string _searchCell = data_matrix[i][_lookupCol];

		if(_regex){
			if(matchRegex(_searchCell, _lookupKey)){
				_matchedRow = i;
			}
		} else{
			if(_searchCell == _lookupKey){
				_matchedRow = i;
			}
		}
		if(_matchedRow>-1) break;
	}
	if(_matchedRow < 0){
		std::cout<<"vLookup failed to "<<(_regex ? "match " : "find ")<<"key "<<_lookupKey<<" in file "<<_inFile<<std::endl;
		return "";
	}

	return data_matrix[_matchedRow][_valCol];
};


Bool_t isROOTfile(std::string _filepath){
	if(!file_exists(_filepath)){
		return 0;
	}
	if(_filepath.substr(_filepath.find_last_of(".") + 1) == "root"){
		if(TFile(_filepath.c_str()).GetListOfKeys() == nullptr){
			return 0;
		}
		return 1;
	}
	return 0;
};


std::vector<Float_t> getXlimits(std::vector<TH1*> _hists, Float_t _binThreshold){
	std::vector<Float_t> _limits;
	for(auto _hist: _hists){
		Double_t first = _hist->GetXaxis()->GetBinLowEdge(_hist->FindFirstBinAbove(_binThreshold));
		Double_t last = _hist->GetXaxis()->GetBinUpEdge(_hist->FindLastBinAbove(_binThreshold));
		_limits.push_back(first);
		_limits.push_back(last);
	}

	std::vector<Float_t> _max_min;
	_max_min.push_back(*std::max_element(_limits.begin(), _limits.end()));
	_max_min.push_back(*std::min_element(_limits.begin(), _limits.end()));
	return _max_min;
};


void clearHeap(){
	TList* obList = gDirectory->GetList();
	for(auto obj : *obList){
		obj->Delete();
	}
};


Double_t weightedYmean(TH1 *_hist){
	Double_t _weightedSumY = 0.;
	Double_t _weightSum = 0.;
	for(Int_t i = 0; i < _hist->GetNbinsX(); i++){
		if(!(_hist->GetBinContent(i) > 0.)) continue;
		// if(!(_hist->GetBinError(i) > 0.)) continue;
		// Double_t _binWeight = 1./(_hist->GetBinError(i) * _hist->GetBinError(i));
		Double_t _binWeight = 1.;
		_weightSum += _binWeight;
		_weightedSumY += _binWeight * _hist->GetBinContent(i);
	}
	return (_weightedSumY/_weightSum);
};


Double_t weightedYspread(TH1 *_hist){
	Double_t _weightedSumDeltaY2 = 0.;
	Double_t _weightSum = 0.;
	Double_t _weightedYmean = weightedYmean(_hist);

	for(Int_t i = 0; i < _hist->GetNbinsX(); i++){
		if(!(_hist->GetBinContent(i) > 0.)) continue;
		// if(!(_hist->GetBinError(i) > 0.)) continue;
		// Double_t _binWeight = 1./(_hist->GetBinError(i) * _hist->GetBinError(i));
		Double_t _binWeight = 1.;
		_weightSum += _binWeight;
		_weightedSumDeltaY2 += _binWeight * (_hist->GetBinContent(i) - _weightedYmean) * (_hist->GetBinContent(i) - _weightedYmean);
	}
	return (_weightedSumDeltaY2/_weightSum);
};


Bool_t branchExists(std::string _branchName, TTree *_tree){
	TBranch* br = (TBranch*) _tree->GetListOfBranches()->FindObject(_branchName.c_str());
	if(br)	return 1;
	else return 0;
};


Float_t getMean(std::vector<Float_t> _set){
	if(_set.empty()) return -9999.;
	Double_t _sum = 0.;
	for(Float_t _num : _set){
		_sum += _num;
	}
	return _sum/((Double_t)_set.size());
};


template <class ObjType>
ObjType copyObjectDeleteSrc(ObjType *_original){
	ObjType _copy(*_original);
	_original->Delete();
	return _copy;
};


template<typename T1, typename T2>
Int_t findSecondaryIndex(T1 searchIndex, std::vector<T2> container){
	T2 secondaryIndex = -999;
	for(UInt_t i = 0; i < container.size(); i++){
		if( container[i] ==  (T2) searchIndex){
			secondaryIndex = i;
			break;
		}
	}
	return secondaryIndex;
};


Short_t findSecondaryIndex(Short_t searchIndex, std::vector<Short_t> container){
	Short_t secondaryIndex = -999;
	for(UInt_t i = 0; i < container.size(); i++){
		if( container[i] ==  searchIndex){
			secondaryIndex = i;
			break;
		}
	}
	return secondaryIndex;
};


std::string getCurrentTime(){
	std::chrono::time_point<std::chrono::system_clock> _now = std::chrono::system_clock::now();
	std::time_t _now_ = std::chrono::system_clock::to_time_t(_now);
	std::string c_time = std::ctime(&_now_);
	c_time.erase(std::remove(c_time.begin(), c_time.end(), '\n'), c_time.end());
	return c_time;
};


template <typename anytype>
void eraseElement(std::vector<anytype> & _vec, UInt_t _Index2Erase){
	_vec.erase(_vec.begin() + _Index2Erase);
};


Double_t getCategoryBoundary(TH1 *_signal, TH1*_background){

	std::cout<<"Signal integral:"<<_signal->Integral()<<std::endl<<"Background integral:"<<_background->Integral()<<std::endl;

	Int_t optimalBoundaryBin = -999;
	Int_t optimalJointSignificance = -999.;

	for(Int_t i = 1; i < _signal ->GetNbinsX() + 1; i++){

		Double_t sIntegralBelow = _signal->Integral(0, i);
		Double_t bIntegralBelow = _background->Integral(0, i);

		Double_t sIntegralAbove = _signal->Integral(i+1, _signal ->GetNbinsX()+1);
		Double_t bIntegralAbove = _background->Integral(i+1, _signal ->GetNbinsX()+1);

		if(sIntegralAbove < 0.1*sIntegralBelow) continue;

		Double_t iJointSignificance = std::sqrt(std::pow(sIntegralBelow/std::sqrt(sIntegralBelow + bIntegralBelow), 2) + std::pow(sIntegralAbove/std::sqrt(sIntegralAbove + bIntegralAbove), 2));

		if(optimalJointSignificance < iJointSignificance){
			optimalJointSignificance = iJointSignificance;
			optimalBoundaryBin = i;
		}
	}

	Double_t sIntegralBelow = _signal->Integral(0, optimalBoundaryBin);
	Double_t bIntegralBelow = _background->Integral(0, optimalBoundaryBin);
	Double_t sIntegralAbove = _signal->Integral(optimalBoundaryBin+1, _signal->GetNbinsX()+1);
	Double_t bIntegralAbove = _background->Integral(optimalBoundaryBin+1, _signal->GetNbinsX()+1);
	Double_t optimalSignificanceBelow = sIntegralBelow/std::sqrt(sIntegralBelow + bIntegralBelow);
	Double_t optimalSignificanceAbove = sIntegralAbove/std::sqrt(sIntegralAbove + bIntegralAbove);
	Double_t uncategorizedSignificance = _signal->Integral(0,_signal->GetNbinsX()+1)/std::sqrt(_signal->Integral(0,_signal->GetNbinsX()+1)+_background->Integral(0,_signal->GetNbinsX()+1));

	std::cout<<"**************************************************************************************************************************"<<std::endl<<
	"Optmal boundary for categorizing in the variable "<<_signal->GetXaxis()->GetTitle()<<std::endl<<
	"\t\tBin # = "<<optimalBoundaryBin<<std::endl<<
	"\t\tUp edge = "<<_signal->GetXaxis()->GetBinUpEdge(optimalBoundaryBin)<<std::endl<<
	"\t\tUncategorized Significance = "<<uncategorizedSignificance<<std::endl<<
	"\t\tOptimized Significance Above = "<<optimalSignificanceAbove<<std::endl<<
	"\t\tOptimized Significance Below = "<<optimalSignificanceBelow<<std::endl<<
	"\t\tJoint Significance = "<<std::sqrt(optimalSignificanceBelow*optimalSignificanceBelow + optimalSignificanceAbove*optimalSignificanceAbove)<<
	"**************************************************************************************************************************"<<std::endl;

	return _signal->GetXaxis()->GetBinUpEdge(optimalBoundaryBin);
};

std::vector<Double_t> vecString2vecDouble(std::vector<std::string> _numStrings){
	std::vector<Double_t> _convertedDoubles;	
	for(UInt_t i = 0; i < _numStrings.size(); i++){
		std::string _iString =  _numStrings[i];
		trim(_iString);
		_convertedDoubles.push_back(std::stod(_iString));
	}
	return _convertedDoubles;
};


Bool_t stringIsNumber(std::string _isThisANumber){
	Bool_t is_a_number = false;
	
	try{
		boost::lexical_cast<double>(_isThisANumber);
		is_a_number = true;
	} catch(boost::bad_lexical_cast &){}

	return is_a_number;
};


TDirectory *mkTFileDir(TFile *_file, std::string _dir){
	if(_file->GetDirectory(_dir.c_str()) == nullptr) _file->mkdir(_dir.c_str());
	return _file->GetDirectory(_dir.c_str());
};


int caselessSubstringPos( std::string str1, std::string str2){
	std::string::const_iterator it = std::search( str1.begin(), str1.end(), 
		str2.begin(), str2.end(), [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2);});
	if ( it != str1.end() ) return it - str1.begin();
	else return -1; // not found
};


Bool_t findStringIC(const std::string & strHaystack, const std::string & strNeedle){
	auto it = std::search(
		strHaystack.begin(), strHaystack.end(),
		strNeedle.begin(),   strNeedle.end(),
		[](char ch1, char ch2){ return std::toupper(ch1) == std::toupper(ch2); }
		);
	return (it != strHaystack.end() );
};


Bool_t stringIsEmpty(std::string str){
	if(str.find_first_not_of(' ') != std::string::npos)	return 0;
	return 1;
};


//line=style,hexcolor-alpha,width
template<typename T>
UChar_t setLineAtts(T* graph, std::string _atts, std::string _delimiter){
	Int_t returnVal = 0;

	indexer<std::string> sequence({"style", "hexcolor-alpha", "width"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	std::cout<<"Line attributes of " + (std::string) graph->GetName() + " set to:";

	for(Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++){
		if((iAtt == sequence["style"]) && stringIsNumber(atts[iAtt])){
			graph->SetLineStyle(std::stoi(atts[iAtt]));
			std::cout<<"\tstyle = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["hexcolor-alpha"]) && !stringIsEmpty(atts[iAtt])){
			Int_t doTransparency = caselessSubstringPos( atts[iAtt], "-");
			if(doTransparency > -1 && stringIsNumber(atts[iAtt].substr(doTransparency+1))){
				Float_t alpha = std::stof(atts[iAtt].substr(doTransparency+1));
				graph->SetLineColorAlpha(TColor::GetColor(atts[iAtt].c_str()), alpha);
				std::cout<<" color = " + atts[iAtt].substr(0,doTransparency-1) + " (alpha=" + atts[iAtt].substr(doTransparency+1) + ") ";
				setBit(returnVal,iAtt,1);
			} else{
				graph->SetLineColor(TColor::GetColor(atts[iAtt].c_str()));
				std::cout<<" color = " + atts[iAtt] + " ";
				setBit(returnVal,iAtt,1);
			}
		} else if((iAtt == sequence["width"]) && stringIsNumber(atts[iAtt])){
			graph->SetLineWidth(std::stof(atts[iAtt]));
			std::cout<<" width = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		}
	}

	std::cout<<std::endl;

	return returnVal;
};


//marker=style,hexcolor-alpha,size
template<typename T>
UChar_t setMarkerAtts(T* graph, std::string _atts, std::string _delimiter){
	Int_t returnVal = 0;

	indexer<std::string> sequence({"style", "hexcolor-alpha", "size"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	std::cout<<"Marker attributes of " + (std::string) graph->GetName() + " set to:";

	for(Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++){
		if((iAtt == sequence["style"]) && stringIsNumber(atts[iAtt])){
			graph->SetMarkerStyle(std::stoi(atts[iAtt]));
			std::cout<<"\tstyle = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["hexcolor-alpha"]) && !stringIsEmpty(atts[iAtt])){
			Int_t doTransparency = caselessSubstringPos( atts[iAtt], "-");
			if(doTransparency > -1 && stringIsNumber(atts[iAtt].substr(doTransparency+1))){
				Float_t alpha = std::stof(atts[iAtt].substr(doTransparency+1));
				graph->SetMarkerColorAlpha(TColor::GetColor(atts[iAtt].c_str()), alpha);
				std::cout<<" color = " + atts[iAtt].substr(0,doTransparency-1) + " (alpha=" + atts[iAtt].substr(doTransparency+1) + ") ";
				setBit(returnVal,iAtt,1);
			} else{
				graph->SetMarkerColor(TColor::GetColor(atts[iAtt].c_str()));
				std::cout<<" color = " + atts[iAtt] + " ";
				setBit(returnVal,iAtt,1);
			}
		} else if((iAtt == sequence["width"]) && stringIsNumber(atts[iAtt])){
			graph->SetMarkerSize(std::stof(atts[iAtt]));
			std::cout<<" size = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		}
	}

	std::cout<<std::endl;

	return returnVal;
};


//fill=style,hexcolor-alpha
template<typename T>
UChar_t setFillAtts(T* graph, std::string _atts, std::string _delimiter){
	Int_t returnVal = 0;

	indexer<std::string> sequence({"style", "hexcolor-alpha"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	std::cout<<"Fill attributes of " + (std::string) graph->GetName() + " set to:";


	for(Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++){
		if((iAtt == sequence["style"]) && stringIsNumber(atts[iAtt])){
			graph->SetFillStyle(std::stoi(atts[iAtt]));
			std::cout<<"\tstyle = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["hexcolor-alpha"]) && !stringIsEmpty(atts[iAtt])){
			Int_t doTransparency = caselessSubstringPos( atts[iAtt], "-");
			if(doTransparency > -1 && stringIsNumber(atts[iAtt].substr(doTransparency+1))){
				Float_t alpha = std::stof(atts[iAtt].substr(doTransparency+1));
				graph->SetFillColorAlpha(TColor::GetColor(atts[iAtt].c_str()), alpha);
				std::cout<<" color = " + atts[iAtt].substr(0,doTransparency-1) + " (alpha=" + atts[iAtt].substr(doTransparency+1) + ") ";
				setBit(returnVal,iAtt,1);
			} else{
				graph->SetFillColor(TColor::GetColor(atts[iAtt].c_str()));
				std::cout<<" color = " + atts[iAtt] + " ";
				setBit(returnVal,iAtt,1);
			}
		}
	}

	std::cout<<std::endl;

	return returnVal;
};


//axis=range(min-max),limits(min-max),title_size,label_size,ndivs(n1-n2-n3-opt),title_offset,label_offset,more_log_labels,maxDigits
UChar_t setAxisAtts(TAxis* _axis, std::string _atts, std::string _delimiter){
	Int_t returnVal = 0;
	
	indexer<std::string> sequence({"range(min-max)", "limits(min-max)", "title_size", "label_size", "ndivs", "title_offset", "label_offset", "more_log_labels", "maxDigits"});
	
	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	
	std::cout<<"TAxis attributes of " + (std::string) _axis->GetName();
	if(_axis->GetParent()) std::cout<< " (" + (std::string) _axis->GetParent()->GetName() + ")";
	std::cout<< " set to:";
	
	for(Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++){
		if((iAtt == sequence["range(min-max)"]) && !stringIsEmpty(atts[iAtt])){
			std::vector<Float_t> rangeAtts = strToFloatList(atts[iAtt], "-");
			if(rangeAtts.size() == 2){
				_axis->SetRangeUser(rangeAtts[0], rangeAtts[1]);
				std::cout<<" range = " + atts[iAtt] + " ";
			}
		} else if((iAtt == sequence["limits(min-max)"]) && !stringIsEmpty(atts[iAtt])){
			std::vector<Float_t> limitsAtts = strToFloatList(atts[iAtt], "-");
			if(limitsAtts.size() == 2) {
				_axis->SetLimits(limitsAtts[0], limitsAtts[1]);
				std::cout<<" limits = " + atts[iAtt] + " ";
			}
		} else if((iAtt == sequence["title_size"]) && stringIsNumber(atts[iAtt])){
			_axis->SetTitleSize(std::stof(atts[iAtt]));
			std::cout<<" title size = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["label_size"]) && stringIsNumber(atts[iAtt])){
			_axis->SetLabelSize(std::stof(atts[iAtt]));
			std::cout<<" label size = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["ndivs"]) && !stringIsEmpty(atts[iAtt])){
			atts[iAtt] = findAndReplaceAll(atts[iAtt], " ", "");
			std::vector<std::string> nDivAtt = split_string(atts[iAtt], "-");
			_axis->SetNdivisions(std::stoi(nDivAtt[0]), std::stoi(nDivAtt[1]), std::stoi(nDivAtt[2]), (nDivAtt.size()>3) ? nDivAtt[3].c_str() : "");
			std::cout<<" nDivisions = " + nDivAtt[0] + " " + nDivAtt[1] + " " + nDivAtt[2] + " " + ((nDivAtt.size()>3) ? nDivAtt[3] : "") + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["title_offset"]) && stringIsNumber(atts[iAtt])){
			_axis->SetTitleOffset(std::stof(atts[iAtt]));
			std::cout<<" title offset = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["label_offset"]) && stringIsNumber(atts[iAtt])){
			_axis->SetLabelOffset(std::stof(atts[iAtt]));
			std::cout<<" label offset = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["more_log_labels"]) && stringIsNumber(atts[iAtt])){
			_axis->SetMoreLogLabels(std::stoi(atts[iAtt]));
			std::cout<<" more log labels = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["maxDigits"]) && stringIsNumber(atts[iAtt])){
			((TGaxis*)_axis)->SetMaxDigits(std::stoi(atts[iAtt]));
			std::cout<<" max digits = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		}
		
	}

	std::cout<<std::endl;

	return returnVal;
};


// pad=coordinates(x1-y1-x2-y2),margins(L-R-B-T),fillStyle-color-alpha,bordersize,grid(x-y),log(x-y)
template<typename T>
UChar_t setPadAtts(T* _pad, std::string _atts, std::string _delimiter){
	Int_t returnVal = 0;

	indexer<std::string> sequence({"coordinates(x1-y1-x2-y2)", "margins(L-R-B-T)", "fillStyle-color-alpha", "bordersize", "grid(x-y)", "log(x-y)"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	std::cout<<"Pad attributes of " + (std::string) _pad->GetName() + " set to:";

	for(Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++){
		if((iAtt == sequence["coordinates(x1-y1-x2-y2)"]) && !stringIsEmpty(atts[iAtt])){
			atts[iAtt] = findAndReplaceAll(atts[iAtt], " ", "");
			std::vector<std::string> coords = split_string(atts[iAtt], "-");
			if(coords.size() == 4 && stringIsNumber(coords[0]) && stringIsNumber(coords[1]) && stringIsNumber(coords[2]) && stringIsNumber(coords[3])){
				_pad->SetPad(std::stof(coords[0]), std::stof(coords[1]), std::stof(coords[2]), std::stof(coords[3]));
				std::cout<<" (x1,y1,x2,y2) = (" + coords[0] + "," + coords[1] + "," + coords[2] + "," + coords[2] + ")";
				setBit(returnVal,iAtt,1);
			}
		} else if((iAtt == sequence["margins(L-R-B-T)"]) && !stringIsEmpty(atts[iAtt])){
			atts[iAtt] = findAndReplaceAll(atts[iAtt], " ", "");
			std::vector<std::string> margins = split_string(atts[iAtt], "-");
			if(stringIsNumber(margins[0]) && stringIsNumber(margins[1]) && stringIsNumber(margins[2]) && stringIsNumber(margins[3])){
				_pad->SetMargin(std::stof(margins[0]), std::stof(margins[1]), std::stof(margins[2]), std::stof(margins[3]));
				std::cout<<" margins(L,R,B,T) = (" + margins[0] + "," + margins[1] + "," + margins[2] + "," + margins[3] + ")";
				setBit(returnVal,iAtt,1);
			}
		} else if((iAtt == sequence["fillStyle-color-alpha"]) && !stringIsEmpty(atts[iAtt])){
			setFillAtts(_pad, atts[iAtt], "-");
		} else if((iAtt == sequence["bordersize"]) && stringIsNumber(atts[iAtt])){
			_pad->SetBorderSize(std::stof(atts[iAtt]));
			std::cout<<" border size = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if((iAtt == sequence["grid(x-y)"]) && !stringIsEmpty(atts[iAtt])){
			atts[iAtt] = findAndReplaceAll(atts[iAtt], " ", "");

			std::vector<std::string> gridAtts = split_string(atts[iAtt], "-");

			if(gridAtts.size() > 0 && stringIsNumber(gridAtts[0])){
				_pad->SetGridx(std::stoi(gridAtts[0]));
				std::cout<<" grid X = " + gridAtts[0] + " ";
				setBit(returnVal,iAtt,1);
			}

			if(gridAtts.size() > 1 && stringIsNumber(gridAtts[1])){
				_pad->SetGridy(std::stoi(gridAtts[1]));
				std::cout<<" grid Y = " + gridAtts[1] + " ";
				setBit(returnVal,iAtt,1);
			}
		} else if((iAtt == sequence["log(x-y)"]) && !stringIsEmpty(atts[iAtt])){
			std::vector<std::string> logXY = split_string(atts[iAtt], "-");
			if(logXY.size()>0 && stringIsNumber(logXY[0])) _pad->SetLogx(std::stoi(logXY[0]));
			if(logXY.size()>1 && stringIsNumber(logXY[1])) _pad->SetLogy(std::stoi(logXY[1]));
			std::cout<<" log scale (X,Y) = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		}
	};

	std::cout<<std::endl;

	return returnVal;
};


// frame=fill=style,color-alpha;line=style,color-alpha,width
template<typename T>
void setFrameAtts(T* _hist, std::string _style, std::string _del1, std::string _del2){
	_style = findAndReplaceAll(_style, " ", "");
	std::vector<std::string> styleOptions = split_string(_style, _del1);
	for(const std::string styleOption : styleOptions){
		cout<<styleOption<<endl;
		std::vector<std::string> varVal = split_string(styleOption, "=");
		if(stringToUpper(varVal[0]) == stringToUpper("line")) setLineAtts(_hist, varVal[1], _del2);
		else if(stringToUpper(varVal[0]) == stringToUpper("fill")) setFillAtts(_hist, varVal[1], _del2);
	};
};


//histstyle==line=style,hexcolor-alpha,width;marker=style,hexcolor-alpha,size;fill=style,hexcolor-alpha
template<typename T>
void setHistStyle(T* _hist, std::string _style, std::string _del1, std::string _del2){
	_style = findAndReplaceAll(_style, " ", "");
	std::vector<std::string> styleOptions = split_string(_style, _del1);
	for(const std::string styleOption : styleOptions){
		std::vector<std::string> varVal = split_string(styleOption, "=");
		if(stringToUpper(varVal[0]) == stringToUpper("marker")) setMarkerAtts(_hist, varVal[1], _del2);
		else if(stringToUpper(varVal[0]) == stringToUpper("line")) setLineAtts(_hist, varVal[1], _del2);
		else if(stringToUpper(varVal[0]) == stringToUpper("fill")) setFillAtts(_hist, varVal[1], _del2);
	};
};


// legendstyle==bounds=border_size,x1-y1-x2-y2;fill=style,hexcolor-alpha;entries=text_size,nCols
void setLegendStyle(TLegend *_legend, std::string _style, std::string _del1=";", std::string _del2=","){
	std::vector<std::string> legendStyleOptions = split_string(_style, _del1);

	std::cout<<"Legend attributes of "<<_legend->GetName()<< " set to ";
	for(const std::string styleOption : legendStyleOptions){
		std::vector<std::string> subOptions = split_string(styleOption, "=");

		if(stringToUpper(subOptions[0]) == stringToUpper("bounds")){
			std::vector<std::string> boundsOpts = split_string(subOptions[1], _del2);
			if(stringIsNumber(boundsOpts[0])){
				_legend->SetBorderSize(std::stof(boundsOpts[0]));
				std::cout<<"border size = "+boundsOpts[0]<<" ";
			}
			std::vector<Float_t> coordinates = (boundsOpts.size()>0) ? strToFloatList(boundsOpts[1],"-") : (std::vector<Float_t>) {};
			if(coordinates.size() == 4){
				_legend->SetX1(coordinates[0]);
				_legend->SetY1(coordinates[1]);
				_legend->SetX2(coordinates[2]);
				_legend->SetY2(coordinates[3]);
				std::cout<<"coordinates = "+boundsOpts[1]<<" ";
			}
		} else if(stringToUpper(subOptions[0]) == stringToUpper("fill")){
			setFillAtts(_legend, subOptions[1], _del2);
		} else if(stringToUpper(subOptions[0]) == stringToUpper("entries")){
			std::vector<std::string> entryOpts = split_string(subOptions[1], _del2);
			if(entryOpts.size()>0 && stringIsNumber(entryOpts[0])) {
				_legend->SetTextSize(std::stof(entryOpts[0]));
				std::cout<<"text size = "+entryOpts[0]<<" ";
			}
			if(entryOpts.size()>1 && stringIsNumber(entryOpts[1])) {
				_legend->SetNColumns(std::stoi(entryOpts[1]));
				std::cout<<"N columns = "+entryOpts[1]<<" ";
			}
		}
	};
	std::cout<<std::endl;
};


struct DinkyHister{
	
	DinkyHister(){};
	
	~DinkyHister(){
		for(TH1 * hist : histograms){
			hist->Delete();
		};
	};

	DinkyHister(std::string _canvasSize="1600x1200", std::string _padAtts="", std::string _legendStyle = "bounds=0.,0.8-0.7-0.97-0.95;fill=4000,;entries=0.1,1",
		std::string _titles=""){
		std::vector<std::string> canvasSize = split_string(_canvasSize,"x");
		canvas.SetCanvasSize(1500, 1500);
		canvas.SetWindowSize(std::stoi(canvasSize[0]), std::stoi(canvasSize[1]));
		setPadAtts(&canvas, _padAtts);
		setLegendStyle(&legend, _legendStyle);
		histStack.SetTitle(_titles.c_str());
	};

	void applyStyleToAll(std::string _style, std::string _del1=";", std::string _del2=","){
		for(TH1* hist : histograms){
			setHistStyle(hist, _style, _del1, _del2);
		}
	};
	//_legend=legend--option
	void add(TH1* _hist, std::string _style, std::string _drawOption, std::string _legend){
		
		std::string newName = (std::string) _hist->GetName() + "_" + std::to_string(reinterpret_cast<ULong_t>(this)) + "_" + std::to_string(reinterpret_cast<ULong_t>(_hist));
		
		TH1* hist = (TH1*) _hist->Clone(newName.c_str());
		setHistStyle(hist, _style);
		histStack.Add(hist, _drawOption.c_str());
		histograms.push_back(hist);
		
		if(!stringIsEmpty(_legend)){
			std::vector<std::string> legendOpt = split_string(_legend, "--");
			TLegendEntry* legendEntry = nullptr;
			if(legendOpt.size()>1) legendEntry = legend.AddEntry(hist, legendOpt[0].c_str(), legendOpt[1].c_str());
			else legendEntry = legend.AddEntry(_hist, legendOpt[0].c_str());
			legendEntry->SetTextColor(hist->GetLineColor());
		}
	};

	void draw(std::string _option){
		canvas.Draw();
		canvas.cd();
		histStack.Draw(_option.c_str());
		legend.Draw();
	};

	void write(std::string _path, std::string _fileName, std::string _formats = "png,pdf"){
		_formats = findAndReplaceAll(_formats, ".", "");
		std::vector<std::string> extensions = split_string(_formats, ",");
		mkdir(_path);
		for(const std::string & iFormat : extensions){
			if(!stringIsEmpty(iFormat)) canvas.SaveAs((_path + "/" + _fileName + "." + iFormat).c_str());
		}

	};

	void disown(){
		histograms.clear();
	};

	void update(){
		gPad->Update();
		gPad->RedrawAxis();
		gPad->Modified();
		
		canvas.Update();
		canvas.RedrawAxis();
		canvas.Modified();
	};

	std::vector<TH1*> histograms;
	TCanvas canvas;
	TLegend legend;
	THStack histStack;
	std::string titles;
};


struct correlationMatix{
	UInt_t 								nVar;
	std::vector<std::vector<Double_t>>	data2p;
	std::vector<Double_t>				data1p;
	Double_t 							sumW;

	correlationMatix(){};

	correlationMatix(Int_t _nVar){
		nVar 			=	_nVar;
		data2p			=	std::vector<std::vector<Double_t>>(nVar, std::vector<Double_t>(nVar, 0.));
		data1p			=	std::vector<Double_t>(nVar, 0.);
		sumW 			= 	0.;
	};

	void addPoint(std::vector<Double_t> _point, Double_t _weight){
		sumW += _weight;

		for(UInt_t i = 0; i < nVar; i++){

			data1p[i] += _point[i]*_weight;

			for(UInt_t j = 0; j < nVar; j++){
				data2p[i][j] += (_point[i]*_point[j]*_weight);
			}
		}
		
	}

	TH2F* getCorrelationHistAbs(){
		TH2F* 						corrHist 		=		new TH2F("corrHist", "", nVar, 0, nVar, nVar, 0, nVar);
		std::vector<Double_t>		means;

		for(UInt_t i = 0; i < data1p.size(); i++){
			means.push_back(data1p[i]/sumW);
		}		

		for(UInt_t i = 0; i < data2p.size(); i++){
			Double_t 	denI 	= std::sqrt(data2p[i][i] - sumW * means[i]* means[i]);

			for(UInt_t j = 0; j < data2p.size(); j++){

				Double_t  Rij	=	data2p[i][j] - sumW * means[i]*means[j];
				Rij 			/=	denI;
				Rij 			/=	std::sqrt(data2p[j][j] - sumW * means[j]* means[j]);

				std::cout<<"r["<<i<<", "<<j<<"]\t=\t"<<Rij<<std::endl;
				corrHist->SetBinContent(1+i, 1+j, std::abs(Rij));
			}
		}

		return corrHist;
	}

	TH2F* getCorrelationHist(){
		TH2F* 						corrHist 		=		new TH2F("corrHist", "", nVar, 0, nVar, nVar, 0, nVar);
		std::vector<Double_t>		means;

		for(UInt_t i = 0; i < data1p.size(); i++){
			means.push_back(data1p[i]/sumW);
		}		

		for(UInt_t i = 0; i < data2p.size(); i++){
			Double_t 	denI 	= std::sqrt(data2p[i][i] - sumW * means[i]* means[i]);
			for(UInt_t j = 0; j < data2p.size(); j++){
				Double_t  Rij	=	data2p[i][j] - sumW * means[i]*means[j];
				Rij 			/=	denI;
				Rij 			/=	std::sqrt(data2p[j][j] - sumW * means[j]* means[j]);
				std::cout<<"r["<<i<<", "<<j<<"]\t=\t"<<Rij<<std::endl;
				corrHist->SetBinContent(1+i, 1+j, Rij);
			}
		}

		return corrHist;
	}

};


std::string stringToUpper(std::string _str){
	boost::to_upper(_str);
	return _str;
};


std::string stringToLower(std::string _str){
	boost::to_lower(_str);
	return _str;
};


template<typename T1, typename T2>
Int_t findIndex(const std::vector<T1> & _haystack, T2 _needle){
	Int_t Index = -999;
	for(UInt_t i = 0; i < _haystack.size(); i++){
		if( _haystack[i] == (T1) _needle){
			Index = i;
			break;
		}
	}
	return Index;
};


Bool_t objectExists(TFile *_File, std::string _objectName){
	TObject *_object = _File->Get(_objectName.c_str());
	if(!_object) return 0;
	_object->Delete();
	return 1;
};


Bool_t objectExists(std::string _FilePath, std::string _objectName){
	if(!file_exists(_FilePath)) return 0;
	TFile _file(_FilePath.c_str(), "READ");
	TObject *_object = _file.Get(_objectName.c_str());
	if(!_object) return 0;
	_object->Delete();
	_file.Close();
	return 1;
};


std::vector<Float_t> strToFloatList(std::string _listString, std::string _delimiter){
	std::vector<std::string> numberStrings = split_string(_listString, _delimiter);
	std::vector<Float_t> convertedNumbers;
	for(const std::string & numberString : numberStrings){
		if(stringIsEmpty(numberString)) convertedNumbers.push_back(std::numeric_limits<Float_t>::max());
		else if(stringIsNumber(numberString)) convertedNumbers.push_back(std::stof(numberString));
	};
	return convertedNumbers;
};

std::vector<Double_t> strToDoubleList(std::string _listString, std::string _delimiter){
	std::vector<std::string> numberStrings = split_string(_listString, _delimiter);
	std::vector<Double_t> convertedNumbers;
	for(const std::string & numberString : numberStrings){
		if(stringIsEmpty(numberString)) convertedNumbers.push_back(std::numeric_limits<Double_t>::max());
		else if(stringIsNumber(numberString)) convertedNumbers.push_back(std::stod(numberString));
	};
	return convertedNumbers;
};


Int_t hex2rootColor(std::string _hexCode){
	trim(_hexCode);
	if(!match("#*",_hexCode)) _hexCode = "#" + _hexCode;
	return TColor::GetColor(_hexCode.c_str());
};


std::string sysexec(std::string cmd) {
	std::array<char, 128> buffer;
	std::string result;
	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
	if (!pipe) {
		throw std::runtime_error("popen() failed!");
	}
	while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
		result += buffer.data();
	}
	return result;
};


std::vector<Double_t> getNpercentMinInterval(TH1 *_hist, Float_t _N){
	if (_hist->IsZombie()) {
		cout << "Error: Null histogram!" << endl;
		exit (EXIT_FAILURE);
	}

	UInt_t Nbins = _hist->GetNbinsX();
	UInt_t firstBin = 1, lastBin = Nbins;
	Float_t P = _N / 100.0;
	TH1D* cum = (TH1D*) _hist->GetCumulative(1, "cum");
	TH1D* rcum = (TH1D*) _hist->GetCumulative(0, "rcum");
	TAxis *xaxis = _hist->GetXaxis();
	Double_t integral = cum->GetBinContent(Nbins);
	Double_t tempInterval = (xaxis->GetXmax()) - (xaxis->GetXmin());
	Double_t partialIntegral = 0;
	Double_t leftEdge, rightEdge;

	for (UInt_t i = 1; i < Nbins + 1; i++) {
		if (P * integral > rcum->GetBinContent(i)) break;
		leftEdge = xaxis->GetBinLowEdge(i);
		for (UInt_t j = i; j < Nbins + 1; j++) {
			partialIntegral = cum->GetBinContent(j) - cum->GetBinContent(i - 1);
			rightEdge = xaxis->GetBinUpEdge(j);
			if ((partialIntegral / integral) >= P && (tempInterval > (rightEdge - leftEdge))) {
				lastBin = j;
				firstBin = i;
				tempInterval = rightEdge - leftEdge;
				break;
			}
		}
	}
	cum->Delete();
	rcum->Delete();

	std::vector<Double_t> minRange = {xaxis->GetBinLowEdge(firstBin), xaxis->GetBinUpEdge(lastBin)};
	
	return minRange;
};



struct logStream{

	logStream(){};

	logStream(std::string _logFilePath){
		init(_logFilePath);
	};

	~logStream(){
		logFile.close();
	};

	void init(std::string _logFilePath){
		logFile.open(_logFilePath);
		std::cout<< "Writing logs to file: "<<_logFilePath<<std::endl;		
	};

	ofstream logFile;

	template<typename T>
	logStream & operator<<(const T& mValue){
		std::cout << mValue;
		logFile << mValue;
		return *this;
	};
};


// std::vector<Double_t> getNpercentLinearInterval(TH1 *_hist, Float_t _N){
// 	return {};
// };


std::vector<std::string> prefixVecString(std::vector<std::string> _vecStr, std::string _prefix, Int_t _insPos){
	std::vector<std::string> prefixedVecStr = _vecStr;
	if(_insPos < 0) {
		std::cout<<"Error! negative insert position passed to prefixVecString()"<<std::endl;
		return {};
	}

	UInt_t tmpInsPos = _insPos;
	if(_insPos < 0){
		std::for_each(prefixedVecStr.begin(), prefixedVecStr.end(), [&](std::string & iString){ 
			iString.append(_prefix);
		});
	} else{
		std::for_each(prefixedVecStr.begin(), prefixedVecStr.end(), [&](std::string & iString){
			if(_insPos > (Int_t) iString.length())  tmpInsPos = iString.length();
			iString.insert(tmpInsPos, _prefix);
		});
	}
	return prefixedVecStr;
};



void normalizeHist(TH1* _hist, Double_t _norm){
	_hist->Scale(_norm/_hist->Integral());
};


void normalizeHist(TH1& _hist, Double_t _norm){
	_hist.Scale(_norm/_hist.Integral());
};


std::vector<Int_t> strToIntList(std::string _listString, std::string _delimiter){
	std::vector<std::string> numberStrings = split_string(_listString, _delimiter);
	std::vector<Int_t> convertedNumbers;
	for(const std::string & numberString : numberStrings){
		if(stringIsEmpty(numberString)) convertedNumbers.push_back(std::numeric_limits<Int_t>::max());
		else if(stringIsNumber(numberString)) convertedNumbers.push_back(std::stoi(numberString));
	};
	return convertedNumbers;
};


Double_t spearmanR(TH2 *_hist, Int_t _Nbins){
	
	TH1D*_xProj = (TH1D*)_hist->ProjectionX();
	normalizeHist(_xProj, 1.);
	TH1D*_xCum = (TH1D*)_xProj->GetCumulative();
	delete _xProj;

	TH1D*_yProj = (TH1D*)_hist->ProjectionY();
	normalizeHist(_yProj, 1.);
	TH1D*_yCum = (TH1D*)_yProj->GetCumulative();
	delete _yProj;

	std::string xyRanksName = (std::string)_hist->GetName() + "_xyRanks";
	TH2D xyRanks(xyRanksName.c_str(), "", _Nbins, 0., 1.000001, _Nbins,  0., 1.000001);

	Int_t _xNbins = _xCum->GetNbinsX() + 1;
	Int_t _yNbins = _yCum->GetNbinsX() + 1;

	for(Int_t iX = 1; iX < _xNbins; iX++){
		Double_t iXrank = _xCum->GetBinContent(iX);
		for(Int_t iY = 1; iY < _yNbins; iY++){
			Double_t iYrank = _yCum->GetBinContent(iY);
			Double_t iXYweight = _hist->GetBinContent(iX, iY);
			xyRanks.Fill(iXrank, iYrank, iXYweight);
		}
	}

	delete _xCum;
	delete _yCum;

	Double_t spearmanRval = xyRanks.GetCorrelationFactor();

	return spearmanRval;
};


Float_t relInvMass(Float_t _pT1, Float_t _eta1, Float_t _phi1, Float_t _pT2, Float_t _eta2, Float_t _phi2){
	Float_t invMassVal = std::sqrt(2.*_pT1*_pT2*(std::cosh(_eta1-_eta2) - std::cos(_phi1-_phi2)));
	return invMassVal;
};


template<typename T>
void remove_intersection(std::vector<T>& a, std::vector<T>& b){
	std::unordered_multiset<T> st;
	st.insert(a.begin(), a.end());
	st.insert(b.begin(), b.end());
	auto predicate = [&st](const T& k){ return st.count(k) > 1; };
	a.erase(std::remove_if(a.begin(), a.end(), predicate), a.end());
	b.erase(std::remove_if(b.begin(), b.end(), predicate), b.end());
};


TH1* getSqrtHist(TH1* _hist){
	std::string sqrtName = _hist->GetName();
	sqrtName += "_sqrt";
	TH1* _sqrtHist = (TH1*) _hist->Clone(sqrtName.c_str());
	_sqrtHist->Reset();

	Int_t nBins = _hist->GetNbinsX();

	for(Int_t iBin = 1; iBin < nBins+1; iBin++){
		_sqrtHist->SetBinContent(iBin, std::sqrt(_hist->GetBinContent(iBin)));
	}

	return _sqrtHist;
};


void removeNegativeBins(TH1* _hist){
	for(Int_t iBin = 1; iBin < _hist->GetNbinsX() +1; iBin++){
		if(_hist->GetBinContent(iBin) < 0.) _hist->SetBinContent(iBin, 0.);
	}
};


bool isInteger(const std::string & s) {
	if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;
	char * p;
	strtol(s.c_str(), &p, 10);
	return (*p == 0);
};


Float_t foldedPhi(Float_t _phi){
	if(std::abs(_phi) <= TMath::PiOver2()) return _phi;
	if(_phi < - TMath::PiOver2()) return -(_phi + TMath::Pi());
	if(_phi > TMath::PiOver2()) return -(_phi - TMath::Pi());
	return _phi;
};


Char_t RGB2ColPlt(std::string _CSVFile, Int_t _firstCol=0){
	CSVReader csvFile(_CSVFile, ",");
	std::vector<std::vector<std::string> > csvDat = csvFile.getData();

	std::vector<Int_t> palette;

	for(UInt_t iRow = 0; iRow < csvDat.size(); iRow++){
		if(!stringIsNumber(csvDat[iRow][_firstCol+0]) || !stringIsNumber(csvDat[iRow][_firstCol+1]) || !stringIsNumber(csvDat[iRow][_firstCol+2])) continue;

		// std::cout<<__LINE__<<std::endl;
		if(isInteger(csvDat[iRow][_firstCol+0]) && isInteger(csvDat[iRow][_firstCol+1]) && isInteger(csvDat[iRow][_firstCol+2])) palette.push_back(TColor::GetColor(std::stoi(csvDat[iRow][_firstCol+0]), std::stoi(csvDat[iRow][_firstCol+1]), std::stoi(csvDat[iRow][_firstCol+2])));
		else palette.push_back(TColor::GetColor(std::stof(csvDat[iRow][_firstCol+0]), std::stof(csvDat[iRow][_firstCol+1]), std::stof(csvDat[iRow][_firstCol+2])));
	}

	gStyle->SetPalette(palette.size(), palette.data());
	
	std::cout<<"Loaded color palette from file "<<_CSVFile<<"\t\tN = "<<palette.size()<<std::endl;

	return 0;
};


TH2* getAbsHist(TH2* _hist){
	TH2* _clone = (TH2*) _hist->Clone("abs_clone");
	_clone->Reset();

	for(Int_t iBinX = 1; iBinX <= _hist->GetNbinsX(); iBinX++){
		for(Int_t iBinY = 1; iBinY <= _hist->GetNbinsY(); iBinY++){
			_clone->SetBinContent(iBinX, iBinY, std::abs(_hist->GetBinContent(iBinX, iBinY)));
		}		
	}
	return _clone;
};

#endif
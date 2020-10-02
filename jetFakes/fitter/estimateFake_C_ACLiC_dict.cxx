// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIdatadIcmszfs1dIuserdIwaduddIaNTGCmetdIaNTGC_analysisdIjetFakesdIfitterdIestimateFake_C_ACLiC_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/jetFakes/fitter/./estimateFake.C"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_estimateFake_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./estimateFake.C",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf7/include",
"/data/cmszfs1/user/wadud/aNTGCmet/CMSSW_10_2_23/src",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_23/src",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/coral/CORAL_2_3_21-gnimlf11/include/LCG",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/mctester/1.25.0a-gnimlf8/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf7/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dd4hep/v01-08x-gnimlf4/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtDesigner",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/herwigpp/7.2.0-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tauolapp/1.1.5-gnimlf6/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/charybdis/1.003-gnimlf4/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sherpa/2.2.8-gnimlf4/include/SHERPA-MC",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtOpenGL",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtGui",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/thepeg/2.2.0-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia8/230-gnimlf6/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/herwig/6.521-gnimlf4/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/rivet/2.5.4-gnimlf8/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/Qt3Support",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lwtnn/2.4-gnimlf4/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/geant4/10.04-gnimlf3/include/Geant4",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/classlib/3.1.3-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf4/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2-gnimlf2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tkonlinesw/4.2.0-1_gcc7-gnimlf8/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/starlight/r193-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/Qt",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtCore",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/qt/4.8.7-omkpbe2/include/QtXml",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/mcdb/1.0.3-omkpbe2/interface",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libungif/4.1.4-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libtiff/4.0.3-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libpng/1.6.16-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/frontier_client/2.9.1-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pcre/8.37-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xrootd/4.8.3-gnimlf/include/xrootd/private",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/vdt/0.4.0-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/valgrind/3.13.0-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/utm/utm_0.7.1-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/toprex/4.23-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tbb/2018_U1-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tauola/27.121.5-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sigcpp/2.6.2-omkpbe2/include/sigc++-2.0",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/sqlite/3.22.0-omkpbe/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/protobuf/3.5.2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pacparser/1.3.5-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/oracle/12.1.0.2.0/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/meschach/1.2.pCMS1-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libuuid/2.22.2-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libhepml/0.2.1-omkpbe2/interface",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/ktjet/1.06-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/jimmy/4.2-gnimlf4/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/jemalloc/5.1.0/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/madgraph5amcatnlo/2.6.0-gnimlf12",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/heppdt/3.03.00-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/hector/1.3.4_patch1-gnimlf8/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gsl/2.2.1-gnimlf2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libjpeg-turbo/1.3.1-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/giflib/4.2.3-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/vecgeom/v00.05.00-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gdbm/1.10-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/freetype/2.5.3-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fftw3/3.3.2-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fftjet/1.5.0-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fastjet/3.3.0-omkpbe/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/expat/2.1.0-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dpm/1.8.0.1-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/hepmc/2.06.07-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xerces-c/3.1.3-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xz/5.2.2-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/dcap/2.47.8-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libxml2/2.9.1-omkpbe2/include/libxml2",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/curl/7.59.0/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cppunit/1.12.1-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/clhep/2.4.0.0-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/openssl/1.0.2d-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia6/426-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/photos/215.5-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/zlib-x86_64/1.2.11-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cascade/2.2.04-gnimlf4/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/bz2lib/1.0.6-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-omkpbe4/include/python2.7",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/tinyxml2/6.2.0/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/ittnotify/16.06.18-gnimlf/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gosamcontrib/2.0-20150803-omkpbe2/include",
"/usr/local/include",
"/usr/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf7/etc",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf7/etc/cling",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf7/include",
"/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_10_2_21-slc7_amd64_gcc700/build/CMSSW_10_2_21-build/slc7_amd64_gcc700/external/zlib-x86_64/1.2.11-omkpbe2/include",
"/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_10_2_21-slc7_amd64_gcc700/build/CMSSW_10_2_21-build/slc7_amd64_gcc700/external/tbb/2018_U1-omkpbe2/include",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf7/include",
"/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/jetFakes/fitter/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "estimateFake_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "estimateFake_C_ACLiC_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "./estimateFake.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"fitNumerator", payloadCode, "@",
"getDataDenominator", payloadCode, "@",
"getDataNumerator", payloadCode, "@",
"getDataSideband", payloadCode, "@",
"getMCSideband", payloadCode, "@",
"getPromptTemplate", payloadCode, "@",
"init", payloadCode, "@",
"options", payloadCode, "@",
"optionsFile", payloadCode, "@",
"writeString", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("estimateFake_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_estimateFake_C_ACLiC_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_estimateFake_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_estimateFake_C_ACLiC_dict() {
  TriggerDictionaryInitialization_estimateFake_C_ACLiC_dict_Impl();
}

///////////////////////////////////////////////////////////////////////////////
// $Id: CaloClusterCheck_module.cc,v 1.3 2014/11/17 21:29:51 murat Exp $
// $Author: murat $
// $Date: 2014/11/17 21:29:51 $
//
//
// 2018-02-06: compare the fast clustering to the regular one
//
// bit  1: print run/event number
// bit 11: print diag if ecl > 0, efcl < 0
//
// .fcl file to use: murat/test/trackRecoCheck.fcl
///////////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Selector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
// #include "Mu2eUtilities/inc/SimParticlesWithHits.hh"

#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "Stntuple/mod/StntupleModule.hh"

#include "TH1F.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class CaloClusterCheck : public StntupleModule {
  private:
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    struct Hist_t {
      TH1F*   ncl;                   // radial distance between the StepPointMC and the corresponding hit
      TH1F*   nfcl;			// hit residual
      TH1F*   emax;
      TH1F*   efmax;
      TH1F*   de;
    } ;

    Hist_t _hist;

    string  _caloClusterCollTag;
    string  _fastCaloClusterCollTag;

    const CaloClusterCollection  *_cccol, *_fcccol;

  public:
    explicit CaloClusterCheck(fhicl::ParameterSet const& pset);
    virtual ~CaloClusterCheck();
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual void     endJob  ();
    virtual bool     filter (art::Event& Evt);
  };


//-----------------------------------------------------------------------------
  CaloClusterCheck::CaloClusterCheck(fhicl::ParameterSet const& pset): 
    StntupleModule(pset,"CaloClusterCheck"),
    _caloClusterCollTag     (pset.get<std::string>("caloClusterCollTag"        )),
    _fastCaloClusterCollTag (pset.get<std::string>("fastCaloClusterCollTag"    ))
  {

  }

//-----------------------------------------------------------------------------
  CaloClusterCheck::~CaloClusterCheck() { 
  }


//-----------------------------------------------------------------------------
  void CaloClusterCheck::endJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }

//-----------------------------------------------------------------------------
  void CaloClusterCheck::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;

    _hist.ncl          = tfs->make<TH1F>("ncl" ,"N(clusters)"     ,100,0,100);
    _hist.nfcl         = tfs->make<TH1F>("nfcl","N(fast clusters)",100,0,100);
    _hist.emax         = tfs->make<TH1F>("emax","Emax"            ,400,0,200);
    _hist.efmax        = tfs->make<TH1F>("efmax","Efmax"          ,400,0,200);
    _hist.de           = tfs->make<TH1F>("de"   ,"emax-efmax"     ,200,-100,100);
//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------
  }


  //-----------------------------------------------------------------------------
  bool CaloClusterCheck::filter(art::Event& Evt) {
    const char* oname = "CaloClusterCheck::filter";

    if (DebugBit(1)) printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());
//-----------------------------------------------------------------------------
// find data
//-----------------------------------------------------------------------------
    auto ccH     = Evt.getValidHandle<CaloClusterCollection>(_caloClusterCollTag);
    _cccol  = ccH.product();
    auto fccH    = Evt.getValidHandle<CaloClusterCollection>(_fastCaloClusterCollTag);
    _fcccol = fccH.product();

    if ((_cccol == NULL) || (_fcccol == NULL)) {
      printf(" ERROR: no cluster collections\n");
    }

    int ncl  =  _cccol->size();
    int nfcl =  _fcccol->size();

    double ecl(-1.), efcl(-1.);
    
    if (ncl > 0) {
      const CaloCluster* cl  = &_cccol->at(0);
      ecl  = cl->energyDep();
    }

    if (nfcl > 0) {
      const CaloCluster* fcl = &_fcccol->at(0);
      efcl = fcl->energyDep();
    }

    _hist.ncl->Fill(ncl);
    _hist.nfcl->Fill(nfcl);

    _hist.emax->Fill(ecl);
    _hist.efmax->Fill(efcl);

    if (DebugBit(11)) {
      if ((ecl > 0) && (efcl < 0)) {
	printf("[%s]:011 RUN: %10i EVENT: %10i ecl = %8.3f efcl = %8.3f\n",oname,Evt.run(),Evt.event(),ecl,efcl);
      }
    }

    if ((ecl > 0) && (efcl > 0)) {
      double de = ecl-efcl;
      _hist.de->Fill(de);

      if (DebugBit(12)) {
	if (de > 20) {
	  printf("[%s]:012 RUN: %10i EVENT: %10i ecl = %8.3f efcl = %8.3f\n",oname,Evt.run(),Evt.event(),ecl,efcl);
	}
      }
    }
    
    return 1;
  }


}

using mu2e::CaloClusterCheck;
DEFINE_ART_MODULE(CaloClusterCheck);

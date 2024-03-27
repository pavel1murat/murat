///////////////////////////////////////////////////////////////////////////////
// $Id: GenpHist_module.cc,v 1.11 2014/10/02 17:15:09 murat Exp $
// $Author: murat $
// $Date: 2014/10/02 17:15:09 $
//
// .fcl file to use: murat/test/genp_hist.fcl
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

#include "Offline/MCDataProducts/inc/GenParticle.hh"

#include "Stntuple/mod/StntupleModule.hh"

#include "TH1F.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class GenpHist : public StntupleModule {
  private:
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    art::InputTag      _genpCollTag;

//-----------------------------------------------------------------------------
// histogramming
//-----------------------------------------------------------------------------
    enum { kMaxGenpHistSets = 100 };

    struct GenpHist_t {
      TH1F*   energy [2];
      TH1F*   mom    [2];//to be filled
      TH1F*   pdgId;
      TH1F*   costh;
      TH1F*   phi;
    } _hist[kMaxGenpHistSets];

  public:
    explicit GenpHist(fhicl::ParameterSet const& pset);
    virtual ~GenpHist();

    void     book_histograms();
    void     bookGenpHistograms(art::TFileDirectory* Dir, GenpHist_t* Hist);

    void     fill_genp_histograms(GenpHist_t* Hist, const GenParticle* Genp);
    // void     fill_histograms();
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob()                      override;
    virtual void     beginRun(const art::Run& )      override;
    virtual void     endJob  ()                      override;
    virtual void     analyze (const art::Event& Evt) override;
  };


//-----------------------------------------------------------------------------
  GenpHist::GenpHist(fhicl::ParameterSet const& pset): StntupleModule (pset,"GenpHist"),
    _genpCollTag              (pset.get<art::InputTag>("genpCollTag"))
  {

  }

//-----------------------------------------------------------------------------
  GenpHist::~GenpHist() { 
  }

//-----------------------------------------------------------------------------
  void GenpHist::endJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }


//-----------------------------------------------------------------------------
  void GenpHist::bookGenpHistograms(art::TFileDirectory* Dir, GenpHist_t* Hist) {

    Hist->energy[0] = Dir->make<TH1F>("e_0"  , "Energy[0]"   , 2400,   0.0,   120);
    Hist->energy[1] = Dir->make<TH1F>("e_1"  , "Energy[1]"   , 2400,   0.0,  1200);
    Hist->mom   [0] = Dir->make<TH1F>("p_0"  , "Momentum"    ,  500,   0.0,  1000);
    Hist->mom   [1] = Dir->make<TH1F>("p_1"  , "Momentum"    ,  500,   0.0,  1000); //want to know z momentum as well
    Hist->pdgId     = Dir->make<TH1F>("pdgId", "PDG ID"      ,  500,  -250,   250);
    Hist->costh     = Dir->make<TH1F>("costh", "cos(th)"     ,  200,  -1  ,   1);
    Hist->phi       = Dir->make<TH1F>("phi"  , "phi"         ,  315,  -3.15,   3.15);
  }

//-----------------------------------------------------------------------------
  void GenpHist::book_histograms() {

    art::ServiceHandle<art::TFileService> tfs;

    int book_genp_histograms[kMaxGenpHistSets];
    for (int i=0; i<kMaxGenpHistSets; i++)  book_genp_histograms[i] = 0;

    book_genp_histograms[0] = 1;     // all
    book_genp_histograms[1] = 1;     // electrons

    for (int i=0; i<kMaxGenpHistSets; i++) {
      if (book_genp_histograms[i] == 0) continue;

      art::TFileDirectory dir = tfs->mkdir(Form("genp_%02i",i));
      bookGenpHistograms(&dir,&_hist[i]);
    }
  }

//-----------------------------------------------------------------------------
  void GenpHist::beginJob() {
    book_histograms();
  }

//-----------------------------------------------------------------------------
  void GenpHist::beginRun(const art::Run& ) {
  }


  //-----------------------------------------------------------------------------
  void GenpHist::fill_genp_histograms(GenpHist_t* Hist, const GenParticle* Genp) {

    Hist->pdgId->Fill(Genp->pdgId());
    
    double e = Genp->momentum().e();
    double p = Genp->momentum().vect().mag();
    
    Hist->energy[0]->Fill(e);
    Hist->energy[1]->Fill(e);
    Hist->mom   [0]->Fill(p);
    Hist->mom   [1]->Fill(p);

    
    double costh(-100.), phi(-100.);

    if (p > 0) { 
      costh = Genp->momentum().z()/p;
      phi   = atan2(Genp->momentum().y(),Genp->momentum().x());
    }

    Hist->costh->Fill(costh);
    Hist->phi->Fill(phi);
  }


//-----------------------------------------------------------------------------
  void GenpHist::analyze(const art::Event& rEvent) {

   auto genColl = rEvent.getValidHandle<GenParticleCollection>( _genpCollTag);

    for (const auto& genp: *genColl) {
      fill_genp_histograms(&_hist[0],&genp);

      if (abs(genp.pdgId()) == 11) fill_genp_histograms(&_hist[1],&genp);
    }
  }
}

using mu2e::GenpHist;
DEFINE_ART_MODULE(GenpHist)

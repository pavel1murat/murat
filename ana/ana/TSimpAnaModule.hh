///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TSimpAnaModule_hh
#define murat_ana_TSimpAnaModule_hh

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TTrackStrawHitBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"
#include "Stntuple/geom/TStnCrystal.hh"

#include "murat/ana/TAnaModule.hh"

namespace murat {
class TSimpAnaModule: public murat::TAnaModule {
public:
//-----------------------------------------------------------------------------
  enum { kNEventHistSets  = 100 };
  enum { kNGenpHistSets   = 100 };
  enum { kNSimpHistSets   = 100 };

  struct Hist_t {
    EventHist_t*           fEvent  [kNEventHistSets  ];
    GenpHist_t*            fGenp   [kNGenpHistSets   ];
    SimpHist_t*            fSimp   [kNSimpHistSets   ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TGenpBlock*       fGenpBlock;
  TSimpBlock*       fSimpBlock;

  EventPar_t        fEvtPar;

  SimPar_t          fSimPar;
					// histograms filled
  Hist_t            fHist;

  TGenParticle*     fParticle;		// electron or muon

  int               fNElectrons;
  int               fNMuons[2];         // [0]: total, [1]: 'trigger', P > 90
  int               fMinNStrawHits;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TSimpAnaModule(const char* name="SimpAna", const char* title="SimpAna");
  ~TSimpAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int     BeginJob();
  // virtual int     BeginRun();   // use TAnaModule::BeginRun()
  virtual int     Event   (int ientry);
  virtual int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookHistograms();
  void    FillHistograms();

  void    Debug();

  ClassDef(murat::TSimpAnaModule,0)
};
}
#endif

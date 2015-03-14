//--------------------------------------------------------------------------
// Class: LLHRFilter (Stntuple-based)
//
// Environment: Software developed for the CDF at FNAL.
//
// Copyright Information: 
//	Copyright (C) 1997		Fermilab
//
//  revision history:
//  -----------------
//------------------------------------------------------------------------
#ifndef murat_mod_LLHRFilter
#define murat_mod_LLHRFilter

#include "TNamed.h"

#include "Stntuple/mod/StntupleModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"

#include "Stntuple/alg/TEmuLogLH.hh"
#include "Stntuple/alg/TStnTrackID.hh"

namespace mu2e {
class LLHRFilter : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:
					// process name, default - PROD
  std::string   fProcessName;
//-----------------------------------------------------------------------------
// module parameters
//-----------------------------------------------------------------------------
  std::string   fG4ModuleLabel;
  std::string   fStrawHitMaker;
  std::string   fTrkPatRecDem;
  std::string   fTrkPatRecUem;
  std::string   fCaloCrystalHitMaker;
  std::string   fCaloClusterMaker;
  std::string   fTrkExtrapol;
  std::string   fTrkCalMatch;
  std::string   fPidDem;

  int           fFilterEp;
  double        fMinEP;
  int           fFilterLLHRCal;
  double        fMaxLLHRCal;
  
  double        fMinTActive  ;  // start of the active window

  TNamed*       fVersion;

  TStnTrackBlock*  fTrackBlock;

  TEmuLogLH*     fLogLH;
  TStnTrackID*   fTrackID;

  int           fPointersInitialized;
//------------------------------------------------------------------------------
// function members
//------------------------------------------------------------------------------
public:
					// constructors and destructor

  LLHRFilter(fhicl::ParameterSet const& pset);

  ~LLHRFilter();
//-----------------------------------------------------------------------------
// functions of the module
//-----------------------------------------------------------------------------
  // static int InitCalDataBlock (TStnDataBlock* Block);
  // static int InitElectronBlock(TStnDataBlock* Block);
  int InitTrackBlock   (TStnDataBlock* Block);
  //  static int InitTriggerBlock (TStnDataBlock* Block);

					// ****** setters

//-----------------------------------------------------------------------------
// overloaded virtual functions of EDFilter
//-----------------------------------------------------------------------------
  virtual bool beginRun(art::Run& ARun);
  virtual bool endRun  (art::Run& ARun);
  virtual void beginJob();
  virtual void endJob  ();
  virtual bool filter  (AbsEvent& event);

  //  ClassDef(LLHRFilter,0)
};
} // end namespace mu2e

#endif

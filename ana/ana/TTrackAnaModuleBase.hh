///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTrackAnaModuleBase_hh
#define murat_ana_TTrackAnaModuleBase_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TVdetDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/obj/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TTrackAnaModuleBase: public TStnModule {
public:
//-----------------------------------------------------------------------------
// track and sim particle additional parameters
//-----------------------------------------------------------------------------
#include "murat/ana/TrackPar_t.hh"
#include "murat/ana/SimPar_t.hh"

  enum { kNDisks    = 2,
	 kMaxNTrkID = 20		// max number of trtck ID objects
  };
//-----------------------------------------------------------------------------
//  data members: no data blocks! 
//-----------------------------------------------------------------------------
public:
					// additional track parameters (assume ntracks < 20)
  TrackPar_t        fTrackPar[20];
  SimPar_t          fSimPar;

  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  int               fCalorimeterType;
  TDiskCalorimeter* fDiskCalorimeter;
					// best track in the event
  TStnTrack*        fTrack;
  TStnCluster*      fCluster;

  int               fNTrkID;            // 0:SetC 1:DaveTrkQual>0.1 2:DaveTrkQual>0.4
  TStnTrackID*      fTrackID[kMaxNTrkID];
  int               fBestID;

  TEmuLogLH*        fLogLH;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackAnaModuleBase(const char* name="TrackAnaBase", const char* title="TrackAnaBase");
  ~TTrackAnaModuleBase();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TStnTrackID*       GetTrackID(int I) { return fTrackID[I];   }
  TEmuLogLH*         GetLogLH       () { return fLogLH;        }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void   SetPdgCode      (int Code) { fPdgCode       = Code; }
  void   SetGeneratorCode(int Code) { fGeneratorCode = Code; }
//-----------------------------------------------------------------------------
// other methods: defaults, could be redefined
//-----------------------------------------------------------------------------
  virtual int     InitTrackPar(TStnTrackBlock*   TrackBlock  , 
			       TStnClusterBlock* ClusterBlock, 
			       TrackPar_t*       TrackPar    );

  virtual void    FillEfficiencyHistograms(TStnTrackBlock* TrackBlock, 
					   TStnTrackID*    TrackID   , 
					   int HistSet);
//-----------------------------------------------------------------------------
// virtual functions
//-----------------------------------------------------------------------------
  virtual void  BookEventHistograms   (HistData_t* Hist, const char* Folder) = 0;
  virtual void  BookTrackHistograms   (HistData_t* Hist, const char* Folder) = 0;

  virtual void  FillEventHistograms   (HistData_t* Hist                    ) = 0;
  virtual void  FillTrackHistograms   (HistData_t* Hist, TStnTrack*  Trk   ) = 0;
  
  virtual void    BookHistograms() = 0;
  virtual void    FillHistograms() = 0;

  ClassDef(TTrackAnaModuleBase,0)
};

#endif
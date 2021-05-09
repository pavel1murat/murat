//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 27 23:03:20 2021 by ROOT version 6.18/04
// from TTree t/combo hits
// found on file: Memory Directory
//////////////////////////////////////////////////////////

#ifndef murat_gui_combo_hits_h
#define murat_gui_combo_hits_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TGraph.h>

#include "murat/gui/TComboHitData.hh"
//-----------------------------------------------------------------------------
class combo_hits : public TNamed {

  struct StationData_t {
  public:
    int       fN;
    TObjArray fHits;

    TComboHitData* GetHit  (int I) { return (TComboHitData*) fHits.UncheckedAt(I); }

    int            GetNHits     () { return fHits.GetEntriesFast(); }
  };

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
//-----------------------------------------------------------------------------
   StationData_t   fStationData [18];
//-----------------------------------------------------------------------------
// Fixed size dimensions of array or collections stored in the TTree if any.
// Declaration of the leaf types
//-----------------------------------------------------------------------------
   Int_t           i;
   Int_t           nsh;
   Int_t           sid;
   Char_t          flags[9];
   Int_t           pln;
   Int_t           pnl;
   Int_t           lay;
   Int_t           str;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         time;
   Float_t         edep;
   Int_t           end;
   Float_t         drtime;
   Float_t         prtime;
   Float_t         tres;
   Float_t         wdist;
   Float_t         wres;
   Int_t           pdg;
   Int_t           pdgm;
   Int_t           gen;
   Int_t           id;
   Float_t         p;
   Float_t         pz;
//-----------------------------------------------------------------------------
// List of branches
//-----------------------------------------------------------------------------
   TBranch        *b_i;   //!
   TBranch        *b_nsh;   //!
   TBranch        *b_sid;   //!
   TBranch        *b_flags;   //!
   TBranch        *b_pln;   //!
   TBranch        *b_pnl;   //!
   TBranch        *b_lay;   //!
   TBranch        *b_str;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_time;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_end;   //!
   TBranch        *b_drtime;   //!
   TBranch        *b_prtime;   //!
   TBranch        *b_tres;   //!
   TBranch        *b_wdist;   //!
   TBranch        *b_wres;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_pdgm;   //!
   TBranch        *b_gen;   //!
   TBranch        *b_id;   //!
   TBranch        *b_p;   //!
   TBranch        *b_pz;   //!
//-----------------------------------------------------------------------------
// fumctions
//-----------------------------------------------------------------------------
   combo_hits(const char* Name, const char* Fn = 0);
   virtual ~combo_hits();

   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
//-----------------------------------------------------------------------------
   virtual void     BookHistograms();
   void             PlotXY  (int Station = -1, float TMin = 0., float TMax = 2000.);
   void             PlotTZ  (float TMin = 0., float TMax = 2000.);
   void             PlotTime(int Station = -1);
//-----------------------------------------------------------------------------
// flag = 0: all
//      = 1: electrons P > 20
//      = 2: electrons P < 20
//      = 3: positrons
//      = 4: mu-
//      = 5: mu+
//      = 6: protons
//      = 7: the rest
//-----------------------------------------------------------------------------
   void             Loop(int Flag, float MinE = -1., float TMin = 0, float TMax = 2000);
};

#endif

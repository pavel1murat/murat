////////////////////////////////////////////////////////////////////////////////
// 2020-08-01 P.Murat
// simulate timing distribution for the out-of-time RPC component 
////////////////////////////////////////////////////////////////////////////////
#ifndef __murat_rpc_timing_hh__
#define __murat_rpc_timing_hh__

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TRandom3.h"
#include "TH1F.h"

class rpc_timing {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   struct Hist_t {
     TH1D*           fTProt;
     TH1D*           fTStop;
     TH1D*           fTStopW;
     TH1D*           fTVertUU;             // vertex time, unwrapped, unweighted
     TH1D*           fTVertWU;
     TH1D*           fTVertWW;
   } fHist;

   TRandom3*       fRn3;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         stops_x;
   Float_t         stops_y;
   Float_t         stops_z;
   Float_t         stops_time;
   Float_t         stops_tauNormalized;

   // List of branches
   TBranch        *b_stops;   //!

   rpc_timing(TTree *tree=0);
   virtual ~rpc_timing();

   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int N=-1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

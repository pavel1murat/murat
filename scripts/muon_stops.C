//////////////////////////////////////////////////////////
// handle new stopped muon file
// This class has been automatically generated on
// Thu Jan 23 13:17:30 2014 by ROOT version 5.34/09
// from TTree stops/Stopped particles ntuple
// found on file: /mu2e/data/tdr/beam/g4s3p3/mergedMuonStops/mustops.1025a_1025a_1316a.14278089.root
//////////////////////////////////////////////////////////

#ifndef stops_h
#define stops_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class stops {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         stops_x;
   Float_t         stops_y;
   Float_t         stops_z;
   Float_t         stops_time;

   // List of branches
   TBranch        *b_stops;   //!

   stops(TTree *tree=0);
   virtual ~stops();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


  TH1F*  h_time;
  TH1F*  h_z;

};

#endif

#ifdef stops_cxx
stops::stops(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mu2e/data/tdr/beam/g4s3p3/mergedMuonStops/mustops.1025a_1025a_1316a.14278089.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/mu2e/data/tdr/beam/g4s3p3/mergedMuonStops/mustops.1025a_1025a_1316a.14278089.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/mu2e/data/tdr/beam/g4s3p3/mergedMuonStops/mustops.1025a_1025a_1316a.14278089.root:/stoppedMuonDumper");
      dir->GetObject("stops",tree);

   }
   Init(tree);

   h_time = new TH1F("h_time","muon time",1000,0,10000);
   h_z    = new TH1F("h_z"   ,"stopped muon Z",1000,5400,6400);
}

stops::~stops()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();

   delete h_time;
   delete h_z;
}

Int_t stops::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t stops::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void stops::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("stops", &stops_x, &b_stops);
   Notify();
}

Bool_t stops::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void stops::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t stops::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef stops_cxx



#define stops_cxx
#include "stops.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void stops::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L stops.C
//      Root > stops t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      h_time->Fill(stops_time);
      h_z->Fill(stops_z);
   }
}


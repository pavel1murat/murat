//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 23 17:02:20 2014 by ROOT version 5.34/09
// from TTree stops/Stopped muons
// found on file: stoppedMuons_02.root
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class muon_stops_2012 {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
  Float_t         x;
  Float_t         y;
  Float_t         z;
  Float_t         time;

   // List of branches
  TBranch        *b_x;   //!
  TBranch        *b_y;   //!
  TBranch        *b_z;   //!
  TBranch        *b_time;   //!
  
  TH1F*  h_time;
  TH1F*  h_z;
  TH2F*  h_xy;
 
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  muon_stops_2012(TTree *tree=0);

  virtual ~muon_stops_2012();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

//-----------------------------------------------------------------------------
muon_stops_2012::muon_stops_2012(const char* FnTree) : fChain(0) {

  TTree* tree;

  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(FnTree);
  if (!f || !f->IsOpen()) {
    f = new TFile(FnTree);
  }
  
  f->GetObject("stops",tree);

  Init(tree);
  
  h_time = new TH1F("h_time","muon time",1000,0,10000);
  h_z    = new TH1F("h_z"   ,"stopped muon Z",500,5400,6400);
  h_xy   = new TH2F("h_xy"  ,"stopped muon X:Y",100,-100,100,100,-100,100);
  
}

muon_stops_2012::~muon_stops_2012()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t muon_stops_2012::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t muon_stops_2012::LoadTree(Long64_t entry)
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

void muon_stops_2012::Init(TTree *tree)
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

   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("time", &time, &b_time);
   Notify();
}

Bool_t muon_stops_2012::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void muon_stops_2012::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t muon_stops_2012::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void muon_stops_2012::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L muon_stops_2012.C
//      Root > muon_stops_2012 t
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
      nb = fChain->GetEntry(jentry);   
      nbytes += nb;


      h_time->Fill(time);
      h_z->Fill(z);
      h_xy->Fill(x+3904.,y);

   }
}

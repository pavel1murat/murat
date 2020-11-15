//-----------------------------------------------------------------------------
//   In a ROOT session, you can do:
//      root> .L stops.C
//      root> stops t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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
//-----------------------------------------------------------------------------
#include "ana/rpc_timing.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//-----------------------------------------------------------------------------
void rpc_timing::Loop(int N) {

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   fHist.fTProt  = new TH1D("tprot" ,"t proton",2000,0,2000);

   fHist.fTStop  = new TH1D("ts" ,"t stop, unweighted",2000,0,2000);
   fHist.fTStopW = new TH1D("tsw","t stop, weighted"  ,2000,0,2000);

   fHist.fTVertUU  = new TH1D("tvuu","T(vertex), unwrapped, unweighted",2000,0,2000);
   fHist.fTVertWU  = new TH1D("tvwu","T(vertex), wrapped  , unweighted",2000,0,2000);
   fHist.fTVertWW  = new TH1D("tvww","T(vertex), wrapped  , weighted"  ,2000,0,2000);

   Long64_t nbytes = 0, nb = 0;

   int nsamples = 1000;

   fRn3 = new TRandom3();

   Long64_t nent = nentries;

   if (N > 0) nent = N;

   for (Long64_t jentry=0; jentry<nent;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      fHist.fTStop->Fill(stops_time);
      fHist.fTStopW->Fill(stops_time,exp(-stops_tauNormalized));

      for (int is=0; is<nsamples; is++) {
	float tprot = fRn3->Rndm()*(1570-125)+125;

	fHist.fTProt->Fill(tprot);

	float t0 = stops_time+tprot;

	fHist.fTVertUU->Fill(t0);

	float tw  = t0;
	if (tw > 1695.) {
	  int ix = (int) (tw/1695.);
	  tw = tw-ix*1695;
	}

	fHist.fTVertWU->Fill(tw);
	fHist.fTVertWW->Fill(tw,exp(-stops_tauNormalized));
      }
   }
}

//-----------------------------------------------------------------------------
rpc_timing::rpc_timing(TTree *tree) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("stops/nts.mu2e.su2020.bpim0s31b0.001000_10000000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("stops/nts.mu2e.su2020.bpim0s31b0.001000_10000000.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("stops/nts.mu2e.su2020.bpim0s31b0.001000_10000000.root:/stoppedPionDumper");
      dir->GetObject("stops",tree);

   }
   Init(tree);
}

rpc_timing::~rpc_timing() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t rpc_timing::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t rpc_timing::LoadTree(Long64_t entry) {
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

//-----------------------------------------------------------------------------
void rpc_timing::Init(TTree *tree) {
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

//-----------------------------------------------------------------------------
Bool_t rpc_timing::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

//-----------------------------------------------------------------------------
void rpc_timing::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

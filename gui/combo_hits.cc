///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TComboHitData.hh"
#include "murat/gui/TComboHitVisNode.hh"
#include "murat/gui/combo_hits.hh"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFolder.h"

using murat::TComboHitVisNode;

//-----------------------------------------------------------------------------
combo_hits::combo_hits(const char* Name, const char* Fn) : TNamed(Name,Name), fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  TTree* tree;

  TFolder* f = (TFolder*) gROOT->GetRootFolder()->FindObject("ROOT Memory");

  TString name = Form("%s_t_combo_hits",GetName());
  tree = (TTree*) f->FindObject(name.Data());

  if (tree == nullptr) {
    TString format = "i/I:nsh/I:sid/I:flags/C:pln/I:pnl/I:lay/I:str/I:x/F:y/F:z/F:time/F:";
    format        += "edep/F:end/I:drtime/F:prtime/F:tres/F:wdist/F:wres/F:pdg/I:pdgm/I:gen/I:id/I:p/F:pz/F";

    tree = new TTree(name.Data(),"Combo hit tree");
    tree->ReadFile(Fn,format.Data(),' ');
  }

  printf(" input tree: Nhits: %lli\n",tree->GetEntries());

  fHitData = new HitData_t();
  fNode    = new TComboHitVisNode("emoe_TComboHitVisNode",fHitData);

  printf("in the combo_hits constructor\n");
  fNode->Print();

  TStnView* view  = new TStnView(TEvdManager::kTZ,-1,"TZView","TZView") ;
  view->AddNode(fNode);

  TEvdManager::Instance()->AddView(view);

  Init(tree);
}

//-----------------------------------------------------------------------------
combo_hits::~combo_hits() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   
   delete fHitData;
}

//-----------------------------------------------------------------------------
void combo_hits::BookHistograms() {
}

//-----------------------------------------------------------------------------
void combo_hits::Loop(int Flag, float EMin, float TMin, float TMax) {
//-----------------------------------------------------------------------------
//   In a ROOT session, you can do:
//      root> .L combo_hits.C
//      root> combo_hits t
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
  if (fChain == 0) return;

  for (int i=0; i<18; i++) {
    fHitData->fStation[i].fHits.Delete();
  }

  Long64_t nentries = fChain->GetEntriesFast();

  printf(" ---- nentries = %lli EMin = %5.1f TMin = %10.3f TMax=%10.3f\n",nentries,EMin,TMin,TMax);
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);   // nbytes += nb;
    printf("ientry, edep,time: %5lli %10.3f %10.3f\n",ientry,edep,time);

    if ((edep < EMin) or (time < TMin) or (time > TMax)) continue; 
//-----------------------------------------------------------------------------
// prepare data to plot Y:X by station
//-----------------------------------------------------------------------------
    int station       = pln/2;
    StationData_t* sd = &fHitData->fStation[station];

    float t = time-drtime-prtime;
    
    TComboHitData* hit = new TComboHitData(x,y,z,t,edep,pdg,p);
    
    hit->fPln  = pln;
    hit->fPnl  = pnl;
    hit->fLay  = lay;
    hit->fStr  = str;
    
    printf("adding hit to station %i\n",station);
    sd->fHits.Add(hit);
  }
//-----------------------------------------------------------------------------
// plot TZ 
//-----------------------------------------------------------------------------
  fNode->InitEvent();
  printf("--- right after InitEvent()\n");
  fNode->Print();

  TEvdManager::Instance()->OpenTrkTZView();
  // PlotTZ(Flag);
}


//-----------------------------------------------------------------------------
Long64_t combo_hits::LoadTree(Long64_t entry) {
//-----------------------------------------------------------------------------
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
void combo_hits::Init(TTree *tree) {
//-----------------------------------------------------------------------------
// The Init() function is called when the selector needs to initialize
// a new tree or chain. Typically here the branch addresses and branch
// pointers of the tree will be set.
// It is normally not necessary to make changes to the generated
// code, but the routine can be extended by the user if needed.
// Init() will be called many times when running on PROOF
// (once per file to be processed).

// Set branch addresses and branch pointers
//-----------------------------------------------------------------------------
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("i", &i, &b_i);
   fChain->SetBranchAddress("nsh", &nsh, &b_nsh);
   fChain->SetBranchAddress("sid", &sid, &b_sid);
   fChain->SetBranchAddress("flags", flags, &b_flags);
   fChain->SetBranchAddress("pln", &pln, &b_pln);
   fChain->SetBranchAddress("pnl", &pnl, &b_pnl);
   fChain->SetBranchAddress("lay", &lay, &b_lay);
   fChain->SetBranchAddress("str", &str, &b_str);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("edep", &edep, &b_edep);
   fChain->SetBranchAddress("end", &end, &b_end);
   fChain->SetBranchAddress("drtime", &drtime, &b_drtime);
   fChain->SetBranchAddress("prtime", &prtime, &b_prtime);
   fChain->SetBranchAddress("tres", &tres, &b_tres);
   fChain->SetBranchAddress("wdist", &wdist, &b_wdist);
   fChain->SetBranchAddress("wres", &wres, &b_wres);
   fChain->SetBranchAddress("pdg", &pdg, &b_pdg);
   fChain->SetBranchAddress("pdgm", &pdgm, &b_pdgm);
   fChain->SetBranchAddress("gen", &gen, &b_gen);
   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   Notify();
}

Bool_t combo_hits::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

//-----------------------------------------------------------------------------
void combo_hits::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}


//-----------------------------------------------------------------------------
Int_t combo_hits::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


//-----------------------------------------------------------------------------
void combo_hits::PlotXY(int Station, float TMin, float TMax) {

  TH2F* h2 = new TH2F("h2_xy",Form("h2 xy station %2i",Station),400,-800,800,800,-800,800);
  h2->SetStats(0);
  h2->Draw();

  for (int ist=0; ist<18; ist++) {
    if ((Station < 0) or (ist == Station)) {
      StationData_t* sd = &fHitData->fStation[ist];
      int n = sd->GetNHits();
      for (int i=0; i<n; i++) {
	TComboHitData* hit = sd->GetHit(i);

	if ((hit->fT >= TMin) and (hit->fT < TMax)) {
	  hit->fXYMarker.Draw();
	}
      }
    }
  }
}

//-----------------------------------------------------------------------------
void combo_hits::PlotTime(int Station) {
  TH1F* h1 = new TH1F("h1_time",Form("h1 time station %2i",Station),300,300,1800);
  h1->SetStats(0);

  for (int ist=0; ist<18; ist++) {
    if ((Station < 0) or (ist == Station)) {
      StationData_t* sd = &fHitData->fStation[ist];
      int n = sd->GetNHits();
      for (int i=0; i<n; i++) {
	TComboHitData* hit = sd->GetHit(i);
	h1->Fill(hit->fT);
	hit->Draw();
      }
    }
  }

  h1->Draw();
}


//-----------------------------------------------------------------------------
void combo_hits::PlotTZ(float TMin, float TMax) {

  TH2F* h2 = new TH2F("h2","combo hit time vs Z",1600,-1600,1600,300,300,1800);
  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("Z, mm");
  h2->GetYaxis()->SetTitle("time, ns");
  h2->Draw();
   
  for (int ist = 0; ist<18; ist++) {
    StationData_t* sd = &fHitData->fStation[ist];
    int nh = sd->GetNHits();
    for (int i=0; i<nh; i++) {
      TComboHitData* hit = sd->GetHit(i);
      if ((hit->fT >= TMin) and (hit->fT < TMax)) {
	hit->fTZMarker.Draw();
      }
    }
  }
}

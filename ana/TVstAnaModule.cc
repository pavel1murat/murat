///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
// bit 0 : all events
// bit 1 : passed events
// bit 2 : rejected events
//
// bit 3 : events with 2 hit faces
// 
// panel=0: face=0
//       1:      1
//       2:      0
//       3:      1
//       4:      0
//       5:      1
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TControlBar.h"
#include "TROOT.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TStrawHit.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include "murat/ana/TVstAnaModule.hh"

ClassImp(TVstAnaModule)
//-----------------------------------------------------------------------------
TVstAnaModule::TVstAnaModule(const char* name, const char* title): TStnModule(name,title)
{
//-----------------------------------------------------------------------------
// not noise, but not saturation as well
//-----------------------------------------------------------------------------
  fMinEDep[0] = 0.0006;
  fMinEDep[1] = 0.0004;
  fMinEDep[2] = 0.0003;
  fMinEDep[3] = 0.0004;
  fMinEDep[4] = 0.0004;
  fMinEDep[5] = 0.0002;

  fMaxEDep[0] = 0.0058;
  fMaxEDep[1] = 0.0056;
  fMaxEDep[2] = 0.0058;
  fMaxEDep[3] = 0.0058;
  fMaxEDep[4] = 0.0055;
  fMaxEDep[5] = 0.0058;

  for (int i=0; i<6; i++) {
    Panel_t* panel = fPanel+i;
    panel->fNHits      = 0;
    panel->fNGoodHits  = 0;
  }
}

//-----------------------------------------------------------------------------
TVstAnaModule::~TVstAnaModule() {
  for (int i=0; i<6; i++) {
    Panel_t* panel = fPanel+i;
    for (int ich=0; ich<100; ich++) {
      panel->fChannel[ich].Delete();
    }
  }
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TVstAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("StrawHitBlock" ,"TStrawHitBlock" ,&fStrawHitBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//_____________________________________________________________________________
int TVstAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//-----------------------------------------------------------------------------
void TVstAnaModule::BookPanelHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  PanelHist_t* hist = (PanelHist_t*) Hist;

  HBook1F(hist->fNHits            ,"nh"      ,Form("%s: N(hits)"                ,Folder),  100,    0, 100,Folder);
  HBook1F(hist->fNGoodHits        ,"nhg"     ,Form("%s: N(good hits)"           ,Folder),  100,    0, 100,Folder);
}

//-----------------------------------------------------------------------------
void TVstAnaModule::BookChannelHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  ChannelHist_t* hist = (ChannelHist_t*) Hist;

  HBook1F(hist->fNSig             ,"nsig"       ,Form("%s: N(signals)"             ,Folder),    5,    0,   5,Folder);
  HBook1F(hist->fDT10             ,"dt10"       ,Form("%s: T(1)-T(0)"              ,Folder),  200, -100, 100,Folder);
  HBook1F(hist->fDTot10           ,"dtot10"     ,Form("%s: TOT(1)-TOT(0)"          ,Folder),  200, -100, 100,Folder);
  HBook1F(hist->fDTNext           ,"dtnext"     ,Form("%s: <T>(i+1)-<T>(i)"        ,Folder),  200, -100, 100,Folder);
}

//-----------------------------------------------------------------------------
void TVstAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fRunNumber        ,"runnum"     ,Form("%s: Run   Number"           ,Folder), 1000, 1e5, 1.e5+1000,Folder);
  HBook1F(hist->fEventNumber      ,"evtnum"     ,Form("%s: Event Number"           ,Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fInstLum          ,"inst_lum"   ,Form("%s: Inst Lum"               ,Folder), 1000, 0,  20e7,Folder);
  HBook1F(hist->fNStrawHits [0]   ,"nsh_0"      ,Form("%s: N(Straw Hits)[0]"       ,Folder),  200, 0,   200,Folder);
  HBook1F(hist->fNStrawHits [1]   ,"nsh_1"      ,Form("%s: N(Straw Hits)[1]"       ,Folder), 1000, 0, 10000,Folder);
  HBook1F(hist->fNGoodStrawHits   ,"nsh_good"   ,Form("%s: N(Good Straw Hits)"     ,Folder),  200, 0,   200,Folder);

  for (int i=0; i<6; i++) {
    TString name  = Form("nshp_%i",i);
    TString title = Form("N(straw hits_%i",i);
    HBook1F(hist->fNHitsPerPanel[i],name.Data(),title.Data(),200, 0,   200,Folder);
  }

  HBook1F(hist->fNHitFaces        ,"nhfaces"     ,Form("%s: N(hit faces)"          ,Folder),  10, 0,     10,Folder);
  HBook1F(hist->fNHitPanels       ,"nhpanels"    ,Form("%s: N(hit panels)"         ,Folder),  10, 0,     10,Folder);

  HBook1F(hist->fNHitsPerFace[0]  ,"nshf_0"      ,Form("%s: N(hits face = 0)"      ,Folder), 100, 0,    100,Folder);
  HBook1F(hist->fNHitsPerFace[1]  ,"nshf_1"      ,Form("%s: N(hits face = 1)"      ,Folder), 100, 0,    100,Folder);

  for (int i=0; i<6; i++) {
    TString name  = Form("nshp_%i_good",i);
    TString title = Form("N(good straw hits panel %i",i);
    HBook1F(hist->fNGoodHitsPerPanel[i],name.Data(),title.Data(),100, 0,   100,Folder);
  }

  HBook1F(hist->fNGoodHitFaces    ,"nhfaces_good" ,Form("%s: N(good hit faces)"    ,Folder),  10, 0,     10,Folder);
  HBook1F(hist->fNGoodHitPanels   ,"nhpanels_good",Form("%s: N(good hit panels)"   ,Folder),  10, 0,     10,Folder);

  HBook1F(hist->fNGoodHitsPerFace[0],"nshf_0_good",Form("%s: N(good hits face = 0)",Folder), 100, 0,    100,Folder);
  HBook1F(hist->fNGoodHitsPerFace[1],"nshf_1_good",Form("%s: N(good hits face = 1)",Folder), 100, 0,    100,Folder);
}

//-----------------------------------------------------------------------------
void TVstAnaModule::BookStrawHitHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  StrawHitHist_t* hist = (StrawHitHist_t*) Hist;

  HBook1F(hist->fEDep         ,"edep0"      ,Form("%s: EDep"                 ,Folder),  200,    0,  0.01,Folder);
  HBook1F(hist->fTime[0]      ,"t0"         ,Form("%s: Time[0]"              ,Folder),  600,    0, 30000,Folder);
  HBook1F(hist->fTime[1]      ,"t1"         ,Form("%s: Time[1]"              ,Folder),  600,    0, 30000,Folder);
  HBook1F(hist->fTOT [0]      ,"tot0"       ,Form("%s: TOT [0]"              ,Folder),  100,    0,   100,Folder);
  HBook1F(hist->fTOT [1]      ,"tot1"       ,Form("%s: TOT [1]"              ,Folder),  100,    0,   100,Folder);
  HBook1F(hist->fDt           ,"dt"         ,Form("%s: DeltaT = T1-T2"       ,Folder),  100,  -50,    50,Folder);
  HBook1F(hist->fStation      ,"station"    ,Form("%s: Station"              ,Folder),   40,    0,    40,Folder);
  HBook1F(hist->fPanel        ,"panel"      ,Form("%s: Panel"                ,Folder),   10,    0,    10,Folder);
  HBook1F(hist->fLayer        ,"layer"      ,Form("%s: Layer"                ,Folder),    5,    0,     5,Folder);
  HBook1F(hist->fStraw        ,"straw"      ,Form("%s: Straw"                ,Folder),  100,    0,   100,Folder);
  HBook1F(hist->fPreamp       ,"preamp"     ,Form("%s: Preamp"               ,Folder),  100,    0,   100,Folder);
  HBook1F(hist->fFace         ,"face"       ,Form("%s: Face"                 ,Folder),   10,    0,    10,Folder);
}


//_____________________________________________________________________________
void TVstAnaModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;		// events with      hits in two faces
  book_event_histset[ 2] = 1;		// events with good hits in two faces

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book straw hit histograms
//-----------------------------------------------------------------------------
  TString*  hit_selection  [kNStrawHitHistSets];
  for (int i=0; i<kNStrawHitHistSets; i++) hit_selection[i] = nullptr;

  hit_selection[  0] = new TString("all      hits");

  // split by panel

  hit_selection[80] = new TString("all hits panel=0"); 
  hit_selection[81] = new TString("all hits panel=1");
  hit_selection[82] = new TString("all hits panel=2");
  hit_selection[83] = new TString("all hits panel=3");
  hit_selection[84] = new TString("all hits panel=4");
  hit_selection[85] = new TString("all hits panel=5");
 
  hit_selection[90] = new TString("good hits panel=0"); 
  hit_selection[91] = new TString("good hits panel=1");
  hit_selection[92] = new TString("good hits panel=2");
  hit_selection[93] = new TString("good hits panel=3");
  hit_selection[94] = new TString("good hits panel=4");
  hit_selection[95] = new TString("good hits panel=5");
 
  // 100-196 : panel=0 plot all channels individually
  for (int i=0; i<96; i++) {
    hit_selection[100+i] = new TString(Form("panel=0, channel=%3i",i));
  }
  
  for (int i=0; i<kNStrawHitHistSets; i++) {
    if (hit_selection[i] != 0) {
      sprintf(folder_name,"strh_%i",i);
      fol          = (TFolder*) hist_folder->FindObject(folder_name);
      const char* title = hit_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,title);
      fHist.fStrawHit[i] = new StrawHitHist_t;
      BookStrawHitHistograms(fHist.fStrawHit[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book channel histograms
//-----------------------------------------------------------------------------
  char panel_folder_name[100], channel_folder_name[100];
  for (int ip=0; ip<6; ip++) {
    sprintf(panel_folder_name,"panel_%i",ip);
    hist_folder->AddFolder(panel_folder_name,panel_folder_name);
    TFolder* panel_fol  = (TFolder*) hist_folder->FindObject(panel_folder_name);
    BookPanelHistograms(&fHist.fPanel[ip],Form("Hist/%s",panel_folder_name));
    for (int ich=0; ich<96; ich++) { 
      sprintf(channel_folder_name,"ch_%02i",ich);
      panel_fol->AddFolder(channel_folder_name,channel_folder_name);
      // TFolder* channel_fol  = (TFolder*) panel_fol->FindObject(channel_folder_name);
      BookChannelHistograms(&fHist.fPanel[ip].fChannel[ich],Form("Hist/%s/%s",panel_folder_name, channel_folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TVstAnaModule::FillChannelHistograms(HistBase_t* Hist, ChannelData_t* ChData) {

  ChannelHist_t* hist = (ChannelHist_t*) Hist;
  
  hist->fNSig->Fill(ChData->fNSig);
  hist->fDT10->Fill(ChData->fDT10);
  hist->fDTot10->Fill(ChData->fDTot10);
  hist->fDTNext->Fill(ChData->fDTNext);
}

//-----------------------------------------------------------------------------
void TVstAnaModule::FillPanelHistograms(HistBase_t* Hist, Panel_t* Panel) {

  PanelHist_t* hist = (PanelHist_t*) Hist;
  
  hist->fNHits->Fill(Panel->fNHits);
  hist->fNGoodHits->Fill(Panel->fNGoodHits);
}

//-----------------------------------------------------------------------------
void TVstAnaModule::FillStrawHitHistograms(HistBase_t* Hist, TStrawHit* Hit, StrawHitPar_t* Shp) {

  StrawHitHist_t* hist = (StrawHitHist_t*) Hist;
  
  hist->fEDep->Fill(Hit->EDep());
  hist->fTime[0]->Fill(Hit->Time(0));
  hist->fTime[1]->Fill(Hit->Time(1));
  hist->fTOT [0]->Fill(Hit->TOT (0));
  hist->fTOT [1]->Fill(Hit->TOT (1));
  float dt  = Hit->Time(1) - Hit->Time(0);
  hist->fDt->Fill(dt);
  hist->fStation->Fill(Hit->Station());
  hist->fPanel->Fill(Hit->Panel());
  hist->fStraw->Fill(Hit->Straw());
  hist->fFace->Fill(Hit->Face());
  hist->fLayer->Fill(Hit->Layer());
  hist->fPreamp->Fill(Hit->Preamp());
}

//-----------------------------------------------------------------------------
void TVstAnaModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int   event_number = GetHeaderBlock()->EventNumber();
  int   run_number   = GetHeaderBlock()->RunNumber();
  float inst_lum     = GetHeaderBlock()->InstLum();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);
  hist->fInstLum->Fill(inst_lum);


  hist->fNStrawHits[0]->Fill(fAllHits.fNHits);
  hist->fNStrawHits[1]->Fill(fAllHits.fNHits);

  hist->fNGoodStrawHits->Fill(fGoodHits.fNHits);

  for (int i=0; i<2; i++) {
    hist->fNHitsPerFace[i]->Fill(fAllHits.fNHitsPerFace[i]);
    hist->fNGoodHitsPerFace[i]->Fill(fGoodHits.fNHitsPerFace[i]);
  }

  for (int i=0; i<6; i++) {
    hist->fNGoodHitsPerPanel[i]->Fill(fGoodHits.fNHitsPerPanel[i]);
  }

  hist->fNGoodHitPanels->Fill(fAllHits.fNHitPanels);
  hist->fNGoodHitFaces ->Fill(fAllHits.fNHitFaces );

}

//_____________________________________________________________________________
void TVstAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
// EVT_1: events with hits in two faces
// EVT_2: events with good hits in both faces
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

  if (fAllHits.fNHitFaces  > 1)  FillEventHistograms(fHist.fEvent[1]);
  if (fGoodHits.fNHitFaces > 1)  FillEventHistograms(fHist.fEvent[2]);
//-----------------------------------------------------------------------------
// straw hit histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNStrawHits; i++) {
    if (i > kMaxNStrawHits) {
      GetHeaderBlock()->Print(Form("ERROR: hit number %i greater than kMaxNStrawHits = %i, SKIP.",i,kMaxNStrawHits));
      break;
    }

    TStrawHit* hit     = fStrawHitBlock->Hit(i);
    StrawHitPar_t* shp = fStrawHitPar+i;

    FillStrawHitHistograms(fHist.fStrawHit[0],hit,shp);

    int panel = hit->Panel();
    FillStrawHitHistograms(fHist.fStrawHit[80+panel],hit,shp);

    if ((hit->EDep() > fMinEDep[panel]) and (hit->EDep() < fMaxEDep[panel])) {
					// good hit
      FillStrawHitHistograms(fHist.fStrawHit[90+panel],hit,shp);
    } 
  }
//-----------------------------------------------------------------------------
// panel and channel histograms
//-----------------------------------------------------------------------------
  for (int ip=0; ip<6; ip++) {
    Panel_t* panel = &fPanel[ip];

    FillPanelHistograms(&fHist.fPanel[ip],panel);

    if (panel->fNGoodHits == 0)                                       continue;

    for (int ich=0; ich<96; ich++) {
      int nh = panel->fChannel[ich].GetEntriesFast();
      if (nh > 0) {
	for (int ih=0; ih<nh; ih++) {
	  ChannelData_t* chd = (ChannelData_t*) panel->fChannel[ich].UncheckedAt(ih);

  
	  int ich_next = ich+1;
	  int nh_next = panel->fChannel[ich_next].GetEntriesFast();

	  if (nh_next > 0) {
	    for (int ih1=0; ih1<nh_next; ih1++) {
	      ChannelData_t* chd_next = (ChannelData_t*) panel->fChannel[ich_next].UncheckedAt(ih1);

	      float t      = (chd->fHit->Time(0)+chd->fHit->Time(1))/2.;
	      float t_next = (chd_next->fHit->Time(0)+chd_next->fHit->Time(1))/2.;
	      
	      chd->fDTNext = t_next-t;
	    }
	  
	    FillChannelHistograms(&fHist.fPanel[ip].fChannel[ich],chd);
	  }
	}
      }
    }
  }
}

//_____________________________________________________________________________
int TVstAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fStrawHitBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// clear the local storage
//-----------------------------------------------------------------------------
  fAllHits.Clear();
  fGoodHits.Clear();

  for (int i=0; i<6; i++) fPanel[i].Clear();

  fNStrawHits        = fStrawHitBlock->NHits();

  TStrawHit* hit;
  for (int i=0; i<fNStrawHits; i++) {
    if (i > kMaxNStrawHits) {
      GetHeaderBlock()->Print(Form("ERROR: hit number %i greater than kMaxNStrawHits = %i, skip.",i,kMaxNStrawHits));
      break;
    }

    hit     = fStrawHitBlock->Hit(i);

    int ip  = hit->Panel();
    int ich = hit->Straw();
//-----------------------------------------------------------------------------
// within the panel, everything is ordered by channel
// for the moment, ignore the case ot 2 hits in the same channel
//-----------------------------------------------------------------------------
    Panel_t* panel = &fPanel[ip];

    ChannelData_t* chd = new ChannelData_t();

					// this belongs to a constructor
    chd->fHit    = hit;
    chd->fDTNext = -999;

    if ((hit->Time(0) > 0) and (hit->Time(1) > 0)) {
      chd->fDT10   = hit->Time(1)-hit->Time(0);
      chd->fNSig   = 2;
    }
    else {
      chd->fNSig   = 1;
      if      (hit->Time(0) > 0) chd->fDT10   = -999;
      else if (hit->Time(1) > 0) chd->fDT10   =  999;
      else {
	chd->fNSig  = 0;
      }
    }

    if ((hit->TOT(0) > 0) and (hit->TOT(1) > 0)) {
      chd->fDTot10 = hit->TOT (1)-hit->TOT (0);
    }
    else { 
      if      (hit->TOT(0) > 0) chd->fDTot10 = -999;
      else if (hit->TOT(1) > 0) chd->fDTot10 =  999;
    }
    
    fAllHits.fNHits++;
    fAllHits.fNHitsPerPanel[ip]++;

    if (GoodHit(hit)) {
      fGoodHits.fNHits++;
      fGoodHits.fNHitsPerPanel[ip]++;
//-----------------------------------------------------------------------------
// for more detailed analysis, store only good hits
//-----------------------------------------------------------------------------
      panel->fNHits++;
      panel->fNGoodHits++;
      panel->fChannel[ich].Add(chd);
    }
  }

  for (int i=0; i<6; i++) {
    int face = kFace[i];
    fAllHits.fNHitsPerFace [face] += fAllHits.fNHitsPerPanel[i];
    fGoodHits.fNHitsPerFace[face] += fGoodHits.fNHitsPerPanel[i];

    if (fAllHits.fNHitsPerPanel [i] > 0) fAllHits.fNHitPanels += 1;
    if (fGoodHits.fNHitsPerPanel[i] > 0) fGoodHits.fNHitPanels += 1;
  }

  for (int i=0; i<2; i++) {
    if (fAllHits.fNHitsPerFace [i] > 0) fAllHits.fNHitFaces += 1;
    if (fGoodHits.fNHitsPerFace[i] > 0) fGoodHits.fNHitFaces += 1;
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TVstAnaModule::Debug() {

//-----------------------------------------------------------------------------  
// bit 3 : events with hits in both faces
//-----------------------------------------------------------------------------
  if (GetDebugBit(3) == 1) {
    if (fAllHits.fNHitFaces > 1) PrintStrawHitBlock();
  }
//-----------------------------------------------------------------------------  
// bit 4 : events with hits in both faces
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
    if (fGoodHits.fNHitFaces > 1) PrintStrawHitBlock();
  }
//-----------------------------------------------------------------------------  
// bit 5 : events with N(good hits) >= 2 in both faces
//-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {
    if ((fGoodHits.fNHitFaces       >  1) and 
	(fGoodHits.fNHitsPerFace[0] >= 2) and 
	(fGoodHits.fNHitsPerFace[1] >= 2)     ) {
      PrintStrawHitBlock();
    }
  }
//-----------------------------------------------------------------------------  
// bit 6 : print hits with 1315 < SppTime < 1360
//-----------------------------------------------------------------------------
  if (GetDebugBit(6) == 1) {
    int banner_printed(0);
    for (int i=0; i<fNStrawHits; i++) {
      TStrawHit* hit = fStrawHitBlock->Hit(i);
      StrawHitPar_t* shp = fStrawHitPar+i;
      if ((shp->fSppTime > 1315) && (shp->fSppTime <= 1360)) {
	if (banner_printed == 0) { 
	  GetHeaderBlock()->Print("bit 006: Event of interest");
	  PrintStrawHit(hit,"banner") ; 
	  banner_printed = 1; 
	}
	PrintStrawHit(hit,"data") ;
      }
    }
  }
//-----------------------------------------------------------------------------  
// bit 7 : print hits with SppTime > 1 ms
//-----------------------------------------------------------------------------
  if (GetDebugBit(7) == 1) {
    int nhits = 0;
    for (int i=0; i<fNStrawHits; i++) {
      // TStrawHitData* hit = fStrawHitDataBlock->Hit(i);
      StrawHitPar_t* shp = fStrawHitPar+i;
      if (shp->fSppTime > 1e6) {
	nhits += 1;
      }
    }
    if (nhits > 0) GetHeaderBlock()->Print(Form("bit 007: hits with spptime > 1e6"));
  }
}

//_____________________________________________________________________________
int TVstAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TVstAnaModule::PrintStrawHitBlock() {
  int banner_printed(0);

  for (int i=0; i<fNStrawHits; i++) {
    TStrawHit* hit = fStrawHitBlock->Hit(i);
    if (banner_printed == 0) {
      PrintStrawHit(hit,"banner");
      banner_printed = 1;
    }
    PrintStrawHit(hit,"data");
  }
}

//_____________________________________________________________________________
void TVstAnaModule::PrintStrawHit(TStrawHit* Hit, const char* Option) {

  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("----------------------------------------------------------------------------------------\n");
    printf("Index Plane Face Panel Straw Layer Preamp    T(0)   T(1) TOT(0) TOT(1)   Dt     Energy  \n");
    printf("----------------------------------------------------------------------------------------\n");
  }

  if ((opt != "") and (opt.Index("data") < 0)) return;

  printf("%5i %4i %5i %5i %5i %5i %6i  %6i %6i %5i %5i %5i %9.5f\n",
	 Hit->StrawID(), 
	 Hit->Station(),
	 Hit->Face (),
	 Hit->Panel(),
	 Hit->Straw(),
	 Hit->Layer(),
	 Hit->Preamp(),
	 Hit->Time(0),
	 Hit->Time(1),
	 Hit->TOT (0),
	 Hit->TOT (1),
	 Hit->Time(1)-Hit->Time(0),
	 Hit->EDep()
	 );
}


//-----------------------------------------------------------------------------
void TVstAnaModule::ControlBar() {

   //Add the tutorials directory to the macro path
   //This is necessary in case this macro is executed from another user directory
   TString dirName = gSystem->UnixPathName(__FILE__);
   dirName.ReplaceAll("demos.C","");
   dirName.ReplaceAll("/./","");
   const char *current = gROOT->GetMacroPath();
   gROOT->SetMacroPath(Form("%s:%s",current,dirName.Data()));

   TControlBar *bar = new TControlBar("vertical", "VST control bar",10,10);
   bar->AddButton("Next Event","g.x->Continue(1)",        "Click Here For Help on Running the Demos");
   // bar->AddButton("browser",   "new TBrowser;",         "Start the ROOT Browser");
   // bar->AddButton("framework", ".x graphics/framework.C","An Example of Object Oriented User Interface");
   // bar->AddButton("first",     ".x graphics/first.C",   "An Example of Slide with Root");
   // bar->AddButton("hsimple",   ".x hsimple.C",          "An Example Creating Histograms/Ntuples on File");
   // bar->AddButton("hsum",      ".x hist/hsum.C",        "Filling Histograms and Some Graphics Options");
   // bar->AddButton("formula1",  ".x graphics/formula1.C","Simple Formula and Functions");
   // bar->AddButton("surfaces",  ".x graphs/surfaces.C",  "Surface Drawing Options");
   // bar->AddButton("fillrandom",".x hist/fillrandom.C",  "Histograms with Random Numbers from a Function");
   // bar->AddButton("fit1",      ".x fit/fit1.C",         "A Simple Fitting Example");
   // bar->AddButton("multifit",  ".x fit/multifit.C",     "Fitting in Subranges of Histograms");
   // bar->AddButton("h1draw",    ".x hist/h1draw.C",      "Drawing Options for 1D Histograms");
   // bar->AddButton("graph",     ".x graphs/graph.C",     "Example of a Simple Graph");
   // bar->AddButton("gerrors",   ".x graphs/gerrors.C",   "Example of a Graph with Error Bars");
   // bar->AddButton("tornado",   ".x graphics/tornado.C", "Examples of 3-D PolyMarkers");
   // bar->AddButton("shapes",    ".x geom/shapes.C",      "The Geometry Shapes");
   // bar->AddButton("geometry",  ".x geom/geometry.C",    "Creation of the NA49 Geometry File");
   // bar->AddButton("na49view",  ".x geom/na49view.C",    "Two Views of the NA49 Detector Geometry");
   // bar->AddButton("file",      ".x io/file.C",          "The ROOT File Format");
   // bar->AddButton("fildir",    ".x io/fildir.C",        "The ROOT File, Directories and Keys");
   // bar->AddButton("tree",      ".x tree/tree.C",        "The Tree Data Structure");
   // bar->AddButton("ntuple1",   ".x tree/ntuple1.C",     "Ntuples and Selections");
   // bar->AddButton("rootmarks", ".x rootmarks.C",        "Prints an Estimated ROOTMARKS for Your Machine");
   bar->SetButtonWidth(90);
   bar->Show();
   gROOT->SaveContext();
}

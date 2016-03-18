///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
// 0  : all events
// 1  : passed events
// 2  : rejected events
// 
// 3  : events with set C tracks and 70mm < |dx|  < 90 mm
// 4  : events with DpF > 1 MeV : obviously, misreconstructed ones
// 5  : events with N(tracks) > 1
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TCalAnaModule.hh"

ClassImp(TCalAnaModule)
//-----------------------------------------------------------------------------
TCalAnaModule::TCalAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{

  fDiskCalorimeter = new TDiskCalorimeter();

  fMinT0 = 0; // do not cut on time by default
}

//-----------------------------------------------------------------------------
TCalAnaModule::~TCalAnaModule() {
}


//-----------------------------------------------------------------------------
void TCalAnaModule::BookCalHitHistograms(CalHitHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------

  HBook1F(Hist->fDiskID,"diskid",Form("%s: Disk ID   ",Folder),  2, 0,   2,Folder);
  HBook1F(Hist->fEnergy,"energy",Form("%s: Hit Energy",Folder),500, 0, 100,Folder);
  HBook1F(Hist->fTime  ,"time"  ,Form("%s: Hit Time  ",Folder),200, 0,2000,Folder);
}


//-----------------------------------------------------------------------------
void TCalAnaModule::BookCaloHistograms(CaloHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------
  HBook1F(Hist->fDiskID ,"disk_id",Form("%s: Disk ID"       ,Folder), 10, 0,  10,Folder);

  for (int i=0; i<kNDisks; i++) {
    HBook1F(Hist->fEnergy  [i],Form("energy_%i",i),Form("%s: Hit Energy[%i]",Folder,i),200, 0, 100,Folder);
    HBook1F(Hist->fTime    [i],Form("time_%i"  ,i),Form("%s: Hit time  [%i]",Folder,i),200, 0,2000,Folder);
    HBook1F(Hist->fNHits   [i],Form("nhits_%i" ,i),Form("%s: NHits     [%i]",Folder,i), 50, 0,  50,Folder);
    HBook1F(Hist->fRadius  [i],Form("r_%i"     ,i),Form("%s: Radius    [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(Hist->fRadiusWE[i],Form("rwe_%i"   ,i),Form("%s: RadiusWE  [%i]",Folder,i),100, 0,1000,Folder);

    HBook1F(Hist->fE700    [i],Form("e700_%i",i),Form("%s: Hit Energy[%i] (T > 700ns)",Folder,i),200, 0, 100,Folder);
    HBook1F(Hist->fT700    [i],Form("t700_%i",i),Form("%s: Hit time  [%i] (T > 700ns)",Folder,i),200, 0,2000,Folder);
    HBook1F(Hist->fN700    [i],Form("n700_%i",i),Form("%s: NHits     [%i] (T > 700ns)",Folder,i), 50, 0,  50,Folder);

    HBook1F(Hist->fR700    [i],Form("r700_%i"  ,i),Form("%s: Radius (T>700) [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(Hist->fRWE700  [i],Form("rwe700_%i",i),Form("%s: Radius*E(T>700)[%i]",Folder,i),100, 0,1000,Folder);
  }
}

//-----------------------------------------------------------------------------
void TCalAnaModule::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fDiskID ,"disk_id",Form("%s: Disk ID"       ,Folder), 10, 0,  10,Folder);
  HBook1F(Hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),500, 0, 250,Folder);
  HBook1F(Hist->fT0     ,"t0"     ,Form("%s: cluster T0"    ,Folder),200, 0,2000,Folder);
  HBook1F(Hist->fRow    ,"row"    ,Form("%s: cluster Row"   ,Folder),200, 0, 200,Folder);
  HBook1F(Hist->fCol    ,"col"    ,Form("%s: cluster column",Folder),200, 0, 200,Folder);
  HBook1F(Hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, -5000,5000,Folder);
  HBook1F(Hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200,-1000,1000,Folder);
  HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 11500,13500,Folder);
  HBook1F(Hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  1000,Folder);
  HBook1F(Hist->fYMean  ,"ymean"  ,Form("%s: cluster YMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fZMean  ,"zmean"  ,Form("%s: cluster ZMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fSigY   ,"sigy"   ,Form("%s: cluster SigY"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigZ   ,"sigz"   ,Form("%s: cluster SigZ"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigR   ,"sigr"   ,Form("%s: cluster SigR"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr0   ,"ncr0"   ,Form("%s: cluster NCR[0]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr1   ,"ncr1"   ,Form("%s: cluster NCR[1]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fSigE1  ,"sige1"   ,Form("%s: SigmaE/Etot"  ,Folder),200, 0, 10,Folder);
  HBook1F(Hist->fSigE2  ,"sige2"   ,Form("%s: SigmaE/Emean" ,Folder),200, 0, 10,Folder);
}

//-----------------------------------------------------------------------------
void TCalAnaModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP      ,"p"       ,Form("%s: Momentum"     ,Folder),1000,     0, 200,Folder);
  HBook1F(Hist->fPdgCode[0],"pdg_code_0",Form("%s: PDG Code[0]"     ,Folder),200, -100, 100,Folder);
  HBook1F(Hist->fPdgCode[1],"pdg_code_1",Form("%s: PDG Code[1]"     ,Folder),500, -2500, 2500,Folder);
  HBook1F(Hist->fGenID  ,"gen_id"  ,Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
  HBook1F(Hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
  HBook1F(Hist->fT0     ,"t0"      ,Form("%s: T0"           ,Folder), 200,     0, 2000,Folder);
  HBook1F(Hist->fR0     ,"r"       ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
  HBook1F(Hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);
}

//-----------------------------------------------------------------------------
void TCalAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 200,Folder);
  HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 400, 5000,7000,Folder);
  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder), 100,-1,    1,Folder);
  HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder), 500, 0, 1000,Folder);
  HBook1F(Hist->fNCl20     ,"ncl20"    ,Form("%s: N(Clusters), E > 20"             ,Folder), 500, 0, 1000,Folder);
  HBook1F(Hist->fNCl50     ,"ncl50"    ,Form("%s: N(Clusters), E > 50"             ,Folder), 500, 0, 1000,Folder);
  HBook1F(Hist->fNCl70     ,"ncl70"    ,Form("%s: N(Clusters), E > 70"             ,Folder), 500, 0, 1000,Folder);
  HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder), 250, 0,  500,Folder);
  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder), 500, 0,10000,Folder);

  for (int i=0; i<kNDisks; i++) {
    HBook1F(Hist->fNHits       [i],Form("nhits_%i"   ,i),Form("%s: N(Hits) [%i]"           ,Folder,i),  500,0,10000,Folder);
    HBook1F(Hist->fETot        [i],Form("etot_%i"    ,i),Form("%s: Etot[%i]"               ,Folder,i), 1000,0,10000,Folder);
    HBook1D(Hist->fECrVsR      [i],Form("ecr_vs_r_%i",i),Form("%s: E Cr Vs R [%i]"         ,Folder,i),  100,0, 1000,Folder);
    HBook1D(Hist->fNHitsVsR    [i],Form("nh_vs_r_%i" ,i),Form("%s: N hits Vs R [%i]"       ,Folder,i),  100,0, 1000,Folder);
    HBook1D(Hist->fNCrystalsVsR[i],Form("ncr_vs_r_%i",i),Form("%s: N Hit Crystals[%i] vs R",Folder,i),  100,0, 1000,Folder);

    HBook2F(Hist->fEHitVsR     [i],Form("eh_vs_r_%i",i),Form("%s: E Hit Vs R [%i]",Folder,i),100,0,1000,250,0,500,Folder);
  }

  HBook1F(Hist->fNHitsTot       ,"nhits_tot",Form("%s: Total number of hits",Folder), 400, 0,20000,Folder);
  HBook1F(Hist->fNHitCrystalsTot,"nhcr_tot" ,Form("%s: NHit Crystals Tot"   ,Folder), 800, 0,4000,Folder);

  HBook1F(Hist->fECaloTot       ,"ecal"     ,Form("%s: Total energy in the calorimeter",Folder), 500,0,250,Folder);
}

// //-----------------------------------------------------------------------------
// void TCalAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
//   //  char name [200];
//   //  char title[200];

//   HBook1F(Hist->fPdgCode   ,"pdg"         ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
//   HBook1F(Hist->fNStrawHits,"nsth"        ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
//   HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
//   HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
// }


//_____________________________________________________________________________
void TCalAnaModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

//-----------------------------------------------------------------------------
// book crystal histograms
//-----------------------------------------------------------------------------
  HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
  HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;		// EMin = 0.1
  book_event_histset[ 2] = 1;		// EMin = 0.5
  book_event_histset[ 3] = 1;		// EMin = 1.0
  book_event_histset[ 4] = 1;		// EMin = 1.5
  book_event_histset[ 5] = 1;		// EMin = 2.0
  book_event_histset[ 6] = 1;		// EMin = 1.0 , TMin = 400
  book_event_histset[ 7] = 1;		// EMin = 1.0 , TMin = 500
  book_event_histset[ 8] = 1;		// EMin = 1.0 , TMin = 600

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
// book simp histograms
//-----------------------------------------------------------------------------
//   int book_simp_histset[kNSimpHistSets];
//   for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

//   book_simp_histset[ 0] = 1;		// all events

//   for (int i=0; i<kNSimpHistSets; i++) {
//     if (book_simp_histset[i] != 0) {
//       sprintf(folder_name,"sim_%i",i);
//       fol = (TFolder*) hist_folder->FindObject(folder_name);
//       if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
//       fHist.fSimp[i] = new SimpHist_t;
//       BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
//     }
//   }
//-----------------------------------------------------------------------------
// book cluster histograms
//-----------------------------------------------------------------------------
  int book_cluster_histset[kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

  book_cluster_histset[0] = 1;		// all clusters
  book_cluster_histset[1] = 1;		// E > 20 MeV
  book_cluster_histset[2] = 1;		// E > 50 MeV
  book_cluster_histset[3] = 1;		// E > 70 MeV

  for (int i=0; i<kNClusterHistSets; i++) {
    if (book_cluster_histset[i] != 0) {
      sprintf(folder_name,"cls_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book cal hit histograms
//-----------------------------------------------------------------------------
  int book_calhit_histset[kNCalHitHistSets];
  for (int i=0; i<kNCalHitHistSets; i++) book_calhit_histset[i] = 0;

  book_calhit_histset[0] = 1;		// all clusters

  for (int i=0; i<kNCalHitHistSets; i++) {
    if (book_calhit_histset[i] != 0) {
      sprintf(folder_name,"clh_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCalHit[i] = new CalHitHist_t;
      BookCalHitHistograms(fHist.fCalHit[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book calorimeter histograms
//-----------------------------------------------------------------------------
  int book_calo_histset[kNCaloHistSets];
  for (int i=0; i<kNCaloHistSets; i++) book_calo_histset[i] = 0;

  book_calo_histset[0] = 1;		// all crystals
  book_calo_histset[1] = 1;		// all crystals, e > 0
  book_calo_histset[2] = 1;		// all crystals, e > 0.1
  book_calo_histset[3] = 1;		// all crystals, e > 0.5
  book_calo_histset[4] = 1;		// all crystals, e > 1.0

  for (int i=0; i<kNCaloHistSets; i++) {
    if (book_calo_histset[i] != 0) {
      sprintf(folder_name,"cal_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCalo[i] = new CaloHist_t;
      BookCaloHistograms(fHist.fCalo[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Genp histograms
//-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles

  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TCalAnaModule::FillCaloHistograms(CaloHist_t* Hist, TStnCrystal* Cr) {

  int                    nhits, idisk;
  float                  t, e, r, e700, n700;
  TCalHitData*           hit;

  TDisk* disk = Cr->Disk();

  idisk = disk->SectionID();

  e     = Cr->Energy();
  r     = Cr->Radius();
  nhits = Cr->NHits();

  Hist->fDiskID->Fill(idisk);

  Hist->fEnergy  [idisk]->Fill(e);
  Hist->fNHits   [idisk]->Fill(nhits);
  //    Hist->fTime    [idisk]->Fill(t);
  Hist->fRadius  [idisk]->Fill(r);
  Hist->fRadiusWE[idisk]->Fill(r,e);
    
  e700 = 0;
  n700 = 0;

  for (int i=0; i<nhits; i++) {
    hit  = Cr->CalHitData(i);
    t   = hit->Time();
    Hist->fTime[idisk]->Fill(t);
    if (t > 700.) {
      n700 += 1;
      e700 += hit->Energy();
      Hist->fT700[idisk]->Fill(t);
    }
  }

  Hist->fE700   [idisk]->Fill(e700);
  Hist->fN700   [idisk]->Fill(n700);

  if (n700 > 0) {
    Hist->fR700  [idisk]->Fill(r);
    Hist->fRWE700[idisk]->Fill(r,e700);
  }
}


//-----------------------------------------------------------------------------
void TCalAnaModule::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {
  int   row, col;
  float  x, y, z, r;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX+3904.;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  Hist->fDiskID->Fill(Cluster->DiskID());
  Hist->fEnergy->Fill(Cluster->Energy());
  Hist->fT0->Fill(Cluster->Time());
  Hist->fRow->Fill(row);
  Hist->fCol->Fill(col);
  Hist->fX->Fill(x);
  Hist->fY->Fill(y);
  Hist->fZ->Fill(z);
  Hist->fR->Fill(r);

  Hist->fYMean->Fill(Cluster->fYMean);
  Hist->fZMean->Fill(Cluster->fZMean);
  Hist->fSigY ->Fill(Cluster->fSigY);
  Hist->fSigZ ->Fill(Cluster->fSigZ);
  Hist->fSigR ->Fill(Cluster->fSigR);
  Hist->fNCr0 ->Fill(Cluster->fNCrystals);
  Hist->fNCr1 ->Fill(Cluster->fNCr1);
  Hist->fFrE1 ->Fill(Cluster->fFrE1);
  Hist->fFrE2 ->Fill(Cluster->fFrE2);
  Hist->fSigE1->Fill(Cluster->fSigE1);
  Hist->fSigE2->Fill(Cluster->fSigE2);
}

//-----------------------------------------------------------------------------
void TCalAnaModule::FillCalHitHistograms(CalHitHist_t* Hist, TCalHitData* Hit) {

  Hist->fDiskID->Fill(-1);
  Hist->fEnergy->Fill(Hit->Energy());
  Hist->fTime  ->Fill(Hit->Time  ());
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TCalAnaModule::FillEventHistograms(EventHist_t* Hist, double EMin, double TMin) {
  double            cos_th, xv, yv, rv, zv, p;
  TLorentzVector    mom;

  fElectron->Momentum(mom);

  p      = mom.P();
  cos_th = mom.Pz()/p;

  xv     = fElectron->Vx()+3904.;
  yv     = fElectron->Vy();
  rv     = sqrt(xv*xv+yv*yv);
  zv     = fElectron->Vz();

  Hist->fEleCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  Hist->fNClusters->Fill(fNClusters);
  Hist->fNCl20->Fill(fNCl20);
  Hist->fNCl50->Fill(fNCl50);
  Hist->fNCl70->Fill(fNCl70);

  double emax = -1;

  TStnCluster* cluster(0);
  if (fNClusters > 0) cluster = fClusterBlock->Cluster(0);

  if (cluster) {
    emax   = cluster->Energy();
  }

  Hist->fEMax->Fill(emax);
  Hist->fNGenp->Fill(fNGenp);
//-----------------------------------------------------------------------------
// calorimeter
//-----------------------------------------------------------------------------
  int      ndisks, nc, nh, nhh, bin;

  int      n_hits_tot, n_hit_crystals_tot;
  int      n_hit_crystals[kNDisks], n_hits[kNDisks];
  int      n_hits_r[kNDisks][100], n_hit_crystals_r[kNDisks][100];

  float    etot[kNDisks];

  TCalHitData* hit;
  TStnCrystal* cr;

  float        ehit, r;

  ndisks = fDiskCalorimeter->NDisks();

  n_hits_tot         = 0;
  n_hit_crystals_tot = 0;

  for (int idisk=0; idisk<kNDisks; idisk++) {

    TDisk* disk = fDiskCalorimeter->Disk(idisk);

    etot          [idisk] = 0;
    n_hits        [idisk] = 0;
    n_hit_crystals[idisk] = 0;

    for (int ib=0; ib<100; ib++) {
      n_hits_r        [idisk][ib] = 0;
      n_hit_crystals_r[idisk][ib] = 0;
    }

    nc = disk->NCrystals();
    for (int ic=0; ic<nc; ic++) {
      cr   = disk->Crystal(ic);
      r    = cr->Radius();
      nh   = cr->NHits ();
      ehit = cr->Energy();
      bin  = (int) (r/10.);
					// total energy deposited in the disk
      etot[idisk] += cr->Energy();

      nhh = 0;
      for (int ih=0; ih<nh; ih++) {
	hit = cr->CalHitData(ih);
	if ((hit->Energy() > EMin) && (hit->Time() > TMin)) {
					// number of hits above the threshold (EMin)

	  nhh                         += 1; // N(hits) per crystal
	  n_hits         [idisk]      += 1; // total N(hits) in the disk
	  n_hits_r       [idisk][bin] += 1; // total N(hits) in a given radial bin
	}
      }
    
      Hist->fECrVsR   [idisk]->Fill(r,ehit);
      Hist->fEHitVsR  [idisk]->Fill(r,ehit);

      if (nhh > 0) {
	n_hit_crystals  [idisk]      += 1;
	n_hit_crystals_r[idisk][bin] += 1;
      }
    }
    n_hit_crystals_tot += n_hit_crystals[idisk];
    n_hits_tot         += n_hits[idisk];
  }

  double ecalo(0);
  for (int idisk=0; idisk<ndisks; idisk++) {
    ecalo += etot[idisk];
//-----------------------------------------------------------------------------
// fill 'per-disk' histograms
//-----------------------------------------------------------------------------
    Hist->fETot [idisk]->Fill(etot  [idisk]);
    Hist->fNHits[idisk]->Fill(n_hits[idisk]);
//-----------------------------------------------------------------------------
// 100 is fixed by the number of bins in the radial distributions
//-----------------------------------------------------------------------------
    for (int ib=0; ib<100; ib++) {
      r = (ib+0.5)*10.;
      Hist->fNHitsVsR    [idisk]->Fill(r,n_hits_r        [idisk][ib]);
      Hist->fNCrystalsVsR[idisk]->Fill(r,n_hit_crystals_r[idisk][ib]);
    }
  }

  Hist->fECaloTot->Fill(ecalo);
  Hist->fNHitCrystalsTot->Fill(n_hit_crystals_tot);
  Hist->fNHitsTot->Fill(n_hits_tot);
}

//-----------------------------------------------------------------------------
void TCalAnaModule::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp) {
  int    gen_id;
  float  p, cos_th, z0, t0, r0, x0, y0;

  TLorentzVector mom, v;

  Genp->Momentum(mom);
  //  Genp->ProductionVertex(v);

  p      = mom.P();
  cos_th = mom.CosTheta();

  x0     = Genp->Vx()+3904.;
  y0     = Genp->Vy();

  z0     = Genp->Vz();
  t0     = Genp->T();
  r0     = sqrt(x0*x0+y0*y0);
  gen_id = Genp->GetStatusCode();

  Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
  Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
  Hist->fGenID->Fill(gen_id);
  Hist->fZ0->Fill(z0);
  Hist->fT0->Fill(t0);
  Hist->fR0->Fill(r0);
  Hist->fP->Fill(p);
  Hist->fCosTh->Fill(cos_th);
}

//-----------------------------------------------------------------------------
// void TCalAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

//   Hist->fPdgCode->Fill(Simp->fPdgCode);
//   Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
//   Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
//   Hist->fNStrawHits->Fill(Simp->fNStrawHits);
// }

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TCalAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
//  RegisterDataBlock("TrackBlock"    ,"TStnTrackBlock"   ,&fTrackBlock  );
  RegisterDataBlock("ClusterBlock"  ,"TStnClusterBlock" ,&fClusterBlock);
  RegisterDataBlock("CalDataBlock"  ,"TCalDataBlock"    ,&fCalDataBlock);
  RegisterDataBlock("GenpBlock"     ,"TGenpBlock"       ,&fGenpBlock);
//   RegisterDataBlock("SimpBlock"     ,"TSimpBlock"       ,&fSimpBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
void TCalAnaModule::FillHistograms() {

  static int first_entry(1);
//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0], 0.0,300);
  FillEventHistograms(fHist.fEvent[1], 0.1,300);
  FillEventHistograms(fHist.fEvent[2], 0.5,300);
  FillEventHistograms(fHist.fEvent[3], 1.0,300);
  FillEventHistograms(fHist.fEvent[4], 1.5,300);
  FillEventHistograms(fHist.fEvent[5], 2.0,300);
//-----------------------------------------------------------------------------
// EVT_6: E> 1 MeV, T>400
// EVT_7: E> 1 MeV, T>500
// EVT_8: E> 1 MeV, T>600
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[6], 1.0,400);
  FillEventHistograms(fHist.fEvent[7], 1.0,500);
  FillEventHistograms(fHist.fEvent[8], 1.0,600);
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
//   if (fSimp) {
//     FillSimpHistograms(fHist.fSimp[0],fSimp);
//   }
//-----------------------------------------------------------------------------
// cluster histograms
//-----------------------------------------------------------------------------
  TStnCluster* cl;
  fNCl20 = 0;
  fNCl50 = 0;
  fNCl70 = 0;
  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);

    FillClusterHistograms(fHist.fCluster[0],cl);

    if (cl->Energy() > 20.) {
      FillClusterHistograms(fHist.fCluster[1],cl);
      fNCl70 += 1;
    }

    if (cl->Energy() > 50.) {
      FillClusterHistograms(fHist.fCluster[2],cl);
      fNCl70 += 1;
    }

    if (cl->Energy() > 70.) {
      FillClusterHistograms(fHist.fCluster[3],cl);
      fNCl70 += 1;
    }
  }
//-----------------------------------------------------------------------------
// calorimeter histograms
//-----------------------------------------------------------------------------
  TDisk*         disk;
  TStnCrystal*   cr;

  for (int i=0; i<kNDisks; i++) {
    disk = fDiskCalorimeter->Disk(i);
    for (int ic=0; ic<disk->NCrystals(); ic++) {
      cr = disk->Crystal(ic);
      FillCaloHistograms(fHist.fCalo[0],cr);
      if (cr->Energy() > 0) {
	FillCaloHistograms(fHist.fCalo[1],cr);
      }
      if (cr->Energy() > 0.1) {
	FillCaloHistograms(fHist.fCalo[2],cr);
      }
      if (cr->Energy() > 0.5) {
	FillCaloHistograms(fHist.fCalo[3],cr);
      }
      if (cr->Energy() > 1.0) {
	FillCaloHistograms(fHist.fCalo[4],cr);
      }
    }
  }
//-----------------------------------------------------------------------------
// hit histograms
//-----------------------------------------------------------------------------
  int          nh;
  TCalHitData* hit;

  for (int i=0; i<kNDisks; i++) {
    disk = fDiskCalorimeter->Disk(i);
    for (int ic=0; ic<disk->NCrystals(); ic++) {
      cr = disk->Crystal(ic);
      nh = cr->NHits();
      for (int ih=0; ih<nh; ih++) {
	hit = cr->CalHitData(ih);
	FillCalHitHistograms(fHist.fCalHit[0],hit);
      }
    }
  }
//-----------------------------------------------------------------------------
// radial distributions for crystals
//-----------------------------------------------------------------------------
  if (first_entry == 1) {
    int nd = fDiskCalorimeter->NDisks();
	
    for (int idisk=0; idisk<nd; idisk++) {
      disk = fDiskCalorimeter->Disk(idisk);
      for (int ic=0; ic<disk->NCrystals(); ic++) {
	cr = disk->Crystal(ic);
	fHist.fCrystalR[idisk]->Fill(cr->Radius());
      }
    }
  }
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }

  first_entry = 0;
}



//_____________________________________________________________________________
int TCalAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TCalAnaModule::Event(int ientry) {

  double                p;
  //  TEmuLogLH::CalData_t  dat;
  //  TStnTrack*            track;
  //  int                   id_word;
  TLorentzVector        mom;

  TDiskCalorimeter::GeomData_t disk_geom;

  fClusterBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  //  fSimpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp    = fGenpBlock->NParticles();
  fElectron = fGenpBlock->Particle(0);
  //  fSimp     = fSimpBlock->Particle(0);

  fElectron->Momentum(mom);
					// this is a kludge, to be removed at the next 
					// ntupling 
  //  fEleE     = fElectron->Energy();
  p         = mom.P();
  fEleE     = sqrt(p*p+0.511*0.511);


  if (fDiskCalorimeter->Initialized() == 0) {
    disk_geom.fNDisks = fCalDataBlock->NDisks();

    for (int i=0; i<disk_geom.fNDisks; i++) {
      //      disk_geom.fNCrystals[i] = fCalDataBlock->fNCrystals[i];
      disk_geom.fRMin[i]      = fCalDataBlock->fRMin[i];
      disk_geom.fRMax[i]      = fCalDataBlock->fRMax[i];
      disk_geom.fZ0  [i]      = fCalDataBlock->fZ0  [i];
    }

    disk_geom.fHexSize          = fCalDataBlock->CrystalSize()*2;
    // kludge , so far
    disk_geom.fMinFraction      = 1.; // fCalDataBlock->MinFraction();
    disk_geom.fWrapperThickness = fCalDataBlock->WrapperThickness();
    disk_geom.fShellThickness   = fCalDataBlock->ShellThickness();

    fDiskCalorimeter->Init(&disk_geom);
  }

  TStnCluster* cl;
  fNClusters  = fClusterBlock->NClusters();

  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);

  fNCl20      = 0;
  fNCl50      = 0;
  fNCl70      = 0;
  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);
    if (cl->Energy() > 20.) fNCl20 += 1;
    if (cl->Energy() > 50.) fNCl50 += 1;
    if (cl->Energy() > 70.) fNCl70 += 1;
  }

  fNCalHits   = fCalDataBlock->NHits();

  fDiskCalorimeter->InitEvent(fCalDataBlock);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCalAnaModule::Debug() {

// //-----------------------------------------------------------------------------
// // bit 5: events with N(tracks) > 1
// //-----------------------------------------------------------------------------
//   if (GetDebugBit(5) == 1) {
//     int ntrk = fTrackBlock->NTracks();
//     if (ntrk > 1) {
//       GetHeaderBlock()->Print(Form("NTracks = %i5",ntrk));
//     }
//   }
}

//_____________________________________________________________________________
int TCalAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TCalAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}


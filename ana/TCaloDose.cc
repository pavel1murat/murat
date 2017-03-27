#define TCaloDose_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TRandom3.h"

#include "murat/ana/TCaloDose.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/geom/TStnCrystal.hh"

ClassImp(TCaloDose)
// //-----------------------------------------------------------------------------
// // z range for different stations
// //-----------------------------------------------------------------------------
// double zplane[36][2] = {
//    8660., 8720.,         // 0
//    8720., 8780.,         // 0
//    8820., 8880.,         // 1
//    8880., 8920.,         // 1
//    8980., 9030.,         // 2
//    9030., 9080.,         // 2
//    9140., 9180.,         // 3
//    9180., 9230.,         // 3
//    9300., 9340.,         // 4
//    9340., 9390.,         // 4
//    9450., 9500.,         // 5
//    9500., 9550.,         // 5

//    9760., 9800.,         // 6
//    9800., 9860.,         // 6
//    9920., 9960.,         // 7
//    9960.,  10010.,       // 7
//    10080., 10120.,       // 8
//    10120., 10180.,       // 8
//    10220., 10280.,       // 9
//    10280., 10320.,       // 9
//    10380., 10440.,       // 10
//    10440., 10480.,       // 10
//    10540., 10580.,       // 11
//    10580., 10640.,       // 11

//    10850., 10900.,       // 12
//    10950., 10950.,       // 12
//    11010., 11050.,       // 13
//    11050., 11100.,       // 13
//    11160., 11210.,       // 14
//    11210., 11260.,       // 14
//    11320., 11360.,       // 15
//    11360., 11420.,       // 15
//    11480., 11520.,       // 16
//    11520., 11580.,       // 16
//    11620., 11680.,       // 17
//    11680., 11720.        // 17
// };
  

// //-----------------------------------------------------------------------------
// int plane_number(float Z) {
//   int ist(-1);

//   for (int i=0; i<36; i++) {
//     if ((Z > zplane[i][0]) && (Z < zplane[i][1])) {
//       ist = i;
//       break;
//     }
//   }

//   return ist;
// }

namespace {
};

//-----------------------------------------------------------------------------
TCaloDose::TCaloDose(const char* Process) : TStnModule("CaloDose","CaloDose") {

  fProcess = Process;
  
  TDiskCalorimeter::GeomData_t data;

  data.fNDisks           = 2;
  data.fNEdges           = 4;
  data.fHexSize          = 34. ; // mm
  data.fMinFraction      = 1.;
  data.fWrapperThickness = 0.150; // 150 um
  data.fShellThickness   = 0;
  
  data.fNCrystals[0]     = 678;  // real number + 4
  data.fRMin     [0]     = 374.;
  data.fRMax     [0]     = 660.;
  data.fZ0       [0]     = 11853.15+100.;

  data.fNCrystals[1]     = 678;  // real number + 4
  data.fRMin     [1]     = 374.;
  data.fRMax     [1]     = 660.;
  data.fZ0       [1]     = 12553.15+100.;

  fCal                   = new TDiskCalorimeter();
  fCal->Init(&data);

  fNPOT = 1.2e20;           // "per year of running"

  TH1::AddDirectory(0);

  InitChain();

  BookHistograms();

  fZMin[0] = 11853.15;
  fZMin[1] = 12553.15;
  fZMax[0] = fZMin[0]+200.;
  fZMax[1] = fZMin[1]+200.;
//-----------------------------------------------------------------------------
// resert energies
//-----------------------------------------------------------------------------
  for (int id=0; id<2; id++) {
    for (int ic=0; ic<678; ic++) {
      fE[id][ic] = 0;
      for (int iz=0; iz<20; iz++) {
	fdEdZ[id][ic][iz] = 0;
      }
    }
  }

  fEDisk[0] = 0;
  fEDisk[1] = 0;
}

//-----------------------------------------------------------------------------
TCaloDose::~TCaloDose() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();

   delete fCal;
}


//-----------------------------------------------------------------------------
int TCaloDose::BookEventHistograms(EventHist_t* Hist, const char* Folder) {

  HBook1F(Hist->fNCrystals ,"ncrystals" ,Form("%s: N(crystals) ",Folder), 2000,   0  ,2000,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int TCaloDose::BookDiskHistograms(DiskHist_t* Hist, const char* Folder) {

  HBook1F(Hist->fZ      ,"dz"     ,Form("%s: DZ"        ,Folder),    20,       0,  200,Folder);
  //  HBook1F(Hist->fEVsZ   ,"e_vs_dz",Form("%s: Edep vs DZ",Folder),    20,       0,  200,Folder);

  HBook1F(Hist->fR      ,"r"      ,Form("%s: R "        ,Folder),   100,     370,  670,Folder);

  HBook1F(Hist->fNVsR   ,"n_vs_r" ,Form("%s: N    vs R ",Folder),    10,     370,  670,Folder);

  HBook2F(Hist->fDzVsIc ,"dz_vs_ic",Form("%s: Dz vs IC ",Folder),   700,      0,   700,20,0,200,Folder);

  HBook2F(Hist->fRVsDz  ,"r_vs_dz", Form("%s: R vs Dz ",Folder)  ,    20,0,200, 10,370,670,Folder);


  for (int ir=0; ir<10; ir++) {
    HBook1F(Hist->fEVsZ[ir],Form("ir_%02i",ir),Form("ir_%02i",ir),20,0,200,Folder);
  }

  for (int iz=0; iz<20; iz++) {
    HBook1F(Hist->fEVsR[iz],Form("iz_%02i",iz),Form("iz_%02i",iz),10,370,670,Folder);
  }

  HBook1F(Hist->fETotVsR,"etot_vs_r","E(tot) vs R",10,370,670,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int TCaloDose::BookHistograms() {

  TFolder*    fol;
  TFolder*    hist_folder;
  char        folder_name[200];
 
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kMaxEventHistSets];
  for (int i=0; i<kMaxEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[  0] = 1;     		// all hits

  for (int i=0; i<kMaxEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// disk histograms
//-----------------------------------------------------------------------------
  int book_disk_histset[kMaxDiskHistSets];
  for (int i=0; i<kMaxDiskHistSets; i++) book_disk_histset[i] = 0;

  book_disk_histset[  0] = 1;     		// first disk
  book_disk_histset[  1] = 1;     		// second disk

  for (int i=0; i<kMaxDiskHistSets; i++) {
    if (book_disk_histset[i] != 0) {
      sprintf(folder_name,"disk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fDisk[i] = new DiskHist_t;
      BookDiskHistograms(fHist.fDisk[i],Form("Hist/%s",folder_name));
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
int TCaloDose::FillEventHistograms(EventHist_t* Hist) {

  Hist->fNCrystals->Fill(nCrystal);
  return 0;
}

//-----------------------------------------------------------------------------
int TCaloDose::FillDiskHistograms(DiskHist_t* Hist) {
  return 0;
}

//-----------------------------------------------------------------------------
int TCaloDose::FillHistograms() {

  FillEventHistograms(fHist.fEvent[0]);

  for (int i=0; i<nCrystal; i++) {
    int idisk = crystalSectionId[i];
    int icr   = crystalId[i];

    TDisk* disk = fCal->Disk(idisk);

    int icc   = icr-disk->FirstChanOffset();  // index within the disk

    float x = crystalPosX[i];
    float y = crystalPosY[i];
    float z = crystalPosZ[i];

    float dz = z-fZMin[idisk];   // ranges from 0 to 200 mm
    //    int   iz = dz/10.;           // in cm
    
    float r = sqrt(x*x+y*y);

    fEDisk[idisk]  +=crystalEdep[i];

    fE[idisk][icc] +=crystalEdep[i];

    disk->Crystal(icc)->AddEnergy(crystalEdep[i]);

    fHist.fDisk[idisk]->fDzVsIc->Fill(icc,dz,crystalEdep[i]);
    fHist.fDisk[idisk]->fRVsDz->Fill(dz  ,r ,crystalEdep[i]);
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TCaloDose::ResetHistograms() {

  return 0;
}

//-----------------------------------------------------------------------------
//   In a ROOT session, you can do:
//      root> .L Calo.C
//      root> Calo t
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
// by  b_branchname->GetEntry(ientry); //read only this branch
//-----------------------------------------------------------------------------
void TCaloDose::Loop(Long64_t NEvents) {

  ResetHistograms();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  Long64_t nent = NEvents;
  if (nent < 0) nent = nentries;
  
  Long64_t nbytes = 0, nb = 0;

  long int nev = 0;

  for (Long64_t jentry=0; jentry<nent;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //    InitEvent();
    
    FillHistograms();

    nev += 1;
  }

  printf("processed %li events, e[0] = %12.3f e[1] = %12.3f\n",nev,fEDisk[0],fEDisk[1]);

  // printf(" fZMin[0], fZMax[0], fZMin[1], fZMax[1] = %10.4f %10.4f %10.4f %10.4f\n",
  // 	 fZMin[0], fZMax[0], fZMin[1], fZMax[1]);

  fNEvents = nev;
//-----------------------------------------------------------------------------
// remaining part : normalize histograms to the total N(POT)
// last term, 1.e6, accounts for energies being measured in MeV
//-----------------------------------------------------------------------------
  // float mev_per_joule = 1.6e-19*1.e6;
  // float sf            = fNPOT*fNPerPOT/(fNEvents+1.e-12)*mev_per_joule;

  int  nbr = fHist.fDisk[0]->fNVsR->GetNbinsX();
  float dr = fHist.fDisk[0]->fNVsR->GetXaxis()->GetBinWidth(1);

  for (int i=0; i<2; i++) {
    for (int ic=0; ic<678; ic++) {

      TStnCrystal* cr = fCal->Disk(i)->Crystal(ic);
      //      float x = cr->Center()->X();
      //      float y = cr->Center()->Y();
      //      float r = sqrt(x*x+y*y);

      for (int ib=1; ib<=nbr; ib++) {
	float r        = fHist.fDisk[i]->fNVsR->GetBinCenter(ib);
	float rmin     = r-dr/2;
	float rmax     = r+dr/2;
	float fraction = insideFraction(cr,rmin,rmax);

	fHist.fDisk[i]->fNVsR->Fill(r,fraction);
      }

    }
//-----------------------------------------------------------------------------
// make projections normalized to dose
//-----------------------------------------------------------------------------
    float density        = 4.5;   // CsI
    float mev_per_joule  = 1.6e-19*1.e6;
    float krad_per_gray  = 10.;
    float crystal_area   = 3.4*3.4; //  cm^2
    float crystal_length = 20.; // cm

    char name[100];

    int nz = fHist.fDisk[i]->fRVsDz->GetNbinsX();
    int nr = fHist.fDisk[i]->fRVsDz->GetNbinsY();

    //    float dr =  fHist.fDisk[i]->fRVsDz->GetYaxis()->GetBinWidth(1)/10.; //  // convert to cm
    float dz =  fHist.fDisk[i]->fRVsDz->GetXaxis()->GetBinWidth(1)/10.; //  // convert to cm

    for (int iz=1; iz<=nz; iz++) {
      sprintf(name,"d%i_iz_%02i_e",i,iz);
      TH1D* hp = fHist.fDisk[i]->fRVsDz->ProjectionY(name,iz,iz);

      for (int ir=1; ir<=nr; ir++) {  // do not include underflows and overflows

	//	float r     = fHist.fDisk[i]->fRVsDz->GetYaxis()->GetBinCenter(ir)/10.; // convert to cm
	float ncr     = fHist.fDisk[i]->fNVsR->GetBinContent(ir);
	float mass  = ncr*crystal_area*dz*density/1.e3; // in kG
	float sf(0);
	if (mass > 0) sf = (fNPOT*fNPerPOT)/(nent+1.e-12)*(nentries/(fNSimulated+1.e-12))*mev_per_joule/mass/krad_per_gray;
      
	float x = hp->GetBinContent(ir);
	float e = hp->GetBinError  (ir);
            
	fHist.fDisk[i]->fEVsR[iz-1]->SetBinContent(ir,x*sf);
	fHist.fDisk[i]->fEVsR[iz-1]->SetBinError  (ir,e*sf);
      }
    }

    for (int ir=1; ir<=nr; ir++) {
      sprintf(name,"d%i_ir_%02i_e",i,ir);
      TH1D* hp = fHist.fDisk[i]->fRVsDz->ProjectionX(name,ir,ir);
      float ncr  = fHist.fDisk[i]->fNVsR->GetBinContent(ir);
      for (int iz=1; iz<=nz; iz++) {
	//	float r     = fHist.fDisk[i]->fRVsDz->GetYaxis()->GetBinCenter(ir)/10.; // convert to cm
	//	float mass  = 2*M_PI*r*dr*dz*density/1.e3; // in kG
	float mass  = ncr*crystal_area*dz*density/1.e3; // in kG
	float sf(0);
	if (mass > 0) sf = (fNPOT*fNPerPOT)/(nent+1.e-12)*(nentries/(fNSimulated+1.e-12))*mev_per_joule/mass/krad_per_gray;
      
	float x = hp->GetBinContent(iz);
	float e = hp->GetBinError  (iz);
            
	fHist.fDisk[i]->fEVsZ[ir-1]->SetBinContent(iz,x*sf);
	fHist.fDisk[i]->fEVsZ[ir-1]->SetBinError  (iz,e*sf);
      }
    }
//-----------------------------------------------------------------------------
// dose vs R, integrated and averaged over the crystal length
//-----------------------------------------------------------------------------
    sprintf(name,"etot_vs_r_proj");
    TH1D* hp = fHist.fDisk[i]->fRVsDz->ProjectionY(name);
    for (int ir=1; ir<=nr; ir++) {
      float ncr  = fHist.fDisk[i]->fNVsR->GetBinContent(ir);
      float mass  = ncr*crystal_area*crystal_length*density/1.e3; // in kG
      float sf(0);
      if (mass > 0) sf = (fNPOT*fNPerPOT)/(nent+1.e-12)*(nentries/(fNSimulated+1.e-12))*mev_per_joule/mass/krad_per_gray;
      
      float x = hp->GetBinContent(ir);
      float e = hp->GetBinError  (ir);
            
      fHist.fDisk[i]->fETotVsR->SetBinContent(ir,x*sf);
      fHist.fDisk[i]->fETotVsR->SetBinError  (ir,e*sf);
    }
  }

}

//-----------------------------------------------------------------------------
Int_t TCaloDose::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


//-----------------------------------------------------------------------------
Long64_t TCaloDose::LoadTree(Long64_t entry) {
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
void TCaloDose::Init(TChain* Chain) {
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!Chain) return;
   fChain = Chain;

   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("caloCrystals", &caloCrystals, &b_caloCrystals);
   fChain->SetBranchAddress("caloDisk0Crystals", &caloDisk0Crystals, &b_caloDisk0Crystals);
   fChain->SetBranchAddress("caloDisk1Crystals", &caloDisk1Crystals, &b_caloDisk1Crystals);
   fChain->SetBranchAddress("caloVolume", &caloVolume, &b_caloVolume);
   fChain->SetBranchAddress("crystalVolume", &crystalVolume, &b_crystalVolume);
   fChain->SetBranchAddress("nGen", &nGen, &b_nGen);
   fChain->SetBranchAddress("genId", genId, &b_genId);
   fChain->SetBranchAddress("genCrCode", genCrCode, &b_genCrCode);
   fChain->SetBranchAddress("genMomX", genMomX, &b_genMomX);
   fChain->SetBranchAddress("genMomY", genMomY, &b_genMomY);
   fChain->SetBranchAddress("genMomZ", genMomZ, &b_genMomZ);
   fChain->SetBranchAddress("genStartX", genStartX, &b_genStartX);
   fChain->SetBranchAddress("genStartY", genStartY, &b_genStartY);
   fChain->SetBranchAddress("genStartZ", genStartZ, &b_genStartZ);
   fChain->SetBranchAddress("genStartT", genStartT, &b_genStartT);
   fChain->SetBranchAddress("nCrystal", &nCrystal, &b_nCrystal);
   fChain->SetBranchAddress("crystalId", &crystalId, &b_crystalId);
   fChain->SetBranchAddress("crystalSectionId", &crystalSectionId, &b_crystalSectionId);
   fChain->SetBranchAddress("crystalPosX", &crystalPosX, &b_crystalPosX);
   fChain->SetBranchAddress("crystalPosY", &crystalPosY, &b_crystalPosY);
   fChain->SetBranchAddress("crystalPosZ", &crystalPosZ, &b_crystalPosZ);
   fChain->SetBranchAddress("crystalEdep", &crystalEdep, &b_crystalEdep);
   fChain->SetBranchAddress("crystalDose", &crystalDose, &b_crystalDose);
   fChain->SetBranchAddress("crystalDose0", &crystalDose0, &b_crystalDose0);
   fChain->SetBranchAddress("crystalDose1", &crystalDose1, &b_crystalDose1);
   fChain->SetBranchAddress("crystalDose2", &crystalDose2, &b_crystalDose2);
   fChain->SetBranchAddress("crystalDose3", &crystalDose3, &b_crystalDose3);
   fChain->SetBranchAddress("crystalDose4", &crystalDose4, &b_crystalDose4);
   fChain->SetBranchAddress("crystalDose5", &crystalDose5, &b_crystalDose5);
   fChain->SetBranchAddress("crystalDose6", &crystalDose6, &b_crystalDose6);
   fChain->SetBranchAddress("crystalDose7", &crystalDose7, &b_crystalDose7);
   fChain->SetBranchAddress("crystalDose8", &crystalDose8, &b_crystalDose8);
   fChain->SetBranchAddress("crystalDose9", &crystalDose9, &b_crystalDose9);
   fChain->SetBranchAddress("crystalDose10", &crystalDose10, &b_crystalDose10);
   fChain->SetBranchAddress("crystalDose11", &crystalDose11, &b_crystalDose11);
   fChain->SetBranchAddress("crystalDose12", &crystalDose12, &b_crystalDose12);
   fChain->SetBranchAddress("crystalDose13", &crystalDose13, &b_crystalDose13);
   fChain->SetBranchAddress("crystalDose14", &crystalDose14, &b_crystalDose14);
   fChain->SetBranchAddress("crystalDose15", &crystalDose15, &b_crystalDose15);
   fChain->SetBranchAddress("crystalDose16", &crystalDose16, &b_crystalDose16);
   fChain->SetBranchAddress("crystalDose17", &crystalDose17, &b_crystalDose17);
   fChain->SetBranchAddress("crystalDose18", &crystalDose18, &b_crystalDose18);
   fChain->SetBranchAddress("crystalDose19", &crystalDose19, &b_crystalDose19);
   fChain->SetBranchAddress("nCrystalRO", &nCrystalRO, &b_nCrystalRO);
   fChain->SetBranchAddress("crystalROSectionId", crystalROSectionId, &b_crystalROSectionId);
   fChain->SetBranchAddress("crystalROCrystalId", crystalROCrystalId, &b_crystalROCrystalId);
   fChain->SetBranchAddress("crystalROEdep", crystalROEdep, &b_crystalROEdep);
   fChain->SetBranchAddress("crystalRODose", crystalRODose, &b_crystalRODose);
   fChain->SetBranchAddress("crystalROX", crystalROX, &b_crystalROX);
   fChain->SetBranchAddress("crystalROY", crystalROY, &b_crystalROY);
   fChain->SetBranchAddress("crystalROZ", crystalROZ, &b_crystalROZ);
   fChain->SetBranchAddress("crystalROR", crystalROR, &b_crystalROR);
   fChain->SetBranchAddress("nCrystalROCard", &nCrystalROCard, &b_nCrystalROCard);
   fChain->SetBranchAddress("crystalROCardSectionId", crystalROCardSectionId, &b_crystalROCardSectionId);
   fChain->SetBranchAddress("crystalROCardCrystalId", crystalROCardCrystalId, &b_crystalROCardCrystalId);
   fChain->SetBranchAddress("crystalROCardEdep", crystalROCardEdep, &b_crystalROCardEdep);
   fChain->SetBranchAddress("crystalROCardDose", crystalROCardDose, &b_crystalROCardDose);
   fChain->SetBranchAddress("crystalROCardX", crystalROCardX, &b_crystalROCardX);
   fChain->SetBranchAddress("crystalROCardY", crystalROCardY, &b_crystalROCardY);
   fChain->SetBranchAddress("crystalROCardZ", crystalROCardZ, &b_crystalROCardZ);
   fChain->SetBranchAddress("crystalROCardR", crystalROCardR, &b_crystalROCardR);
   fChain->SetBranchAddress("nCrateHits", &nCrateHits, &b_nCrateHits);
   fChain->SetBranchAddress("crateEdep", crateEdep, &b_crateEdep);
   //   fChain->SetBranchAddress("crateDose", crateDose, &b_crateDose);
   fChain->SetBranchAddress("crateX", crateX, &b_crateX);
   fChain->SetBranchAddress("crateY", crateY, &b_crateY);
   fChain->SetBranchAddress("crateZ", crateZ, &b_crateZ);
   fChain->SetBranchAddress("crateR", crateR, &b_crateR);
   fChain->SetBranchAddress("cratePdgId", cratePdgId, &b_cratePdgId);
   fChain->SetBranchAddress("crateE", crateE, &b_crateE);
   fChain->SetBranchAddress("crateEkin", crateEkin, &b_crateEkin);
   fChain->SetBranchAddress("crateMass", crateMass, &b_crateMass);
   fChain->SetBranchAddress("crateL", crateL, &b_crateL);
   fChain->SetBranchAddress("vNHits", &vNHits, &b_vNHits);
   fChain->SetBranchAddress("vId", vId, &b_vId);
   fChain->SetBranchAddress("vPdgId", vPdgId, &b_vPdgId);
   fChain->SetBranchAddress("vP", vP, &b_vP);
   fChain->SetBranchAddress("vPx", vPx, &b_vPx);
   fChain->SetBranchAddress("vPy", vPy, &b_vPy);
   fChain->SetBranchAddress("vPz", vPz, &b_vPz);
   fChain->SetBranchAddress("vE", vE, &b_vE);
   fChain->SetBranchAddress("vEKin", vEKin, &b_vEKin);
   fChain->SetBranchAddress("vM", vM, &b_vM);
   fChain->SetBranchAddress("vT", vT, &b_vT);
   fChain->SetBranchAddress("vX", vX, &b_vX);
   fChain->SetBranchAddress("vY", vY, &b_vY);
   fChain->SetBranchAddress("vZ", vZ, &b_vZ);
   fChain->SetBranchAddress("vCosth", vCosth, &b_vCosth);
   fChain->SetBranchAddress("vRadius", vRadius, &b_vRadius);

   Notify();
}

//-----------------------------------------------------------------------------
Bool_t TCaloDose::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

//-----------------------------------------------------------------------------
void TCaloDose::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

//-----------------------------------------------------------------------------
// fNPetPOT for non-flash backgrounds is a product of the number of stopped 
// muons per POT and average yiled of a given particle species per muon stop
//-----------------------------------------------------------------------------
int TCaloDose::InitChain() {

  TChain* chain = new TChain("//calorimeterDose/Calo","Calo");

  if (fProcess == "FLASH") {
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_0.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_1.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_2.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_3.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_4.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_5.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_6.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_7.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_8.root");
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloFLASH_9.root");

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9;
  }
  else if (fProcess == "DIO") {
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloDIO.root");
    fNPerPOT    = 7.27e-4;
    fNSimulated = 1.0e7;
  }
  else if (fProcess == "OOT") {
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloOOT.root");
    fNPerPOT    = 3.97e-3;
    fNSimulated = 3.0e7;
  }
  else if (fProcess == "PHOTON") {
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloPHOTON.root");
    fNPerPOT    = 2.28e-3;
    fNSimulated = 1.0e8;
  }
  else if (fProcess == "PROTON") {
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloPROTON.root");
    fNPerPOT    = 5.69e-5;
    fNSimulated = 1.0e8;
  }
  else if (fProcess == "DEUTERON") {
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloDEUTERON.root");
    fNPerPOT    = 2.84e-5;
    fNSimulated = 1.0e8;
  }
  else if (fProcess == "NEUTRON") {
    chain->Add("/mu2e/data/users/gianipez/hist/treeCaloNEUTRON.root");
    fNPerPOT    = 1.37e-3;
    fNSimulated = 1.0e8;
  }

  fChain = chain;

  Init(fChain);
//-----------------------------------------------------------------------------
// leave only tracker branches active
//-----------------------------------------------------------------------------
  fChain->SetBranchStatus("*",1);

  // fChain->SetBranchStatus("nTtsUpHits",1);
  // fChain->SetBranchStatus("ttsUpEdep" ,1);
  // fChain->SetBranchStatus("ttsUpX"    ,1);
  // fChain->SetBranchStatus("ttsUpY"    ,1);
  // fChain->SetBranchStatus("ttsUpZ"    ,1);
  // fChain->SetBranchStatus("ttsUpR"    ,1);
  // fChain->SetBranchStatus("ttsUpPdgId",1);
  // fChain->SetBranchStatus("ttsUpEkin" ,1);
  // fChain->SetBranchStatus("ttsUpMass" ,1);
  // fChain->SetBranchStatus("ttsUpL"    ,1);
  
  // fChain->SetBranchStatus("nTtsDwHits",1);
  // fChain->SetBranchStatus("ttsDwEdep" ,1);
  // fChain->SetBranchStatus("ttsDwX"    ,1);
  // fChain->SetBranchStatus("ttsDwY"    ,1);
  // fChain->SetBranchStatus("ttsDwZ"    ,1);
  // fChain->SetBranchStatus("ttsDwR"    ,1);
  // fChain->SetBranchStatus("ttsDwPdgId",1);
  // fChain->SetBranchStatus("ttsDwEkin" ,1);
  // fChain->SetBranchStatus("ttsDwMass" ,1);
  // fChain->SetBranchStatus("ttsDwL"    ,1);
  
  return 0;
}

//-----------------------------------------------------------------------------
// brute force: MC with 10K throws per crystal
//-----------------------------------------------------------------------------
float TCaloDose::insideFraction(TStnCrystal* Cr, float RMin, float RMax) {

  float x, y, x0, y0, r, size, fraction(0.);
  int   nin(0), nmax(10000);
  TRandom3 rn;

  int sx[4] = { -1, 1,  1, -1};
  int sy[4] = {  1, 1, -1, -1};

  x0   = Cr->Center()->X();
  y0   = Cr->Center()->Y();
  size = Cr->Size();

  int inside(0);

  for (int i=0; i<4; i++) {
    x = x0+sx[i]*size/2;
    y = y0+sy[i]*size/2;
    r = sqrt(x*x+y*y);

    if ((r > RMin) && (r < RMax)) inside++;
  }


  if (inside > 0) {
    
    for (int i=0; i<nmax; i++) {
      x = x0+size/2*(2*rn.Rndm(i)-1);
      y = y0+size/2*(2*rn.Rndm(i)-1);
      r = sqrt(x*x+y*y);
      if ((r > RMin) && (r < RMax)) nin++;
    }

    fraction = float(nin)/float(nmax);
  }

  return fraction;

}

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 17 16:24:24 2017 by ROOT version 6.08/04
// from TTree TTrackerDose/TTrackerDose
// found on file: /mu2e/data/users/gianipez/hist/treeTrackerFLASH_0.root
//
// use inheritance from TStnModule to book histograms
//////////////////////////////////////////////////////////

#ifndef TTrackerDose_h
#define TTrackerDose_h

#include "TNamed.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "Stntuple/loop/TStnModule.hh"

// Header file for the classes stored in the TTree if any.

class TTrackerDose : public TStnModule {
public :
					// stored in ntuple , but in an inconvenient way...
  struct StepData_t {
    float fEDep;
    float fXt;
    float fYt;
    float fZt;
    float fR;
    float fPdgID;
    float fMass;
    float fE;
    float fEKin;
    float fStep;
  };

  enum { kMaxTrackHistSets = 200 };
  enum { kMaxEventHistSets =  10 };
  enum { kMaxDiskHistSets  =  10 };

  struct EventHist_t {
    TH1F*     fNumber;
  };

  struct TrackHist_t {
    TH1F*     fNHits;
    TH1F*     fPdgID[2];
    TH1F*     fR;
    TH1F*     fZ;
    TH1F*     fEKin[5];
    TH2F*     fEDepVsPlane[2];
  };

  struct DiskHist_t {
    TH1F*     fNHits;
    TH1F*     fPdgID[2];
    TH1F*     fR;
    TH1F*     fZ;
    TH1F*     fEKin[5];
    TH2F*     fEDepVsZ[2];
  };

  TChain*     fChain;   //! pointer to the analyzed TTree or TChain
  Int_t       fCurrent; //! current Tree number in a TChain

  struct Hist_t {
    EventHist_t*    fEvt [kMaxEventHistSets];
    TrackHist_t*    fUp  [kMaxTrackHistSets];
    TrackHist_t*    fDn  [kMaxTrackHistSets];
    DiskHist_t*     fDisk[kMaxDiskHistSets ];
  };

  Hist_t    fHist;

  float     fNPerPOT;   // efficiency for a given bgr source
  float     fNPOT   ;   // number of protons on target (expected)
  long int  fNEvents;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           evt;
   Int_t           run;
   Int_t           caloCrystals;
   Int_t           caloDisk0Crystals;
   Int_t           caloDisk1Crystals;
   Float_t         caloVolume;
   Float_t         crystalVolume;
   Int_t           nGen;
   Int_t           genId[1];   //[nGen]
   Int_t           genCrCode[1];   //[nGen]
   Float_t         genMomX[1];   //[nGen]
   Float_t         genMomY[1];   //[nGen]
   Float_t         genMomZ[1];   //[nGen]
   Float_t         genStartX[1];   //[nGen]
   Float_t         genStartY[1];   //[nGen]
   Float_t         genStartZ[1];   //[nGen]
   Float_t         genStartT[1];   //[nGen]
   Int_t           nCrystal;
   Int_t           crystalId[1];   //[nCrystal]
   Int_t           crystalSectionId[1];   //[nCrystal]
   Float_t         crystalPosX[1];   //[nCrystal]
   Float_t         crystalPosY[1];   //[nCrystal]
   Float_t         crystalPosZ[1];   //[nCrystal]
   Float_t         crystalR[1];   //[nCrystal]
   Float_t         crystalEdep[1];   //[nCrystal]
   Float_t         crystalDose[1];   //[nCrystal]
   Float_t         crystalDose0[1];   //[nCrystal]
   Float_t         crystalDose1[1];   //[nCrystal]
   Float_t         crystalDose2[1];   //[nCrystal]
   Float_t         crystalDose3[1];   //[nCrystal]
   Float_t         crystalDose4[1];   //[nCrystal]
   Float_t         crystalDose5[1];   //[nCrystal]
   Float_t         crystalDose6[1];   //[nCrystal]
   Float_t         crystalDose7[1];   //[nCrystal]
   Float_t         crystalDose8[1];   //[nCrystal]
   Float_t         crystalDose9[1];   //[nCrystal]
   Float_t         crystalDose10[1];   //[nCrystal]
   Float_t         crystalDose11[1];   //[nCrystal]
   Float_t         crystalDose12[1];   //[nCrystal]
   Float_t         crystalDose13[1];   //[nCrystal]
   Float_t         crystalDose14[1];   //[nCrystal]
   Float_t         crystalDose15[1];   //[nCrystal]
   Float_t         crystalDose16[1];   //[nCrystal]
   Float_t         crystalDose17[1];   //[nCrystal]
   Float_t         crystalDose18[1];   //[nCrystal]
   Float_t         crystalDose19[1];   //[nCrystal]
   Int_t           nCrystalRO;
   Int_t           crystalROSectionId[1];   //[nCrystalRO]
   Int_t           crystalROCrystalId[1];   //[nCrystalRO]
   Float_t         crystalROEdep[1];   //[nCrystalRO]
   Float_t         crystalRODose[1];   //[nCrystalRO]
   Float_t         crystalROX[1];   //[nCrystalRO]
   Float_t         crystalROY[1];   //[nCrystalRO]
   Float_t         crystalROZ[1];   //[nCrystalRO]
   Float_t         crystalROR[1];   //[nCrystalRO]
   Int_t           nCrystalROCard;
   Int_t           crystalROCardSectionId[1];   //[nCrystalROCard]
   Int_t           crystalROCardCrystalId[1];   //[nCrystalROCard]
   Float_t         crystalROCardEdep[1];   //[nCrystalROCard]
   Float_t         crystalROCardDose[1];   //[nCrystalROCard]
   Float_t         crystalROCardX[1];   //[nCrystalROCard]
   Float_t         crystalROCardY[1];   //[nCrystalROCard]
   Float_t         crystalROCardZ[1];   //[nCrystalROCard]
   Float_t         crystalROCardR[1];   //[nCrystalROCard]
   Int_t           nCrateHits;
   Float_t         crateEdep[1];   //[nCrateHits]
   Float_t         crateX[1];   //[nCrateHits]
   Float_t         crateY[1];   //[nCrateHits]
   Float_t         crateZ[1];   //[nCrateHits]
   Float_t         crateR[1];   //[nCrateHits]
   Int_t           cratePdgId[1];   //[nCrateHits]
   Float_t         crateE[1];   //[nCrateHits]
   Float_t         crateEkin[1];   //[nCrateHits]
   Float_t         crateMass[1];   //[nCrateHits]
   Float_t         crateL[1];   //[nCrateHits]
   Int_t           vNHits;
   Int_t           vId[82];   //[vNHits]
   Int_t           vPdgId[82];   //[vNHits]
   Float_t         vP[82];   //[vNHits]
   Float_t         vPx[82];   //[vNHits]
   Float_t         vPy[82];   //[vNHits]
   Float_t         vPz[82];   //[vNHits]
   Float_t         vE[82];   //[vNHits]
   Float_t         vEKin[82];   //[vNHits]
   Float_t         vM[82];   //[vNHits]
   Float_t         vT[82];   //[vNHits]
   Float_t         vX[82];   //[vNHits]
   Float_t         vY[82];   //[vNHits]
   Float_t         vZ[82];   //[vNHits]
   Float_t         vCosth[82];   //[vNHits]
   Float_t         vRadius[82];   //[vNHits]
   Int_t           nTtsUpHits;
   Float_t         ttsUpEdep[167];   //[nTtsUpHits]
   Float_t         ttsUpX[167];   //[nTtsUpHits]
   Float_t         ttsUpY[167];   //[nTtsUpHits]
   Float_t         ttsUpZ[167];   //[nTtsUpHits]
   Float_t         ttsUpR[167];   //[nTtsUpHits]
   Int_t           ttsUpPdgId[167];   //[nTtsUpHits]
   Float_t         ttsUpEkin[167];   //[nTtsUpHits]
   Float_t         ttsUpMass[167];   //[nTtsUpHits]
   Float_t         ttsUpL[167];   //[nTtsUpHits]
   Int_t           nTtsDwHits;
   Float_t         ttsDwEdep[186];   //[nTtsDwHits]
   Float_t         ttsDwX[186];   //[nTtsDwHits]
   Float_t         ttsDwY[186];   //[nTtsDwHits]
   Float_t         ttsDwZ[186];   //[nTtsDwHits]
   Float_t         ttsDwR[186];   //[nTtsDwHits]
   Int_t           ttsDwPdgId[186];   //[nTtsDwHits]
   Float_t         ttsDwEkin[186];   //[nTtsDwHits]
   Float_t         ttsDwMass[186];   //[nTtsDwHits]
   Float_t         ttsDwL[186];   //[nTtsDwHits]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_run;   //!
   TBranch        *b_caloCrystals;   //!
   TBranch        *b_caloDisk0Crystals;   //!
   TBranch        *b_caloDisk1Crystals;   //!
   TBranch        *b_caloVolume;   //!
   TBranch        *b_crystalVolume;   //!
   TBranch        *b_nGen;   //!
   TBranch        *b_genId;   //!
   TBranch        *b_genCrCode;   //!
   TBranch        *b_genMomX;   //!
   TBranch        *b_genMomY;   //!
   TBranch        *b_genMomZ;   //!
   TBranch        *b_genStartX;   //!
   TBranch        *b_genStartY;   //!
   TBranch        *b_genStartZ;   //!
   TBranch        *b_genStartT;   //!
   TBranch        *b_nCrystal;   //!
   TBranch        *b_crystalId;   //!
   TBranch        *b_crystalSectionId;   //!
   TBranch        *b_crystalPosX;   //!
   TBranch        *b_crystalPosY;   //!
   TBranch        *b_crystalPosZ;   //!
   TBranch        *b_crystalR;   //!
   TBranch        *b_crystalEdep;   //!
   TBranch        *b_crystalDose;   //!
   TBranch        *b_crystalDose0;   //!
   TBranch        *b_crystalDose1;   //!
   TBranch        *b_crystalDose2;   //!
   TBranch        *b_crystalDose3;   //!
   TBranch        *b_crystalDose4;   //!
   TBranch        *b_crystalDose5;   //!
   TBranch        *b_crystalDose6;   //!
   TBranch        *b_crystalDose7;   //!
   TBranch        *b_crystalDose8;   //!
   TBranch        *b_crystalDose9;   //!
   TBranch        *b_crystalDose10;   //!
   TBranch        *b_crystalDose11;   //!
   TBranch        *b_crystalDose12;   //!
   TBranch        *b_crystalDose13;   //!
   TBranch        *b_crystalDose14;   //!
   TBranch        *b_crystalDose15;   //!
   TBranch        *b_crystalDose16;   //!
   TBranch        *b_crystalDose17;   //!
   TBranch        *b_crystalDose18;   //!
   TBranch        *b_crystalDose19;   //!
   TBranch        *b_nCrystalRO;   //!
   TBranch        *b_crystalROSectionId;   //!
   TBranch        *b_crystalROCrystalId;   //!
   TBranch        *b_crystalROEdep;   //!
   TBranch        *b_crystalRODose;   //!
   TBranch        *b_crystalROX;   //!
   TBranch        *b_crystalROY;   //!
   TBranch        *b_crystalROZ;   //!
   TBranch        *b_crystalROR;   //!
   TBranch        *b_nCrystalROCard;   //!
   TBranch        *b_crystalROCardSectionId;   //!
   TBranch        *b_crystalROCardCrystalId;   //!
   TBranch        *b_crystalROCardEdep;   //!
   TBranch        *b_crystalROCardDose;   //!
   TBranch        *b_crystalROCardX;   //!
   TBranch        *b_crystalROCardY;   //!
   TBranch        *b_crystalROCardZ;   //!
   TBranch        *b_crystalROCardR;   //!
   TBranch        *b_nCrateHits;   //!
   TBranch        *b_crateEdep;   //!
   TBranch        *b_crateX;   //!
   TBranch        *b_crateY;   //!
   TBranch        *b_crateZ;   //!
   TBranch        *b_crateR;   //!
   TBranch        *b_cratePdgId;   //!
   TBranch        *b_crateE;   //!
   TBranch        *b_crateEkin;   //!
   TBranch        *b_crateMass;   //!
   TBranch        *b_crateL;   //!
   TBranch        *b_vNHits;   //!
   TBranch        *b_vId;   //!
   TBranch        *b_vPdgId;   //!
   TBranch        *b_vP;   //!
   TBranch        *b_vPx;   //!
   TBranch        *b_vPy;   //!
   TBranch        *b_vPz;   //!
   TBranch        *b_vE;   //!
   TBranch        *b_vEKin;   //!
   TBranch        *b_vM;   //!
   TBranch        *b_vT;   //!
   TBranch        *b_vX;   //!
   TBranch        *b_vY;   //!
   TBranch        *b_vZ;   //!
   TBranch        *b_vCosth;   //!
   TBranch        *b_vRadius;   //!
   TBranch        *b_nTtsUpHits;   //!
   TBranch        *b_ttsUpEdep;   //!
   TBranch        *b_ttsUpX;   //!
   TBranch        *b_ttsUpY;   //!
   TBranch        *b_ttsUpZ;   //!
   TBranch        *b_ttsUpR;   //!
   TBranch        *b_ttsUpPdgId;   //!
   TBranch        *b_ttsUpEkin;   //!
   TBranch        *b_ttsUpMass;   //!
   TBranch        *b_ttsUpL;   //!
   TBranch        *b_nTtsDwHits;   //!
   TBranch        *b_ttsDwEdep;   //!
   TBranch        *b_ttsDwX;   //!
   TBranch        *b_ttsDwY;   //!
   TBranch        *b_ttsDwZ;   //!
   TBranch        *b_ttsDwR;   //!
   TBranch        *b_ttsDwPdgId;   //!
   TBranch        *b_ttsDwEkin;   //!
   TBranch        *b_ttsDwMass;   //!
   TBranch        *b_ttsDwL;   //!

  TTrackerDose(const char* Name);
  virtual ~TTrackerDose();

  int     InitChain      ();

  int     BookTrackerHistograms(TrackHist_t* Hist, const char* Folder);
  int     BookHistograms ();

  int     FillTrackerHistograms (TrackHist_t* Up, TrackHist_t* Dn);
  int     FillHistograms ();

  int     ResetTrackerHistograms (TrackHist_t* Up, TrackHist_t* Dn);
  int     ResetHistograms();

  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TChain* Chain);
  virtual void     Loop(Long64_t NEvents = -1);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif


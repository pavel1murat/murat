//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 17 16:24:24 2017 by ROOT version 6.08/04
// from TTree TCaloDose/TCaloDose
// found on file: /mu2e/data/users/gianipez/hist/treeTrackerFLASH_0.root
//
// use inheritance from TStnModule to book histograms
//////////////////////////////////////////////////////////

#ifndef TCaloDose_h
#define TCaloDose_h

#include "TNamed.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/geom/TDiskCalorimeter.hh"
#include "Stntuple/geom/TStnCrystal.hh"

// Header file for the classes stored in the TTree if any.

class TCaloDose : public TStnModule {
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

  enum { kMaxEventHistSets =  10 };
  enum { kMaxDiskHistSets  =  10 };

  struct EventHist_t {
    TH1F*     fNCrystals;
  };

  struct TrackHist_t {
    TH1F*     fZ;
  };

  struct DiskHist_t {
    TH1F*     fR;
    TH1F*     fZ;
    TH2F*     fDzVsIc;
    
    TH2F*     fRVsDz;

    TH1F*     fNVsR;

    TH1F*     fEVsZ[10]; // projections for different radii
    TH1F*     fEVsR[20]; // projections for different Z

    TH1F*     fETotVsR;
  };

  TChain*     fChain;   //! pointer to the analyzed TTree or TChain
  Int_t       fCurrent; //! current Tree number in a TChain

  struct Hist_t {
    EventHist_t*    fEvent[kMaxEventHistSets];
    DiskHist_t*     fDisk [kMaxDiskHistSets ];
  };

  Hist_t    fHist;

  TString   fProcess;

  float     fNPerPOT;   // efficiency for a given bgr source
  float     fNPOT   ;   // number of protons on target (expected)
  long int  fNEvents;
  long int  fNSimulated;

  float     fZMin[2];
  float     fZMax[2];

  TDiskCalorimeter* fCal;

  double    fE   [2][678];
  double    fdEdZ[2][678][20];
  double    fEDisk[2];

// Fixed size dimensions of array or collections stored in the TTree if any.


   // Declaration of leaf types
   Int_t          _evt;
   Int_t          _run;
   Int_t          _caloCrystals;
   Int_t          _caloDisk0Crystals;
   Int_t          _caloDisk1Crystals;
   Float_t        _caloVolume;
   Float_t        _crystalVolume;
   Int_t          _nGen;
   Int_t          _genId[1000];   //[nGen]
   Int_t          _genCrCode[1000];   //[nGen]
   Float_t        _genMomX[1000];   //[nGen]
   Float_t        _genMomY[1000];   //[nGen]
   Float_t        _genMomZ[1000];   //[nGen]
   Float_t        _genStartX[1000];   //[nGen]
   Float_t        _genStartY[1000];   //[nGen]
   Float_t        _genStartZ[1000];   //[nGen]
   Float_t        _genStartT[1000];   //[nGen]
   Int_t          _nCrystal;
   Int_t          _crystalId[10000];   //[nCrystal]
   Int_t          _crystalSectionId[10000];   //[nCrystal]
   Float_t        _crystalPosX[10000];   //[nCrystal]
   Float_t        _crystalPosY[10000];   //[nCrystal]
  Float_t         _crystalPosZ[10000];   //[nCrystal]
   Float_t        _crystalEdep[10000];   //[nCrystal]
   Float_t        _crystalDose[10000];   //[nCrystal]
   Float_t        _crystalDose0[10000];   //[nCrystal]
   Float_t        _crystalDose1[10000];   //[nCrystal]
   Float_t        _crystalDose2[10000];   //[nCrystal]
   Float_t        _crystalDose3[10000];   //[nCrystal]
   Float_t        _crystalDose4[10000];   //[nCrystal]
   Float_t        _crystalDose5[10000];   //[nCrystal]
   Float_t        _crystalDose6[10000];   //[nCrystal]
   Float_t        _crystalDose7[10000];   //[nCrystal]
   Float_t        _crystalDose8[10000];   //[nCrystal]
   Float_t        _crystalDose9[10000];   //[nCrystal]
   Float_t        _crystalDose10[10000];   //[nCrystal]
   Float_t        _crystalDose11[10000];   //[nCrystal]
   Float_t        _crystalDose12[10000];   //[nCrystal]
   Float_t        _crystalDose13[10000];   //[nCrystal]
   Float_t        _crystalDose14[10000];   //[nCrystal]
   Float_t        _crystalDose15[10000];   //[nCrystal]
   Float_t        _crystalDose16[10000];   //[nCrystal]
   Float_t        _crystalDose17[10000];   //[nCrystal]
   Float_t        _crystalDose18[10000];   //[nCrystal]
   Float_t        _crystalDose19[10000];   //[nCrystal]
   Int_t          _nCrystalRO;
   Int_t          _crystalROSectionId[100000];   //[nCrystalRO]
   Int_t          _crystalROCrystalId[100000];   //[nCrystalRO]
   Float_t        _crystalROEdep[100000];   //[nCrystalRO]
   Float_t        _crystalRODose[100000];   //[nCrystalRO]
   Float_t        _crystalROX[100000];   //[nCrystalRO]
   Float_t        _crystalROY[100000];   //[nCrystalRO]
   Float_t        _crystalROZ[100000];   //[nCrystalRO]
   Float_t        _crystalROR[100000];   //[nCrystalRO]
   Int_t          _nCrystalROCard;
   Int_t          _crystalROCardSectionId[100000];   //[nCrystalROCard]
   Int_t          _crystalROCardCrystalId[100000];   //[nCrystalROCard]
   Float_t        _crystalROCardEdep[100000];   //[nCrystalROCard]
   Float_t        _crystalROCardDose[100000];   //[nCrystalROCard]
   Float_t        _crystalROCardX[100000];   //[nCrystalROCard]
   Float_t        _crystalROCardY[100000];   //[nCrystalROCard]
   Float_t        _crystalROCardZ[100000];   //[nCrystalROCard]
   Float_t        _crystalROCardR[100000];   //[nCrystalROCard]
   Int_t          _nCrateHits;
   Float_t        _crateEdep[100000];   //[nCrateHits]
   //   Float_t         crateDose[100000];   //[nCrateHits]
   Float_t        _crateX[100000];   //[nCrateHits]
   Float_t        _crateY[100000];   //[nCrateHits]
   Float_t        _crateZ[100000];   //[nCrateHits]
   Float_t        _crateR[100000];   //[nCrateHits]
   Int_t          _cratePdgId[100000];   //[nCrateHits]
   Float_t        _crateE[100000];   //[nCrateHits]
   Float_t        _crateEkin[100000];   //[nCrateHits]
   Float_t        _crateMass[100000];   //[nCrateHits]
   Float_t        _crateL[100000];   //[nCrateHits]
   Int_t          _vNHits;
   Int_t          _vId[100000];   //[vNHits]
   Int_t          _vPdgId[100000];   //[vNHits]
   Float_t        _vP[100000];   //[vNHits]
   Float_t        _vPx[100000];   //[vNHits]
   Float_t        _vPy[100000];   //[vNHits]
   Float_t        _vPz[100000];   //[vNHits]
   Float_t        _vE[100000];   //[vNHits]
   Float_t        _vEKin[100000];   //[vNHits]
   Float_t        _vM[100000];   //[vNHits]
   Float_t        _vT[100000];   //[vNHits]
   Float_t        _vX[100000];   //[vNHits]
   Float_t        _vY[100000];   //[vNHits]
   Float_t        _vZ[100000];   //[vNHits]
   Float_t        _vCosth[100000];   //[vNHits]
   Float_t        _vRadius[100000];   //[vNHits]

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
   //   TBranch        *b_crateDose;   //!
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




  TCaloDose(const char* Name);
  virtual ~TCaloDose();

  int     InitChain      ();

  int     BookDiskHistograms   (DiskHist_t* Hist , const char* Folder);
  int     BookEventHistograms  (EventHist_t* Hist, const char* Folder);
  int     BookHistograms ();

  int     FillDiskHistograms    (DiskHist_t* Hist );
  int     FillEventHistograms   (EventHist_t* Hist);
  int     FillHistograms ();

  int     ResetHistograms();
					// calculate fraction of the crystal area inside the two radii

  float   insideFraction(TStnCrystal* Crystal, float RMin, float RMax);

  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TChain* Chain);
  virtual void     Loop(Long64_t NEvents = -1);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  ClassDef(TCaloDose,0)
};

#endif


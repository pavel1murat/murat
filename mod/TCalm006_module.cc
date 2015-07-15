//-----------------------------------------------------------------------------
// 2013-10-19 P.Murat: analysis of the energy deposition in a crystal block
//                     optimize the crystal size
// bit usage: 
//
//-----------------------------------------------------------------------------

// #include "Stntuple/obj/AbsEvent.hh"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"

#include "murat/mod/TCalm006_module.hh"
#include "Stntuple/mod/TAnaDump.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "TTrackerGeom/inc/TTracker.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  TCalm006::TCalm006(fhicl::ParameterSet const& pset): 
    THistModule(pset,"TCalm006"),
    fG4RunModuleLabel(pset.get<std::string> ("g4ModuleLabel"  , "g4run"        )),
    fHistFileName    (pset.get<std::string> ("histFileName"   , "tcalm006.hist")),
    fProductName     (pset.get<std::string> ("productName"    , "calorimeter"  )),
    fMinTActive      (pset.get<double>      ("minTActive"     ,   710.         )),
    fCrystalHitMaker (pset.get<std::string> ("crystalHitMaker", "CaloCrystalHitsMaker"))
  {
					// reset all histogram pointers

    for (int i=0; i<kNEventHistSets; i++) {
      fHist.fEvent[i] = 0;
    }

    fProcessName = "";
    fHexMap      = new HexMap();
  }


//-----------------------------------------------------------------------------
  TCalm006::~TCalm006() {

    //    delete fHexMap;
  }

//-----------------------------------------------------------------------------
  void TCalm006::BookEventHistograms(EventHist_t* Hist, const char* Folder) {

    char   name[100], title[100];

    HBook1F(Hist->fEventNumber,"event" ,Form("%s: Event Number"          ,Folder),1000,0,10000,Folder);
    HBook1F(Hist->fRunNumber  ,"run"   ,Form("%s: Run  Number"           ,Folder),1000,0,10000,Folder);
    HBook1F(Hist->fEEle       ,"eele"  ,Form("%s: Inital Electron Energy",Folder),110,110,110,Folder);
    
    for (int i=0; i<10; i++) {
      sprintf(name,"e%i",i);
      sprintf(title,"%s: energy in the ring number %i",Folder,i);
      HBook1F(Hist->fEnergy[i] ,name, title,110,0,110,Folder);

      sprintf(name,"esum%i",i);
      sprintf(title,"%s: Sum energy in the rings up to %i",Folder,i);
      HBook1F(Hist->fESum  [i] ,name, title,100,0,110,Folder);
    }
  }

//-----------------------------------------------------------------------------
  void TCalm006::BookHistograms() {
    TFolder  *hist_folder, *fol;
    char     folder_name[200];

    hist_folder = (TFolder*) fFolder->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[0] = 1;		// all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
	sprintf(folder_name,"evt_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fEvent[i] = new EventHist_t;
	BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
      }
    }
  }

//-----------------------------------------------------------------------------
// begin job - book histograms etc
//-----------------------------------------------------------------------------
  void TCalm006::beginJob() {

//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
    BookHistograms();
  }

//-----------------------------------------------------------------------------
// begin run 
//-----------------------------------------------------------------------------
  bool TCalm006::beginRun(art::Run& R) {

//     GeomHandle<TTracker> tt;
//     fTracker = (TTracker*) tt.operator -> ();

    //    GeomHandle<VaneCalorimeter> cg;
    //    fCal = &(*cg);
    //    fCal = cg.operator -> ();

    return true;
  }


  //-----------------------------------------------------------------------------
  void TCalm006::endJob() {
    TDirectory* old_dir = gDirectory;
    
    TFile* f = TFile::Open(fHistFileName.data(),"recreate");

    SaveFolder(fFolder,f);
    f->Write();
    f->Close();
    
    old_dir->cd();
  }

//-----------------------------------------------------------------------------
  void TCalm006::FillEventHistograms(EventHist_t* Hist) {
    
    Hist->fEventNumber->Fill(fEventNumber);
    Hist->fRunNumber->  Fill(fRunNumber);

    for (int i=0; i<10; i++) {
      Hist->fEnergy[i]->Fill(fEnergy[i]);
      Hist->fESum  [i]->Fill(fESum  [i]);
    }
  }

//-----------------------------------------------------------------------------
  void TCalm006::FillHistograms() {
    //    const char oname[] = "TCalm006::FillHistograms";


    //    mu2e::StepPointMC*   step;
    //    int                  ndisks;
    //    double               rin[2], rout[2];
    //    int                  inside, pdg_code;
    //    double               x0, y0, z0, r0;

    art::ServiceHandle<mu2e::GeometryService> geom;
    //    const mu2e::DiskCalorimeter* cal;
    //    const mu2e::Disk* disk;
    
    // if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
      // mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
      // cal  = dc.operator->();
      // ndisks = cal->nDisk();
      // for (int i=0; i<ndisks; i++) {
      // 	rin[i]  = cal->disk(i).innerRadius();
      // 	rout[i] = cal->disk(i).outerRadius();
      // }
    // }
    // else {
    //   printf(">>> ERROR: disk calorimeter not found.\n");
    // }
//-----------------------------------------------------------------------------
// event histograms, everything is in mm!!!
//-----------------------------------------------------------------------------
    FillEventHistograms(fHist.fEvent[0]);
  }

//-----------------------------------------------------------------------------
  void TCalm006::getData(art::Event* Evt) {
    const char* oname = "TCalm006::getData";
//-----------------------------------------------------------------------------
// CaloCrystalHitCollection
//-----------------------------------------------------------------------------
    art::Handle<CaloCrystalHitCollection> ccHandle;
    Evt->getByLabel(fCrystalHitMaker.data(),ccHandle);
    if (ccHandle.isValid()) {
      fListOfCaloCrystalHits = (CaloCrystalHitCollection*) ccHandle.product();
      fNCaloCrystalHits      = fListOfCaloCrystalHits->size();
    }
    else {
      fListOfCaloCrystalHits = NULL;
      fNCaloCrystalHits      = 0;
      printf(">>> ERROR in %s: failed to locate CaloCrystalHitCollection, BAIL OUT.\n",oname);
    }

  }



//-----------------------------------------------------------------------------
  void TCalm006::Init(art::Event* Evt) {
    //    const char oname[] = "TCalm006::Init";

    int             id, ring, l, k;
    CaloCrystalHit* hit;

    for (int i=0; i<10; i++) {
      fEnergy[i] = 0;
    }

    int nhits = fListOfCaloCrystalHits->size();

    for (int i=0; i<nhits; i++) {
      hit = &fListOfCaloCrystalHits->at(i);
      id  = hit->id();

      l   = fHexMap->l(id);
      k   = fHexMap->k(id);

      HexLK lk(l,k);

      ring = fHexMap->ring(lk);

      fEnergy[ring] += hit->energyDep();
    }

    fESum[0] = fEnergy[0];
    for (int i=1; i<10; i++) {
      fESum[i] = fESum[i-1]+fEnergy[i];
    }

  }



//-----------------------------------------------------------------------------
  void TCalm006::Debug(art::Event* Evt) {
    //    const char* oname = "TCalm006::Debug";
//-----------------------------------------------------------------------------
// bit 000: calorimeter step points
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[0] != 0) {
      TModule::fDump->printStepPointMCCollection("g4run","calorimeter");
    }
  }

//-----------------------------------------------------------------------------
  bool TCalm006::filter(art::Event& Evt) {
    const char* oname = "TCalm006::filter";

    bool rc(true);

    fEventNumber = Evt.event();
    fRunNumber   = Evt.run();

    printf(" >>>>>>> [%s] RUN: %10i  EVENT: %10i\n",oname,fRunNumber,fEventNumber);

    getData(&Evt);
    Init   (&Evt);
    Debug  (&Evt);

    FillHistograms();
					// so far - fake
    if (TModule::fDebugBit[51] != 0) {
      rc = 1;
    }

    TModule::filter(Evt);
    
    return rc;
  }
  
}  // end namespace mu2e

using mu2e::TCalm006;

DEFINE_ART_MODULE(TCalm006);

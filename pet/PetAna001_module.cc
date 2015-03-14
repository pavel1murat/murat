//-----------------------------------------------------------------------------
// PetAna001: analyse output of the PetParticleGun
//
// debug bits:
// 000:           event with >= 2 calorimeter hits
//-----------------------------------------------------------------------------
#include "TEnv.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"

#include "Stntuple/mod/TAnaDump.hh"

#include "murat/pet/PetAna001_module.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include "murat/pet/PetGeometryService.hh"
#include "murat/pet/PetGeomHandle.hh"
#include "murat/pet/BrainImager.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  PetAna001::PetAna001(fhicl::ParameterSet const& pset): 
    THistModule(pset,"PetAna001"),
    fHistFileName       (pset.get<std::string>("histFileName"    , "petana001.hist")),
    fGenParticleMaker   (pset.get<std::string>("genParticleMaker", "petgen")),
    fCaloCrystalHitMaker(pset.get<std::string>("crystalHitMaker" , "petMakeCaloCrystalHits")),
    fMinPhotopeakE      (pset.get<double>     ("minPhotopeakE"   , 0.420)),
    fMaxPhotopeakE      (pset.get<double>     ("maxPhotopeakE"   , 0.600)),
    fTriggerDt          (pset.get<double>     ("triggerDt"       , 2.0  )),
    fPrintFrequency     (pset.get<int>        ("printFrequency"  , 10   ))
  {

    fFolder = new TFolder("PetAna001","PetAna001 Folder");
    fFolder->AddFolder("Hist","Hist");
					// reset all histogram pointers

    for (int i=0; i<kNEventHistSets; i++) {
      fHist.fEvent[i] = 0;
    }

    // fMinEPhotopeak = 420.;
    // fMaxEPhotopeak = 600.;
  }



//-----------------------------------------------------------------------------
  PetAna001::~PetAna001() {
    printf(" >>> PetAna001::~PetAna001 : deleting PetAna001\n");

//     fFolder->Delete();
//     delete fFolder;

  }

//-----------------------------------------------------------------------------
  void PetAna001::BookPhotonHistograms(PhotonHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fTime ,"time" ,Form("%s: Photon Time"      ,Folder),1024,   0,1024,Folder);
    HBook1F(Hist->fCosTh,"costh",Form("%s: Photon Cos(Theta)",Folder), 100,  -1,   1,Folder);
    HBook1F(Hist->fPhi  ,"phi"  ,Form("%s: Photon Phi"       ,Folder), 500,  -5,   5,Folder);
    HBook1F(Hist->fZ    ,"z"    ,Form("%s: Photon Z"         ,Folder),1000, -900,100,Folder);
    HBook1F(Hist->fR    ,"r"    ,Form("%s: Photon R"         ,Folder), 250,   0, 250,Folder);
    HBook2F(Hist->fRVsZ ,"r_vs_z",Form("%s: Photon R_Vs_Z"   ,Folder), 1000,-900,100,250, 0, 250,Folder);
  }

//-----------------------------------------------------------------------------
  void PetAna001::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fNPhotons     ,"nph"      ,Form("%s: Number of Photons"       ,Folder),250,0,250,Folder);
    HBook1F(Hist->fNHitCrystals ,"nhcr"     ,Form("%s: Number of Hit Crystals"  ,Folder),50,0,50,Folder);
    HBook1F(Hist->fNHitModules  ,"nhmod"    ,Form("%s: Number of Hit Modules"   ,Folder), 40,0, 40,Folder);
    HBook1F(Hist->fEMax         ,"emax"     ,Form("%s: E(max)"                  ,Folder),200,0,1.,Folder);
    HBook1F(Hist->fETot         ,"etot"     ,Form("%s: E(TOT) module 0"         ,Folder),200,0,1.,Folder);
    HBook2F(Hist->fE2VsE1Cr     ,"e2_vs_e1_cr" ,Form("%s: E2 vs E1 Crystals"    ,Folder),200,0,1.,200,0,1.,Folder);
    HBook2F(Hist->fE2VsE1Wed    ,"e2_vs_e1_wed",Form("%s: E2 vs E1 Wedge"       ,Folder),200,0,1.,200,0,1.,Folder);
    HBook2F(Hist->fW2VsW1       ,"w2_vs_w1"    ,Form("%s: W2 vs W1"             ,Folder), 40,0,40,40,0,40,Folder);
    HBook1F(Hist->fHitTime      ,"hit_time" ,Form("%s: Hit Time"                ,Folder),1200,-100,1100,Folder);
    HBook1F(Hist->fDtMin        ,"dt_min"   ,Form("%s: Dt(Min)"                 ,Folder),500,-50, 50,Folder);
    HBook1F(Hist->fDist         ,"dist"     ,Form("%s: Dist"                    ,Folder),200,  0,200,Folder);
    HBook1F(Hist->fCType        ,"ctype"    ,Form("%s: Coincidence Type"        ,Folder), 10,  0, 10,Folder);
  }

//-----------------------------------------------------------------------------
  void PetAna001::BookHistograms() {
    TFolder  *hist_folder, *fol;
    char     folder_name[200];

    hist_folder = (TFolder*) fFolder->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[0] = 1;		// all events
    book_event_histset[1] = 1;		// events with one hit
    book_event_histset[2] = 1;		// events with two or more hits
    book_event_histset[3] = 1;		// events with two or more hits in the photopeak
    book_event_histset[4] = 1;		// trigger events

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
// book photon histograms
//-----------------------------------------------------------------------------
    int book_photon_histset[kNPhotonHistSets];
    for (int i=0; i<kNPhotonHistSets; i++) book_photon_histset[i] = 0;

    book_photon_histset[0] = 1;		// all photons
    book_photon_histset[1] = 1;		// photons in events with 1 hit module
    book_photon_histset[2] = 1;		// photons in events with >=2 hit modules
    book_photon_histset[3] = 1;		// photons in events with 2 photons in the photopeak
    book_photon_histset[4] = 1;		// photons in "trigger" events 

    for (int i=0; i<kNPhotonHistSets; i++) {
      if (book_photon_histset[i] != 0) {
	sprintf(folder_name,"pho_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fPhoton[i] = new PhotonHist_t;
	BookPhotonHistograms(fHist.fPhoton[i],Form("Hist/%s",folder_name));
      }
    }
  }

//-----------------------------------------------------------------------------
// begin job - book histograms etc
//-----------------------------------------------------------------------------
  void PetAna001::beginJob() {
    BookHistograms();
  }


//-----------------------------------------------------------------------------
// begin run - 
//-----------------------------------------------------------------------------
  bool PetAna001::beginRun(art::Run& aRun) {
    return true;
  }


//-----------------------------------------------------------------------------
  void PetAna001::endJob() {
    TDirectory* old_dir = gDirectory;
    
    TFile* f = TFile::Open(fHistFileName.data(),"recreate");

    SaveFolder(fFolder,f);
    f->Write();
    f->Close();
    
    old_dir->cd();
  }

//-----------------------------------------------------------------------------
// histograms for the photon emission point
//-----------------------------------------------------------------------------
  void PetAna001::FillPhotonHistograms(PhotonHist_t* Hist, mu2e::GenParticle* Photon) {

    double   cos_th, phi, r, z, t;

    cos_th = Photon->momentum().pz()/Photon->momentum().vect().mag();
    phi    = Photon->momentum().phi();
    r      = Photon->position().perp();
    z      = Photon->position().z();
    t      = Photon->time();

    Hist->fTime->Fill(t);
    Hist->fCosTh->Fill(cos_th);
    Hist->fPhi->Fill(phi);
    Hist->fZ->Fill(z);
    Hist->fR->Fill(r);
    Hist->fRVsZ->Fill(z,r);
  }

//-----------------------------------------------------------------------------
  void PetAna001::FillEventHistograms(EventHist_t* Hist) {

    double  e0_max, e1_max, e0_tot, e1_tot, emax, etot;
    int     l0, l1;

    Hist->fNPhotons->Fill(fNPhotons);
    Hist->fNHitCrystals->Fill(fNHitCrystals);
    Hist->fNHitModules->Fill(fNHitModules);

    emax = -1;
    etot = -1;
    if (fNHitModules > 0) {
      l0   = fModIndex[0];
      emax = fEMaxModule[l0];
      etot = fETotModule[l0];
    }

    Hist->fEMax->Fill(emax);
    Hist->fETot->Fill(etot);

    if (fNHitModules > 1) {
      l0 = fModIndex[0];
      l1 = fModIndex[1];

      e0_tot = fETotModule[l0];
      e1_tot = fETotModule[l1];

      e0_max = fEMaxModule[l0];
      e1_max = fEMaxModule[l1];
    }
    else {
      l0 = -1;
      l1 = -1;
      e0_tot = -1.;
      e1_tot = -1.;

      e0_max = -1.;
      e1_max = -1.;
    }

    Hist->fE2VsE1Cr->Fill(e0_max,e1_max);
    Hist->fE2VsE1Wed->Fill(e0_tot,e1_tot);
    Hist->fW2VsW1->Fill(l0,l1);

    Hist->fDtMin->Fill(fDtMin);
    Hist->fDist->Fill(fDistance);
    Hist->fCType->Fill(fCoincidenceType);
  }

//-----------------------------------------------------------------------------
  const mu2e::SimParticle* PetAna001::FindPhoton(int CrystalID) {
    const char* oname = "PetAna001::FindPhoton";
    
    const mu2e::SimParticle* sim (0);

    int nsteps = fStepPointMCCollection->size();

    const mu2e::StepPointMC* step;

    //    int banner_printed = 0;
    for (int i=0; i<nsteps; i++) {
      step = &fStepPointMCCollection->at(i);
      if (step->volumeId() == CrystalID) {
//-----------------------------------------------------------------------------
// hit in the right crystal
//-----------------------------------------------------------------------------
	art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle();
	sim  = simptr.operator ->();
	break;
      }
    }

    return sim;
  }


//-----------------------------------------------------------------------------
// this analysis code is intended to run on the obsolete dataset, so it can 
// contain various kludges
//-----------------------------------------------------------------------------
  void PetAna001::FillHistograms() {
//-----------------------------------------------------------------------------
// EVT[0]: all events
// EVT[1]: events with one hit module
// EVT[2]: events with two or more hit modules
//-----------------------------------------------------------------------------
    FillEventHistograms(fHist.fEvent[0]);

    if (fNHitModules == 1) FillEventHistograms(fHist.fEvent[1]);
    if (fNHitModules >= 2) FillEventHistograms(fHist.fEvent[2]);
//-----------------------------------------------------------------------------
// EVT[3]: events with both photons in the photopeak region
//-----------------------------------------------------------------------------
    int l0 = fModIndex[0];
    int l1 = fModIndex[1];

    bool photopeak = ((fEMaxModule[l0] > fMinPhotopeakE) && 
		      (fEMaxModule[l0] < fMaxPhotopeakE) && 
		      (fEMaxModule[l1] > fMinPhotopeakE) && 
		      (fEMaxModule[l1] < fMaxPhotopeakE)    );

    if (photopeak) FillEventHistograms(fHist.fEvent[3]);

    bool trigger = photopeak && (fabs(fDtMin) < fTriggerDt);

    if (trigger)  FillEventHistograms(fHist.fEvent[4]);
//-----------------------------------------------------------------------------
// fill photon histograms
//-----------------------------------------------------------------------------
    GenParticle* photon;
    int          nph;

    nph = fListOfGenParticles->size();
    for (int i=0; i<nph; i++) {
      photon =  (GenParticle*) &fListOfGenParticles->at(i);
      FillPhotonHistograms(fHist.fPhoton[0],photon);

      if (fNHitModules == 1) FillPhotonHistograms(fHist.fPhoton[1],photon);
      if (fNHitModules >= 2) FillPhotonHistograms(fHist.fPhoton[2],photon);
      if (photopeak        ) FillPhotonHistograms(fHist.fPhoton[3],photon);
      if (trigger          ) FillPhotonHistograms(fHist.fPhoton[4],photon);
    }
  }

//-----------------------------------------------------------------------------
  int PetAna001::getData(const art::Event& Evt) {
    const char* oname = "PetAna001::getData";
//-----------------------------------------------------------------------------
// generated particles
//-----------------------------------------------------------------------------
    int rc = 0;

    art::Handle<GenParticleCollection> genpHandle;
    Evt.getByLabel(fGenParticleMaker.data(),genpHandle);
    if (genpHandle.isValid()) {
      fListOfGenParticles = (GenParticleCollection*) &(*genpHandle);
    }
    else {
      printf(">>> PetAna001::getData ERROR: genpHandle for %s invalid. BAIL OUT\n", fGenParticleMaker.data());
      return -1;
    }

    fNPhotons = fListOfGenParticles->size();


    art::Handle<mu2e::StepPointMCCollection> handle;
    art::Selector  selector(art::ProductInstanceNameSelector("calorimeter") &&
			    art::ProcessNameSelector        ("")            && 
			    art::ModuleLabelSelector        ("PetG4")          );
    Evt.get(selector, handle);

    if (handle.isValid()) fStepPointMCCollection = handle.product();
    else {
      printf(">>> ERROR in %s: failed to locate StepPointMCCollection",oname);
      printf(". BAIL OUT. \n");
      return -1;
    }

    art::Selector  s2(art::ProductInstanceNameSelector("") &&
		      art::ProcessNameSelector("")         && 
		      art::ModuleLabelSelector(fCaloCrystalHitMaker.data()));

    art::Handle<mu2e::CaloCrystalHitCollection> caloCrystalHitsHandle;
    
    Evt.get(s2,caloCrystalHitsHandle);
    
    fListOfCaloCrystalHits = caloCrystalHitsHandle.product();

    fNHitCrystals = fListOfCaloCrystalHits->size();

    const mu2e::CaloCrystalHit* hit;
//-----------------------------------------------------------------------------
// preset times and energies to something unreasonable
//-----------------------------------------------------------------------------
    for (int i=0; i<fNWedges; i++) {
      fNHitsPerModule[i] = 0;
      fETotModule    [i] = 0;
      fEMaxModule    [i] = 0;
      fTMaxModule    [i] = -1.e6;
      fModIndex      [i] = i;
    }

    fNHitModules = 0;

    // printf("----------------------------------------------------------------\n");
    // printf("CrystalID      Time   Energy    EnergyTot  NRoids               \n");
    // printf("----------------------------------------------------------------\n");
    
    art::ServiceHandle<PetGeometryService> geom;    
    if( !(geom->hasElement<BrainImager>()) ) {
      printf(">>> ERROR in PetAna001::getData: BrainImager not defined. RETURN\n");
      return -1;
    }

    mu2e::PetGeomHandle<BrainImager> imager;

    int nr = imager->nCrystalR();
    int nz = imager->nCrystalZ();

    int nc_per_wedge = nr*nz;

    int                      iw, loc;
    double                   ehit;
    const mu2e::Vane         *wedge;
    const mu2e::Crystal      *cr;
    CLHEP::Hep3Vector        pos;

    for (int ic=0; ic<fNHitCrystals; ic++) {
      hit  = &fListOfCaloCrystalHits->at(ic);

      iw   = hit->id()/nc_per_wedge;
      loc  = hit->id()-iw*nc_per_wedge;

      // define pointer to the crystal

      wedge = imager->Vane(iw);
      cr    = &wedge->crystal(loc);
      pos   = imager->fromVaneFrame(iw,cr->position());
      
      fNHitsPerModule[iw] += 1;
      if (fNHitsPerModule[iw] == 1) {
	fNHitModules += 1;
      }

      ehit = hit->energyDep();
      fETotModule[iw] += ehit;
      if (ehit > fEMaxModule[iw]) {
	fEMaxModule[iw] = ehit;
	fTMaxModule[iw] = hit->time();
	fHitPos    [iw] = pos;
	fHit       [iw] = hit;
      }

    //   printf("%7i  %10.3f %10.3f %10.3f %5i\n",
    // 	     hit->id(),
    // 	     hit->time             (),
    // 	     hit->energyDep        (),
    // 	     hit->energyDepTotal   (),
    // 	     hit->numberOfROIdsUsed());
    }
//-----------------------------------------------------------------------------
// finally, need to determine two modules with highest energy depositions
// because of that, module #2 sometimes may have a crystal with higher 
// energy than module #1 ...
//-----------------------------------------------------------------------------
    int    li, lj;
    double ei, ej;

    for (int i=0; i<fNWedges-1; i++) {
      li = fModIndex  [i];
      ei = fEMaxModule[li];

      for (int j=i+1; j<fNWedges; j++) {
	lj = fModIndex  [j];
	ej = fEMaxModule[lj];

	if (ej > ei) {
	  fModIndex[i] = lj;
	  fModIndex[j] = li;
	  li           = lj;
	  ei           = ej;
	}
      }
    }
    // printf(" i fModIndex[i],fEMaxModule[i],fETotModule[i],fTMaxModule[i]\n");
    // for (int i=0; i<fNWedges; i++) {
    //   printf("%3i %3i %10.3f %10.3f %10.3f\n",
    // 	     i,fModIndex[i],fEMaxModule[i],fETotModule[i],fTMaxModule[i]);
    // }
//-----------------------------------------------------------------------------
// trigger emulation, use modules ordered in energy
//-----------------------------------------------------------------------------
    double   dt, dt_min, e1, e2;
    int      l1, l2;

    dt_min = 1.e6;

    fTrigIndex[0] = -1;
    fTrigIndex[1] = -1;

    for (int i1=0; i1<fNWedges-1; i1++)  {
      l1 = fModIndex  [i1];
      e1 = fEMaxModule[l1];

      if ((e1 > fMinPhotopeakE) && (e1 < fMaxPhotopeakE)) {

	for (int i2=i1+1; i2<fNWedges; i2++)  {
	  l2 = fModIndex  [i2];
	  e2 = fEMaxModule[l2];

	  if ((e2 > fMinPhotopeakE) && (e2 < fMaxPhotopeakE)) {
	    dt = fTMaxModule[l1]- fTMaxModule[l2];
	    if (fabs(dt) < fabs(dt_min)) {
//-----------------------------------------------------------------------------
// new "best" trigger pair
//-----------------------------------------------------------------------------
	      dt_min        = dt;
	      fTrigIndex[0] = l1;
	      fTrigIndex[1] = l2;
//
	    }
	  }
	}
      }
    }

    fDtMin = dt_min;
//-----------------------------------------------------------------------------
// finally, determine whether this is a true, scatter, or random coincidence
//-----------------------------------------------------------------------------
    const mu2e::SimParticle  *ph1, *ph2;
    const CLHEP::Hep3Vector  *pos1, *pos2, *v;
    CLHEP::Hep3Vector        dp, n, d;
    int                      iv1, iv2;

    fDistance        = -1.;
    fCoincidenceType = -1;

    if (fTrigIndex[1] >= 0) {
      l1 = fTrigIndex[0];
      l2 = fTrigIndex[1];
//-----------------------------------------------------------------------------
// find photon which generated that hit
//-----------------------------------------------------------------------------
      ph1 = FindPhoton(fHit[l1]->id());
      ph2 = FindPhoton(fHit[l2]->id());

      iv1 = (ph1->id().asInt()-1)/2;
      iv2 = (ph2->id().asInt()-1)/2;

      if (iv1 == iv2) {
//-----------------------------------------------------------------------------
// photons are produced in the same annihilation, 
// check whether this is a true or scattered event - try to use the reconstructed
// quantities only - find crystals and calculate the distance between the vertex 
// and the line, connecting the 2 crystals
//-----------------------------------------------------------------------------
	v    = &ph1->startPosition();
	
	pos1 = &fHitPos[l1];
	pos2 = &fHitPos[l2];

	dp   = *pos2-*pos1;
	n    = dp.unit();

	d         = (*pos1-*v) - n*((*pos1-*v).dot(n));
	fDistance = d.mag();

	if (fDistance < 8.) {
					// "true" coincidence
	  fCoincidenceType = 0;
	}
	else {
					// "scatter"
	  fCoincidenceType = 1;
	}
      }
      else {
//-----------------------------------------------------------------------------
// random coincidence
//-----------------------------------------------------------------------------
	fCoincidenceType = 2;
      }
    }

    return rc;
  }



//-----------------------------------------------------------------------------
  void PetAna001::Init(art::Event& Evt) {
    //
    static int first_call(1);

    if (first_call) {
					// initialize geometry
      PetGeomHandle<BrainImager> h;
      fBrainImager = (BrainImager*) h.get();

      fNWedges = fBrainImager->NWedges();
    }
  }


//-----------------------------------------------------------------------------
  void PetAna001::Debug(const art::Event& Evt) {
    const char* oname = "PetAna001::Debug";
 //-----------------------------------------------------------------------------
// bit 000: 
//-----------------------------------------------------------------------------
   if (TModule::fDebugBit[0] != 0) {
     if (fNHitModules >= 2) {
       TModule::fDump->printStepPointMCCollection("PetG4","calorimeter");
       TModule::fDump->printCaloCrystalHits      ("PetMakeCaloCrystalHits");
     }
   }
  }

//-----------------------------------------------------------------------------
  bool PetAna001::filter(art::Event& Evt) {
    const char* oname = "PetAna001::filter";

    bool rc(true);

    if ((Evt.event() % fPrintFrequency) == 0) {
      printf(" >>>>>>> [%s] Run:EVENT : %10i:%10i\n",oname,Evt.run(),Evt.event());
    }

    TModule::fDump->SetEvent(Evt);

    rc = getData(Evt);
    if (rc < 0) return false;

    Init(Evt);

    Debug(Evt);

    FillHistograms();

//     if (TModule::fDebugBit[51] != 0) {
//       rc = (fNMatchedTracks > 0);
//     }

    TModule::filter(Evt);

    return rc;
  }

}  // end namespace mu2e

using mu2e::PetAna001;

DEFINE_ART_MODULE(PetAna001);

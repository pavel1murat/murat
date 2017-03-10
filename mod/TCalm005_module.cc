//-----------------------------------------------------------------------------
// 2013-08-29 P.Murat: study origin of the particles hitting the calorimeter
//                     assume disk-based geometry
// bit usage: 
//
//-----------------------------------------------------------------------------

// #include "Stntuple/obj/AbsEvent.hh"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"

#include "murat/mod/TCalm005_module.hh"
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
  TCalm005::TCalm005(fhicl::ParameterSet const& pset): 
    THistModule(pset,"TCalm005"),
    fG4RunModuleLabel(pset.get<std::string> ("g4ModuleLabel", "g4run"        )),
    fHistFileName    (pset.get<std::string> ("histFileName" , "tcalm005.hist")),
    fProductName     (pset.get<std::string> ("productName"  , "calorimeter"  )),
    fMinTActive      (pset.get<double>      ("minTActive"   ,   710.         ))
  {
					// reset all histogram pointers

    for (int i=0; i<kNEventHistSets; i++) {
      fHist.fEvent[i] = 0;
    }

    for (int i=0; i<kNStepPointMCHistSets; i++) {
      fHist.fStepPointMC[i] = 0;
    }

    fProcessName = "";
  }

//-----------------------------------------------------------------------------
  TCalm005::~TCalm005() {
  }

//-----------------------------------------------------------------------------
  void TCalm005::BookStepPointMCHistograms(StepPointMCHist_t* Hist, const char* Folder) {

    HBook2F(Hist->fRVsZ      ,"rz"   ,Form("%s: Step Point MC R vs Z"       ,Folder),1500,0,15000,200,0,2000,Folder);
    HBook2F(Hist->fR0VsZ0    ,"r0z0" ,Form("%s: Step Point MC R0 vs Z0"     ,Folder),1500,0,15000,200,0,2000,Folder);
    HBook1F(Hist->fPdgCode[0],"pdg_0",Form("%s: Step Point MC PDG Code"     ,Folder), 200,-100, 100,Folder);
    HBook1F(Hist->fPdgCode[1],"pdg_1",Form("%s: Step Point MC PDG Code"     ,Folder),1000,  0, 10000,Folder);
    HBook1F(Hist->fTime      ,"time" ,Form("%s: Step Point MC Time"         ,Folder),1000,0,10000,Folder);
    HBook1F(Hist->fLength    ,"len"  ,Form("%s: Step Point MC Length (mm)"  ,Folder), 200,0,200,Folder);
    HBook1F(Hist->fEDep      ,"edep" ,Form("%s: Step Point MC EDep (MeV)"   ,Folder), 200,0,100,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm005::BookEventHistograms(EventHist_t* Hist, const char* Folder) {

    HBook1F(Hist->fEventNumber,"event" ,Form("%s: Event Number"                ,Folder),1000,0,10000,Folder);
    HBook1F(Hist->fRunNumber  ,"run"   ,Form("%s: Run  Number"                 ,Folder),1000,0,10000,Folder);
    HBook1F(Hist->fNSteps     ,"nsteps",Form("%s: Number of CAL StepPointMC's" ,Folder), 200,0,  200,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm005::BookHistograms() {
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
//-----------------------------------------------------------------------------
// book StepPointMC histograms
//-----------------------------------------------------------------------------
    int book_step_point_mc_histset[kNStepPointMCHistSets];
    for (int i=0; i<kNStepPointMCHistSets; i++) book_step_point_mc_histset[i] = 0;

    book_step_point_mc_histset[0] = 1;		// all step points
    book_step_point_mc_histset[1] = 1;		// abs(PDG) < 100: electrons, photons 
    book_step_point_mc_histset[2] = 1;		// neutrons : 100  < abs(PDG) 10000
    book_step_point_mc_histset[3] = 1;		// nuclei   : abs(PDG) > 10000

    book_step_point_mc_histset[10] = 1;		// all step points particles outside 
    book_step_point_mc_histset[11] = 1;		// abs(PDG) < 100: electrons, photons 
    book_step_point_mc_histset[12] = 1;		// neutrons : 100  < abs(PDG) 10000
    book_step_point_mc_histset[13] = 1;		// nuclei   : abs(PDG) > 10000

    book_step_point_mc_histset[20] = 1;		// all step points particles inside
    book_step_point_mc_histset[21] = 1;		// abs(PDG) < 100: electrons, photons 
    book_step_point_mc_histset[22] = 1;		// neutrons : 100  < abs(PDG) 10000
    book_step_point_mc_histset[23] = 1;		// nuclei   : abs(PDG) > 10000

    for (int i=0; i<kNStepPointMCHistSets; i++) {
      if (book_step_point_mc_histset[i] != 0) {
	sprintf(folder_name,"stp_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fStepPointMC[i] = new StepPointMCHist_t;
	BookStepPointMCHistograms(fHist.fStepPointMC[i],Form("Hist/%s",folder_name));
      }
    }
  }

//-----------------------------------------------------------------------------
// begin job - book histograms etc
//-----------------------------------------------------------------------------
  void TCalm005::beginJob() {

//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
    BookHistograms();
  }

//-----------------------------------------------------------------------------
// begin run 
//-----------------------------------------------------------------------------
  bool TCalm005::beginRun(art::Run& R) {

//     GeomHandle<TTracker> tt;
//     fTracker = (TTracker*) tt.operator -> ();

    //    GeomHandle<VaneCalorimeter> cg;
    //    fCal = &(*cg);
    //    fCal = cg.operator -> ();

    return true;
  }


  //-----------------------------------------------------------------------------
  void TCalm005::endJob() {
    TDirectory* old_dir = gDirectory;
    
    TFile* f = TFile::Open(fHistFileName.data(),"recreate");

    SaveFolder(fFolder,f);
    f->Write();
    f->Close();
    
    old_dir->cd();
  }

//-----------------------------------------------------------------------------
  void TCalm005::FillStepPointMCHistograms(StepPointMCHist_t* Hist, StepPointMC* Step) {

    static const char        oname[] = "TCalm005::FillStepPointMCHistograms";
    int                      pdg_code;
    double                   x, y, z, r, x0, y0, z0, r0;
    const mu2e::SimParticle* sim; 

    x = Step->position().x()+3904.;
    y = Step->position().y();
    z = Step->position().z();
    r = sqrt(x*x+y*y);

    //    art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle();

    sim  = Step->simParticle().operator ->();

    if (sim == NULL) {
      printf(">>> ERROR: %s sim == NULL\n",oname);
      pdg_code =  0;
      x0       =  0;
      y0       =  0;
      z0       =  0;
      r0       = -1;
    }
    else {
      x0 = sim->startPosition().x()+3904.;
      y0 = sim->startPosition().y();
      z0 = sim->startPosition().z();
      r0 = sqrt(x0*x0+y0*y0);
      pdg_code = sim->pdgId();
    }

    Hist->fRVsZ->Fill(z,r);
    Hist->fR0VsZ0->Fill(z0,r0);
    Hist->fPdgCode[0]->Fill(pdg_code);
    Hist->fPdgCode[1]->Fill(pdg_code);
    Hist->fEDep->Fill(Step->eDep());
    Hist->fTime->Fill(Step->time());
    Hist->fLength->Fill(Step->stepLength());
    
  }

//-----------------------------------------------------------------------------
  void TCalm005::FillEventHistograms(EventHist_t* Hist) {
    
    Hist->fEventNumber->Fill(fEventNumber);
    Hist->fRunNumber->  Fill(fRunNumber);
    Hist->fNSteps->Fill(fNSteps);
    
  }

//-----------------------------------------------------------------------------
  void TCalm005::FillHistograms() {
    const char oname[] = "TCalm005::FillHistograms";


    mu2e::StepPointMC*   step;
    int                  abs_pdg;
    double               dz, clen, rin[2], rout[2];
    int                  inside, pdg_code;
    double               x0, y0, z0, r0;

    art::ServiceHandle<mu2e::GeometryService> geom;
    const mu2e::DiskCalorimeter* cal(NULL);
    const mu2e::Disk* disk;
    
    if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
      mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
      cal  = dc.operator->();

      for (int i=0; i<2; i++) {
	rin[i]  = cal->disk(i).innerRadius();
	rout[i] = cal->disk(i).outerRadius();
      }
    }
    else {
      printf(">>> ERROR: disk calorimeter not found.\n");
    }
//-----------------------------------------------------------------------------
// event histograms, everything is in mm!!!
//-----------------------------------------------------------------------------
    FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// StepPointMC histograms: consider only StepPoints corresponding to particles
// originating outside the calorimeter
//-----------------------------------------------------------------------------
    for (int i=0; i<fNSteps; i++ ) {
      step = (mu2e::StepPointMC*) &fStepPointMCColl->at(i);
      //      art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle();
      const mu2e::SimParticle* sim  = step->simParticle().operator ->();

      if (sim == NULL) {
	printf(">>> ERROR: %s sim == NULL\n",oname);
	pdg_code =  0;
	x0       =  0;
	y0       =  0;
	z0       =  0;
	r0       = -1;
      }
      else {
	x0 = sim->startPosition().x()+3904.;
	y0 = sim->startPosition().y();
	z0 = sim->startPosition().z();
	r0 = sqrt(x0*x0+y0*y0);
	pdg_code = sim->pdgId();
      }

      inside = -1;
      for (int idisk=0; idisk<2; idisk++) {
	disk = &cal->disk(idisk);
	clen = cal->caloInfo().crystalHalfLength();
	dz   = z0-disk->geomInfo().origin().z();
					// add 1mm safety margins 
	if ((fabs(dz) < clen+1.) && (r0 > rin[idisk]-1.) && (r0 < rout[idisk]+1.)) {
	  inside = idisk;
	  break;
	}
      }
      
      abs_pdg = abs(pdg_code);

      FillStepPointMCHistograms(fHist.fStepPointMC[0],step);
      if       (abs_pdg <   100) FillStepPointMCHistograms(fHist.fStepPointMC[1],step);
      else if  (abs_pdg < 10000) FillStepPointMCHistograms(fHist.fStepPointMC[2],step);
      else                       FillStepPointMCHistograms(fHist.fStepPointMC[3],step);

//-----------------------------------------------------------------------------
// particle comes from outside the calorimeter
//----------------------------------------------------------------------------- 
      if (inside < 0) {
	FillStepPointMCHistograms(fHist.fStepPointMC[10],step);
	if       (abs_pdg <   100) FillStepPointMCHistograms(fHist.fStepPointMC[11],step);
	else if  (abs_pdg < 10000) FillStepPointMCHistograms(fHist.fStepPointMC[12],step);
	else                       FillStepPointMCHistograms(fHist.fStepPointMC[13],step);
      }
//-----------------------------------------------------------------------------
// particle comes from inside the calorimeter
//----------------------------------------------------------------------------- 
      else {
	FillStepPointMCHistograms(fHist.fStepPointMC[20],step);
	if       (abs_pdg <   100) FillStepPointMCHistograms(fHist.fStepPointMC[21],step);
	else if  (abs_pdg < 10000) FillStepPointMCHistograms(fHist.fStepPointMC[22],step);
	else                       FillStepPointMCHistograms(fHist.fStepPointMC[23],step);
      }
    }
  }

//-----------------------------------------------------------------------------
  void TCalm005::getData(art::Event* Evt) {
    const char* oname = "TCalm005::getData";
//-----------------------------------------------------------------------------
// generated particles
//-----------------------------------------------------------------------------
    art::Handle<mu2e::StepPointMCCollection> handle;
    //    const mu2e::StepPointMCCollection*       coll(0);

    if (fProductName[0] != 0) {
      art::Selector  selector(art::ProductInstanceNameSelector(fProductName) &&
			      art::ProcessNameSelector        (fProcessName) && 
			      art::ModuleLabelSelector        (fG4RunModuleLabel)    );
      Evt->get(selector, handle);
    }
    else {
      art::Selector  selector(art::ProcessNameSelector(fProcessName)         && 
			      art::ModuleLabelSelector(fG4RunModuleLabel)            );
      Evt->get(selector, handle);
    }
    
    if (handle.isValid()) {
      fStepPointMCColl = handle.product();
      fNSteps          = fStepPointMCColl->size();
    }
    else {
      fStepPointMCColl = 0;
      fNSteps          = 0;
      printf(">>> ERROR in %s: failed to locate collection, BAIL OUT.\n",oname);
    }
  }



//-----------------------------------------------------------------------------
  void TCalm005::Init(art::Event* Evt) {
    //  const char oname[] = "TCalm005::Init";

  }



//-----------------------------------------------------------------------------
  void TCalm005::Debug(art::Event* Evt) {
    //    const char* oname = "TCalm005::Debug";
//-----------------------------------------------------------------------------
// bit 000: calorimeter step points
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[0] != 0) {
      TModule::fDump->printStepPointMCCollection("g4run","calorimeter");
    }
  }

//-----------------------------------------------------------------------------
  bool TCalm005::filter(art::Event& Evt) {
    const char* oname = "TCalm005::filter";

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

using mu2e::TCalm005;

DEFINE_ART_MODULE(TCalm005);

#define TTrackerDose_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "murat/ana/TTrackerDose.hh"
// //-----------------------------------------------------------------------------
// // z range for different stations
// //-----------------------------------------------------------------------------
// double zplane[40][2] = {
//    8660. , 8720. ,       // 0
//    8720. , 8780. ,       // 0
//    8820. , 8880. ,       // 1
//    8880. , 8920. ,       // 1
//    8980. , 9030. ,       // 2
//    9030. , 9080. ,       // 2
//    9140. , 9180. ,       // 3
//    9180. , 9230. ,       // 3
//    9300. , 9340. ,       // 4
//    9340. , 9390. ,       // 4
//    9450. , 9500. ,       // 5
//    9500. , 9550. ,       // 5
//    0.    , -1.   ,       // 6 missing
//    0.    , -1.   ,       // 6 missing
//    9760. , 9800. ,       // 7
//    9800. , 9860. ,       // 7
//    9920. , 9960. ,       // 8
//    9960. , 10010.,       // 8
//    10080., 10120.,       // 9
//    10120., 10180.,       // 9
//    10220., 10280.,       // 10
//    10280., 10320.,       // 10
//    10380., 10440.,       // 11
//    10440., 10480.,       // 11
//    10540., 10580.,       // 12
//    10580., 10640.,       // 12
//    0.    , -1.   ,       // 13 missing
//    0.    , -1.   ,       // 13 missing
//    10850., 10900.,       // 14
//    10900., 10950.,       // 14
//    11010., 11050.,       // 15
//    11050., 11100.,       // 15
//    11160., 11210.,       // 16
//    11210., 11260.,       // 16
//    11320., 11360.,       // 17
//    11360., 11420.,       // 17
//    11480., 11520.,       // 18
//    11520., 11580.,       // 18
//    11620., 11680.,       // 19
//    11680., 11720.        // 19
// };
  

//-----------------------------------------------------------------------------
// 2017-06-14: update for nw configuration: z range for different stations
//-----------------------------------------------------------------------------
double zplane[40][2] = {
   8640. , 8660. ,       // 0
   8700. , 8710. ,       // 0
   8820. , 8830. ,       // 1
   8870. , 8880. ,       // 1
   8990. , 9010. ,       // 2
   9040. , 9060. ,       // 2
   9160. , 9180. ,       // 3  
   9220. , 9240. ,       // 3
   9340. , 9360. ,       // 4
   9390. , 9410. ,       // 4
   9520. , 9524. ,       // 5
   9571. , 9575. ,       // 5
   9694. , 9698. ,       // 6 
   9745. , 9749. ,       // 6 
   9868. , 9872. ,       // 7
   9919. , 9923. ,       // 7
   10042., 10046.,       // 8
   10093., 10097.,       // 8
   10216., 10220.,       // 9
   10267., 10271.,       // 9
   10390., 10394.,       // 10
   10441., 10445.,       // 10
   10564., 10568.,       // 11
   10615., 10619.,       // 11
   10738., 10742.,       // 12
   10789., 10793.,       // 12
   10912., 10916.,       // 13
   10963., 10967.,       // 13
   11086., 11090.,       // 14
   11137., 11141.,       // 14
   11260., 11264.,       // 15
   11311., 11315.,       // 15
   11434., 11348.,       // 16
   11485., 11489.,       // 16
   11608., 11612.,       // 17
   11659., 11663.,       // 17
   -1.,    -1.,
   -1.,    -1.,
   -1.,    -1.,
   -1.,    -1.
};
  
//-----------------------------------------------------------------------------
int plane_number(float Z) {
  int ist(-1);

  for (int i=0; i<40; i++) {
    if ((Z > zplane[i][0]) && (Z < zplane[i][1])) {
      ist = i;
      break;
    }
  }

  return ist;
}

//-----------------------------------------------------------------------------
// estimate radiation dose absorbed by the tracker electronics due to a given 'Process'
// Process = "FLASH","DIO" etc
//-----------------------------------------------------------------------------
TTrackerDose::TTrackerDose(const char* Process) : TStnModule("TrackerDose","TrackerDose") {

  fProcess = Process;

  fNPOT = 1.2e20;           // "per year of running"
  fDensity[1] = 1.8;        // G10
  fDensity[0] = 0  ;        // make sure the code crashes 

  TH1::AddDirectory(0);

  InitChain();

  BookHistograms();
}

//-----------------------------------------------------------------------------
int TTrackerDose::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fNumber ,"evnumber",Form("%s: nhits "    ,Folder),  200, 0, 1e8,Folder);
  HBook1F(Hist->fNVDHits,"nvdhits" ,Form("%s: N(VD hits)",Folder),  250, 0, 500,Folder);
  HBook1F(Hist->fEDepTot,"edeptot" ,Form("%s: E(Dep)tot ",Folder), 1000, 0, 10.,Folder);
  HBook1F(Hist->fETotVDet13,"etot_vdet13" ,Form("%s: ETot VDet13", Folder), 1000, 0, 50.,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int TTrackerDose::BookTrackerHistograms(TrackHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fNHits   ,"nhits" ,Form("%s: nhits ",Folder), 200,   0  ,200,Folder);
  HBook1F(Hist->fPdgID[0],"pdg_0" ,Form("%s: PDG ID",Folder), 200, -5000, 5000,Folder);
  HBook1F(Hist->fPdgID[1],"pdg_1" ,Form("%s: PDG ID",Folder), 500, -250,  250,Folder);
  HBook1F(Hist->fZ       ,"z"     ,Form("%s: Z"     ,Folder), 4000,  8000, 12000,Folder);
  HBook1F(Hist->fR       ,"radius",Form("%s: radius",Folder), 200,   650,   850,Folder);

  HBook1F(Hist->fEKin[0] ,"ekin_e",Form("%s: E(kin)[e]",Folder), 200,     0,   100,Folder);
  HBook1F(Hist->fEKin[1] ,"ekin_g",Form("%s: E(kin)[g]",Folder), 200,     0,   100,Folder);
  HBook1F(Hist->fEKin[2] ,"ekin_p",Form("%s: E(kin)[p]",Folder), 200,     0,   100,Folder);
  HBook1F(Hist->fEKin[3] ,"ekin_n",Form("%s: E(kin)[n]",Folder), 200,     0,   100,Folder);
  HBook1F(Hist->fEKin[4] ,"ekin_m",Form("%s: E(kin)[m]",Folder), 200,     0,   100,Folder);

  HBook2F(Hist->fEDepVsPlane[0],"edep_vs_plane_0",Form("%s: E(dep) vs Plane Raw ",Folder), 40,0,40,20,650,850,Folder);
  HBook2F(Hist->fEDepVsPlane[1],"edep_vs_plane_1",Form("%s: E(dep) vs Plane Norm",Folder), 40,0,40,20,650,850,Folder);
  HBook2F(Hist->fEDepVsPlane[2],"edep_vs_plane_2",Form("%s: E(dep) vs Plane ele" ,Folder), 40,0,40,20,650,850,Folder);
  HBook2F(Hist->fEDepVsPlane[3],"edep_vs_plane_3",Form("%s: E(dep) vs Plane gam" ,Folder), 40,0,40,20,650,850,Folder);

  HBook1D(Hist->fMeanDoseVsPlane,"dmean_vs_plane",Form("%s: Mean Dose Vs Plane"   ,Folder), 40,0,40,Folder);
  HBook1D(Hist->fMaxDoseVsPlane ,"dmax_vs_plane" ,Form("%s: Dose Vs Plane 71<r<72",Folder), 40,0,40,Folder);
  HBook1D(Hist->fDose73VsPlane  ,"d73_vs_plane"  ,Form("%s: Dose Vs Plane 73<r<74",Folder), 40,0,40,Folder);
  HBook1D(Hist->fDose75VsPlane  ,"d75_vs_plane"  ,Form("%s: Dose Vs Plane 75<r<76",Folder), 40,0,40,Folder);

  for (int i=0;i<40; i++) {
    int ipl = i+1;
    HBook1D(Hist->fDoseVsR[i],Form("dose_vs_r_%02i",ipl),Form("%s: Dose Vs R Plane %02i",Folder,ipl), 20,650,850,Folder); // 
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TTrackerDose::BookVDetHistograms(VDetHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fPdgId   ,"pdg"   ,Form("%s: PDG ID",Folder), 500, -5000, 5000,Folder);
  HBook1F(Hist->fE[0]    ,"e_0"   ,Form("%s: e[0]"  ,Folder), 200,     0,   100,Folder);
  HBook1F(Hist->fE[1]    ,"e_1"   ,Form("%s: e[1]"  ,Folder), 200,     0,     2,Folder);
  HBook1F(Hist->fR       ,"radius",Form("%s: radius",Folder), 850,     0,   850,Folder);
  HBook1F(Hist->fCosth   ,"costh" ,Form("%s: costh" ,Folder), 200,    -1,     1,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int TTrackerDose::BookHistograms() {

  TFolder*    fol;
  TFolder*    hist_folder;
  char        folder_name[200];
 
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

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

  int book_tracker_histset[kMaxTrackHistSets];
  for (int i=0; i<kMaxTrackHistSets; i++) book_tracker_histset[i] = 0;

  book_tracker_histset[  0] = 1;     		// all hits

  for (int i=0; i<kMaxTrackHistSets; i++) {
    if (book_tracker_histset[i] != 0) {
      sprintf(folder_name,"up_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fUp[i] = new TrackHist_t;
      BookTrackerHistograms(fHist.fUp[i],Form("Hist/%s",folder_name));

      sprintf(folder_name,"dn_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fDn[i] = new TrackHist_t;
      BookTrackerHistograms(fHist.fDn[i],Form("Hist/%s",folder_name));
    }
  }

  int book_vdet_histset[kMaxVDetHistSets];
  for (int i=0; i<kMaxVDetHistSets; i++) book_vdet_histset[i] = 0;

  book_vdet_histset[  0] = 1;     		// all hits
  book_vdet_histset[  1] = 1;     		// all    hits in the tracker front detector 71 < R < 80
  book_vdet_histset[  2] = 1;     		// electr hits in the tracker front detector 71 < R < 80
  book_vdet_histset[  3] = 1;     		// photon hits in the tracker front detector 71 < R < 80
  book_vdet_histset[  4] = 1;     		// other  hits in the tracker front detector 71 < R < 80
  book_vdet_histset[  5] = 1;     		// all    hits in the tracker front detector 40 < R < 80

  book_vdet_histset[100] = 1;     		// all hits - tracker mid

  book_vdet_histset[200] = 1;     		// all hits - tracker back 

  for (int i=0; i<kMaxVDetHistSets; i++) {
    if (book_vdet_histset[i] != 0) {
      sprintf(folder_name,"vdet_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fVDet[i] = new VDetHist_t;
      BookVDetHistograms(fHist.fVDet[i],Form("Hist/%s",folder_name));
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
int TTrackerDose::InitEvent() {

  fEDepTot = 0;

  for (int i=0; i<nTtsUpHits; i++) {
    fEDepTot += ttsUpEdep[i];
  }

  for (int i=0; i<nTtsDwHits; i++) {
    fEDepTot += ttsDwEdep[i];
  }

  fETotVDet13  = 0;
  
  for (int i=0; i<vNHits; i++) {
    if (vId[i] == 13) {
      if ((vRadius[i] > 400.) && (vRadius[i] < 800.)) {
	fETotVDet13 += vEKin[i];
      }
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TTrackerDose::FillEventHistograms(EventHist_t* Hist) {

  Hist->fNumber->Fill(evt);
  Hist->fNVDHits->Fill(vNHits);
  Hist->fEDepTot->Fill(fEDepTot);
  Hist->fETotVDet13->Fill(fETotVDet13);

  return 0;
}

//-----------------------------------------------------------------------------
// by placing fEDepTot calculation here, implicitly assume that 
// FillTrackerHistograms is called just once. 
// in the calculation, however, do not assume that 
//-----------------------------------------------------------------------------
int TTrackerDose::FillTrackerHistograms(TrackHist_t* Hup, TrackHist_t* Hdn) {

  Hup->fNHits->Fill(nTtsUpHits);

  for (int i=0; i<nTtsUpHits; i++) {

    int pdg_id = ttsUpPdgId[i];
    Hup->fPdgID[0]->Fill(pdg_id);
    Hup->fPdgID[1]->Fill(pdg_id);

    float z = ttsUpZ[i];
    Hup->fZ    ->Fill(z);

    float x = ttsUpX[i];
    float y = ttsUpY[i];
    float r = sqrt(x*x+y*y);

    Hup->fR    ->Fill(r);

    float ekin = ttsUpEkin[i];
    if (abs(pdg_id) ==   11) Hup->fEKin[0] ->Fill(ekin);
    if (abs(pdg_id) ==   22) Hup->fEKin[1] ->Fill(ekin);
    if (abs(pdg_id) == 2212) Hup->fEKin[2] ->Fill(ekin);
    if (abs(pdg_id) == 2112) Hup->fEKin[3] ->Fill(ekin);
    if (abs(pdg_id) ==   13) Hup->fEKin[4] ->Fill(ekin);

    int plane  = plane_number(z);
    float edep = ttsUpEdep[i];

    Hup->fEDepVsPlane[0] ->Fill(plane,r,edep);

    if (abs(pdg_id) ==   11) Hup->fEDepVsPlane[2] ->Fill(plane,r,edep);
    if (abs(pdg_id) ==   22) Hup->fEDepVsPlane[3] ->Fill(plane,r,edep);
  }

  Hdn->fNHits->Fill(nTtsDwHits);

  for (int i=0; i<nTtsDwHits; i++) {
    int pdg_id = ttsDwPdgId[i];

    Hdn->fPdgID[0]->Fill(pdg_id);
    Hdn->fPdgID[1]->Fill(pdg_id);

    float z = ttsDwZ[i];
    Hdn->fZ    ->Fill(z);

    float x = ttsDwX[i];
    float y = ttsDwY[i];
    float r = sqrt(x*x+y*y);

    Hdn->fR    ->Fill(r);

    float ekin = ttsDwEkin[i];
    if (abs(pdg_id) ==   11) Hdn->fEKin[0] ->Fill(ekin);
    if (abs(pdg_id) ==   22) Hdn->fEKin[1] ->Fill(ekin);
    if (abs(pdg_id) == 2212) Hdn->fEKin[2] ->Fill(ekin);
    if (abs(pdg_id) == 2112) Hdn->fEKin[3] ->Fill(ekin);
    if (abs(pdg_id) ==   13) Hdn->fEKin[4] ->Fill(ekin);

    int plane  = plane_number(z);
    float edep = ttsDwEdep[i];

    Hdn->fEDepVsPlane[0] ->Fill(plane,r,edep);

    if (abs(pdg_id) ==   11)  Hdn->fEDepVsPlane[2] ->Fill(plane,r,edep);
    if (abs(pdg_id) ==   22)  Hdn->fEDepVsPlane[3] ->Fill(plane,r,edep);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int TTrackerDose::FillVDetHistograms(VDetHist_t* Hist, VDetData_t* VDet) {

  Hist->fPdgId->Fill(VDet->PdgId);
  Hist->fE[0] ->Fill(VDet->E);
  Hist->fE[1] ->Fill(VDet->E);
  Hist->fR    ->Fill(VDet->Radius);
  Hist->fCosth->Fill(VDet->Costh);

  return 0;
}
//-----------------------------------------------------------------------------
int TTrackerDose::FillHistograms() {

//-----------------------------------------------------------------------------
// 1. fill event histograms 
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// 2. fill tracker histograms 
//-----------------------------------------------------------------------------
  FillTrackerHistograms(fHist.fUp[0], fHist.fDn[0]);
//-----------------------------------------------------------------------------
// 3. fill virtual detecotr histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<vNHits; i++) {
    fVDet.Id     = vId[i];
    fVDet.PdgId  = vPdgId[i];
    fVDet.P      = vP[i];   
    fVDet.Px     = vPx[i];  
    fVDet.Py     = vPy[i];  
    fVDet.Pz     = vPz[i];  
    fVDet.E      = vE[i];  
    fVDet.EKin   = vEKin[i];
    fVDet.M      = vM[i];   
    fVDet.T      = vT[i];   
    fVDet.X      = vX[i];   
    fVDet.Y      = vY[i];   
    fVDet.Z      = vZ[i];   
    fVDet.Costh  = vCosth[i];   
    fVDet.Radius = vRadius[i];   

    if (fVDet.Id == 13) {
      FillVDetHistograms(fHist.fVDet[0],&fVDet);
      
      if ((fVDet.Radius > 710) && (fVDet.Radius < 800)) {
	FillVDetHistograms(fHist.fVDet[1],&fVDet);

	if      (abs(fVDet.PdgId) == 11) FillVDetHistograms(fHist.fVDet[2],&fVDet);
	else if (abs(fVDet.PdgId) == 22) FillVDetHistograms(fHist.fVDet[3],&fVDet);
	else                             FillVDetHistograms(fHist.fVDet[4],&fVDet);

      }
      if ((fVDet.Radius > 400) && (fVDet.Radius < 800)) {
	FillVDetHistograms(fHist.fVDet[5],&fVDet);
      }
    }

    if (fVDet.Id == 11) {
      FillVDetHistograms(fHist.fVDet[100],&fVDet);
    }

    if (fVDet.Id == 15) {
      FillVDetHistograms(fHist.fVDet[200],&fVDet);
    }
  }


  return 0;
}

//-----------------------------------------------------------------------------
int TTrackerDose::ResetTrackerHistograms(TrackHist_t* Hup, TrackHist_t* Hdn) {
  Hup->fNHits   ->Reset();
  Hup->fPdgID[0]->Reset();
  Hup->fPdgID[1]->Reset();
  Hup->fZ       ->Reset();
  Hup->fR       ->Reset();
  Hup->fEKin[0] ->Reset();
  Hup->fEKin[1] ->Reset();
  Hup->fEKin[2] ->Reset();
  Hup->fEKin[3] ->Reset();
  Hup->fEKin[4] ->Reset();

  Hup->fEDepVsPlane[0]->Reset();
  Hup->fEDepVsPlane[1]->Reset();

  Hdn->fNHits   ->Reset();
  Hdn->fPdgID[0]->Reset();
  Hdn->fPdgID[1]->Reset();
  Hdn->fZ       ->Reset();
  Hdn->fR       ->Reset();
  Hdn->fEKin[0] ->Reset();
  Hdn->fEKin[1] ->Reset();
  Hdn->fEKin[2] ->Reset();
  Hdn->fEKin[3] ->Reset();
  Hdn->fEKin[4] ->Reset();

  Hdn->fEDepVsPlane[0]->Reset();
  Hdn->fEDepVsPlane[1]->Reset();

  return 0;
}
//-----------------------------------------------------------------------------
int TTrackerDose::ResetHistograms() {

  ResetTrackerHistograms(fHist.fUp[0], fHist.fDn[0]);

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
void TTrackerDose::Loop(Long64_t NEvents) {

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
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    InitEvent();

    FillHistograms();

    nev += 1;
  }

  printf("processed %li events\n",nev);

  fNEvents = nev;
//-----------------------------------------------------------------------------
// remaining part : normalize histograms to the total N(POT)
// last term, 1.e6, accounts for energies being measured in MeV
// material used to calculate losses: 3mm thick G10 disks 
//-----------------------------------------------------------------------------
  float thickness     = 0.3;
  float mev_per_joule = 1.6e-19*1.e6;
  float krad_per_gray = 10.;

  int nx   =  fHist.fUp[0]->fEDepVsPlane[1]->GetNbinsX();
  int ny   =  fHist.fUp[0]->fEDepVsPlane[1]->GetNbinsY();
  float dr =  fHist.fUp[0]->fEDepVsPlane[1]->GetYaxis()->GetBinWidth(1)/10.; //  // convert to cm

//-----------------------------------------------------------------------------
// "upstream" layer of "electronics"
//-----------------------------------------------------------------------------
  for (int i=0; i<nx; i++) {
    int binx = i+1;
    for (int ir=0; ir<ny; ir++) {

      int biny = ir+1;

      float r     = fHist.fUp[0]->fEDepVsPlane[1]->GetYaxis()->GetBinCenter(biny)/10.; // convert to cm
      float mass  = 2*M_PI*r*dr*thickness*fDensity[0]/1.e3; // in kG

      // the first and the last 1cm bins are effectively more narrow, not 1 cm, but 1-0.27 cm wide
      if ((r > 71.) && (r < 72.)) mass = mass*(1-0.27);
      if ((r > 79.) && (r < 80.)) mass = mass*(1-0.27);

      float sf    = (fNPOT*fNPerPOT)/(nent+1.e-12)*(nentries/(fNSimulated+1.e-12))*mev_per_joule/mass/krad_per_gray;
     
      double x    = fHist.fUp[0]->fEDepVsPlane[0]->GetBinContent(binx,biny);
      double e    = fHist.fUp[0]->fEDepVsPlane[0]->GetBinError  (binx,biny);

      fHist.fUp[0]->fEDepVsPlane[1]->SetBinContent(binx,biny,x*sf);
      fHist.fUp[0]->fEDepVsPlane[1]->SetBinError  (binx,biny,e*sf);
    }
  }

  float rmax = 79.73;
  float rmin = 71.27;

  float mass_up = M_PI*(rmax*rmax-rmin*rmin)*thickness*fDensity[0]/1.e3; // in kG
  float sf_up   = (fNPOT*fNPerPOT)/(nent+1.e-12)*(nentries/(fNSimulated+1.e-12))*mev_per_joule/mass_up/krad_per_gray;

  fHist.fUp[0]->fMeanDoseVsPlane->Add(fHist.fUp[0]->fEDepVsPlane[0]->ProjectionX());
  fHist.fUp[0]->fMeanDoseVsPlane->Scale(sf_up); 

  fHist.fUp[0]->fMaxDoseVsPlane->Add(fHist.fUp[0]->fEDepVsPlane[1]->ProjectionX("px07_up", 7, 7));
  fHist.fUp[0]->fDose73VsPlane ->Add(fHist.fUp[0]->fEDepVsPlane[1]->ProjectionX("px09_up", 9, 9));
  fHist.fUp[0]->fDose75VsPlane ->Add(fHist.fUp[0]->fEDepVsPlane[1]->ProjectionX("px11_up",11,11));

  for (int i=0; i<40; i++) {
    int ipl = i+1;
    fHist.fUp[0]->fDoseVsR[i]->Add(fHist.fUp[0]->fEDepVsPlane[1]->ProjectionY(Form("py%2i_up",ipl),ipl,ipl));
  }
//-----------------------------------------------------------------------------
// "downstream" layer of "electronics"
//-----------------------------------------------------------------------------
  for (int i=0; i<nx; i++) {
    int binx = i+1;
    for (int ir=0; ir<ny; ir++) {
      int biny    = ir+1;
      float r     = fHist.fUp[0]->fEDepVsPlane[1]->GetYaxis()->GetBinCenter(biny)/10.;  // convert to cm
      float mass  = 2*M_PI*r*dr*thickness*fDensity[1]/1.e3; // in kG

      // the first and the last 1cm bins are effectively more narrow, not 1 cm, but 1-0.27 cm wide
      if ((r > 71.) && (r < 72.)) mass = mass*(1-0.27);
      if ((r > 79.) && (r < 80.)) mass = mass*(1-0.27);


      float sf    = (fNPOT*fNPerPOT)/(nent+1.e-12)*(nentries/(fNSimulated+1.e-12))*mev_per_joule/mass/krad_per_gray;
      
      double x    = fHist.fDn[0]->fEDepVsPlane[0]->GetBinContent(binx,biny);
      double e    = fHist.fDn[0]->fEDepVsPlane[0]->GetBinError  (binx,biny);

      fHist.fDn[0]->fEDepVsPlane[1]->SetBinContent(binx,biny,x*sf);
      fHist.fDn[0]->fEDepVsPlane[1]->SetBinError  (binx,biny,e*sf);
    }
  }

  fHist.fDn[0]->fMeanDoseVsPlane->Add(fHist.fDn[0]->fEDepVsPlane[0]->ProjectionX());

  float mass_dn = M_PI*(rmax*rmax-rmin*rmin)*thickness*fDensity[1]/1.e3; // in kG
  float sf_dn   = (fNPOT*fNPerPOT)/(nent+1.e-12)*(nentries/(fNSimulated+1.e-12))*mev_per_joule/mass_dn/krad_per_gray;

  fHist.fDn[0]->fMeanDoseVsPlane->Scale(sf_dn); // full range: (71.3 - 79.7 cm)

  fHist.fDn[0]->fMaxDoseVsPlane->Add(fHist.fDn[0]->fEDepVsPlane[1]->ProjectionX("px07_dn", 7, 7));
  fHist.fDn[0]->fDose73VsPlane ->Add(fHist.fDn[0]->fEDepVsPlane[1]->ProjectionX("px09_dn", 9, 9));
  fHist.fDn[0]->fDose75VsPlane ->Add(fHist.fDn[0]->fEDepVsPlane[1]->ProjectionX("px11_dn",11,11));

  for (int i=0; i<40; i++) {
    int ipl = i+1;
    fHist.fDn[0]->fDoseVsR[i]->Add(fHist.fDn[0]->fEDepVsPlane[1]->ProjectionY(Form("py%02i_dn",ipl),ipl,ipl));
  }
}

//-----------------------------------------------------------------------------
TTrackerDose::~TTrackerDose() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


//-----------------------------------------------------------------------------
Int_t TTrackerDose::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


//-----------------------------------------------------------------------------
Long64_t TTrackerDose::LoadTree(Long64_t entry) {
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
void TTrackerDose::Init(TChain* Chain) {
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
   fChain->SetBranchAddress("crystalR", &crystalR, &b_crystalR);
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
   fChain->SetBranchAddress("crystalROSectionId", &crystalROSectionId, &b_crystalROSectionId);
   fChain->SetBranchAddress("crystalROCrystalId", &crystalROCrystalId, &b_crystalROCrystalId);
   fChain->SetBranchAddress("crystalROEdep", &crystalROEdep, &b_crystalROEdep);
   fChain->SetBranchAddress("crystalRODose", &crystalRODose, &b_crystalRODose);
   fChain->SetBranchAddress("crystalROX", &crystalROX, &b_crystalROX);
   fChain->SetBranchAddress("crystalROY", &crystalROY, &b_crystalROY);
   fChain->SetBranchAddress("crystalROZ", &crystalROZ, &b_crystalROZ);
   fChain->SetBranchAddress("crystalROR", &crystalROR, &b_crystalROR);
   fChain->SetBranchAddress("nCrystalROCard", &nCrystalROCard, &b_nCrystalROCard);
   fChain->SetBranchAddress("crystalROCardSectionId", &crystalROCardSectionId, &b_crystalROCardSectionId);
   fChain->SetBranchAddress("crystalROCardCrystalId", &crystalROCardCrystalId, &b_crystalROCardCrystalId);
   fChain->SetBranchAddress("crystalROCardEdep", &crystalROCardEdep, &b_crystalROCardEdep);
   fChain->SetBranchAddress("crystalROCardDose", &crystalROCardDose, &b_crystalROCardDose);
   fChain->SetBranchAddress("crystalROCardX", &crystalROCardX, &b_crystalROCardX);
   fChain->SetBranchAddress("crystalROCardY", &crystalROCardY, &b_crystalROCardY);
   fChain->SetBranchAddress("crystalROCardZ", &crystalROCardZ, &b_crystalROCardZ);
   fChain->SetBranchAddress("crystalROCardR", &crystalROCardR, &b_crystalROCardR);
   fChain->SetBranchAddress("nCrateHits", &nCrateHits, &b_nCrateHits);
   fChain->SetBranchAddress("crateEdep", &crateEdep, &b_crateEdep);
   fChain->SetBranchAddress("crateX", &crateX, &b_crateX);
   fChain->SetBranchAddress("crateY", &crateY, &b_crateY);
   fChain->SetBranchAddress("crateZ", &crateZ, &b_crateZ);
   fChain->SetBranchAddress("crateR", &crateR, &b_crateR);
   fChain->SetBranchAddress("cratePdgId", &cratePdgId, &b_cratePdgId);
   fChain->SetBranchAddress("crateE", &crateE, &b_crateE);
   fChain->SetBranchAddress("crateEkin", &crateEkin, &b_crateEkin);
   fChain->SetBranchAddress("crateMass", &crateMass, &b_crateMass);
   fChain->SetBranchAddress("crateL", &crateL, &b_crateL);
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
   fChain->SetBranchAddress("nTtsUpHits", &nTtsUpHits, &b_nTtsUpHits);
   fChain->SetBranchAddress("ttsUpEdep", ttsUpEdep, &b_ttsUpEdep);
   fChain->SetBranchAddress("ttsUpX", ttsUpX, &b_ttsUpX);
   fChain->SetBranchAddress("ttsUpY", ttsUpY, &b_ttsUpY);
   fChain->SetBranchAddress("ttsUpZ", ttsUpZ, &b_ttsUpZ);
   fChain->SetBranchAddress("ttsUpR", ttsUpR, &b_ttsUpR);
   fChain->SetBranchAddress("ttsUpPdgId", ttsUpPdgId, &b_ttsUpPdgId);
   fChain->SetBranchAddress("ttsUpEkin", ttsUpEkin, &b_ttsUpEkin);
   fChain->SetBranchAddress("ttsUpMass", ttsUpMass, &b_ttsUpMass);
   fChain->SetBranchAddress("ttsUpL", ttsUpL, &b_ttsUpL);
   fChain->SetBranchAddress("nTtsDwHits", &nTtsDwHits, &b_nTtsDwHits);
   fChain->SetBranchAddress("ttsDwEdep", ttsDwEdep, &b_ttsDwEdep);
   fChain->SetBranchAddress("ttsDwX", ttsDwX, &b_ttsDwX);
   fChain->SetBranchAddress("ttsDwY", ttsDwY, &b_ttsDwY);
   fChain->SetBranchAddress("ttsDwZ", ttsDwZ, &b_ttsDwZ);
   fChain->SetBranchAddress("ttsDwR", ttsDwR, &b_ttsDwR);
   fChain->SetBranchAddress("ttsDwPdgId", ttsDwPdgId, &b_ttsDwPdgId);
   fChain->SetBranchAddress("ttsDwEkin", ttsDwEkin, &b_ttsDwEkin);
   fChain->SetBranchAddress("ttsDwMass", ttsDwMass, &b_ttsDwMass);
   fChain->SetBranchAddress("ttsDwL", ttsDwL, &b_ttsDwL);
   Notify();
}

//-----------------------------------------------------------------------------
Bool_t TTrackerDose::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

//-----------------------------------------------------------------------------
void TTrackerDose::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

//-----------------------------------------------------------------------------
int TTrackerDose::AddFiles(const char* Fn) {

  char   name[200], c[1000];

  FILE* f = fopen(Fn,"r");
  if (f) {
    while ( (c[0]=getc(f)) != EOF) {
//-----------------------------------------------------------------------------
// skip comment lines starting from '#'
//-----------------------------------------------------------------------------
      if ((c[0] != '#') && (c[0] != '\n')) {
	ungetc(c[0],f);
	fscanf(f,"%s",name);
	printf("file : %s\n",name);
	fChain->Add(name);
      }
      fgets(c,100,f);
    }
    fclose(f);
    return 0;
  }
  else {				// can't open the file
    printf(" CANT open %s, BAIL OUT\n",Fn);
    return -1;
  }
  
}

//-----------------------------------------------------------------------------
// fNPetPOT for non-flash backgrounds is a product of the number of stopped 
// muons per POT and average yiled of a given particle species per muon stop
//-----------------------------------------------------------------------------
int TTrackerDose::InitChain() {

  fChain = new TChain("//calorimeterDose/Calo","Calo");

  if (fProcess == "FLASH_C360BRASS_01") { // read catalog file
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-c360brass-01-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9 * 5./51.;
    fDensity[0] = 1.7;
  }
  if (fProcess == "622_0000") {
    fDataset = Form("%s/ts3_tooth/g4s4_flash_%s/ts3_tooth.g4s4_flash_%s.cal_dose",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 15.6e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "g4s3_622_0001") {
    fDataset = Form("%s/ts3_tooth/%s/ts3_tooth.%s.mothers.g4s4_flash.cal_dose",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 21.2e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "622_0002") {
    fDataset = Form("%s/ts3_tooth/cal_dose_%s/ts3_tooth.cal_dose_%s",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 50.e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "622_0003") {
    fDataset = Form("%s/ts3_tooth/cal_dose_%s/ts3_tooth.cal_dose_%s",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 50.e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "622_0004") {
    fDataset = Form("%s/ts3_tooth/cal_dose_%s/ts3_tooth.cal_dose_%s",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 1.e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "622_0005") {
    fDataset = Form("%s/ts3_tooth/%s/cal_dose/ts3_tooth.%s_cal_dose",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 63.2e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "622_0006") {
    fDataset = Form("%s/ts3_tooth/%s/cal_dose/ts3_tooth.%s_cal_dose",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 99.6e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "622_0007") {
    fDataset = Form("%s/ts3_tooth/%s/cal_dose/ts3_tooth.%s_cal_dose",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 100e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "622_0008") {
    fDataset = Form("%s/ts3_tooth/%s/cal_dose/ts3_tooth.%s_cal_dose",getenv("CATALOG_DIR"),fProcess.Data(),fProcess.Data());
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 100e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "FLASH_C360BRASS_02") { // read catalog file
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-c360brass-02-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9 * 5./51.;
    fDensity[0] = 1.7;
  }
  if (fProcess == "FLASH_C360BRASS_03") { // read catalog file
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-c360brass-03-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9 * 5./51.;
    fDensity[0] = 1.7;
  }
  if (fProcess == "FLASH_C360BRASS_04") { // read catalog file
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-c360brass-04-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9 * 5./51.;
    fDensity[0] = 1.7;
  }
  if (fProcess == "FLASH_C360BRASS_05") { // read catalog file
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-c360brass-05-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9 * 5./51.;
    fDensity[0] = 1.7;
  }
  if (fProcess == "FLASH_C360BRASS_06") { // read catalog file
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-c360brass-06-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9 * 5./51.;
    fDensity[0] = 1.7;
  }
  if (fProcess == "FLASH") {
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-default-giani-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9;
    fDensity[0] = 1.7;
  }
  if (fProcess == "FLASH_V2") {
    fDataset = Form("%s/datasets/ts3-tooth/flash-v2-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 8e6;
    fDensity[0] = 1.7;
  }
  if (fProcess == "FLASH_PASHA") {
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-default-pasha-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_G4V10") {
    fDataset = Form("%s/datasets/tracker-rad-dose/flash-g4v10-giani-nt",getenv("MU2E_BASE_RELEASE"));
    AddFiles(fDataset.Data());

    fNPerPOT    = 1.;
    fNSimulated = 5.1e9;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_CU050") { // updated 2017-04-11
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748704/00/00000/nts.gianipez.bbb.g4v10.001002_00170003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748704/00/00001/nts.gianipez.bbb.g4v10.001002_00210003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748704/00/00002/nts.gianipez.bbb.g4v10.001002_00230003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748704/00/00003/nts.gianipez.bbb.g4v10.001002_00380101.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748704/00/00004/nts.gianipez.bbb.g4v10.001002_00410125.root");
    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 8.96*0.50;
  }
  else if (fProcess == "FLASH_CU100") {
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18671676/00/00000/nts.gianipez.bbb.g4v10.001002_00170003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18671676/00/00001/nts.gianipez.bbb.g4v10.001002_00210003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18671676/00/00002/nts.gianipez.bbb.g4v10.001002_00230003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18671676/00/00003/nts.gianipez.bbb.g4v10.001002_00380101.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18671676/00/00004/nts.gianipez.bbb.g4v10.001002_00410125.root");
    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 8.96*1.00;
  }
  else if (fProcess == "FLASH_CU100_MURAT") {
    // stiffener ring out of Cu (x1), the rest - default
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/cu100/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00170003.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/cu100/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00210003.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/cu100/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00230003.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/cu100/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00380101.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/cu100/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00410125.root");
    fNPerPOT    = 1.;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_CU215") {  // updated 2017-04-11
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748681/00/00000/nts.gianipez.bbb.g4v10.001002_00170003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748681/00/00001/nts.gianipez.bbb.g4v10.001002_00210003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748681/00/00002/nts.gianipez.bbb.g4v10.001002_00230003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748681/00/00003/nts.gianipez.bbb.g4v10.001002_00380101.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18748681/00/00004/nts.gianipez.bbb.g4v10.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 8.96*2.15;
  }
  else if (fProcess == "FLASH_CU100_680") {  // updated 2017-04-11
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18817868/00/00000/nts.gianipez.bbb.g4v10innerRing680.001002_00170003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18817868/00/00001/nts.gianipez.bbb.g4v10innerRing680.001002_00210003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18817868/00/00002/nts.gianipez.bbb.g4v10innerRing680.001002_00230003.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18817868/00/00003/nts.gianipez.bbb.g4v10innerRing680.001002_00380101.root");
    fChain->Add("/pnfs/mu2e/scratch/users/gianipez/workflow/tracker-dose-ana-g4v10/outstage/18817868/00/00004/nts.gianipez.bbb.g4v10innerRing680.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 8.96*1.00;
  }
  else if (fProcess == "FLASH_CU100_RING_ONLY") {  // updated 2017-04-21
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10/nts.gianipez.bbb.g4v10innerRingCu10.001002_00170003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10/nts.gianipez.bbb.g4v10innerRingCu10.001002_00210003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10/nts.gianipez.bbb.g4v10innerRingCu10.001002_00230003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10/nts.gianipez.bbb.g4v10innerRingCu10.001002_00380101.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10/nts.gianipez.bbb.g4v10innerRingCu10.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_CU100_IR_CU") {  // updated 2017-04-23
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10-ir-cu/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00170003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10-ir-cu/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00210003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10-ir-cu/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00230003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10-ir-cu/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00380101.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/sr-cu100-700-ele-g10-ir-cu/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_SMCO_X1") {  // updated 2017-04-27
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_1/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00170003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_1/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00210003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_1/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00230003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_1/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00380101.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_1/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_SMCO_X2") {  // updated 2017-04-27
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_2/step2/nts.MU2EGRIDDSOWNER.flash-sf-2.MU2EGRIDDSCONF.001002_00170003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_2/step2/nts.MU2EGRIDDSOWNER.flash-sf-2.MU2EGRIDDSCONF.001002_00210003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_2/step2/nts.MU2EGRIDDSOWNER.flash-sf-2.MU2EGRIDDSCONF.001002_00230003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_2/step2/nts.MU2EGRIDDSOWNER.flash-sf-2.MU2EGRIDDSCONF.001002_00380101.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_2/step2/nts.MU2EGRIDDSOWNER.flash-sf-2.MU2EGRIDDSCONF.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_SMCO_X3") {  // updated 2017-04-27
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_3/step2/nts.MU2EGRIDDSOWNER.flash-sf-3.MU2EGRIDDSCONF.001002_00170003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_3/step2/nts.MU2EGRIDDSOWNER.flash-sf-3.MU2EGRIDDSCONF.001002_00210003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_3/step2/nts.MU2EGRIDDSOWNER.flash-sf-3.MU2EGRIDDSCONF.001002_00230003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_3/step2/nts.MU2EGRIDDSOWNER.flash-sf-3.MU2EGRIDDSCONF.001002_00380101.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/smart_copper_sf_3/step2/nts.MU2EGRIDDSOWNER.flash-sf-3.MU2EGRIDDSCONF.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_DEFAULT") {  // updated 2017-05-03
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/default/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00170003.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/default/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00210003.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/default/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00230003.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/default/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00380101.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/default/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_LEAD_2MM") {  // updated 2017-04-29
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/lead_2mm/step2/nts.MU2EGRIDDSOWNER.flash_lead_2mm.MU2EGRIDDSCONF.001002_00170003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/lead_2mm/step2/nts.MU2EGRIDDSOWNER.flash_lead_2mm.MU2EGRIDDSCONF.001002_00210003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/lead_2mm/step2/nts.MU2EGRIDDSOWNER.flash_lead_2mm.MU2EGRIDDSCONF.001002_00230003.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/lead_2mm/step2/nts.MU2EGRIDDSOWNER.flash_lead_2mm.MU2EGRIDDSCONF.001002_00380101.root");
    fChain->Add("/mu2e/data/users/gianipez/g4v10-flash-tracker/lead_2mm/step2/nts.MU2EGRIDDSOWNER.flash_lead_2mm.MU2EGRIDDSCONF.001002_00410125.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*5/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "FLASH_LEAD_2MM_V1") {  // updated 2017-05-03
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/flash-lead-2mm-v1/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00230003.root");
    fChain->Add("/mu2e/data/users/murat/datasets/tracker-rad-dose/flash-lead-2mm-v1/step2/nts.MU2EGRIDDSOWNER.bbb.MU2EGRIDDSCONF.001002_00380101.root");

    fNPerPOT    = 1.;
    //    fNSimulated = 5.1e9;
    fNSimulated = 5.1e9*2/51.;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "DIO") {
    fChain->Add("/mu2e/data/users/gianipez/hist/treeTrackerDIO.root");
    fNPerPOT    = 7.27e-4;
    fNSimulated = 1.0e7;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "OOT") {
    fChain->Add("/mu2e/data/users/gianipez/hist/treeTrackerOOT.root");
    fNPerPOT    = 3.97e-3;
    fNSimulated = 3.0e7;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "PHOTON") {
    fChain->Add("/mu2e/data/users/gianipez/hist/treeTrackerPHOTON.root");
    fNPerPOT    = 2.28e-3;
    fNSimulated = 1.0e8;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "PROTON") {
    fChain->Add("/mu2e/data/users/gianipez/hist/treeTrackerPROTON.root");
    fNPerPOT    = 5.69e-5;
    fNSimulated = 1.0e8;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "DEUTERON") {
    fChain->Add("/mu2e/data/users/gianipez/hist/treeTrackerDEUTERON.root");
    fNPerPOT    = 2.84e-5;
    fNSimulated = 1.0e8;
    fDensity[0] = 1.7;
  }
  else if (fProcess == "NEUTRON") {
    fChain->Add("/mu2e/data/users/gianipez/hist/treeTrackerNEUTRON.root");
    fNPerPOT    = 1.37e-3;
    fNSimulated = 1.0e8;
    fDensity[0] = 1.7;
  }

  Init(fChain);
//-----------------------------------------------------------------------------
// leave only tracker branches active
//-----------------------------------------------------------------------------
  fChain->SetBranchStatus("*",0);

  fChain->SetBranchStatus("nTtsUpHits",1);
  fChain->SetBranchStatus("ttsUpEdep" ,1);
  fChain->SetBranchStatus("ttsUpX"    ,1);
  fChain->SetBranchStatus("ttsUpY"    ,1);
  fChain->SetBranchStatus("ttsUpZ"    ,1);
  fChain->SetBranchStatus("ttsUpR"    ,1);
  fChain->SetBranchStatus("ttsUpPdgId",1);
  fChain->SetBranchStatus("ttsUpEkin" ,1);
  fChain->SetBranchStatus("ttsUpMass" ,1);
  fChain->SetBranchStatus("ttsUpL"    ,1);
  
  fChain->SetBranchStatus("nTtsDwHits",1);
  fChain->SetBranchStatus("ttsDwEdep" ,1);
  fChain->SetBranchStatus("ttsDwX"    ,1);
  fChain->SetBranchStatus("ttsDwY"    ,1);
  fChain->SetBranchStatus("ttsDwZ"    ,1);
  fChain->SetBranchStatus("ttsDwR"    ,1);
  fChain->SetBranchStatus("ttsDwPdgId",1);
  fChain->SetBranchStatus("ttsDwEkin" ,1);
  fChain->SetBranchStatus("ttsDwMass" ,1);
  fChain->SetBranchStatus("ttsDwL"    ,1);
  
  fChain->SetBranchStatus("vNHits"    ,1);
  fChain->SetBranchStatus("vId"       ,1);
  fChain->SetBranchStatus("vPdgId"    ,1);
  fChain->SetBranchStatus("vP"        ,1);
  fChain->SetBranchStatus("vPx"       ,1);
  fChain->SetBranchStatus("vPy"       ,1);
  fChain->SetBranchStatus("vPz"       ,1);
  fChain->SetBranchStatus("vE"        ,1);
  fChain->SetBranchStatus("vEKin"     ,1);
  fChain->SetBranchStatus("vM"        ,1);
  fChain->SetBranchStatus("vT"        ,1);
  fChain->SetBranchStatus("vX"        ,1);
  fChain->SetBranchStatus("vY"        ,1);
  fChain->SetBranchStatus("vZ"        ,1);
  fChain->SetBranchStatus("vCosth"    ,1);
  fChain->SetBranchStatus("vRadius"   ,1);
  
  return 0;
}

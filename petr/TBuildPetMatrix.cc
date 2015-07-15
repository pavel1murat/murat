///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
#include "petr/TBuildPetMatrix.hh"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"

ClassImp(TBuildPetMatrix)

TBuildPetMatrix::TBuildPetMatrix() {
//-----------------------------------------------------------------------------
// describe imager, crystals approx 5x5
//-----------------------------------------------------------------------------
  fZImager     =  IMAGER_ZMAX;		// in mm
  fRImager     =  IMAGER_R   ;		// in mm

  fCrystalDPhi = 2*TMath::Pi()/IMAGER_PHI_SEG;
  fCrystalDz   = 2*fZImager/IMAGER_Z_SEG;
//-----------------------------------------------------------------------------
// describe phantom, voxel 5x5x5mm, 
//-----------------------------------------------------------------------------
  fPhantomR    =  PHANTOM_R;
  fPhantomMaxZ =  PHANTOM_ZMAX;

  fWriteVoxelData = 0;

  fVoxelBuffer = new TVoxelData();

  fVoxelDx     =  VOXEL_DX;
  fVoxelDy     =  VOXEL_DY;
  fVoxelDz     =  VOXEL_DZ;

  fNVoxelsX    = 2*fPhantomR   /fVoxelDx;
  fNVoxelsY    = 2*fPhantomR   /fVoxelDy;
  fNVoxelsZ    = 2*fPhantomMaxZ/fVoxelDz;

  fRn          = new TRandom3();
//-----------------------------------------------------------------------------
// create voxel map
//-----------------------------------------------------------------------------
  for (int ix=0; ix<fNVoxelsX; ix++) {
    for (int iy=0; iy<fNVoxelsY; iy++) {
      for (int iz=0; iz<fNVoxelsZ; iz++) {
	fVoxel[ix][iy][iz] = new TVoxelData(ix,iy,iz);
      }
    }
  }

  fHist       = new TH2F("h2","h2",IMAGER_PHI_SEG,0,IMAGER_PHI_SEG,IMAGER_Z_SEG,0,IMAGER_Z_SEG);

  BookHist(&fDiagHist[0],0);
  BookHist(&fDiagHist[1],1);

  for (int i=0;i<100; i++) {
    fDebugFlag[i] = 0;
  }
}

//-----------------------------------------------------------------------------
TBuildPetMatrix::~TBuildPetMatrix() {
}



//-----------------------------------------------------------------------------
// generate occupancy matrix corresponding to one voxel
//-----------------------------------------------------------------------------
int TBuildPetMatrix::BookHist(Hist_t* Hist, int ISet) {
  Hist->fX     = new TH1F(Form("x_%i",ISet),Form("x_%i",ISet),200,-100,100);
  Hist->fY     = new TH1F(Form("y_%i",ISet),Form("y_%i",ISet),200,-100,100);
  Hist->fZ     = new TH1F(Form("z_%i",ISet),Form("z_%i",ISet),200,-100,100);
  Hist->fTheta = new TH1F(Form("theta_%i",ISet),Form("theta_%i",ISet),200,-2,2);
  Hist->fPhi1  = new TH1F(Form("phi1_%i" ,ISet),Form("phi1_%i" ,ISet),250, 0,6.25);
  Hist->fIp2VsIp1 = new TH2F(Form("ip2_vs_ip1_%i" ,ISet),Form("ip2_vs_ip1_%i" ,ISet),
			     IMAGER_PHI_SEG,0,IMAGER_PHI_SEG,
			     IMAGER_PHI_SEG,0,IMAGER_PHI_SEG);
  Hist->fIz2VsIz1 = new TH2F(Form("iz2_vs_iz1_%i" ,ISet),Form("iz2_vs_iz1_%i" ,ISet),
			     IMAGER_Z_SEG,0,IMAGER_Z_SEG,
			     IMAGER_Z_SEG,0,IMAGER_Z_SEG);

  return 0;
}

//-----------------------------------------------------------------------------
// reset histogram set
//-----------------------------------------------------------------------------
int TBuildPetMatrix::ResetDiagHist(Hist_t* Hist) {
  Hist->fX->Reset();
  Hist->fY->Reset();
  Hist->fZ->Reset();
  Hist->fTheta->Reset();
  Hist->fPhi1->Reset();

  return 0;
}


//-----------------------------------------------------------------------------
// generate occupancy matrix corresponding to one voxel
//-----------------------------------------------------------------------------
int TBuildPetMatrix::GenerateVoxelData(int Ix, int Iy, int Iz, int NEvents) {

  int                              rc (0), iz1, iz2, iphi1, iphi2, nc, found ;
  double                           x, y, z,  r2, nx, ny, nz, nxy, xnxy;
  double                           theta, phi1, phi2, t1, t2, z1, z2;
  TVoxelData::Occ_t                data;
  std::vector<TVoxelData::Occ_t>*  v;

  if (fVoxel[Ix][Iy][Iz] == 0) fVoxel[Ix][Iy][Iz] = new TVoxelData(Ix,Iy,Iz);
  else                            fVoxel[Ix][Iy][Iz]->Clear();

  TVoxelData* vd = fVoxel[Ix][Iy][Iz];

  vd->fNEvents   = NEvents;
  vd->fNin       = 0;
  vd->fNCellsTot = 0;

  r2 = fRImager*fRImager;

  ResetDiagHist(&fDiagHist[0]);
  ResetDiagHist(&fDiagHist[1]);

  for (int ievent=0; ievent<NEvents; ievent++) {

				        // generate point within the voxel

    x = (Ix+fRn->Rndm(ievent))*fVoxelDx-fPhantomR;
    y = (Iy+fRn->Rndm(ievent))*fVoxelDy-fPhantomR;

    // this check is already done at a level of voxels
    // perhaps the edge effects need to be handled more accurately and for  voxels 
    // on the edge the check needs to be done for each point
    // NOT NOW! 

    //    if (x*x+y*y > r2)                                       goto NEXT_EVENT;

    z = (Iz+fRn->Rndm(ievent))*fVoxelDz-fPhantomMaxZ;
//-----------------------------------------------------------------------------
// generate direction - polar angle, then - azimuthal angle
// 0 < phi < 2pi
//-----------------------------------------------------------------------------
    theta = TMath::Pi()/2*(2*fRn->Rndm(ievent)-1);
    phi1  = 2*TMath::Pi()*fRn->Rndm(ievent);

    fDiagHist[0].fX->Fill(x);
    fDiagHist[0].fY->Fill(y);
    fDiagHist[0].fZ->Fill(z);
    fDiagHist[0].fTheta->Fill(theta);
    fDiagHist[0].fPhi1->Fill(phi1);

    nx    = TMath::Cos(theta)*TMath::Cos(phi1);
    ny    = TMath::Cos(theta)*TMath::Sin(phi1);
    nz    = TMath::Sin(theta);

    nxy   = sqrt(nx*nx+ny*ny);

    xnxy  = x*nx+y*ny;

    t1    = (-xnxy+sqrt(xnxy*xnxy+r2-x*x-y*y))/(nxy+1.e-16);
    t2    = (-xnxy-sqrt(xnxy*xnxy+r2-x*x-y*y))/(nxy+1.e-16);

    z1    = z+t1*nz;
    z2    = z+t2*nz;

    if (fabs(z1) > fZImager)                                goto NEXT_EVENT;
//-----------------------------------------------------------------------------
// one more first photon within the acceptance for this voxel
// look at the second one
//-----------------------------------------------------------------------------
    iz1                      = (z1+fZImager)/fCrystalDz;
    iphi1                    = phi1/fCrystalDPhi;
    vd->fMatrix[iphi1][iz1] += 1;

    if (fabs(z2) > fZImager)                                   goto NEXT_EVENT;
//-----------------------------------------------------------------------------
// the second photon is also within the acceptance
//-----------------------------------------------------------------------------
    vd->fNin  += 1;
//-----------------------------------------------------------------------------
// calculate "hit crystal" numbers, fill diagnostics histograms
//-----------------------------------------------------------------------------
    iz2   = (z2+fZImager)/fCrystalDz;

    if   (phi1 < TMath::Pi()) phi2 = phi1+TMath::Pi();
    else                      phi2 = phi1-TMath::Pi();
	
    iphi2 = phi2/fCrystalDPhi;
    v     = vd->fCellData[iphi1][iz1];
    nc    = v->size();
//-----------------------------------------------------------------------------
// the vector exists, first check 
//-----------------------------------------------------------------------------
    found = 0;
    for (int i=0; i<nc; i++) {
      TVoxelData::Occ_t* o = &v->at(i);
      if ((o->fIPhi == iphi2) && (o->fIz == iz2)) {
	found = 1;
	o->fNEntries++;
	break;
      }
    }
    if (found == 0) {
      data.fIPhi     = iphi2;
      data.fIz       = iz2;
      data.fNEntries = 1;
      v->push_back(data);		// new cell

      vd->fNCells[iphi1][iz1] += 1;
      vd->fNCellsTot          += 1;
    }
//-----------------------------------------------------------------------------
// fill the second set of diagnostics histograms
//-----------------------------------------------------------------------------
    fDiagHist[1].fX->Fill(x);
    fDiagHist[1].fY->Fill(y);
    fDiagHist[1].fZ->Fill(z);
    fDiagHist[1].fTheta->Fill(theta);
    fDiagHist[1].fPhi1->Fill(phi1);
    fDiagHist[1].fIp2VsIp1->Fill(iphi1,iphi2);
    fDiagHist[1].fIz2VsIz1->Fill(iz1,iz2);
//-----------------------------------------------------------------------------
// next event
//-----------------------------------------------------------------------------
  NEXT_EVENT:;
  }

  if (fDebugFlag[1] > 0) {
    //    printf(">>>  TBuildPetMatrix::GenerateVoxelData print fVoxelBuffer \n");
    vd->Print();
  }

  return rc;
}

//-----------------------------------------------------------------------------
int TBuildPetMatrix::Generate(int NEventsPerVoxel) {

  int rc = 0;

  //  double x0, y0; // , z0; // , r2;

  if (fWriteVoxelData) {
    fFile      = new TFile("build_pet_matrix.root","recreate");
    fTree      = new TTree("VoxelDataTree","VoxelDataTree");

    fVoxelDataBranch = fTree->Branch("TVoxelData","TVoxelData",&fVoxelBuffer,64000,0);
    fVoxelDataBranch->SetCompressionLevel(1);
  }

  fNEventsPerVoxel = NEventsPerVoxel;

  //  r2 = fPhantomR*fPhantomR;

  for (int iz=0; iz<fNVoxelsZ; iz++) {
    //    z0 = -fPhantomMaxZ+(iz+0.5)*fVoxelDz;

    for (int ix=0; ix<fNVoxelsX; ix++) {
      //      x0 = -fPhantomR+(ix+0.5)*fVoxelDx;

      for (int iy=0; iy<fNVoxelsY; iy++) {
	//	y0 = -fPhantomR+(iy+0.5)*fVoxelDy;

	//	double rho2 = x0*x0+y0*y0;
	//	if (rho2 < r2) {

	GenerateVoxelData(ix,iy,iz,NEventsPerVoxel);
//-----------------------------------------------------------------------------
// and write the event out
// start from printing the vector (or plotting it...)
// voxel inside the phantom; generate events
// generation also writes the output ntuple
//-----------------------------------------------------------------------------
	WriteVoxelData(fVoxel[ix][iy][iz]);
	//	}
      }
    }
  }

  if (fWriteVoxelData) {
    fTree->Write();
    fFile->Close();
  }
  return rc;
}


//-----------------------------------------------------------------------------
int TBuildPetMatrix::WriteVoxelData(TVoxelData* Voxel) {
  int rc (0);
  if (fWriteVoxelData) {
					// fill voxel data
    Voxel->Copy(fVoxelBuffer);

    if (fDebugFlag[0] == 1) {
      //      printf(">>>  TBuildPetMatrix::WriteVoxelData print fVoxelBuffer \n");
      fVoxelBuffer->Print();
    }
    else if (fDebugFlag[0] == 2) {
      fVoxelBuffer->Print("full");
    }
					// and fill the tree
    fTree->Fill();
  }
  return rc;
}

//-----------------------------------------------------------------------------
// results are in fVector
//-----------------------------------------------------------------------------
int TBuildPetMatrix::ReadVoxelMap(const char* Filename) {
  int ix, iy, iz;
  TFile* f = TFile::Open(Filename);

  TTree* t = (TTree*) f->Get("VoxelDataTree");
  
  TBranch* b = t->GetBranch("TVoxelData");
  b->SetAddress(&fVoxelBuffer);
  b->SetAutoDelete(0);

  int nv = t->GetEntries();

  for (int iv=0; iv<nv; iv++) {
    t->GetEntry(iv);

    ix = fVoxelBuffer->fIx;
    iy = fVoxelBuffer->fIy;
    iz = fVoxelBuffer->fIz;

    if (fVoxel[ix][iy][iz] == 0) fVoxel[ix][iy][iz] = new TVoxelData();

    if (fDebugFlag[1] == 1) {
      //      printf(">>>  TBuildPetMatrix::ReadVoxelMap print fVoxelBuffer \n");
      // fVoxelBuffer->Print();
    }

    fVoxelBuffer->Copy(fVoxel[ix][iy][iz]);

    if (fDebugFlag[1] == 1) {
      //      printf(">>>  TBuildPetMatrix::ReadVoxelMap print fVoxel[ix][iy][iz] \n");
      fVoxel[ix][iy][iz]->Print();
    }
  }
 
  f->Close();
  
  return 0;
}

//-----------------------------------------------------------------------------
// results are in fVector
//-----------------------------------------------------------------------------
int TBuildPetMatrix::PlotResults(int VoxelIX, int VoxelIY, int VoxelIZ, int IPhi, int Iz) {
  double x0, y0, z0;

  fHist->Reset();

  TVoxelData* vd = fVoxel[VoxelIX][VoxelIY][VoxelIZ]; 

  x0 = -fPhantomR   +(VoxelIX+0.5)*fVoxelDx;
  y0 = -fPhantomR   +(VoxelIY+0.5)*fVoxelDy;
  z0 = -fPhantomMaxZ+(VoxelIZ+0.5)*fVoxelDz;

  printf(" VoxelIX, VoxelIY, VoxelIZ, x0, y0, z0 = %4i %4i %4i %10.3f %10.3f %10.3f\n",
	 VoxelIX,VoxelIY,VoxelIZ,x0,y0,z0);

  if (vd == 0) {
    printf(">>> ERROR in TBuildPetMatrix::PlotResults: voxel Ix=%4i, Iy=%4i, Iz=%4i not defined\n",
	   VoxelIX, VoxelIY, VoxelIZ);
    return -1;
  }

  int n = vd->fMatrix[IPhi][Iz];

  printf("n = %5i\n",n);

  std::vector<TVoxelData::Occ_t>* v = vd->fCellData[IPhi][Iz];

  if (v == 0) {
    printf(">>> ERROR in TBuildPetMatrix::PlotResults: vector not defined for IPhi=%3i Iz=%3i\n",IPhi,Iz);
    return -1;
  }

  int ncells = v->size();

  printf("n = %5i  nin  = %5i  ncells = %5i\n",n, vd->fNin, ncells);

  for (int ic=0;ic<ncells; ic++) {
    TVoxelData::Occ_t* cell = &v->at(ic);
    
    for (int i=0; i<cell->fNEntries; i++) {
      fHist->Fill(cell->fIPhi,cell->fIz);
    } 
    printf("[cell data] IPhi = %3i Iz = %3i NEntries = %6i\n",cell->fIPhi,cell->fIz,cell->fNEntries);
  }
  
  fHist->Draw("box");

  return 0;
}


//-----------------------------------------------------------------------------
int TBuildPetMatrix::PlotIPhiIZ(int VoxelIX, int VoxelIY, int VoxelIZ, int IPhi, int IZ) {

  int                             nc, iphi1, iphi2, iz1, iz2;
  std::vector<TVoxelData::Occ_t>* cd;
  TVoxelData::Occ_t*              cell;
  
  TH2F   *h_ip2_vs_ip1, *h_iz2_vs_iz1;
  
  h_ip2_vs_ip1 = new TH2F("h_ip2_vs_ip1","h_ip2_vs_ip1",
			  IMAGER_PHI_SEG,0,IMAGER_PHI_SEG,
			  IMAGER_PHI_SEG,0,IMAGER_PHI_SEG);
  
  h_iz2_vs_iz1 = new TH2F("iz2_vs_iz1_%i","iz2_vs_iz1",
			  IMAGER_Z_SEG,0,IMAGER_Z_SEG,
			  IMAGER_Z_SEG,0,IMAGER_Z_SEG);

  TVoxelData* vd = GetVoxel(VoxelIX,VoxelIY,VoxelIZ);

  if (IPhi >= 0) { iphi1 = IPhi; iphi2 = iphi1+1;       }
  else           { iphi1 = 0   ; iphi2 = IMAGER_PHI_SEG;}

  if (IZ >= 0) { iz1 = IZ; iz2 = iz1+1;     }
  else         { iz1 = 0 ; iz2 = IMAGER_Z_SEG; }

  for (int iphi=iphi1; iphi<iphi2; iphi++) {
    for (int iz=iz1; iz<iz2; iz++) {
      nc = vd->GetNCells(iphi,iz);

      if (nc > 0) {
	cd = vd->GetCellData(iphi,iz);
	for (int ic=0; ic<nc; ic++) {
	  cell = &cd->at(ic);
	  for (int i=0; i<cell->fNEntries; i++) {
	    h_ip2_vs_ip1->Fill(iphi,cell->fIPhi);
	    h_iz2_vs_iz1->Fill(iz  ,cell->fIz);
	  }
	}
      }

    }
  }
 
  TCanvas* c = new TCanvas("c_PlotIPhiIZ","c_PlotIPhiIZ",1200,600);
  c->Divide(2,1);

  c->cd(1);
  h_ip2_vs_ip1->Draw("box");

  c->cd(2);
  h_iz2_vs_iz1->Draw("box");
  
  return 0;
}


//-----------------------------------------------------------------------------
void TBuildPetMatrix::PrintVoxel(int IX, int IY, int IZ) {
  GetVoxel(IX,IY,IZ)->Print("full");
}

//-----------------------------------------------------------------------------
// results are in fVector
//-----------------------------------------------------------------------------
void TBuildPetMatrix::Validate(int Ix, int Iy, int Iz, int CPhi, int Cz) {

  TVoxelData* vd;

  int   nc1, nc2;

  int ix1, ix2, iy1, iy2, iz1, iz2, cphi1, cphi2, cz1, cz2;

  if (Ix >= 0) { ix1 = Ix; ix2 = ix1+1;     }
  else         { ix1 = 0 ; ix2 = fNVoxelsX; }

  if (Iy >= 0) { iy1 = Iy; iy2 = iy1+1;   }
  else         { iy1 = 0 ; iy2 = fNVoxelsY; }

  if (Iz >= 0) { iz1 = Iz; iz2 = iz1+1;     }
  else         { iz1 = 0 ; iz2 = fNVoxelsZ; }

  if (CPhi >= 0) { cphi1 = CPhi; cphi2 = cphi1+1;       }
  else           { cphi1 = 0   ; cphi2 = fNCrystalsPhi; }

  if (Cz >= 0) { cz1 = Cz; cz2 = cz1+1;     }
  else         { cz1 = 0 ; cz2 = fNVoxelsZ; }

  for (int ix=ix1; ix<ix2; ix++) {
    for (int iy=iy1; iy<iy2; iy++) {
      for (int iz=iz1; iz<iz2; iz++) {

	vd = fVoxel[ix][iy][iz]; 	

	for (int cphi=cphi1; cphi<cphi2; cphi++) {
	  for (int cz=cz1; cz<cz2; cz++) {

	    nc1 = vd->fNCells  [cphi][cz];
	    nc2 = vd->fCellData[cphi][cz]->size();

	    if (nc1 != nc2) {
	      printf(">>> ERROR in TBuildPetMatrix::Validate: (ix,iy,iz):(%3i,%3i,%3i)",
		     ix,iy,iz);
	      printf(" [cphi][cz]:[%4i][%4i]", cphi,cz);
	      printf(" nc1 = %i nc2 = %10ii\n",nc1,nc2);
	    }
	  }
	}
      }
    }
  }
 
  return;
}


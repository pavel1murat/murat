///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "murat/petr/TZubalPhantom.hh"

#include "TCanvas.h"
#include "TH2.h"

ClassImp(TZubalPhantom)

using std::abs;

//-----------------------------------------------------------------------------
TZubalPhantom::TZubalPhantom(const char* Filename) {

  fIXMin = 0;
  fIXMax = 256;
  fIYMin = 0;
  fIYMax = 256;
  fIZMin = 0;
  fIZMax = 128;

  FILE* f = fopen(Filename,"r");

  int nb = fread(fData,sizeof(char),128*256*256,f);

  printf("nb = %i\n",nb);

  for (int iz=0; iz<128; iz++) {
    for (int ix=0; ix<256; ix++) {
      for (int iy=0; iy<256; iy++) {
	VolID(ix,iy,iz) = -1;
      }
    }
  }

  fNVolumes      = 0;
  fListOfVolumes = new TObjArray(100);
}


//-----------------------------------------------------------------------------
TZubalPhantom::~TZubalPhantom() {
}


//-----------------------------------------------------------------------------
// Mode = 0: check neigbours
// Mode = 1: add neighbours
//-----------------------------------------------------------------------------
int TZubalPhantom::CheckSurroundingVoxels(int Ix, int Iy, int Iz, int Color, TZubalVolume* Vol) {

  // so far - within one slice

  int dist, col, id;

  TZubalVolume* vol;

  int  rc = -1;

  for (int iz=Iz-1; iz<Iz+2; iz++) {
    if ((iz < fIZMin) || (iz >= fIZMax)) goto NEXT_IZ;
    for (int ix=Ix-1; ix<Ix+2; ix++) {
      if ((ix < fIXMin) || (ix >= fIXMax)) goto NEXT_IX;
      for (int iy=Iy-1; iy<Iy+2; iy++) {
	if ((iy < fIYMin) || (iy >= fIYMax))  goto NEXT_IY;
//-----------------------------------------------------------------------------
// the data go by Z slice
// skip central voxel itself
//-----------------------------------------------------------------------------
	dist = abs(ix-Ix)+abs(iy-Iy)+abs(iz-Iz);

	if (dist == 1) {
	  col = Data (ix,iy,iz);
	  id  = VolID(ix,iy,iz);
	
	  if (col == Color) {
	    if (Vol == NULL) {
//-----------------------------------------------------------------------------
// check if this voxel need to be added
//-----------------------------------------------------------------------------
	      if (id >= 0) {
//-----------------------------------------------------------------------------
// neigbor is a part of a volume, add this one to the same volume
//-----------------------------------------------------------------------------
		vol = Volume(id);
		rc = vol->AddVoxel(Ix,Iy,Iz);
		if (rc == 0) {
		  VolID(Ix,Iy,Iz) = id;
		  rc = 0;
	                                                    goto RETURN;
		}
		else {
		  printf("ERROR : 001 voxel re-added \n");
		}
	      }
	    }
	    else {
//-----------------------------------------------------------------------------
// see if neighbour needs to be added to this volume
//-----------------------------------------------------------------------------
	      if (id < 0) {
					// add
		rc = Vol->AddVoxel(ix,iy,iz);
		if (rc == 0) {
		  VolID(ix,iy,iz) = Vol->ID();
		}
		else {
		  printf("ERROR : 002 voxel re-added \n");
		}
	      }
	    }
					//  neighbor is not associated with any volume yet, continue
	  }
					// neighbor has a different color, continue
	}
      NEXT_IY:;
      }
    NEXT_IX:;
    }
  NEXT_IZ:;
  }

 RETURN:
  return rc;
}

//-----------------------------------------------------------------------------
int TZubalPhantom::ExpandVolume(TZubalVolume* Vol) {

  int   ix, iy, iz, n, loc;
					// to be deleted in the end
  short* list;

  Vol->ResetListOfAddedVoxels(n,list);

  for (int i=0; i<n; i++) {
    loc = 3*i;
    ix = list[loc  ];
    iy = list[loc+1];
    iz = list[loc+2];
    CheckSurroundingVoxels(ix,iy,iz,Vol->fColor,Vol);
  }
  
  delete list;

  return Vol->NAddedVoxels();
}

//-----------------------------------------------------------------------------
int TZubalPhantom::Parse3DMap(int IZMin, int IZMax) {

  int color = 0;
  int nv;

  // fIXMin=120;
  // fIXMax=122;
  // fIYMin=120;
  // fIYMax=122;
  fIZMin=IZMin;
  fIZMax=IZMax;

  for (int iz=IZMin; iz<IZMax; iz++) {
    //    for (int ix=0; ix<256; ix++) {
    for (int ix=fIXMin; ix<fIXMax; ix++) {
      //      for (int iy=0; iy<256; iy++) {
      for (int iy=fIYMin; iy<fIYMax; iy++) {
	color = Data(ix,iy,iz);

	if ((color != 0) && (fVolID[iz][ix][iy] == -1)) {
//-----------------------------------------------------------------------------
// check surrounding voxels - start from the slice
// if rc == 0, the voxel has been added to existing volume
//-----------------------------------------------------------------------------
	  int rc = CheckSurroundingVoxels(ix,iy,iz,color,0);
	  
	  if (rc < 0) {
//-----------------------------------------------------------------------------
// start new volume and continue checking surrounding voxels
//-----------------------------------------------------------------------------
	    printf(" >>> New Volume: ix,iy,iz,color = %3i %3i %3i %3i\n",
		   ix,iy,iz,color);
	    fVol = new TZubalVolume(fNVolumes,ix,iy,iz,color);

	    VolID(ix,iy,iz) = fNVolumes;
	    fListOfVolumes->Add(fVol);
	    fNVolumes += 1;

	    while (fVol->fNAddedVoxels > 0) {
	      nv = ExpandVolume(fVol);
	      if (nv < 0) printf("ERROR: nv = %i",nv);
	    }
	  }
	}
      }
    }
  }
  return 0;
}


//-----------------------------------------------------------------------------
// tests
//-----------------------------------------------------------------------------
void TZubalPhantom::Test01(const char* Filename, int IZ) {
  
  TZubalPhantom* z = new TZubalPhantom(Filename);

  TH2F* h2 = new TH2F(Form("h%03i",IZ),Form("h2 %03i",IZ),256,0,256,256,0,256);

  for (int i=0; i<256; i++) {
    for (int j=0; j<256; j++) {
      h2->SetBinContent(i,j,z->fData[IZ][i][j]);
    }
  }

  h2->Draw("col");

}


//-----------------------------------------------------------------------------
void TZubalPhantom::Test02(int IZMin, int IZMax, int IVol) {

  TZubalPhantom* z = new TZubalPhantom("/home/murat/Downloads/det_head_u2.med");

  z->Parse3DMap(IZMin,IZMax);

  // at this point parsed, print results

  printf("NVolumes = %i\n",z->fNVolumes);
//-----------------------------------------------------------------------------
// count N(voxels) in volumes
//-----------------------------------------------------------------------------
  int nvol = z->NVolumes();
  int nvv = 0;
  for (int i=0; i<nvol; i++) {
    TZubalVolume* v = z->Volume(i);
    printf("volID = %5i color = %3i NVoxels = %8i\n",i,v->Color(),v->NVoxels());
    nvv += v->NVoxels();
  }
//-----------------------------------------------------------------------------
// count total number N(voxels) in the data
//-----------------------------------------------------------------------------
  int nvd = 0;
  for (int iz=IZMin; iz<IZMax; iz++) {
    for (int ix=z->fIXMin; ix<z->fIXMax; ix++) {
      for (int iy=z->fIYMin; iy<z->fIYMax; iy++) {
	if (z->Data(ix,iy,iz) > 0) nvd++;
      }
    }
  }

  printf("nv(data) = %6i nv(volumes) = %6i\n",nvd,nvv);


  return;

  int has_voxel;
  TH2C* h3;

  TZubalVolume* v;
  TCanvas*      c;

  for (int iv=0; iv<nvol; iv++){
    v = z->Volume(iv);

    if (iv == IVol) {
      for (int iz=IZMin; iz<IZMax; iz++) {
	c = new TCanvas(Form("c_%03i_%03i",iv,iz),Form("c_%03i_%03i",iv,iz),800,800);
	c->cd();
	h3 = new TH2C(Form("h3_%03i_%03i",iv,iz),Form("h3 %03i %03i",iv,iz),256,0,256,256,0,256);

	for (int ix=0; ix<256; ix++) {
	  for (int iy=0; iy<256; iy++) {
	    has_voxel = v->HasVoxel(ix,iy,iz);
	    
	    if (has_voxel != 0) {
		h3->SetBinContent(ix,iy,iz,v->Color());
	    }
	  }
	}

	h3->Draw("col");
      }
    }

    // for (int ix=0; ix<256; ix++) {
    //   for (int iy=0; iy<256; iy++) {
    // 	color = z->fVolID[IZMin][ix][iy];

    // 	if (color >= 0) {
    // 	  h2->SetBinContent(ix,iy,color+1);
    // 	}
    //   }
    // }
  }



}

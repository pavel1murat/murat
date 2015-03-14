///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
#include "murat/petr/TVoxelData.hh"
#include "TMath.h"

ClassImp(TVoxelData)

//-----------------------------------------------------------------------------
TVoxelData::TVoxelData() {

  fIx = -1;
  fIy = -1;
  fIz = -1;

  for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
    for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
      fMatrix[iphi][iz] = 0;
      fNCells[iphi][iz] = 0;
      fCellData[iphi][iz] = new std::vector<Occ_t> ;
    }
  }

  fNEvents   = 0;
  fNin       = 0;
  fNCellsTot = 0;
}


//-----------------------------------------------------------------------------
TVoxelData::TVoxelData(int Ix, int Iy, int Iz) {
  fIx = Ix;
  fIy = Iy;
  fIz = Iz;

  for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
    for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
      fMatrix[iphi][iz] = 0;
      fNCells[iphi][iz] = 0;
      fCellData[iphi][iz] = new std::vector<Occ_t> ;
    }
  }

  fNEvents   = 0;
  fNin       = 0;
  fNCellsTot = 0;

  
}

//-----------------------------------------------------------------------------
TVoxelData::~TVoxelData() {

   for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
    for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
      delete fCellData[iphi][iz];
    }
  }
 
}


//______________________________________________________________________________
void TVoxelData::Streamer(TBuffer &R__b) {
  // Stream an object of class TVoxelData.

  int   nc;
  int   cd[3];
  Occ_t cell;


  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
//-----------------------------------------------------------------------------
// read in
//-----------------------------------------------------------------------------
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
    if (R__v) { }
    TObject::Streamer(R__b);

    R__b >> fIx;
    R__b >> fIy;
    R__b >> fIz;

    R__b >> fNEvents;
    R__b >> fNin;
    R__b >> fNCellsTot;

    R__b.ReadFastArray((int*) fMatrix, IMAGER_PHI_SEG*IMAGER_Z_SEG);
    R__b.ReadFastArray((int*) fNCells, IMAGER_PHI_SEG*IMAGER_Z_SEG);

    for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
      for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
	nc = fNCells[iphi][iz];
	for (int i=0; i<nc; i++) {
	  R__b.ReadFastArray(cd,3);
	  cell.fIPhi     = cd[0];
	  cell.fIz       = cd[1];
	  cell.fNEntries = cd[2];
	  fCellData[iphi][iz]->push_back(cell);
	}
      }
    }

    R__b.CheckByteCount(R__s, R__c, TVoxelData::IsA());
  } 
  else {
//-----------------------------------------------------------------------------
// write out
//-----------------------------------------------------------------------------
    R__c = R__b.WriteVersion(TVoxelData::IsA(), kTRUE);
    TObject::Streamer(R__b);

    R__b << fIx;
    R__b << fIy;
    R__b << fIz;

    R__b << fNEvents;
    R__b << fNin;
    R__b << fNCellsTot;

    R__b.WriteFastArray((int*) fMatrix,IMAGER_PHI_SEG*IMAGER_Z_SEG);
    R__b.WriteFastArray((int*) fNCells,IMAGER_PHI_SEG*IMAGER_Z_SEG);
//-----------------------------------------------------------------------------
// now write out cells
//-----------------------------------------------------------------------------
    for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
      for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
	nc = fNCells[iphi][iz];
	for (int i=0;  i<nc; i++) { 
	  Occ_t* cell = &fCellData[iphi][iz]->at(i);
	  cd[0] = cell->fIPhi;
	  cd[1] = cell->fIz;
	  cd[2] = cell->fNEntries;
	  R__b.WriteFastArray(cd,3);
	}
      }
    }
    
    R__b.SetByteCount(R__c, kTRUE);
  }
}

//-----------------------------------------------------------------------------
void TVoxelData::Copy(TVoxelData* Voxel) const {
  
  int nc; 

  Voxel->fIx        = fIx;
  Voxel->fIy        = fIy;
  Voxel->fIz        = fIz;

  Voxel->fNEvents   = fNEvents;
  Voxel->fNin       = fNin;
  Voxel->fNCellsTot = fNCellsTot;

  memcpy(Voxel->fMatrix,fMatrix,IMAGER_PHI_SEG*IMAGER_Z_SEG*sizeof(int));
  memcpy(Voxel->fNCells,fNCells,IMAGER_PHI_SEG*IMAGER_Z_SEG*sizeof(int));

  for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
    for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
      if (Voxel->fCellData[iphi][iz]->size() > 0) Voxel->fCellData[iphi][iz]->clear();
      nc = fNCells[iphi][iz];
      if (nc > 0) {
	for (int ic=0; ic<nc; ic++) {
	  const Occ_t* occ = &fCellData[iphi][iz]->at(ic);
	  Voxel->fCellData[iphi][iz]->push_back(*occ); 
	}
      }
    }
  }
}

//-----------------------------------------------------------------------------
void TVoxelData::Clear(Option_t* Opt) {
  
  fNEvents = 0;
  fNin     = 0;

  memset((int*) fMatrix,0,IMAGER_PHI_SEG*IMAGER_Z_SEG*sizeof(int));
  memset((int*) fNCells,0,IMAGER_PHI_SEG*IMAGER_Z_SEG*sizeof(int));

  for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
    for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
      fCellData[iphi][iz]->clear();
    }
  }
}

//-----------------------------------------------------------------------------
// print only non-empty cells 
void TVoxelData::Print(Option_t* Opt) const {
  
  int  nc1, nc2, nc, nct(0), nent(0);

  for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
    for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
      nc = fCellData[iphi][iz]->size();
      nct += nc;
      for (int i=0; i<nc; i++) {
	Occ_t* cell = &fCellData[iphi][iz]->at(i);
	nent += cell->fNEntries;
      }
    }
  }
  
  printf (">>> IX,IY,IZ = %5i %5i %5i", fIx, fIy, fIz);
  printf (" NEvents: %8i Nin: %8i NCellsTot: %5i nct: %5i nent: %5i\n",
	  fNEvents, fNin, fNCellsTot, nct, nent);

  if (strcmp(Opt,"full") == 0) {

    for (int iphi=0; iphi<IMAGER_PHI_SEG; iphi++) {
      for (int iz=0; iz<IMAGER_Z_SEG; iz++) {
	if (fMatrix[iphi][iz] != 0) {
	  printf (">>> IX,IY,IZ = %5i %5i %5i", fIx, fIy, fIz);
	  printf("    ... iphi1 = %3i iz1 = %3i fMatrix[iphi][iz] = %10i fNCells[iphi][iz] = %3i\n",
		 iphi,iz,fMatrix[iphi][iz],fNCells[iphi][iz]);
	  
	  nc1 = fNCells[iphi][iz];
	  nc2 = fCellData[iphi][iz]->size();
	  if (nc1 != nc2) {
	    printf(" TVoxelData::Print ERROR: iphi=%3i iz=%3i nc1=%5i nc2=%5i\n",
		   iphi,iz,nc1,nc2);
	    return;
	  }
	  for (int ic=0; ic<nc1; ic++) {
	    const Occ_t* cell = & fCellData[iphi][iz]->at(ic);
	    printf (">>> IX,IY,IZ = %5i %5i %5i", fIx, fIy, fIz);
	    printf("        iphi2 = %3i iz2 = %3i  N(entries)=%10i\n",
		   cell->fIPhi,cell->fIz,cell->fNEntries);
	  }
	}
      }
    }
  }
}

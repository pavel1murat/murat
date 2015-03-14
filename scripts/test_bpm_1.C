/*
.L lib/libmu2e_Mu2eG4.so
.L lib/libmurat_pet.so 
.L murat/scripts/test_bpm_1.C

*/

TBuildPetMatrix* bpm(0);

//-----------------------------------------------------------------------------
// test 1: build and visualize the PET reco matrix for one voxel
//-----------------------------------------------------------------------------
int test1(int ix=10, int iy=10, int iz=10, int N = 10000) {

  bpm = new TBuildPetMatrix();

  bpm->SetDebugFlag(0,1);

  bpm->GenerateVoxelData(ix,iy,iz,N);

  //  bpm->PlotResults(0,1);
}

//-----------------------------------------------------------------------------
// test 2: build PET reco matrix for all voxels and write it out
//-----------------------------------------------------------------------------
int test2(int N = 10000, int Debug = 0) {

  bpm = new TBuildPetMatrix();

  bpm->SetDebugFlag(0,Debug);
  bpm->SetWriteVoxelData(1);
  bpm->Generate(N);

  //  bpm->PlotResults(0,1);
}

//-----------------------------------------------------------------------------
// test 3: read BPM matrix back and plot data for one of the voxels
//-----------------------------------------------------------------------------
int test3(const char* Filename, int Debug = 0) {

  bpm = new TBuildPetMatrix();

  bpm->SetDebugFlag(1,Debug);
  bpm->ReadVoxelMap(Filename);

  printf(">>> test3: validating\n");
  bpm->Validate();
}
//-----------------------------------------------------------------------------
// test 3: read BPM matrix back and plot data for one of the voxels
//-----------------------------------------------------------------------------
int test_print_all() {

  printf(">>> test3: validating\n");
  bpm->Validate();

  for (int ix=0; ix<bpm->fNVoxelsX; ix++) {
    for (int iy=0; iy<bpm->fNVoxelsY; iy++) {
      for (int iz=0; iz<bpm->fNVoxelsZ; iz++) {
	bpm->fVoxelMap[ix][iy][iz]->Print();
      }
    }
  }

}


//-----------------------------------------------------------------------------
// 2013-01-24 P.M: make an ntuple starting from the stopped muon text file 
// in a format used by the Mu2e offline :
// .../mu2e/DataFiles/ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt
//-----------------------------------------------------------------------------
void make_ntuple(const char* FnData, const char* FnNtuple) {

  TFile* f = new TFile(FnNtuple,"recreate");
  TNtuple *ntuple = new TNtuple("stops","Stopped muons","x:y:z:time");

  FILE* fdata = fopen(FnData,"r");

  int    done = 0;
  char   c[1000];

  // read the very first line

  fgets (c, 1000, fdata);

  float  x,y,z,time;
  
  while ( ((c[0]=getc(fdata)) != EOF) && !done) {

					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],fdata);
      // read channel number
      fscanf(fdata,"%f" ,&x   );
      fscanf(fdata,"%f" ,&y    );
      fscanf(fdata,"%f" ,&z );
      fscanf(fdata,"%f" ,&time );

      //      printf("x,y,z,time: %15.8f %15.8f %15.8f %15.8f\n",x,y,z,time);
      //      break;
      
      ntuple->Fill(x,y,z,time);
    }

  }

  fclose(fdata);
  
  
  ntuple->Write();
  
  f->Write();
  
  delete f;

}

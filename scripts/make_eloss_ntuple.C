//
// FnData: name of the text file with the data
//
void make_eloss_ntuple(const char* FnData, const char* FnNtuple) {

  TFile* f = new TFile(FnNtuple,"recreate");
  TNtuple *ntuple = new TNtuple("ntuple","Demo ntuple","event:pcode:z:step:eloss");


  FILE* fdata = fopen(FnData,"r");

  int    done = 0;
  char   c[1000];

  int    event, pcode;
  float  z, eloss,step;
  
  while ( ((c[0]=getc(fdata)) != EOF) && !done) {

					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],fdata);
      // read channel number
      fscanf(fdata,"%i" ,&event   );
      fscanf(fdata,"%i" ,&pcode    );
      fscanf(fdata,"%f" ,&z );
      fscanf(fdata,"%f" ,&step );
      fscanf(fdata,"%f",&eloss    );
      
      ntuple->Fill(event,pcode,z,step,eloss);
    }

  }

  fclose(fdata);
  
  
  ntuple->Write();
  
  f->Write();
  
  delete f;
}

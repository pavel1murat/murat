//


//-----------------------------------------------------------------------------
void make_my_example() {

  struct data_t {
    float x1;
    float x2;
    float w;
  } data;

  TFile* f = new TFile("my_tmva_training_example.root","RECREATE");

  TTree* ts = new TTree("signal"    ,"S");
  TTree* tb = new TTree("background","B");

  TBranch* b_sig_x1 = ts->Branch("x1",&data.x1,"F");
  TBranch* b_sig_x2 = ts->Branch("x2",&data.x2,"F");

  TBranch* b_bgr_x1 = tb->Branch("x1",&data.x1,"F");
  TBranch* b_bgr_x2 = tb->Branch("x2",&data.x2,"F");
  TBranch* b_bgr_wt = tb->Branch("weight",&data.w,"F");

  // TBranch* b_sig_x3 = ts->Branch("x3",&data.x3,"F");
  // TBranch* b_sig_x4 = ts->Branch("x4",&data.x4,"F");

  ts->Print();

  TRandom3 rn3;

  int ntest (10000);

  for (int i=0; i<ntest; i++) {
    double phi = TMath::Pi()*2*rn3.Rndm(i);
    data.x1    = TMath::Sin(phi);
    data.x2    = TMath::Cos(phi);

    // data.x3    = 10*TMath::Sin(phi);
    // data.x4    = 10*TMath::Cos(phi);

    ts->Fill();
  }

  for (int i=0; i<ntest; i++) {
    double phi = TMath::Pi()*2*rn3.Rndm(i);
    data.x1    = 0.5+TMath::Sin(phi);
    data.x2    = 0.5+TMath::Cos(phi);
    data.w     = 0.5+TMath::Cos(phi);

    // data.x3    = 10*TMath::Sin(phi);
    // data.x4    = 10*TMath::Cos(phi);

    tb->Fill();
  }

  f->Write();

  f->Close();

  delete f;
}
//-----------------------------------------------------------------------------
void make_my_example_2() {

  struct data_t {
    float x1;
    float x2;
    float w;
  } data;

  TFile* f = new TFile("my_tmva_training_example_2.root","RECREATE");

  TTree* ts = new TTree("signal"    ,"S");
  TTree* tb = new TTree("background","B");

  TBranch* b_sig_x1 = ts->Branch("x1",&data.x1,"F");
  TBranch* b_sig_x2 = ts->Branch("x2",&data.x2,"F");

  TBranch* b_bgr_x1 = tb->Branch("x1",&data.x1,"F");
  TBranch* b_bgr_x2 = tb->Branch("x2",&data.x2,"F");
  TBranch* b_bgr_wt = tb->Branch("weight",&data.w,"F");

  // TBranch* b_sig_x3 = ts->Branch("x3",&data.x3,"F");
  // TBranch* b_sig_x4 = ts->Branch("x4",&data.x4,"F");

  ts->Print();

  TRandom3 rn3;

  int ntest (10000);

  double r0(1.);

  for (int i=0; i<ntest; i++) {

    double phi = TMath::Pi()*2*rn3.Rndm(i);
    double r   = r0*rn3.Rndm(i);

    data.x1    = r*TMath::Sin(phi);
    data.x2    = r*TMath::Cos(phi);

    // data.x3    = 10*TMath::Sin(phi);
    // data.x4    = 10*TMath::Cos(phi);

    ts->Fill();
  }


  double xb(1.5), yb(1.5);

  for (int i=0; i<ntest; i++) {
    double phi = TMath::Pi()*2*rn3.Rndm(i);
    double r   = r0*rn3.Rndm(i);

    data.x1    = xb+r*TMath::Sin(phi);
    data.x2    = yb+r*TMath::Cos(phi);

    data.w     = 1;

    // data.x3    = 10*TMath::Sin(phi);
    // data.x4    = 10*TMath::Cos(phi);

    tb->Fill();
  }

  f->Write();

  f->Close();

  delete f;
}

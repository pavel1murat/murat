//
//-----------------------------------------------------------------------------
void print_canvas_with_date(TCanvas* C, const char* Name) {
  TDatime d;
  C->Print(Form("%i-%02i-%02i-%s.eps",d.GetYear(),d.GetMonth(),d.GetDay(),Name));
}

//-----------------------------------------------------------------------------
// calculate eigenvectors
//-----------------------------------------------------------------------------
void test_eigenv(int NEvents = 1000) {
  TRandom3 rn;

  TMatrixDSym      m(2);
  TVectorD         ei_val(2);
  TMatrixD         ei_vec(2,2);
					// eigen vectors
  TVectorD         v1(2), v2(2);

  TH2D* h2 = new TH2D("h2","h2",1000,-5,5,1000,-5,5);

  TH1D* h11 = new TH1D("h11","h11",1000,-0.0005,0.0005);
  TH1D* h12 = new TH1D("h12","h12",200,-1,1);
  TH1D* h13 = new TH1D("h13","tan-tan2",1000,-0.0005,0.0005);

  int     ncrystals;
  double  x,y;

  for (int ievent=0; ievent<NEvents; ievent++) {

    m(0,0) = 0;
    m(1,1) = 0;
    m(0,1) = 0;
    m(1,0) = 0;
//-----------------------------------------------------------------------------
// simulate a single "event"
//-----------------------------------------------------------------------------
    ncrystals = 10;

    for (int i=0; i<ncrystals; i++) {

      x = rn.Gaus(0,1.0);
      y = rn.Gaus(0,0.2);
    
      m(0,0) += x*x;
      m(0,1) += x*y;
      m(1,0) += x*y;
      m(1,1) += y*y;
      
      h2->Fill(x,y);
    }

    //    m.Print();

    TMatrixDSymEigen md(m);
    
    ei_val = md.GetEigenValues();
    
    //    ei_val.Print();
    ei_vec = md.GetEigenVectors();
    
    //    ei_vec.Print();

    v1(0) = ei_vec(0,0);
    v1(1) = ei_vec(1,0);

    //    v1.Print();
    //    printf("v1.norm2sqr = %12.5f\n",v1.Norm2Sqr());

    v2(0) = ei_vec(0,1);
    v2(1) = ei_vec(1,1);

    //    v2.Print();
    //    printf("v2.norm2sqr = %12.5f\n",v2.Norm2Sqr());

    double dot = v1(0)*v2(0)+v1(1)*v2(1);

    double tan = v1(1)/v1(0);

    double phi2 = 0.5*atan2(2*m(0,1),m(0,0)-m(1,1));

    double tan2 = tan(phi2);

    h11->Fill(dot);

    h12->Fill(tan);
    h13->Fill(tan-tan2);
  
  }

  TCanvas* c = new TCanvas("c","c",0,0,1400,800);

  c->Divide(2,2);

  c->cd(1);
  h2->Draw();

  c->cd(2);
  h11->Draw();

  c->cd(3);
  h12->Draw();

  c->cd(4);
  h13->Draw();


  print_canvas_with_date(c,"test_eigen");				// 
}

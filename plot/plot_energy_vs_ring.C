//

void plot_energy_vs_ring(const char* FnLyso="hist/egun_lyso_30_110_tcalm006.hist", 
			 const char* FnBaF2="hist/egun_baf2_30_180_tcalm006.hist") {

  double   x[10], e1[10], esum1[10], e2[10], esum2[10];

  char  name[100];

  TH1F  *h1, *h2;

  for (int i=0; i<10; i++) {
//-----------------------------------------------------------------------------
// energy per ring
//-----------------------------------------------------------------------------
    sprintf(name,"evt_0/e%i",i);

    h1    = gh1(FnLyso,"TCalm006",name);
    e1[i] = h1->GetMean();

    h2 = gh1(FnBaF2,"TCalm006",name);
    e2[i] = h2->GetMean();
//-----------------------------------------------------------------------------
// sum energy
//-----------------------------------------------------------------------------
    sprintf(name,"evt_0/esum%i",i);

    h1 = gh1(FnLyso,"TCalm006",name);
    esum1[i] = h1->GetMean();

    h2 = gh1(FnBaF2,"TCalm006",name);
    esum2[i] = h2->GetMean();

    x[i] = i;
  }

  TCanvas* c1 = new_slide("c1","c1",2,1,1100,600);

  TPad* p1 = (TPad*) c1->GetPrimitive("p1");

  p1->cd(1);

  TGraph* gr1 = new TGraph(10,x,e1);
  gr1->SetName("lyso");
  TGraph* gr2 = new TGraph(10,x,e2);
  gr2->SetName("baf2");


  TH2F* h3 = new TH2F("h3","h3",110,0,11,1100,0,110);
  h3->GetYaxis()->SetRangeUser(0.1,110);
  //  gPad->SetLogy(1);


  h3->Draw();

  gr1->Draw("LP");
  gr2->Draw("LP");

  p1->cd(2);

  TGraph* gr3 = new TGraph(10,x,esum1);
  gr3->SetName("lyso_sum");

  TGraph* gr4 = new TGraph(10,x,esum2);
  gr4->SetName("baf2_sum");

  TH2F* h4 = new TH2F("h3","h3",110,0,11,1100,0,110);
  h4->GetYaxis()->SetRangeUser(60.,109.99);
  h4->Draw();
  
  gr3->Draw("lp");
  gr4->Draw("lp");
}

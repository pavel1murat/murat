//

void plot_esum_vs_length_50() {

  double   x[10], esum_lyso[10], esum_11[10], esum_15[10], esum_18[10], esum_20[10];

  char  name[100];

  TH1F  *h1, *h2;

    char fnlyso       [] = "hist/egun_lyso_30_110_tcalm006.hist";

    char fnbaf2_50_110[] = "hist/egun_baf2_50_110_tcalm006.hist";
    char fnbaf2_50_150[] = "hist/egun_baf2_50_150_tcalm006.hist";
    char fnbaf2_50_180[] = "hist/egun_baf2_50_180_tcalm006.hist";
    char fnbaf2_50_200[] = "hist/egun_baf2_50_200_tcalm006.hist";

  for (int i=0; i<10; i++) {
//-----------------------------------------------------------------------------
// energy per ring
//-----------------------------------------------------------------------------
    sprintf(name,"evt_0/esum%i",i);

    h1    = gh1(fnlyso,"TCalm006",name);
    esum_lyso[i] = h1->GetMean();

    h1 = gh1(fnbaf2_50_110,"TCalm006",name);
    esum_11[i] = h1->GetMean();
    h1 = gh1(fnbaf2_50_150,"TCalm006",name);
    esum_15[i] = h1->GetMean();
    h1 = gh1(fnbaf2_50_180,"TCalm006",name);
    esum_18[i] = h1->GetMean();
    h1 = gh1(fnbaf2_50_200,"TCalm006",name);
    esum_20[i] = h1->GetMean();

    x[i] = i;
  }

  TGraph* gr_lyso = new TGraph(10,x,esum_lyso);
  gr_lyso->SetName("lyso");
  TGraph* gr_11 = new TGraph(10,x,esum_11);
  gr_11->SetName("baf2_11");
  TGraph* gr_15 = new TGraph(10,x,esum_15);
  gr_15->SetName("baf2_15");
  TGraph* gr_18 = new TGraph(10,x,esum_18);
  gr_18->SetName("baf2_18");
  TGraph* gr_20 = new TGraph(10,x,esum_20);
  gr_20->SetName("baf2_20");


  TH2F* h3 = new TH2F("h3","h3",110,0,11,1100,0,110);
  h3->GetYaxis()->SetRangeUser(0.1,110);
  h3->SetTitle("Collected Energy vs number of used rings, BaF2: 50mm HEX");
  //  gPad->SetLogy(1);

  h3->GetXaxis()->SetRangeUser(0,9);
  h3->GetYaxis()->SetRangeUser(60,104.99);
  h3->GetXaxis()->SetTitle("Number of used rings");
  h3->GetYaxis()->SetTitle("Energy, MeV");
  h3->SetStats(0);
  h3->Draw();

  gr_lyso->SetLineColor(2);
  gr_lyso->SetLineWidth(2);
  gr_lyso->Draw("LP");

  gr_11->Draw("LP");
  gr_15->SetLineColor(3);
  gr_15->Draw("LP");
  gr_18->SetLineColor(4);
  gr_18->Draw("LP");
  gr_20->SetLineColor(6);
  gr_20->Draw("LP");

  TLegend* leg = new TLegend(0.5,0.15,0.8,0.35);
  leg->SetFillStyle(0);

  leg->AddEntry(gr_lyso,"LYSO 30x110 mm","l");
  leg->AddEntry(gr_11  ,"BaF2 50x110 mm","l");
  leg->AddEntry(gr_15  ,"BaF2 50x150 mm","l");
  leg->AddEntry(gr_18  ,"BaF2 50x180 mm","l");
  leg->AddEntry(gr_20  ,"BaF2 50x200 mm","l");

  leg->Draw();

}


//-----------------------------------------------------------------------------
void plot_esum_vs_length_30() {

  double   x[10], esum_lyso[10], esum_11[10], esum_15[10], esum_18[10], esum_20[10];

  char  name[100];

  TH1F  *h1, *h2;

    char fnlyso       [] = "hist/egun_lyso_30_110_tcalm006.hist";

    char fnbaf2_50_110[] = "hist/egun_baf2_30_110_tcalm006.hist";
    char fnbaf2_50_150[] = "hist/egun_baf2_30_150_tcalm006.hist";
    char fnbaf2_50_180[] = "hist/egun_baf2_30_180_tcalm006.hist";
    char fnbaf2_50_200[] = "hist/egun_baf2_30_200_tcalm006.hist";

  for (int i=0; i<10; i++) {
//-----------------------------------------------------------------------------
// energy per ring
//-----------------------------------------------------------------------------
    sprintf(name,"evt_0/esum%i",i);

    h1    = gh1(fnlyso,"TCalm006",name);
    esum_lyso[i] = h1->GetMean();

    h1 = gh1(fnbaf2_50_110,"TCalm006",name);
    esum_11[i] = h1->GetMean();
    h1 = gh1(fnbaf2_50_150,"TCalm006",name);
    esum_15[i] = h1->GetMean();
    h1 = gh1(fnbaf2_50_180,"TCalm006",name);
    esum_18[i] = h1->GetMean();
    h1 = gh1(fnbaf2_50_200,"TCalm006",name);
    esum_20[i] = h1->GetMean();

    x[i] = i;
  }

  TGraph* gr_lyso = new TGraph(10,x,esum_lyso);
  gr_lyso->SetName("lyso");
  TGraph* gr_11 = new TGraph(10,x,esum_11);
  gr_11->SetName("baf2_11");
  TGraph* gr_15 = new TGraph(10,x,esum_15);
  gr_15->SetName("baf2_15");
  TGraph* gr_18 = new TGraph(10,x,esum_18);
  gr_18->SetName("baf2_18");
  TGraph* gr_20 = new TGraph(10,x,esum_20);
  gr_20->SetName("baf2_20");


  TH2F* h3 = new TH2F("h3","h3",110,0,11,1100,0,110);
  h3->GetYaxis()->SetRangeUser(0.1,110);
  h3->SetTitle("Collected Energy vs number of used rings, BaF2: 30mm HEX");
  //  gPad->SetLogy(1);

  h3->GetXaxis()->SetRangeUser(0,9);
  h3->GetYaxis()->SetRangeUser(60,104.99);
  h3->GetXaxis()->SetTitle("Number of used rings");
  h3->GetYaxis()->SetTitle("Energy, MeV");

  h3->SetStats(0);
  h3->Draw();

  gr_lyso->SetLineColor(2);
  gr_lyso->SetLineWidth(2);
  gr_lyso->Draw("LP");

  gr_11->Draw("LP");
  gr_15->SetLineColor(3);
  gr_15->Draw("LP");
  gr_18->SetLineColor(4);
  gr_18->Draw("LP");
  gr_20->SetLineColor(6);
  gr_20->Draw("LP");

  TLegend* leg = new TLegend(0.5,0.15,0.8,0.35);
  leg->SetFillStyle(0);

  leg->AddEntry(gr_lyso,"LYSO 30x110 mm","l");
  leg->AddEntry(gr_11  ,"BaF2 30x110 mm","l");
  leg->AddEntry(gr_15  ,"BaF2 30x150 mm","l");
  leg->AddEntry(gr_18  ,"BaF2 30x180 mm","l");
  leg->AddEntry(gr_20  ,"BaF2 30x200 mm","l");

  leg->Draw();

}



//-----------------------------------------------------------------------------
void plot_esum_vs_length() {

  TCanvas* c1 = new_slide("c1","c1",2,1,1100,600);

  TPad* p1 = (TPad*) c1->GetPrimitive("p1");
  
  p1->cd(1);
  plot_esum_vs_length_50();

  p1->cd(2);
  plot_esum_vs_length_30();
  
}

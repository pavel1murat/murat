//

void plot_esum9_vs_length_30() {

  double   x[10], esum_lyso[10], esum_11[10], esum_15[10], esum_18[10], esum_20[10];

  char  name[100];

  TH1F  *h_lyso, *h_11, *h_15, *h_18, *h_20;

  char fnlyso       [] = "hist/egun_lyso_30_110_tcalm006.hist";
  
  char fnbaf2_30_110[] = "hist/egun_baf2_30_110_tcalm006.hist";
  char fnbaf2_30_150[] = "hist/egun_baf2_30_150_tcalm006.hist";
  char fnbaf2_30_180[] = "hist/egun_baf2_30_180_tcalm006.hist";
  char fnbaf2_30_200[] = "hist/egun_baf2_30_200_tcalm006.hist";

//-----------------------------------------------------------------------------
// energy per ring
//-----------------------------------------------------------------------------
  sprintf(name,"evt_0/esum9");

  h_lyso    = gh1(fnlyso,"TCalm006",name);

  h_11 = gh1(fnbaf2_30_110,"TCalm006",name);
  h_15 = gh1(fnbaf2_30_150,"TCalm006",name);
  h_18 = gh1(fnbaf2_30_180,"TCalm006",name);
  h_20 = gh1(fnbaf2_30_200,"TCalm006",name);

  h_20->SetStats(0);
  h_20->SetMaximum(350);
  h_20->Draw();

  h_lyso->SetLineColor(2);
  h_lyso->SetLineWidth(2);
  h_lyso->Draw("same");

  h_18->SetLineColor(3);
  h_18->Draw("same");
  h_15->SetLineColor(4);
  h_15->Draw("same");
  h_11->SetLineColor(6);
  h_11->Draw("same");

  TLegend* leg = new TLegend(0.15,0.15,0.45,0.35);
  leg->SetFillStyle(0);

  leg->AddEntry(h_lyso,"LYSO 30x110 mm","l");
  leg->AddEntry(h_11  ,"BaF2 30x110 mm","l");
  leg->AddEntry(h_15  ,"BaF2 30x150 mm","l");
  leg->AddEntry(h_18  ,"BaF2 30x180 mm","l");
  leg->AddEntry(h_20  ,"BaF2 30x200 mm","l");

  leg->Draw();

}
void plot_esum9_vs_length_50() {

  double   x[10], esum_lyso[10], esum_11[10], esum_15[10], esum_18[10], esum_20[10];

  char  name[100];

  TH1F  *h_lyso, *h_11, *h_15, *h_18, *h_20;

  char fnlyso       [] = "hist/egun_lyso_30_110_tcalm006.hist";
  
  char fnbaf2_50_110[] = "hist/egun_baf2_50_110_tcalm006.hist";
  char fnbaf2_50_150[] = "hist/egun_baf2_50_150_tcalm006.hist";
  char fnbaf2_50_180[] = "hist/egun_baf2_50_180_tcalm006.hist";
  char fnbaf2_50_200[] = "hist/egun_baf2_50_200_tcalm006.hist";

//-----------------------------------------------------------------------------
// energy per ring
//-----------------------------------------------------------------------------
  sprintf(name,"evt_0/esum9");

  h_lyso    = gh1(fnlyso,"TCalm006",name);

  h_11 = gh1(fnbaf2_50_110,"TCalm006",name);
  h_15 = gh1(fnbaf2_50_150,"TCalm006",name);
  h_18 = gh1(fnbaf2_50_180,"TCalm006",name);
  h_20 = gh1(fnbaf2_50_200,"TCalm006",name);

  h_20->SetStats(0);
  h_20->SetMaximum(350);
  h_20->Draw();

  h_lyso->SetLineColor(2);
  h_lyso->SetLineWidth(2);
  h_lyso->Draw("same");

  h_18->SetLineColor(3);
  h_18->Draw("same");
  h_15->SetLineColor(4);
  h_15->Draw("same");
  h_11->SetLineColor(6);
  h_11->Draw("same");

  TLegend* leg = new TLegend(0.15,0.15,0.45,0.35);
  leg->SetFillStyle(0);

  leg->AddEntry(h_lyso,"LYSO 30x110 mm","l");
  leg->AddEntry(h_11  ,"BaF2 50x110 mm","l");
  leg->AddEntry(h_15  ,"BaF2 50x150 mm","l");
  leg->AddEntry(h_18  ,"BaF2 50x180 mm","l");
  leg->AddEntry(h_20  ,"BaF2 50x200 mm","l");

  leg->Draw();

}





//-----------------------------------------------------------------------------
void plot_esum9_vs_length() {

  TCanvas* c1 = new_slide("c1","c1",2,1,1200,700);

  TPad* p1 = (TPad*) c1->GetPrimitive("p1");
  
  p1->cd(1);
  plot_esum9_vs_length_50();

  p1->cd(2);
  plot_esum9_vs_length_30();
  
}

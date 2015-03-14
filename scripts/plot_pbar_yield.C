//


//-----------------------------------------------------------------------------
void plot_pbar_yield(float P0, float Theta) {

  float  p[1000], xs[1000];
  double e;

  int np = 160;

  double const MP = 0.93825; // in GeV

  for (int i=0; i<np; i++) {

    p[i]  = i*0.05+1.e-12; // to avoid zeroes

    // calculate correction factor

    e = sqrt(p[i]*p[i]+MP*MP);

    double cf = 1.539e6*e/(p[i]*p[i])/(2*TMath::Pi()*TMath::Sin(Theta));

    xs[i] = get_pbar_yield(P0,p[i],Theta)*cf;

    printf(" i, p, xs : %3i   %10.3f  %12.5e\n",i,p[i],xs[i]);
  }


  TGraph* g = new TGraph(np,p,xs);

  g->SetName(Form("gr_%.0f_%.0f",P0,Theta*180/TMath::Pi()));
  g->SetTitle(Form("MARS #bar{p} yield, P0 =%5.1f GeV/c, #theta_{lab} =%6.1f deg",
		   P0,Theta*180/TMath::Pi()));
  g->GetXaxis()->SetRangeUser(0,7.999);
  g->GetYaxis()->SetRangeUser(1.e-15,5.e2);
  g->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}, #mub/(GeV/c)^{2}");
  g->GetYaxis()->SetTitleOffset(1.3);

  gPad->SetLogy(1);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  g->Draw("ALP");
}

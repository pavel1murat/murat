//-----------------------------------------------------------------------------
/*

****************************************
Minimizer is Linear
Chi2                      =     0.533647
NDf                       =            3
p0                        =      85.2736   +/-   12.8564     
p1                        =     -1069.48   +/-   282.23      
p2                        =      6412.45   +/-   2287.49     
p3                        =       -21557   +/-   9184.41     
p4                        =      42503.2   +/-   20199.8     
p5                        =     -48701.2   +/-   24751       
p6                        =      29955.5   +/-   15843       
p7                        =     -7637.13   +/-   4126.39     

*/
void plot_muon_decay_rad_corr() {

  int const npt = 11;

  float eta   [npt] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.96, 0.98  };

  float corr  [npt] = {24.68, 9.69, 5.54, 3.433, 2.01, 0.85, -0.23, -1.42, -3.06, -4.16, -6.45 };


  TGraph* gr = new TGraph(npt,eta,corr);

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize (1);

  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("rad corr");

  //  gr->GetYaxis()->SetRangeUser(0.05,0.12);
  gr->Draw("alp");

  // TText* txt = new TText(0.03,0.115,"combined TrkPatRec+CalPatRec reconstruction efficiency for CE");
  // txt->SetTextSize(0.04);
  // txt->SetTextFont(52);

  // txt->Draw();
}

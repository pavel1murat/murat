//

namespace {
  int cnv00202_mergepatrec[10] = { 
    1400, 686, 646, 646, 631, 461, 245, 220, 177, 177
  };

  int cnv00302_mergepatrec[10] = {
    1000, 452, 423, 423, 393, 278, 155, 143,  125,  125
  };

  int cnv00402_mergepatrec[10] = {
    1000, 465, 439, 439, 328, 238, 157, 141,  107,  107
  };

  char* cut_label[10] = {
    "Total"   , "N(MC)>20", "Pf>100", "CE Pitch", "Track Fit", 
    "Fit Qual", "T0>700"  , "Reco Pitch", "cosmics", "P > 103.2"
  };

};

//-----------------------------------------------------------------------------
a1_mpr() {

  TH1F* h_cnv00202_mergepatrec = new TH1F("h_cnv00202_mergepatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00202_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    float val = cnv00202_mergepatrec[i-1]/(cnv00202_mergepatrec[0]+0.);

    h_cnv00202_mergepatrec->SetBinContent(i,val);
  }

  h_cnv00202_mergepatrec->SetTitle("Mergepatrec+CalPatRec efficiency");
  h_cnv00202_mergepatrec->SetStats(0);
  h_cnv00202_mergepatrec->SetMinimum(0.);
  h_cnv00202_mergepatrec->Draw();
  h_cnv00202_mergepatrec->Draw("same,text45");


  TH1F* h_cnv00302_mergepatrec = new TH1F("h_cnv00302_mergepatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00302_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    float val = cnv00302_mergepatrec[i-1]/(cnv00302_mergepatrec[0]+0.);

    h_cnv00302_mergepatrec->SetBinContent(i,val);
  }

  h_cnv00302_mergepatrec->SetFillStyle(3003);
  h_cnv00302_mergepatrec->SetFillColor(kBlue-7);
  h_cnv00302_mergepatrec->Draw("same");


  TH1F* h_cnv00402_mergepatrec = new TH1F("h_cnv00402_mergepatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00402_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    float val = cnv00402_mergepatrec[i-1]/(cnv00402_mergepatrec[0]+0.);

    h_cnv00402_mergepatrec->SetBinContent(i,val);
  }
  h_cnv00402_mergepatrec->SetFillStyle(3001);
  h_cnv00402_mergepatrec->SetFillColor(kBlue-7);
  h_cnv00402_mergepatrec->Draw("same");


  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_cnv00202_mergepatrec,"CE only","f");
  leg->AddEntry(h_cnv00302_mergepatrec,"CE+MIXP2","f");
  leg->AddEntry(h_cnv00402_mergepatrec,"CE+MIXP3","f");

  leg->Draw();
}


//-----------------------------------------------------------------------------
a2_mpr() {

  TH1F* h_cnv00202_mergepatrec = new TH1F("h_cnv00202_mergepatrec","Efficiency",10,0,10);


    float val;
    for (int i=1; i<=10; i++) {
    h_cnv00202_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00202_mergepatrec[i-1]/(cnv00202_mergepatrec[i-2]+0.);

    h_cnv00202_mergepatrec->SetBinContent(i,val);
  }

    h_cnv00202_mergepatrec->SetMaximum(1.5);
  h_cnv00202_mergepatrec->SetTitle("Mergepatrec+CalPatRec efficiency");
  h_cnv00202_mergepatrec->SetStats(0);
  h_cnv00202_mergepatrec->SetMinimum(0.);
  h_cnv00202_mergepatrec->Draw();
  h_cnv00202_mergepatrec->Draw("same,text45");


  TH1F* h_cnv00302_mergepatrec = new TH1F("h_cnv00302_mergepatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00302_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00302_mergepatrec[i-1]/(cnv00302_mergepatrec[i-2]+0.);

    h_cnv00302_mergepatrec->SetBinContent(i,val);
  }

  h_cnv00302_mergepatrec->SetFillStyle(3003);
  h_cnv00302_mergepatrec->SetFillColor(kBlue-7);
  h_cnv00302_mergepatrec->Draw("same");


  TH1F* h_cnv00402_mergepatrec = new TH1F("h_cnv00402_mergepatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00402_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00402_mergepatrec[i-1]/(cnv00402_mergepatrec[i-2]+0.);

    h_cnv00402_mergepatrec->SetBinContent(i,val);
  }
  h_cnv00402_mergepatrec->SetFillStyle(3001);
  h_cnv00402_mergepatrec->SetFillColor(kBlue-7);
  h_cnv00402_mergepatrec->Draw("same");


  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_cnv00202_mergepatrec,"CE only","f");
  leg->AddEntry(h_cnv00302_mergepatrec,"CE+MIXP2","f");
  leg->AddEntry(h_cnv00402_mergepatrec,"CE+MIXP3","f");

  leg->Draw();
}

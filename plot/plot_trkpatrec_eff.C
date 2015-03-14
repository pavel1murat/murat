//

namespace {
  int cnv00202_trkpatrec[10] = { 
    1400, 686, 646, 646, 628, 473, 253, 226, 177, 177
  };

  int cnv00302_trkpatrec[10] = {
    1000, 452, 423, 423, 335, 238, 135, 127,  111,  111
  };

  int cnv00402_trkpatrec[10] = {
    1000, 472, 439, 439, 204, 159, 118, 110,  81,  81
  };

  char* cut_label[10] = {
    "Total"   , "N(MC)>20", "Pf>100", "CE Pitch", "Track Fit", 
    "Fit Qual", "T0>700"  , "Reco Pitch", "cosmics", "P > 103.2"
  };

};

//-----------------------------------------------------------------------------
a1() {

  TH1F* h_cnv00202_trkpatrec = new TH1F("h_cnv00202_trkpatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00202_trkpatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    float val = cnv00202_trkpatrec[i-1]/(cnv00202_trkpatrec[0]+0.);

    h_cnv00202_trkpatrec->SetBinContent(i,val);
  }

  h_cnv00202_trkpatrec->SetTitle("TrkPatRec efficiency");
  h_cnv00202_trkpatrec->SetStats(0);
  h_cnv00202_trkpatrec->SetMinimum(0.);
  h_cnv00202_trkpatrec->Draw();
  h_cnv00202_trkpatrec->Draw("same,text45");


  TH1F* h_cnv00302_trkpatrec = new TH1F("h_cnv00302_trkpatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00302_trkpatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    float val = cnv00302_trkpatrec[i-1]/(cnv00302_trkpatrec[0]+0.);

    h_cnv00302_trkpatrec->SetBinContent(i,val);
  }

  h_cnv00302_trkpatrec->SetFillStyle(3003);
  h_cnv00302_trkpatrec->SetFillColor(kBlue-7);
  h_cnv00302_trkpatrec->Draw("same");


  TH1F* h_cnv00402_trkpatrec = new TH1F("h_cnv00402_trkpatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00402_trkpatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    float val = cnv00402_trkpatrec[i-1]/(cnv00402_trkpatrec[0]+0.);

    h_cnv00402_trkpatrec->SetBinContent(i,val);
  }
  h_cnv00402_trkpatrec->SetFillStyle(3001);
  h_cnv00402_trkpatrec->SetFillColor(kBlue-7);
  h_cnv00402_trkpatrec->Draw("same");


  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_cnv00202_trkpatrec,"CE only","f");
  leg->AddEntry(h_cnv00302_trkpatrec,"CE+MIXP2","f");
  leg->AddEntry(h_cnv00402_trkpatrec,"CE+MIXP3","f");

  leg->Draw();
}


//-----------------------------------------------------------------------------
a2() {

  TH1F* h_cnv00202_trkpatrec = new TH1F("h_cnv00202_trkpatrec","Efficiency",10,0,10);


    float val;
    for (int i=1; i<=10; i++) {
    h_cnv00202_trkpatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00202_trkpatrec[i-1]/(cnv00202_trkpatrec[i-2]+0.);

    h_cnv00202_trkpatrec->SetBinContent(i,val);
  }

    h_cnv00202_trkpatrec->SetMaximum(1.5);
  h_cnv00202_trkpatrec->SetTitle("TrkPatRec efficiency");
  h_cnv00202_trkpatrec->SetStats(0);
  h_cnv00202_trkpatrec->SetMinimum(0.);
  h_cnv00202_trkpatrec->Draw();
  h_cnv00202_trkpatrec->Draw("same,text45");


  TH1F* h_cnv00302_trkpatrec = new TH1F("h_cnv00302_trkpatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00302_trkpatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00302_trkpatrec[i-1]/(cnv00302_trkpatrec[i-2]+0.);

    h_cnv00302_trkpatrec->SetBinContent(i,val);
  }

  h_cnv00302_trkpatrec->SetFillStyle(3003);
  h_cnv00302_trkpatrec->SetFillColor(kBlue-7);
  h_cnv00302_trkpatrec->Draw("same");


  TH1F* h_cnv00402_trkpatrec = new TH1F("h_cnv00402_trkpatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00402_trkpatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00402_trkpatrec[i-1]/(cnv00402_trkpatrec[i-2]+0.);

    h_cnv00402_trkpatrec->SetBinContent(i,val);
  }
  h_cnv00402_trkpatrec->SetFillStyle(3001);
  h_cnv00402_trkpatrec->SetFillColor(kBlue-7);
  h_cnv00402_trkpatrec->Draw("same");


  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_cnv00202_trkpatrec,"CE only","f");
  leg->AddEntry(h_cnv00302_trkpatrec,"CE+MIXP2","f");
  leg->AddEntry(h_cnv00402_trkpatrec,"CE+MIXP3","f");

  leg->Draw();
}

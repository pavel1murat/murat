//


//-----------------------------------------------------------------------------
void compare_stnana_hist(const char* HistName = "trk_1/p",double XMin = 1, double XMax = -1) {

  const char* fn_tcalm002 = "hist/conversion_01_tcalm002.hist.1000";
  const char* fn_stnana   = "hist/conversion_01_stnana.hist.1000";


  TString name = HistName;

  name.ReplaceAll("/","_");

  TCanvas* c = new TCanvas(Form("c_%s",name.Data()),name.Data(),900,600);

  c->Draw();
  
  TH1F  *h_tcalm, *h_stn;

  h_tcalm = (TH1F*) gh1(fn_tcalm002,"TCalm002",HistName)->Clone(Form("%s_tcalm002",name.Data()));
  h_stn   = (TH1F*) gh1(fn_stnana  ,"TrackAna",HistName)->Clone(Form("%s_stnana",name.Data()));

  h_tcalm->SetMarkerStyle(20);
  h_tcalm->SetMarkerSize (0.8);
  h_tcalm->SetLineColor(1);

  if (XMax > XMin) h_tcalm->GetXaxis()->SetRangeUser(XMin,XMax);

  h_tcalm->Draw("pe");
  gPad->Update();

  TPaveStats *ps1 = (TPaveStats*) h_tcalm->FindObject("stats");
  ps1->SetX1NDC(0.12);
  ps1->SetX2NDC(0.3);
  ps1->SetY1NDC(0.7);
  ps1->SetY2NDC(0.9);


  h_stn->SetLineColor(kRed-3);
  h_stn->SetFillColor(kRed-3);
  h_stn->SetFillStyle(3002);

  if (XMax > XMin) h_stn->GetXaxis()->SetRangeUser(XMin,XMax);

  h_stn->Draw("sames");
  gPad->Update();

  TPaveStats *ps2 = (TPaveStats*) h_stn->FindObject("stats");
  ps2->SetX1NDC(0.7);
  ps2->SetX2NDC(0.88);
  ps2->SetY1NDC(0.7);
  ps2->SetY2NDC(0.9);


  TLegend* leg = new TLegend(0.35,0.77,0.5,0.85,"");
  leg->AddEntry(h_tcalm,"TCalm002","pe");
  leg->AddEntry(h_stn  ,"Stnana"  ,"f");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);

  leg->Draw();

  c->Modified();
  c->Update();


  c->Print(Form("validate_stnana_%s.eps",name.Data()));
  c->Print(Form("validate_stnana_%s.png",name.Data()));
}


//-----------------------------------------------------------------------------
void validate_stnana(const char* Folder) {

  TString fol = Folder;

  if (fol.Index("trk") == 0) {
    compare_stnana_hist(Form("%s/p2"      ,Folder),90.,110.);
    compare_stnana_hist(Form("%s/chi2"    ,Folder));
    compare_stnana_hist(Form("%s/ep"      ,Folder));
    compare_stnana_hist(Form("%s/tdip"    ,Folder));
    compare_stnana_hist(Form("%s/dx"      ,Folder));
    compare_stnana_hist(Form("%s/dy"      ,Folder));
    compare_stnana_hist(Form("%s/dt"      ,Folder),-20,15);
    compare_stnana_hist(Form("%s/t0"      ,Folder));
    compare_stnana_hist(Form("%s/z0"      ,-3000,2000));
    compare_stnana_hist(Form("%s/ncl"     ,Folder));
    compare_stnana_hist(Form("%s/llhr_cal",Folder));
  }
  else if (fol.Index("cls") == 0) {
//-----------------------------------------------------------------------------
// fol = "cls_0" or such
//-----------------------------------------------------------------------------
    compare_stnana_hist(Form("%s/vane_id"   ,Folder),-4,6);
    compare_stnana_hist(Form("%s/energy"    ,Folder));
    compare_stnana_hist(Form("%s/t0"        ,Folder));
    compare_stnana_hist(Form("%s/x"         ,Folder));
    compare_stnana_hist(Form("%s/y"         ,Folder));
    compare_stnana_hist(Form("%s/z"         ,Folder));
    compare_stnana_hist(Form("%s/ymean"     ,Folder));
    compare_stnana_hist(Form("%s/ncr0"      ,Folder));
    compare_stnana_hist(Form("%s/ncr1"      ,Folder));
    compare_stnana_hist(Form("%s/fre1"      ,Folder),0,2);
    compare_stnana_hist(Form("%s/fre2"      ,Folder),0,2);
    compare_stnana_hist(Form("%s/sige1"     ,Folder),0,2);
    compare_stnana_hist(Form("%s/sige2"     ,Folder),0,2);
  }
  else if (fol.Index("rc") == 0) {
//-----------------------------------------------------------------------------
// "rc", no folder, just histograms "rc_0" or "rc_1"
//-----------------------------------------------------------------------------
    compare_stnana_hist("rc_0");
    compare_stnana_hist("rc_1");
  }
  else if (fol.Index("cal") == 0) {
//-----------------------------------------------------------------------------
// fol = "cal_0" or such
//-----------------------------------------------------------------------------
    compare_stnana_hist(Form("%s/vane_id"   ,Folder));
    compare_stnana_hist(Form("%s/energy_0"  ,Folder));
    compare_stnana_hist(Form("%s/time_0"    ,Folder));
    compare_stnana_hist(Form("%s/r_0"       ,Folder));
    compare_stnana_hist(Form("%s/rwe_0"     ,Folder));
    compare_stnana_hist(Form("%s/energy_1"  ,Folder));
    compare_stnana_hist(Form("%s/time_1"    ,Folder));
    compare_stnana_hist(Form("%s/r_1"       ,Folder));
    compare_stnana_hist(Form("%s/rwe_1"     ,Folder));
  }
  

  //  compare_stnana_hist("");
}

//-----------------------------------------------------------------------------
// count the number of wires with hits in per-panel occupancy histograms
// filled by murat/ana/TPanelDtAnaModule.cc
//-----------------------------------------------------------------------------
int count_active_wires(const char* FileName) {
  int rc(0);

  TFile* f = TFile::Open(FileName);

  int nwires(0), nwtot(0);

  for (int slot=0; slot<18; slot++) {
    for (int panel=0; panel<12; panel++) {
      std::string hname = std::format("//Ana/PanelDtAna/sel_01/slot_{:02d}/pnl_{:03}/occup",slot,panel);
      // std::cout << std::format("hname:{}\n",hname);
      TH1* h = (TH1*) f->Get(hname.data());
      h->Draw();
      
      for (int bin=1; bin<97; bin++) {
        int qn = h->GetBinContent(bin);
        nwtot++;
        if (qn > 0) {
          nwires++;
        }
      }
      // return 0;
    }
  }

  std::cout << std::format("nwtot:{} nwires:{} fraction:{}\n",nwtot,nwires,nwires/(nwtot+1.e-12));
  return rc;
}

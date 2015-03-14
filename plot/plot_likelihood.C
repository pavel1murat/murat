//
plot_likelihood() {

  gInterpreter->LoadMacro("murat/scripts/plot_ele_lh_1.C");

  plot_ele_lh_1(100000);

  Hist.fDelLH[0]->Draw();

  c1->Print("log_lhr_electrons.eps");

  Hist.fDelLH[1]->Draw();

  c1->Print("log_lhr_muons.eps");
}

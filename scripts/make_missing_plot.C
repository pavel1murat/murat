//

void make_missing_plot() {

  TCanvas* c = new TCanvas("c_missing_plot", "Missing Plot");
  TLatex* txt = new TLatex(0.25,0.45,"missing plot");
  txt->SetTextSize(0.15);
  txt->Draw();
}

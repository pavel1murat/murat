//
void plot_beam_timing() {
  TTree* t = new TTree("t","t");
  t->ReadFile("ConditionsService/data/potTimingDistribution_20160511.txt","t/F:w/F",' ');
  t->Draw("t>>h(300,-300,300)","w","l");
}

// test_panel.C
// Panel view: the panel envelope is the TOP volume of the geometry,
// and straws are daughters of the panel (same nesting as
// ConstructTrackerDetail5.cc: straws placed in StrawPanelEnvelope).
// Run with:  root -l test_panel.C

void test_panel()
{
  TGeoManager *geom = new TGeoManager("panel", "Straw panel as top volume");
  TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
  TGeoMedium   *med = new TGeoMedium("Vacuum", 1, mat);

  // Top volume = panel envelope (box), 2 layers x 8 straws per layer
  const double strawR   = 2.5;    // straw outer radius
  const double strawHL  = 50.;    // straw half-length (along z)
  const double pitchX   = 6.;     // straw spacing in x (within layer)
  const double pitchY   = 6.;     // layer spacing in y
  const int    nPerLay  = 8;
  const int    nLayers  = 2;

  const double panHX = nPerLay * pitchX / 2. + strawR;
  const double panHY = nLayers * pitchY / 2. + strawR;
  const double panHZ = strawHL + 2.;

  TGeoVolume *panel = geom->MakeBox("PANEL", med, panHX, panHY, panHZ);
  panel->SetLineColor(kGray);
  panel->SetTransparency(60);
  geom->SetTopVolume(panel);

  // Straw: gas tube, placed as daughter of the panel
  TGeoTube   *strawShape = new TGeoTube("StrawGas", 0., strawR, strawHL);
  TGeoVolume *straw      = new TGeoVolume("STRAW", strawShape, med);
  straw->SetLineColor(kRed);
  straw->SetFillColor(kRed);
  straw->SetTransparency(0);

  int copy = 0;
  for (int l = 0; l < nLayers; ++l) {
    for (int s = 0; s < nPerLay; ++s) {
      double x = -panHX + strawR + s * pitchX + (l % 2) * pitchX / 2.;
      double y = (l - (nLayers - 1) / 2.) * pitchY;   // centered in panel
      panel->AddNode(straw, copy++, new TGeoTranslation(x, y, 0.));
    }
  }

  geom->CloseGeometry();

  gROOT->SetBatch(1);
  TCanvas *canvas = new TCanvas("c_panel", "Straw panel", 800, 800);
  gROOT->SetBatch(0);
  canvas->cd();

  // Warm-up draw: first GL draw in a session does not paint; draw twice.
  panel->Draw("ogl");
  canvas->Update();
  gSystem->ProcessEvents();
  panel->Draw("ogl");
  canvas->Update();
}

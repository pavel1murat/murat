// test_tubes.C
// Two concentric TGeoTube volumes rendered in the ROOT OpenGL viewer.
// Run with:  root -l test_tubes.C

void test_tubes()
{

  TGeoManager *geom = new TGeoManager("tubes", "Two concentric tubes");
   TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
   TGeoMedium   *med = new TGeoMedium("Vacuum", 1, mat);

   // Top (world) volume: a box that contains both tubes
   TGeoVolume *top = geom->MakeBox("TOP", med, 40, 40, 60);
   geom->SetTopVolume(top);

  //geom->SetTopVisible(kTRUE);   // do not draw the world box
   // Inner tube: rmin = 5, rmax = 10, half-length dz = 50 (red)
   TGeoTube   *innerShape = new TGeoTube("InnerTube", 5., 10., 50.);
   TGeoVolume *innerVol   = new TGeoVolume("INNER", innerShape, med);
   innerVol->SetLineColor(kRed);
   innerVol->SetFillColor(kRed);
   innerVol->SetTransparency(0);

   // Outer tube: rmin = 15, rmax = 25, same axis (blue)
   TGeoTube   *outerShape = new TGeoTube("OuterTube", 15., 25., 50.);
   TGeoVolume *outerVol   = new TGeoVolume("OUTER", outerShape, med);
   outerVol->SetLineColor(kBlue);
   outerVol->SetFillColor(kBlue);
   outerVol->SetTransparency(0);

   top->AddNode(innerVol, 1);
   top->AddNode(outerVol, 1);
   geom->CloseGeometry();

//-----------------------------------------------------------------------------
  gROOT->SetBatch(1);
  TCanvas *canvas = new TCanvas("c_tubes", "Concentric tubes", 800, 800);
  gROOT->SetBatch(0);
  canvas->cd();
//-----------------------------------------------------------------------------
   // Warm-up draw: the first GL draw in a session does not paint.
   // A second draw into the same canvas paints correctly, so do both.
   top->Draw("ogl");
   canvas->Update();
   gSystem->ProcessEvents();

}

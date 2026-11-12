// test_panel.C
// Panel view: the panel envelope is the TOP volume of the geometry,
// and straws are daughters of the panel (same nesting as
// ConstructTrackerDetail5.cc: straws placed in StrawPanelEnvelope).
// Straw geometry follows Offline/GeometryService/src/TrackerMaker.cc with the
// parameters of Offline/Mu2eG4/geom/tracker_v7.txt:
//   - local panel frame: X along the straw, Y radially outwards, Z perpendicular to the panel
//   - straws are chords of the inner support cylinder: straw at radial distance r
//     has half-length sqrt(innerSupportRadius^2 - r^2)
//   - r positions (along Y): rInnermostWire + 2*manifoldHalfLength[0]*iman, manifolds 0..47,
//     layer 1 shifted by 2*layerHalfShift
//   - the two layers are offset along Z by 2*layerHalfSpacing
// Run with:  root -l test_panel.C

void test_panel()
{
  TGeoManager *geom = new TGeoManager("panel", "Straw panel as top volume");
  TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
  TGeoMedium   *med = new TGeoMedium("Vacuum", 1, mat);

  // straw wall and gas: material names as in tracker.strawMaterials (tracker_v7.txt);
  // the density is only a placeholder, the viewer does not use it
  TGeoMaterial *mylarMat = new TGeoMaterial("G4_MYLAR", 6.0, 3.0, 1.39);
  TGeoMedium   *mylar    = new TGeoMedium("G4_MYLAR", 2, mylarMat);
  TGeoMaterial *gasMat   = new TGeoMaterial("StrawGas", 6.0, 3.0, 0.0018);
  TGeoMedium   *gasMed   = new TGeoMedium("StrawGas", 3, gasMat);

  // tracker_v7.txt
  const double strawR          = 2.5235;   // tracker.strawOuterRadius
  const double strawWallThickness = 0.015; // tracker.strawWallThickness
  const double strawGap        = 1.203;    // tracker.strawGap
  const double rInnermostWire  = 380.;     // tracker.rInnermostWire
  const double innerSupportR   = 700.;     // tracker.innerSupportRadius
  const double manifoldHL0     = 3.125;    // tracker.manifoldHalfLengths[0]
  const int    nManifolds      = 48;       // tracker.manifoldsPerEnd
  const int    nLayers         = 2;

  // TrackerMaker::computeLayerSpacingAndShift
  const double rSpacing        = 2.*strawR + strawGap;
  const double layerHalfSpacing = 0.25*sqrt(3.0)*rSpacing;
  const double layerHalfShift   = 0.25*rSpacing;

  // panel envelope in the local panel frame: straws run along X, spread along Y, layers along Z
  const double yMin  = rInnermostWire - strawR;
  const double yMax  = rInnermostWire + 2.*manifoldHL0*(nManifolds-1) + 2.*layerHalfShift + strawR;
  const double yMid  = 0.5*(yMin + yMax);
  const double panHX = sqrt(innerSupportR*innerSupportR - rInnermostWire*rInnermostWire) + 2.;   // longest straw
  const double panHY = 0.5*(yMax - yMin);
  const double panHZ = layerHalfSpacing + strawR + 1.;

  TGeoVolume *panel = geom->MakeBox("PANEL", med, panHX, panHY, panHZ);
  panel->SetLineColor(kGray);
  panel->SetTransparency(60);
  // the panel has daughters, so it is a container; TGeoVolume::Draw hides containers
  // unless their visibility is set explicitly, leaving only the straws drawn
  panel->SetVisContainers(kTRUE);
  geom->SetTopVolume(panel);

  // straw axis along X: rotate the tube (axis z) by 90 degrees around y, z -> x
  TGeoRotation *alongX = new TGeoRotation("alongX", 0., 0., 0.);
  alongX->RotateY(90.);

  // each straw has its own length, so each straw gets its own shape
  int copy = 0;
  for (int l = 0; l < nLayers; ++l) {
    for (int iman = 0; iman < nManifolds; ++iman) {
      double rstraw = rInnermostWire + 2.*manifoldHL0*iman + l*2.*layerHalfShift;
      double halfLength = sqrt(innerSupportR*innerSupportR - rstraw*rstraw);
      double z = (l - 0.5*(nLayers-1))*2.*layerHalfSpacing;

      // straw wall: a thin tube, inner radius = outer radius - wall thickness (TrackerMaker.cc)
      TGeoTube   *wallShape = new TGeoTube(Form("StrawWallShape_%d",copy),
                                           strawR - strawWallThickness, strawR, halfLength);
      TGeoVolume *wall      = new TGeoVolume(Form("StrawWall_%d",copy), wallShape, mylar);
      wall->SetLineColor(kGray);
      wall->SetFillColor(kGray);
      wall->SetTransparency(0);

      // straw gas: the inside of the wall
      TGeoTube   *gasShape = new TGeoTube(Form("StrawGasShape_%d",copy), 0.,
                                          strawR - strawWallThickness, halfLength);
      TGeoVolume *gas      = new TGeoVolume(Form("StrawGas_%d",copy), gasShape, gasMed);
      gas->SetLineColor(kGray);
      gas->SetFillColor(kGray);
      gas->SetTransparency(50);
      wall->AddNode(gas, 1);
      wall->SetVisContainers(kTRUE);             // it has a daughter, see the panel above

      panel->AddNode(wall, copy, new TGeoCombiTrans(0., rstraw - yMid, z, new TGeoRotation(*alongX)));
      ++copy;
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

// draw Mu2e Calorimeter : v4_0_6

namespace {
  TGeoManager* gm(0);
};

void init_gm(const char* Fn) {
  const char default_fn[] = "/home/murat/figures/mu2e/gdml/mu2e_geometry_parsed_v5_7_6.gdml";
  const char* fn = Fn;

  if (fn == NULL) fn = Fn;

  if (gm == 0) {
    TGeoManager::Import(fn);
    gm = gGeoManager;
  }
}

//-----------------------------------------------------------------------------
void set_colors() {

  TGeoVolume  *pa1, *pa2, *pa3, *v;
  
  pa1 = gm->GetVolume("protonabs1");
  pa1->SetTransparency(1);
  pa1->SetLineColorAlpha(808,0.1);
  //  gm->GetVolume("protonabs2")->SetLineColor(804);

  v = gm->GetVolume("protonabs3");
  v->SetTransparency(1);
  v->SetLineColorAlpha(825,0.03);

  for (int i=0; i<34; i++) {
    gm->GetVolume(Form("Foil_%02i",i))->SetLineColor(20);
  }

  int pipe_color = 807;

  gm->GetVolume("CaloPipe0_0")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_1")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_2")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_3")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_4")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_5")->SetLineColor(pipe_color);

  gm->GetVolume("CaloPipe1_0")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_1")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_2")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_3")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_4")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_5")->SetLineColor(pipe_color);

}


//-----------------------------------------------------------------------------
void draw_mu2e_calorimeter() {

  init_gm();

  gm->GetVolume("VirtualDetector_EMC_Disk_0_SurfIn")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_0_SurfOut")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_0_EdgeOut")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_0_EdgeIn")->SetVisibility(0);

  gm->GetVolume("VirtualDetector_EMC_Disk_1_SurfIn")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_1_SurfOut")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_1_EdgeOut")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_1_EdgeIn")->SetVisibility(0);

  set_colors();

  gm->GetVolume("CalorimeterMother")->Draw("ogl");
}


//-----------------------------------------------------------------------------
void draw_mu2e_disk() {

  init_gm();

  gm->GetVolume("VirtualDetector_EMC_Disk_0_SurfIn")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_0_SurfOut")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_0_EdgeOut")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_0_EdgeIn")->SetVisibility(0);

  gm->GetVolume("VirtualDetector_EMC_Disk_1_SurfIn")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_1_SurfOut")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_1_EdgeOut")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_EMC_Disk_1_EdgeIn")->SetVisibility(0);

  int pipe_color = 807;

  gm->GetVolume("CaloPipe0_0")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_1")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_2")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_3")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_4")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe0_5")->SetLineColor(pipe_color);

  gm->GetVolume("CaloPipe1_0")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_1")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_2")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_3")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_4")->SetLineColor(pipe_color);
  gm->GetVolume("CaloPipe1_5")->SetLineColor(pipe_color);
  

  gm->GetVolume("DiskCalorimeter_0")->Draw("ogl");
}

//-----------------------------------------------------------------------------
// clipping: xc = -387.93
//           zlen = 1429
//-----------------------------------------------------------------------------
void draw_detector_solenoid(const char* Fn = "/home/murat/figures/mu2e/gdml/mu2e_geometry_parsed_v5_7_6.gdml") {

  init_gm(Fn);
				// turn off dirt

  // gm->GetVolume("worldDirtBottom")->SetVisibility(0);
  // gm->GetVolume("worldDirtNW")->SetVisibility(0);
  // gm->GetVolume("worldDirtSW")->SetVisibility(0);
  // gm->GetVolume("worldDirtSE")->SetVisibility(0);
  // gm->GetVolume("worldDirtNE")->SetVisibility(0);

  TGeoVolume* hall = gm->GetVolume("HallAir");

  TObjArray* list_of_nodes = hall->GetNodes();

  int n_nodes = list_of_nodes->GetEntries();

  TGeoNode    *node, *node2;
  const char  *name, *name2;

  for (int i=0; i<n_nodes; i++) {
    node  = hall->GetNode(i);

    name = node->GetName();

    printf(" -- node name: %s\n",name);

    if ((strstr(name,"DS2Vacuum") != 0) ||
	(strstr(name,"DS3Vacuum") != 0)    ) {

      node->SetVisibility(0);
      node->SetVisDaughters(1);
      node->SetVisLeaves(1);


      TObjArray* list_of_nodes_2 = node->GetNodes();
      int n_nodes_2 = list_of_nodes_2->GetEntries();

      for (int i2=0; i2<n_nodes_2; i2++) {
	node2  = (TGeoNode*) node->GetVolume()->GetNode(i2);
	name2 = node2->GetName();
	printf("            -- node2 name: %s\n",name2);
//------------------------------------------------------------------------------
// do not show virtual detectors
//------------------------------------------------------------------------------
	if (strstr(name2,"VirtualDetector") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
//------------------------------------------------------------------------------
// do not show muon beam stop
//------------------------------------------------------------------------------
	else if (strstr(name2,"MBSMother") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
	else if (strstr(name2,"VPSP") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
	else if (strstr(name2,"IFB") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
      }
    }
    else {
//------------------------------------------------------------------------------
// everything outside DS2 and DS3 is invisible
//-----------------------------------------------------------------------------
      node->SetVisibility(0);
      node->SetVisDaughters(0);
      node->SetVisLeaves(0);
    }
  }


  gm->GetVolume("TSdA4")->SetVisibility(0);
  gm->GetVolume("DiskFEB_0")->SetVisibility(0);
  gm->GetVolume("DiskFEB_1")->SetVisibility(0);

  set_colors();

  hall->Draw("ogl");
}


//-----------------------------------------------------------------------------
// the names are this is for v4_0_6
// clipping: xc   = -387.93
//           zlen = 1429
//-----------------------------------------------------------------------------
void turn_off_everything(const char* Fn = "mu2e_geometry_v5_7_6") {

  const char* fn (0);

  char fnn[100];

  const char dir[] = "/home/murat/figures/mu2e/gdml";

  sprintf(fnn,"%s.root",Fn);
  
  fn = gSystem->FindFile(dir,fnn);

  printf(" searching for fnn=%s in dir=%s, fn = .%s.\n", fnn, dir, fn);

  if (gSystem->FindFile(dir,fnn) == NULL) {
    sprintf(fnn,"%s.gdml",Fn);

    if (gSystem->FindFile(dir,fnn) == NULL) {
      printf (" ERROR: can\'t find geometry file for %s",fn);
      return;
    }
  }

  char filename[200];

  sprintf(filename,"%s/%s",dir,fnn);

  printf("filename : %s\n",filename);
  
  init_gm(filename);
				// turn off dirt 

  // gm->GetVolume("worldDirtBottom")->SetVisibility(0);
  // gm->GetVolume("worldDirtNW")->SetVisibility(0);
  // gm->GetVolume("worldDirtSW")->SetVisibility(0);
  // gm->GetVolume("worldDirtSE")->SetVisibility(0);
  // gm->GetVolume("worldDirtNE")->SetVisibility(0);

  TGeoVolume* master = gm->GetVolume("Worldd460");

  TGeoVolume* hall = gm->GetVolume("HallAire550");

  TObjArray* list_of_nodes = master->GetNodes();

  int n_nodes = list_of_nodes->GetEntries();

  TGeoNode    *node, *node2;
  const char  *name, *name2, *vname;

  for (int i=0; i<n_nodes; i++) {
    node  = master->GetNode(i);

    name  = node->GetName();
    vname = node->GetVolume()->GetName();

    printf(" -- node name: %s, volume name: %s\n",name,vname);

    node->GetVolume()->SetVisibility(0);
    node->GetVolume()->SetVisDaughters(0);
    node->GetVolume()->SetVisLeaves(0);
  }
  
  //  hall->Draw("ogl");
  master->Draw("ogl");
}

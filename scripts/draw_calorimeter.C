// draw Mu2e Calorimeter : v4_0_6

namespace {
  TGeoManager* gm(0);
};

void init_gm(const char* Fn = "/home/murat/figures/mu2e/gdml/mu2e_v4_0_6.gdml") {
  if (gm == 0) {
    gm = new TGeoManager();
    gm->Import(Fn);
  }
}

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
  

  gm->GetVolume("CalorimeterMother")->Draw("ogl");
}


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
// the names are this is for v4_0_6
// clipping: xc = -387.93
//           zlen = 1429
//-----------------------------------------------------------------------------
void draw_detector_solenoid(const char* Fn) {

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

  TGeoVolume* node;

  for (int i=0; i<n_nodes; i++) {
    node  = (TGeoVolume*) list_of_nodes->At(i);

    const char* name = node->GetName();

    //    printf(" -- node name: %s\n",name);

    if ((strstr(name,"DS2Vacuum") != 0) ||
	(strstr(name,"DS3Vacuum") != 0)    ) {
      node->SetVisibility(1);
      node->SetVisDaughters(1);
      node->SetVisLeaves(1);
    }
    else {
      node->SetVisibility(0);
      node->SetVisDaughters(0);
      node->SetVisLeaves(0);
    }
  }

  gm->GetVolume("DS1Vacuum")->SetVisibility(0);

  gm->GetVolume("DS2Vacuum")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_ST_In")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_ST_Out")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_TT_Mid")->SetVisibility(0);


  gm->GetVolume("DS3Vacuum")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_TT_MidInner")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_TT_FrontHollow")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_TT_Back")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_TT_Back")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_TT_OutSurf")->SetVisibility(0);
  gm->GetVolume("VirtualDetector_TT_InSurf")->SetVisibility(0);

  gm->GetVolume("protonabs1")->SetLineColor(807);


  TGeoVolume* mbs = gm->GetVolume("MBSMother");
  mbs->SetVisibility(0);
  mbs->SetVisDaughters(0);
  mbs->SetVisLeaves(0);

  gm->GetVolume("protonabs1")->SetLineColor(804);
  gm->GetVolume("protonabs3")->SetLineColor(808);

  gm->GetVolume("InternalNeutronAbsorber1")->SetLineColor(900); // default = 920

  gm->GetVolume("InternalNeutronAbsorber2")->SetLineColor(850); // default = 920
  gm->GetVolume("InternalNeutronAbsorber3a")->SetLineColor(860); // default = 920

  gm->GetVolume("Foil")->SetLineColor(20);

  hall->Draw("ogl");
}


//-----------------------------------------------------------------------------
// the names are this is for v4_0_6
// clipping: xc = -387.93
//           zlen = 1429
//-----------------------------------------------------------------------------
void draw_detector_solenoid_dev2(const char* Fn = "/home/murat/figures/mu2e/gdml/mu2e_dev2.gdml") {

  init_gm(Fn);
				// turn off dirt

  // gm->GetVolume("worldDirtBottom")->SetVisibility(0);
  // gm->GetVolume("worldDirtNW")->SetVisibility(0);
  // gm->GetVolume("worldDirtSW")->SetVisibility(0);
  // gm->GetVolume("worldDirtSE")->SetVisibility(0);
  // gm->GetVolume("worldDirtNE")->SetVisibility(0);

  TGeoVolume* hall = gm->GetVolume("HallAir39e0");

  TObjArray* list_of_nodes = hall->GetNodes();

  int n_nodes = list_of_nodes->GetEntries();

  TGeoNode    *node, *node2;
  const char  *name, *name2;

  for (int i=0; i<n_nodes; i++) {
    node  = hall->GetNode(i);

    name = node->GetName();

    //    printf(" -- node name: %s\n",name);

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
	//	printf("            -- node2 name: %s\n",name2);
	if (strstr(name2,"VirtualDetector") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
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
      node->SetVisibility(0);
      node->SetVisDaughters(0);
      node->SetVisLeaves(0);
    }
  }

  gm->GetVolume("protonabs184c0")->SetLineColor(804);
  gm->GetVolume("protonabs38770")->SetLineColor(808);

  gm->GetVolume("InternalNeutronAbsorber18440")->SetLineColor(900); // default = 920
  gm->GetVolume("InternalNeutronAbsorber281b0")->SetLineColor(850); // default = 920


  //  gm->GetVolume("InternalNeutronAbsorber3a")->SetLineColor(860); // default = 920

  //  gm->GetVolume("Foil")->SetLineColor(20);

  hall->Draw("ogl");
}

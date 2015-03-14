//


#include <cstdio>

//-----------------------------------------------------------------------------
// initialization
//-----------------------------------------------------------------------------
void init() {
  gSystem->Load("lib/libmu2e_GeometryService.so");
  gSystem->Load("lib/libmu2e_Mu2eG4.so");
  gSystem->Load("lib/libmurat_pet.so");

  // then .L murat/scripts/read_zubal.C and find the phantom data location

  // TZubalPhantom::Test01("/home/murat/Downloads/det_head_u2.med",40);
  // TZubalPhantom::Test02(0,128,0);
}

//-----------------------------------------------------------------------------
void read_zubal(const char* Filename, int IZ=40) {
  
  TZubalPhantom* z = new TZubalPhantom(Filename);

  TH2F* h2 = new TH2F(Form("h%03i",IZ),Form("h2 %03i",IZ),256,0,256,256,0,256);

  for (int i=0; i<256; i++) {
    for (int j=0; j<256; j++) {
      h2->SetBinContent(i,j,z->fData[IZ][i][j]);
    }
  }

  h2->Draw("col");

}

//-----------------------------------------------------------------------------
void read_zubal_3d(const char* Filename) {
  

  char zubal[128][256][256];


  FILE* f = fopen(Filename,"r");

  size_t nb = fread(zubal,sizeof(char),128*256*256,f);

  printf("nb = %i\n",nb);

  TGeoManager *geom    = new TGeoManager("geom","TOF-PET geometry");
  
  TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
  TGeoMaterial *Fe     = new TGeoMaterial("Fe",55.845,26,7.87);
  TGeoMaterial *Cu     = new TGeoMaterial("Cu",63.549,29,8.92);
  
  TGeoMedium *Air      = new TGeoMedium("Vacuum",0,vacuum);
  TGeoMedium *Iron     = new TGeoMedium("Iron",1,Fe);
  TGeoMedium *Copper   = new TGeoMedium("Copper",2,Cu);

  TGeoVolume *top      = geom->MakeBox("top",Air,100,100,100);

  geom->SetTopVolume(top);
  geom->SetTopVisible(0);
  
  double dx(0.11), dy(0.11), dz(0.14);
 
  TGeoRotation* rot = new TGeoRotation("rot",0,0,0);
  TGeoRotation* rot90 = new TGeoRotation("rot90",0,0,90);

  // TGeoVolume *box      = geom->MakeBox("box",Air,40,40,40);
  // box->SetLineColor(2);
  // box->SetVisibility(1);

  // top->AddNode(box,1,new TGeoTranslation(0.,0.,0.));

  TGeoVolume    *voxel[256], *layer[128];

  int izmin = 40; 
  int izmax = 51; 

  for (int iz=izmin; iz<izmax; iz++) {
    layer[iz] = geom->MakeBox(Form("Layer%03i",iz),Air,256*dx/2., 256*dy/2., dz/2);
    layer[iz]->SetVisibility(1);
  }

  printf(">>> created layers\n");

  TGeoVolume  *linex[128][256];

  for (int iz=izmin; iz<izmax; iz++) {
    for (int ix=0; ix<256; ix++) {
      linex[iz][ix] = geom->MakeBox(Form("LineX_%03i_%03i",iz,ix),Air,dx/2., 256*dy/2., dz/2);
      // linex[iz][ix]->SetVisibility(1);
      // linex[iz][ix]->SetVisDaughters(1);
    }
  }

  printf(">>> created lines\n");

  int copy_number[256];

  for (int ic=0; ic<256; ic++) {
    voxel[ic] = geom->MakeBox(Form("Voxel%03i",ic),Iron,dx/2., dy/2., dz/2);
    voxel[ic]->SetFillColor(ic);
    voxel[ic]->SetLineColor(ic);
    voxel[ic]->SetFillStyle(3001);
    copy_number[ic] = 0;
    //    voxel[i]->SetVisibility(1);
  }

  TGeoTranslation   *trans_y[256], *trans_x[256];

  for (int ix=0; ix<256; ix++) {
    trans_x[ix] = new TGeoTranslation((ix-127.5)*dx,0.,0.);
  }

  for (int iy=0; iy<256; iy++) {
    trans_y[iy] = new TGeoTranslation(0,(iy-127.5)*dy,0.);
  }

  printf(">>> created voxels\n");




  for (int iz=izmin; iz<izmax; iz++) {
					// create layer volume

    printf(">>>>>> iz = %3i\n",iz);

    top->AddNode(layer[iz],iz+1,new TGeoTranslation(0.,0.,(iz-63.5)*dz));

    //    for (int ix=0; ix<256; ix++) {


    int ixmin = 0;
    int ixmax = 256;

    int iymin = 0;
    int iymax = 256;

    for (int ix=ixmin; ix<ixmax; ix++) {
      printf(">>> ix = %3i\n",ix);
      
      linex[iz][ix]->AddNode(linex[iz][ix],256*iz+ix+1, trans_x[ix]);

      for (int iy=iymin; iy<iymax; iy++) {
	int ic = zubal[iz][iy][ix];

	//	printf(">>> ix, iy, ic: = %3i %3i %3i\n",ix,iy,ic);
	
	if (ic != 0) {
	  linex[iz][ix]->AddNode(voxel[ic], copy_number[ic], trans_y[iy]);
	  copy_number[ic] += 1;
	}
      }
    }
  }

  //  top->SetVisibility(0);
  geom->CloseGeometry();

  geom->SetVisLevel(4);
  //  geom->SetVisOption(0);
  
  top->Draw("ogl");
  
}

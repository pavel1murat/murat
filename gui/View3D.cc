///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TEnv.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TGeoPgon.h"
#include "TGeoTube.h"
#include "TGeoCompositeShape.h"
#include "TVirtualPad.h"

#include "TGeoManager.h"

#include "Stntuple/gui/TEvdCrvSection.hh"
#include "murat/gui/View3D.hh"

ClassImp(View3D);

//-----------------------------------------------------------------------------
View3D::View3D() : TObject() {
  const char* mu2e_subproject_dir = "/mu2e/app/users/murat/subprojects";
  InitGeometry();

/*  mu2egpvm* :
  /mu2e/app/users/murat/subprojects/su2020/crv_satellite_peak/cry32s41b0_trajectories:
  total used in directory 13408 available 5142859264
  drwxr-sr-x  2 murat mu2e 34816 May 22  2021 .
  drwxr-sr-x 13 murat mu2e  4096 Jun 27  2021 ..
  -rw-r--r--  1 murat mu2e  5368 May 22  2021 1002_000331_369087.txt
  -rw-r--r--  1 murat mu2e 37820 May 22  2021 1002_000404_096896.txt
*/
  const char* sp_dir = gEnv->GetValue("muse.SubprojectDir",mu2e_subproject_dir);
  TString dir = Form("%s/su2020/crv_satellite_peak/cry32s41b0_trajectories",sp_dir);
  ReadTrajectories(dir.Data());
  
  double x0[3] = { -390.4, 0., 500.};
  double v [3] = { 200,  800,  250};
  double w     = 30.;
  double zrange[2] = { 500, 1200};

  fTrack = new TEvdTrack("track_001",x0,v,w,zrange);
  fTrack->SetLineColor(kRed);
}

//-----------------------------------------------------------------------------
View3D::~View3D() {
   // Clear out fShapes
  delete fTrack;
  delete fListOfTrajectories;
}

//-----------------------------------------------------------------------------
int View3D::InitGeometry() {
  fGeoManager = TStnGeoManager::Instance("Mu2eGM",nullptr);

  //--- define some materials
  fMaterials.fVacuum = new TGeoMaterial("Vacuum", 0,  0,   0     );
  fMaterials.fVac1   = new TGeoMaterial("Vac1"  , 0,  0,   1.e-15);
  fMaterials.fAl     = new TGeoMaterial("Al"    , 26.98,13,2.7   );

  //--- define some media
  fMedia.fVacuum = new TGeoMedium("Vacuum"  ,1, fMaterials.fVacuum);
  fMedia.fAl     = new TGeoMedium("Aluminum",2, fMaterials.fAl    );
  fMedia.fVac1   = new TGeoMedium("Vac1"    ,3, fMaterials.fVac1  );

  //--- make one panel: first - a sector

  fMu2eGeom      = gGeoManager->MakeBox("mu2e", fMedia.fVacuum, 2500., 2500., 1500.);
  gGeoManager->SetTopVolume(fMu2eGeom);

  InitTrackerGeometry();
  // InitCalorimeterGeometry();
  InitCrvGeometry();

  //--- close the geometry
  // MT: request node-ids to be calculated.
  gGeoManager->CloseGeometry("i");

  return 0;
}

//-----------------------------------------------------------------------------
int View3D::InitTrackerGeometry() {

  // define shape components with names
  
  TGeoTubeSeg* tubs = new TGeoTubeSeg("tubs", 0., 70, 1,0,120);
  TGeoPgon*    pgon = new TGeoPgon   ("pg",-0.05,120.1,1,2);
  
  pgon->DefineSection(0, -1.001, 0, 35.001);
  pgon->DefineSection(1,  1.001, 0, 35.001);

  TGeoCompositeShape* pshape = new TGeoCompositeShape("panel", Form("%s-%s",tubs->GetName(),pgon->GetName()));
  TGeoVolume*  panel         = new TGeoVolume        ("PANEL", pshape);

  panel->SetLineColor(kBlue);
  panel->SetTransparency(90);

  TGeoVolume *panel2 = new TGeoVolume("PANEL2",pgon);
  panel2->SetLineColor(kRed);

  // TGeoTranslation tra1(0., he, he*gr);
  // TGeoTranslation tra2(0.,-he, he*gr);
  // TGeoTranslation tra3(he*gr, 0., he);
  // TGeoTranslation tra4(-he*gr,0., he);
  // TGeoTranslation tra5(he,-he*gr, 0.);
  // TGeoTranslation tra6(-he,-he*gr,0.);

  TGeoRotation rot0("rot0",  0,0,0); 
  TGeoRotation rot1("rot1",120,0,0);
  TGeoRotation rot2("rot2",240,0,0);

  TGeoRotation rot3("rot3", 60,0,0);
  TGeoRotation rot4("rot4",180,0,0);
  TGeoRotation rot5("rot5",300,0,0);

  // TGeoRotation rot2("rot2"); 
  // rot2.RotateZ(180.); // x -> y
  // rot2.RotateX(ATan(gr/(1.+gr))*180./Pi()); // y -> z

  // TGeoRotation rot3("rot3"); 
  // rot3.RotateZ(90.); // x -> y
  // rot3.RotateY(ATan((1.+gr)/gr)*180./Pi()); // z -> x

  // TGeoRotation rot4("rot4"); 
  // rot4.RotateZ(-90.); // y -> x
  // rot4.RotateY(-ATan((1.+gr)/gr)*180./Pi()); // x -> z

  // TGeoRotation rot5("rot5"); 
  // rot5.RotateY(90.); // z -> x
  // rot5.RotateZ(-ATan((1.+gr)/gr)*180./Pi()); // y -> x

  // TGeoRotation rot6("rot6"); 
  // rot6.RotateY(-90.); // x -> z
  // rot6.RotateZ(ATan((1.+gr)/gr)*180./Pi()); // x -> y

  
  TGeoVolume* tracker = gGeoManager->MakeTube("tracker", fMedia.fAl, 30., 71., 160.);
  tracker->SetFillColor(35);
  tracker->SetLineColor(35);
  tracker->SetTransparency(90.);

  fMu2eGeom->AddNode(tracker,1, new TGeoTranslation(-390.4, 0., 1020.)); // this is approximate
 
  TGeoTranslation tra0(0., 0,  0);
  TGeoCombiTrans *com1 = new TGeoCombiTrans(tra0, rot0);
  TGeoCombiTrans *com2 = new TGeoCombiTrans(tra0, rot1);
  TGeoCombiTrans *com3 = new TGeoCombiTrans(tra0, rot2);

  tracker->AddNode(panel, 1, com1);
  tracker->AddNode(panel, 2, com2);
  tracker->AddNode(panel, 3, com3);

  TGeoTranslation tra1(0., 0,  3);
  TGeoCombiTrans *com4 = new TGeoCombiTrans(tra1, rot3);
  TGeoCombiTrans *com5 = new TGeoCombiTrans(tra1, rot4);
  TGeoCombiTrans *com6 = new TGeoCombiTrans(tra1, rot5);

  tracker->AddNode(panel, 4, com4);
  tracker->AddNode(panel, 5, com5);
  tracker->AddNode(panel, 6, com6);

  //  tracker->AddNode(panel2, 2, new TGeoTranslation(0., 0., 30.));
  return 0;
}

//-----------------------------------------------------------------------------
int View3D::InitCalorimeterGeometry() {
  return 0;
}

//-----------------------------------------------------------------------------
int View3D::InitCrvGeometry() {
//   2880   7   0  0   0 1 2 0  -3904.000   2663.210  -2125.450   3000.000      9.900     25.650
//-----------------------------------------------------------------------------
// initializa sectors
//-----------------------------------------------------------------------------
  char c[10000];
  
  const char* fn = "crv_sectors.txt";
  
  FILE* f  = fopen(fn,"r");
  if (f == 0) {
    TString msg = Form("missing file %s",fn);
    Error("Init",(const char*) msg.Data());
    return -2;
  }

  int    is, nm, nl, nbars, first_bar, iwx, iwy, iwz;
  float  dx, dy, dz;
  char   name[100];
  int    done(0);
  double dxx[3];
  
  while ( ((c[0]=getc(f)) != EOF) && !done) {
					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
					// parse line
      fscanf(f,"%i" ,&is);
      fscanf(f,"%s" ,name      );
      fscanf(f,"%i" ,&nm       );
      fscanf(f,"%i" ,&nl       );
      fscanf(f,"%i" ,&nbars    );
      fscanf(f,"%i" ,&first_bar);
      fscanf(f,"%i" ,&iwx      );
      fscanf(f,"%i" ,&iwy      );
      fscanf(f,"%i" ,&iwz      );
      fscanf(f,"%f" ,&dy       );
      fscanf(f,"%f" ,&dx       );
      fscanf(f,"%f" ,&dz       );
					// translate from mm to cm - do we really need to do that ?
      dxx[0]    = dx/10;
      dxx[1]    = dy/10;
      dxx[2]    = dz/10;

      TEvdCrvSection* crvs = new TEvdCrvSection(is,name);

      crvs->fNModules = nm;
      crvs->fNLayers  = nl;
      crvs->fNBars    = nbars;
      crvs->fFirstBar = first_bar;

      printf("is: %2i dx,dy,dz: %10.3f %10.3f %10.3f\n",is,dxx[0],dxx[1],dxx[2]);
      
      crvs->fBarShape = new TGeoBBox(Form("bar_shape_%02i",is), dxx[0], dxx[1], dxx[2]);

      fGeoManager->fCrvSection[is] = crvs;
    }
					// skip line
    fgets(c,1000,f);
  }

  fclose(f);
//-----------------------------------------------------------------------------
// now read the counters file
//-----------------------------------------------------------------------------
  const char* fn2 = "crv_counter_geom.txt";
  
  f  = fopen(fn2,"r");
  if (f == 0) {
    Error("Init",Form("missing file %s",fn2));
    return -2;
  }

  int    bar_index, im, il, ib;
  float  x0, y0, z0;
  double origin[3];

  while ( ((c[0]=getc(f)) != EOF) && !done) {
					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
					// parse line
      fscanf(f,"%i" ,&bar_index);
      fscanf(f,"%i" ,&is       );
      fscanf(f,"%i" ,&im       );
      fscanf(f,"%i" ,&il       );
      fscanf(f,"%i" ,&ib       );
      fscanf(f,"%i" ,&iwy      );
      fscanf(f,"%i" ,&iwx      );
      fscanf(f,"%i" ,&iwz      );
      fscanf(f,"%f" ,&x0       );
      fscanf(f,"%f" ,&y0       );
      fscanf(f,"%f" ,&z0       );
      fscanf(f,"%f" ,&dy       );        // thickness
      fscanf(f,"%f" ,&dx       );        // width
      fscanf(f,"%f" ,&dz       );
					// translate from mm to cm - do we really need to do that ?
      origin[0] = x0/10;
      origin[1] = y0/10;
      origin[2] = z0/10;

      TGeoBBox* bar_shape = fGeoManager->CrvSection(is)->BarShape();

      TGeoVolume*  bar    = new TGeoVolume (Form("%s_bar_%05i",fGeoManager->CrvSection(is)->GetName(),bar_index),bar_shape);
      bar->SetMedium(fMedia.fAl);
      bar->SetLineColor(kGreen-2);
      bar->SetTransparency(80.);
      
      TString rname;
      
      double theta1(0), phi1(0), theta2(0), phi2(0), theta3(0), phi3(0);

      if ((iwy == 1) and (iwx == 2) and (iwz == 0)) {
	rname = "rot_00";
	theta1 =   0;
	phi1   =   0;
	theta2 =  90;
	phi2   =  90;
	theta3 =  90;
	phi3   =   0;
      }
      else if ((iwy == 0) and (iwx == 2) and (iwz == 1)) {
	rname = "rot_01";
	theta1 =   0;
	phi1   =   0;
	theta2 =  90;
	phi2   =   0;
	theta3 =  90;
	phi3   =  90;
      }
      else if ((iwy == 1) and (iwx == 0) and (iwz == 2)) {
	rname = "rot_02";
	theta1 =  90;
	phi1   =   0;
	theta2 =  90;
	phi2   =  90;
	theta3 =   0;
	phi3   =   0;
      }
      else if ((iwy == 2) and (iwx == 1) and (iwz == 0)) {
	rname = "rot_03";
	theta1 =  90;
	phi1   =  90;
	theta2 =   0;
	phi2   =   0;
	theta3 =  90;
	phi3   =   0;
      }
      else if ((iwy == 1) and (iwx == 0) and (iwz == 2)) {
	rname = "rot_03";
	theta1 =  90;
	phi1   =   0;
	theta2 =  90;
	phi2   =  90;
	theta3 =   0;
	phi3   =   0;
      }

      TGeoRotation    rot0(rname.Data(),theta1,phi1,theta2,phi2,theta3,phi3);
      TGeoTranslation tr0 (origin[0], origin[1], origin[2]);
      
      TGeoCombiTrans* com0 = new TGeoCombiTrans(tr0, rot0);

      fMu2eGeom->AddNode(bar,bar_index,com0);
    }
					// skip line
    fgets(c,1000,f);
  }

  fclose(f);

  return 0;
}

//-----------------------------------------------------------------------------
void View3D::Draw(Option_t *option) {
   TObject::Draw(option);
					// Ask pad to create 3D viewer of type 'option'
   gPad->GetViewer3D(option);
}

//-----------------------------------------------------------------------------
void View3D::Paint(Option_t* Option /*option*/) {
  //   TVirtualViewer3D * viewer = gPad->GetViewer3D();

   // If View3D derives from TAtt3D then pad will recognise
   // that the object it is asking to paint is 3D, and open/close
   // the scene for us. If not Open/Close are required
   //viewer->BeginScene();

   // We are working in the master frame - so we don't bother
   // to ask the viewer if it prefers local. Viewer's must
   // always support master frame as minimum. c.f. with
   // viewer3DLocal.C

   fMu2eGeom->Paint();

   fTrack->Paint();

   fListOfTrajectories->Paint();

   // Not required as we are TAtt3D subclass
   // viewer->EndScene();
}

//-----------------------------------------------------------------------------
int View3D::ReadTrajectories(const char* DirName) {

  fListOfTrajectories = new TObjArray();
  fListOfTrajectories->SetOwner(kTRUE);

  char buf[1000], c[1000];

  int index;
  float x[10000], y[10000], z[10000], t[10000], ekin[10000];

  TString cmd = Form("ls %s/*.txt",DirName);

  TObjArray list_of_files;

  FILE* pipe = gSystem->OpenPipe(cmd.Data(),"r");

  while (fgets(buf,1000,pipe)) {
    char fn[1000];
    sscanf(buf,"%s",fn);
    TObjString* s = new TObjString(fn);
    list_of_files.Add(s);
  }
  
  gSystem->ClosePipe(pipe);

  int nfiles = list_of_files.GetEntriesFast();
  for (int i=0; i<nfiles; i++) {
    TObjString* os = (TObjString*) list_of_files.At(i);
    const char* fn = os->String().Strip(TString::kTrailing,'\n').Data();
    FILE* f = fopen(fn,"r");
    // printf("--------------------------- opened %s.\n",fn);
    int npt = 0;
    while ((c[0]=getc(f)) != EOF)  {
      // check if it is a comment line
      // printf("read: %s",c);
      if (c[0] != '#') {
	ungetc(c[0],f);
	// parse line

	float xx, yy, zz;
	
	fscanf(f,"%i" ,&index);
	fscanf(f,"%f" ,&xx  );
	fscanf(f,"%f" ,&yy );
	fscanf(f,"%f" ,&zz  );        // thickness
	fscanf(f,"%f" ,&t[npt]  );        // width
	fscanf(f,"%f" ,&ekin[npt]);

	x[npt] = xx/10.;
	y[npt] = yy/10.;
	z[npt] = zz/10.;
	
	npt++;
      }
					// skip end of line
      fgets(c,1000,f);
    }

    fclose(f);

    TString sfn = fn;
    TObjArray* afn = sfn.Tokenize("/");
    int nel = afn->GetEntries();
    TObjString* ss = (TObjString*) (*afn)[nel-1];

    TEvdTrajectory* trj = new TEvdTrajectory(ss->String().Data(),npt,x,y,z);
    
    trj->SetLineColor(kRed+1);

    fListOfTrajectories->Add(trj);
  }

  return 0;
}

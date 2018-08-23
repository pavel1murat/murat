//
#include "murat/ana/TTsMisalignment.hh"

//-----------------------------------------------------------------------------
//  measurements in the local, for each TS piece, coordinate system
//  TSUL : TS1  collimator
//  TSUR : TS3u collimator
//  TSDL : TS5  collimator
//  TSDR : TS3d collimator
//-----------------------------------------------------------------------------
namespace {

  struct Data_t {
    TVector3 tsul [2];			// nominal 
    TVector3 tsulm[2];			// measured
    TVector3 tsur [2];			// nominal 
    TVector3 tsurm[2];			// measured

    TVector3 tsdl [2];			// nominal 
    TVector3 tsdlm[2];			// measured
    TVector3 tsdr [2];			// nominal 
    TVector3 tsdrm[2];			// measured
  };

//  positions of the outer and inner points on the cylinder axis
//
  Data_t data = {
//            (x,y,z) outer                   (x,y,z) inner
    {{  0.    ,   0.    ,  0.    }, {  0.0000,  65.5034,  0.0000}},  // TSUL nominal
    {{  0.4997,  -0.2817, -0.3827}, {  0.273 ,  65.5034, -0.2076}},  // TSUL measured

    {{405.0   , 417.5   ,  0.    }, {296.9065, 417.5   ,  0.0000}},  // nominal
    {{403.0654, 417.0060,  0.3523}, {296.9065, 417.3526,  0.1170}},  // measured

    {{  0.    ,   0.    ,  0.0000}, {  0.0000,  54.7010,  0.0000}},  // nominal
    {{ -0.0267,  -1.1981, -0.5653}, { -0.1483,  54.7010, -0.2530}},  // measured
    
    {{405.0000, 435.0000,  0.0000}, {324.3169, 435.0000,  0.0000}},  // nominal
    {{403.5686, 434.9822,  0.0944}, {324.3169, 435.0884,  0.1290}}   // measured
  };

  Data_t ds;   // data with scaled displacements
  
};

//-----------------------------------------------------------------------------
TTsMisalignment::TTsMisalignment() {
  fR    = 15;  // cm, this is just an example

  fTsuL     = NULL;
  fTsuR     = NULL;
  fTsuLMeas = NULL;
  fTsuRMeas = NULL;

  fTsdL     = NULL;
  fTsdR     = NULL;
  fTsdLMeas = NULL;
  fTsdRMeas = NULL;

  fScale = 1.;
}

//-----------------------------------------------------------------------------
TTsMisalignment::~TTsMisalignment() {
  if (fTsuL) {
    delete fTsuL;
    delete fTsuLMeas;
    delete fTsuR;
    delete fTsuRMeas;

    delete fTsdL;
    delete fTsdLMeas;
    delete fTsdR;
    delete fTsdRMeas;
  }
}

//-----------------------------------------------------------------------------
// so far, no rotation in Z
//-----------------------------------------------------------------------------
void TTsMisalignment::InitPolyLine(TVector3* Pos, TPolyLine*& Pl) {

  double dx,dy,dz, len, x[5], y[5];
					// Pos[1] - inner end, Pos[0] - outer end
  dx = Pos[1].X()-Pos[0].X();
  dy = Pos[1].Y()-Pos[0].Y();
  dz = Pos[1].Z()-Pos[0].Z();

  len  = sqrt(dx*dx+dy*dy+dz*dz); // length of the cylinder

  double nx1, ny1, nx2, ny2;
  
  nx1 =  dx/len;
  ny1 =  dy/len;

  nx2 =  ny1;
  ny2 = -nx1;
					// polyline points
  x[0] = Pos[0].X()+nx2*fR;
  y[0] = Pos[0].Y()+ny2*fR;

  x[1] = x[0]+len*nx1;
  y[1] = y[0]+len*ny1;
  
  x[2] = x[1]-2*fR*nx2;
  y[2] = y[1]-2*fR*ny2;
  
  x[3] = x[2]-len*nx1;
  y[3] = y[2]-len*ny1;
  
  x[4] = x[0];
  y[4] = y[0];

  if (Pl) delete Pl;

  Pl = new TPolyLine(5,x,y);
 }


//-----------------------------------------------------------------------------
void TTsMisalignment::Transform(double X0, double Y0, double Phi,
				TVector3* X, TVector3* Xr) {

  // rotate 2 points X= {x1,y1,x2,y2}

  // double nx = cos(Phi);
  // double ny = sin(Phi);
					// i=0: outer end, i=1: inner end of the measured cylinder
  for (int i=0; i<2; i++) {
    double dx = X[i].X()-X0;
    double dy = X[i].Y()-Y0;

    Xr[i].SetX(dx*cos(Phi)-dy*sin(Phi));
    Xr[i].SetY(dx*sin(Phi)+dy*cos(Phi));

    printf("i:%3i dx = %12.5f dy = %12.5f\n",i,dx,dy);
    printf("i:%3i X,Xr: %10.4f %10.4f %10.4f %10.4f \n",i,X[i].X(),X[i].Y(),Xr[i].X(),Xr[i].Y());
  }
}


//-----------------------------------------------------------------------------
void TTsMisalignment::Init(double Scale) {
//-----------------------------------------------------------------------------
// scale real displacements
//-----------------------------------------------------------------------------
  fScale = Scale;

  for (int i=0; i<2; i++) {
    ds.tsul [i] = data.tsul[i];
    ds.tsulm[i] = data.tsul[i] + (data.tsulm[i]-data.tsul[i])*fScale;
    ds.tsur [i] = data.tsur[i];
    ds.tsurm[i] = data.tsur[i] + (data.tsurm[i]-data.tsur[i])*fScale;

    ds.tsdl [i] = data.tsdl[i];
    ds.tsdlm[i] = data.tsdl[i] + (data.tsdlm[i]-data.tsdl[i])*fScale;
    ds.tsdr [i] = data.tsdr[i];
    ds.tsdrm[i] = data.tsdr[i] + (data.tsdrm[i]-data.tsdr[i])*fScale;
  }
//-----------------------------------------------------------------------------
// TSU transformation: tsur nominal at (0,0); rotation by -90 deg
//-----------------------------------------------------------------------------
					// nominal position of the outer TSUR cylinder end
  double x0 = ds.tsur[0].X();
  double y0 = ds.tsur[0].Y();
  double z0 = ds.tsur[0].Z();

  double phi = -TMath::Pi()/2;
					// measured position of the outer TSUR cylinder end
  double x0m  = ds.tsurm[0].X();
  double y0m  = ds.tsurm[0].Y();
  double z0m  = ds.tsurm[0].Z();
//-----------------------------------------------------------------------------
// rotate bore pieces as measured to align TS3 end perfectly
//-----------------------------------------------------------------------------
  TVector3 tsurm_dr = ds.tsurm[1] - ds.tsurm[0];
  
  double dx = tsurm_dr.X();
  double dy = tsurm_dr.Y();
  double dz = tsurm_dr.Z();

  double phiy = -atan(dy/dx)-TMath::Pi()/2;
  double phiz = -atan(dz/dx);

  printf("dx =  %12.5f dy = %12.5f dz = %12.5f z0 = %12.5f phi = %12.5f\n",dx,dy,dz,z0,phi);
  
  printf("x0m = %12.5f y0m = %12.5f z0m = %12.5f phiy = %12.5f phiz = %12.5f \n", x0m, y0m, z0m, phiy, phiz);

  //  double x[5], y[5], xr[6];

  TVector3 xr[2];
					// TSuL nominal
  printf("TSuL nominal\n");
  Transform(x0,y0,phi,ds.tsul,xr);
  InitPolyLine(xr,fTsuL);
					// TSuR nominal
  printf("TSuR nominal\n");
  Transform(x0,y0,phi,ds.tsur,xr);
  InitPolyLine(xr,fTsuR);

  // uncomment to see measurements
  // x0m  = x0;
  // y0m  = y0;
  // phim = phi;
					// TSuL measured
  printf("TSuL measured\n");
  Transform(x0m,y0m,phiy,ds.tsulm,xr);
  InitPolyLine(xr,fTsuLMeas);
  fTsuLMeas->SetLineColor(2);
					// TSuR measured
  printf("TSuR measured\n");
  Transform(x0m,y0m,phiy,ds.tsurm,xr);
  InitPolyLine(xr,fTsuRMeas);
  fTsuRMeas->SetLineColor(2);
//-----------------------------------------------------------------------------
// TSd transformation: tsdr nominal at (0,0); rotation by 90 deg
//-----------------------------------------------------------------------------
  double x0d  = ds.tsdr[0].X();
  double y0d  = ds.tsdr[0].Y();
  double z0d  = ds.tsdr[0].Z();
  
  double phid = TMath::Pi()/2;

  printf("x0d = %12.5f y0d = %12.5f z0d = %12.5f phid = %12.5f\n",
	 x0d, y0d, z0d, phid);
					// TSdL nominal
  printf("TSdL nominal\n");
  Transform(x0d,y0d,phid,ds.tsdl,xr);
  InitPolyLine(xr,fTsdL);
					// TSdR nominal
  printf("TSdR nominal\n");
  Transform(x0d,y0d,phid,ds.tsdr,xr);
  InitPolyLine(xr,fTsdR);
					// measured position of the TSd right cylinder
  double x0d_m  = ds.tsdrm[0].X();
  double y0d_m  = ds.tsdrm[0].Y();
  double z0d_m  = ds.tsdrm[0].Z();

  TVector3 drd  = ds.tsdrm[1]-ds.tsdrm[0];
  
  double phiyd_m = -atan(drd.Y()/drd.X())+TMath::Pi()/2;
  double phizd_m = -atan(drd.Z()/drd.X());

  printf("x0d_m = %12.5f y0d_m = %12.5f z0d_m = %12.5f dxd = %12.5f dyd = %12.5f dzd = %12.5f phiyd_m = %12.5f phizd_m = %12.5f\n",
	 x0d_m, y0d_m, z0d_m,
	 drd.X(), drd.Y(), drd.Z(),
	 phiyd_m, phizd_m);

  // uncomment to see the measurements
  // x0d_m   = x0d;
  // y0d_m   = y0d;
  // phid_m = phid;
					// TSdL measured
  printf("TSdL measured\n");
  Transform(x0d_m,y0d_m,phiyd_m,ds.tsdlm,xr);
  InitPolyLine(xr,fTsdLMeas);
  fTsdLMeas->SetLineColor(2);
					// TSdR measured
  printf("TSdR measured\n");
  Transform(x0d_m,y0d_m,phiyd_m,ds.tsdrm,xr);
  InitPolyLine(xr,fTsdRMeas);
  fTsdRMeas->SetLineColor(2);

}

//-----------------------------------------------------------------------------
void TTsMisalignment::Paint(Option_t* Opt) {
  fTsuL->Paint();
  fTsuR->Paint();
  fTsuLMeas->Paint();
  fTsuRMeas->Paint();

  fTsdL->Paint();
  fTsdR->Paint();
  fTsdLMeas->Paint();
  fTsdRMeas->Paint();
}

//-----------------------------------------------------------------------------
// there is a class named Disk in Mu2e offline
//-----------------------------------------------------------------------------
#include "tdr/Disk.hh"


Disk* Disk::gDisk(0);

//-----------------------------------------------------------------------------
Disk::Disk() {
  //  fListOfPolyLines = 0;
}


//-----------------------------------------------------------------------------
Disk::Disk(double RMin, double RMax, double HexSize, double MinFraction) {

  fRMin        = RMin;
  fRMax        = RMax;
  fHexSize     = HexSize;
  fMinFraction = MinFraction;

  fNInsideTot  = 0;

  fPos[0]      = HexIndex( 1, 0);  // up right
  fPos[1]      = HexIndex( 1,-1);  // up 
  fPos[2]      = HexIndex( 0,-1);  // up left
  fPos[3]      = HexIndex(-1, 0);  // down left
  fPos[4]      = HexIndex(-1, 1);  // down
  fPos[5]      = HexIndex( 0, 1);  // down right

  fNRings      = (int(fRMax/fHexSize-0.2)+1)*1.2+2;
  fNCrystals   = 1 + 3*fNRings*(fNRings-1);

  fFirst            = new int[fNRings];
  fNCrystalsPerRing = new int[fNRings];
  fNInside          = new int[fNRings];

  fFirst           [0] = 0;
  fNCrystalsPerRing[0] = 1;
  fNInside         [0] = 0;

  for (int i=1; i<fNRings; i++) {
    fFirst           [i] = 1 + 3*i*(i-1);
    fNCrystalsPerRing[i] = 6*i;
    fNInside         [i] = 0;
  }

  int       loc, ring, inside, first; 
  double    fraction;
  HexIndex  hex_index;
  TVector2  pos;

  for (int ir=0; ir<fNRings; ir++) {
    first = fFirst[ir];
    ring  = -1;
    for (int ic=0; ic<fNCrystalsPerRing[ir]; ic++) {
      loc       = first+ic;
      hex_index = GetHexIndex(loc);

      ring      = GetRing  (&hex_index);
      GetPosition(&hex_index,&pos);

      inside    = IsInside(&hex_index, &fraction);

      if (inside) {
	fNInside[ring] += 1;
      }

      //      printf("loc=%3i ring=%2i x,y=%10.5f %10.5f y=inside=%2i fraction=%7.3f minfr=%7.3f\n",
      //	     loc,ring,pos.X(),pos.Y(),inside,fraction,fMinFraction);
    }

    fNInsideTot += fNInside[ring];
  }

  //  printf(">>> constructor: fMinFractin = %7.3f\n",fMinFraction);
}

//-----------------------------------------------------------------------------
Disk::~Disk() {
  delete [] fFirst;
  delete [] fNCrystalsPerRing;
  delete [] fNInside;

  // if (fListOfPolyLines) {
  //   fListOfPolyLines->Delete();
  //   delete fListOfPolyLines;
  // }
}


//-----------------------------------------------------------------------------
HexIndex Disk::GetHexIndex(int I) {

  HexIndex hex_index(0,0);
  int      ring, loc, l, k, n1, n2;

  if (I > 0) {
				// find the ring number
    for (int i=1; i>0; i++) {
      n1 = 3*(i-1)*i+1;
      n2 = 3*i*(i+1)+1;
      if (I < n2) {
	ring = i;
	loc  = I-n1;
	break;
      }
    }

    int seg1 = loc / ring; 
    int pos  = loc % ring;

    int seg2 = seg1+1;
    if (seg2 == 6) seg2 = 0;

    l = fPos[seg1].fL*ring + (fPos[seg2].fL-fPos[seg1].fL)*pos;
    k = fPos[seg1].fK*ring + (fPos[seg2].fK-fPos[seg1].fK)*pos;
  
    hex_index.Set(l,k);
  }
	 	 	 
  return hex_index;
}

//-----------------------------------------------------------------------------
int Disk::GetRing(HexIndex* Index) {

  int  ring;

  if (Index->fL*Index->fK > 0) {
    ring = std::abs(Index->fL+Index->fK);
  }
  else if ( std::abs(Index->fL) > std::abs(Index->fK) ) {
    ring = std::abs(Index->fL);
  }
  else {
    ring = std::abs(Index->fK);
  }

  return ring;
}

//-----------------------------------------------------------------------------
int Disk::GetRing(int I) {
  HexIndex  hex_index;
  
  hex_index = GetHexIndex(I);

  return GetRing(&hex_index);
}

//-----------------------------------------------------------------------------
void Disk::GetPosition(HexIndex* Index, TVector2* Pos) {
  double x, y;

  x = fHexSize*(Index->fL+Index->fK)*sqrt(3.)/2.;
  y = fHexSize*(Index->fL-Index->fK)/2.;
  Pos->Set(x,y);

}

//-----------------------------------------------------------------------------
void Disk::GetPosition(int I, TVector2* Pos) {
  HexIndex  hex_index;
  
  hex_index = GetHexIndex(I);

  GetPosition(&hex_index,Pos);
}

//-----------------------------------------------------------------------------
int Disk::IsInside(HexIndex* Index, double* Fraction) {

  int       inside, nvin(0), nbelow(0), nabove(0);
  double    x0, y0, r, x, y, phi, s, s0, /*s1,*/ r0, dr, adr;
  double    fr;

  x0 = fHexSize*(Index->fL+Index->fK)*sqrt(3.)/2.;
  y0 = fHexSize*(Index->fL-Index->fK)/2.;

  // loop over 6 vertices and check if all of them are inside

  double rho = fHexSize/sqrt(3.);

  for (int i=0; i<6; i++) {
    phi = i*TMath::Pi()/3;
    x = x0+rho*TMath::Cos(phi);
    y = y0+rho*TMath::Sin(phi);

    r = sqrt(x*x+y*y);

    if (r < fRMin) {
      nbelow += 1;
    }
    else if (r > fRMax) {
      nabove += 1;
    }
    else {
      nvin += 1;
    }
  }

  if (nvin == 0) {
    fr = 0;
  }
  else if (nvin == 6) {
    fr = 1.;
  }
  else { 
				// at this point can calculate fraction of the crystal area inside the ring
    r0  = sqrt(x0*x0+y0*y0);
    dr  = r0-fRMin;

    if (nabove > 0) {
//-----------------------------------------------------------------------------
// crystal crosses the outer ring
//-----------------------------------------------------------------------------
      dr  = r0 -fRMax;
      adr = fabs(dr);
      if (adr > fHexSize/2) adr = fHexSize/2.;

      s   = fHexSize*fHexSize*sqrt(3)/2;
      //      s1  = (2*fHexSize-adr)*adr/sqrt(3.);
      s0  = (3*fHexSize-2*adr)*(fHexSize-2*adr)/4/sqrt(3);

      if (dr <= 0) {
	fr = 1-s0/s;
      }
      else {
        fr = s0/s;
      }
    }
    else {
//-----------------------------------------------------------------------------
// crystal crosses the inner ring
//-----------------------------------------------------------------------------
      adr = fabs(dr);
      if (adr > fHexSize/2) adr = fHexSize/2.;

      s   = fHexSize*fHexSize*sqrt(3)/2;
      //      s1  = (2*fHexSize-adr)*adr/sqrt(3.);
      s0  = (3*fHexSize-2*adr)*(fHexSize-2*adr)/4/sqrt(3);
      if (dr > 0) {
	fr = 1-s0/s;
      }
      else {
        fr = s0/s;
      }
    }
  }

  if (fr >= fMinFraction) inside = 1;
  else                    inside = 0;

  *Fraction = fr;

  return inside;
}

//-----------------------------------------------------------------------------
void Disk::Test1(double RMin, double RMax, double HexSize, double Fraction, const char* Opt) {

  Disk* disk;

  double rmin, rmax;

  disk = new Disk(RMin, RMax, HexSize, Fraction);

  rmin = disk->GetRMin();
  rmax = disk->GetRMax();

  double s0 = TMath::Pi()*(rmax*rmax-rmin*rmin);

  int nin = disk->GetNInsideTot();

  double area = disk->GetCrystalArea();

  double s1 = nin*area;

  double r  = s1/s0;

  int ntotal = disk->GetNCrystals();

  // for (int i=0;  i<disk->GetNRings(); i++) {
  //   printf(" i, NTot NIn : %5i %5i %5i\n",i,disk->GetNCrystalsPerRing(i),disk->GetNInside(i));
  // }

  printf("rin,rout,hsize,nr = %7.2f %7.2f %6.2f %3i ",RMin,RMax,HexSize,disk->GetNRings());

  printf(" ntot,nin,area: %4i %4i %8.5f",ntotal, nin, area);

  printf(" s0,s1,r: %10.4f %10.4f %10.4f \n",s0, s1, r);

  if (index(Opt,'g') != NULL) {
    disk->Draw("");
  }
}

//-----------------------------------------------------------------------------
void Disk::OptimizeGeometry(double RMin, double RMax, double HexCrystalSize) {
  
  
  //  Disk* disk;

  gDisk    = new Disk(RMin,RMax,HexCrystalSize);

  //  int ncr = gDisk->GetNCrystals();
  //  int nr  = gDisk->GetNRings();

  // look at the first and the last rings



}


//-----------------------------------------------------------------------------
void Disk::Draw(Option_t* Opt) {

  TVector2 pos;

  HexIndex  hex_index;

  TPolyLine p;
  TEllipse  e;

  double x[7], y[7];

  TCanvas* c = new TCanvas("c","c",800,800);

  double sf = 1.3;
  TH2F* h2 = new TH2F("h2","h2",100,-sf*fRMax,sf*fRMax,100,-sf*fRMax,sf*fRMax);

  h2->SetStats(0);
  h2->Draw();

  // fListOfPolyLines->Delete();

  for (int i=0; i<fNCrystals; i++) {

    hex_index = GetHexIndex(i);

    GetPosition(i,&pos);

    //    int ring = GetRing(i);

    //    printf(" i, r, l,k, x,y : %5i %5i %5i %5i %10.4lf %10.4lf\n",i,ring, hex_index.fL,hex_index.fK,pos.X(),pos.Y());

    double rho = fHexSize/sqrt(3.);

    for (int i=0; i<7; i++) {
      int l = i % 6;
      x[i] = pos.X() + rho*cos(TMath::Pi()*l/3);
      y[i] = pos.Y() + rho*sin(TMath::Pi()*l/3);
    }

    p.SetLineColor(1);
    p.SetLineWidth(1);

    double fr;

    if (IsInside(&hex_index,&fr)) {

      if (fr == 1) {
	p.SetFillColor(2);                     // fully inside
	p.SetFillStyle(3001);
      }
      else {
	p.SetFillColor(kBlue);             // partially inside
	p.SetFillStyle(1001);
      }

      //      p.SetFillStyle(3001);
      p.DrawPolyLine(7,x,y,"F");
    }
    else {
				// fully outside
      p.SetFillStyle(0);
      p.SetFillColor(0);
    }

    p.DrawPolyLine(7,x,y,"");
  } 

  e.SetFillStyle(0);
  e.SetFillColor(0);
  e.DrawEllipse(0,0,fRMin,fRMin,0,360,0,"");
  e.DrawEllipse(0,0,fRMax,fRMax,0,360,0,"");
  

  c->Update();
  c->Draw();
}

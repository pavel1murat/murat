///////////////////////////////////////////////////////////////////////////////
// use:
// mumem* ana = new mumem(1);
// ana->
///////////////////////////////////////////////////////////////////////////////
const char* FiguresDir         = "/projects/mu2e/talks/2021-03-11-su2020-backgrounds-collab/figures";

#include "murat/limits/mumem.hh"

namespace murat {
//-----------------------------------------------------------------------------
int mumem::scan_pmin(double PMin, double PMax, double TMin, double TMax) {
  //  init_channels();

  TFeldmanCousinsA fc("mumem::Print",-1);
  fc.SetNExp(100000); // enough for defining 5

  double smin(1.), smax(11.);  // range of tested signal values
  int  npt(20);	       // number of points
  double s[npt], msign[npt];

  double qbgr[100];

  int nbgr = GetNBgrChannels();

  double pmax = PMax;
  
  printf("      pmin       pmax       tmin       tmax    %-10s     %-10s %-10s      %-10s Total        Signal         s5          SES         Rmue\n",
	 GetBgrChannel(0)->GetName(),
	 GetBgrChannel(1)->GetName(),
	 GetBgrChannel(2)->GetName(),
	 GetBgrChannel(3)->GetName()
	 );
    
  printf("--------------------------------------------------------------------------------");
  printf("-------------------------------------------------------------------------------\n");

  for (int i=0; i<20; i++) {
    double pmin = PMin-i*0.050;
    double sig = GetSigIntegral(pmin,pmax,TMin,TMax);
    double ses = GetSES        (pmin,pmax,TMin,TMax);

    double bgr_tot = 0;
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      su2020::channel* bgr = GetBgrChannel(ibgr);
      qbgr[ibgr] = bgr->GetIntegral(pmin,pmax,TMin,TMax);
      bgr_tot   += qbgr[ibgr] ;
    }
//-----------------------------------------------------------------------------
// for given bgr_tot, calculate the expected mean 5-sigma sensitivity
// assume bgr_tot < 0.5, [1,11] covers the range in all cases
//-----------------------------------------------------------------------------
    fc.DiscoveryProbMean(bgr_tot, smin, smax, npt, s, msign);

    double sigd(5.), s5;                    // signal coresponding to 1-sided "5-sigma"
    fc.SolveFor(sigd,s,msign,npt,&s5);
    printf("%10.3f %10.3f %10.3f %10.3f",pmin,PMax,TMin,TMax);
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      printf(" %12.5e", qbgr[ibgr]);
    }
    printf(" %12.5e %12.5e %12.5e %12.5e %12.5e\n",bgr_tot,sig,s5,ses,s5*ses);
  }
  return 0;
}


//-----------------------------------------------------------------------------
int mumem::scan_tmin(double PMin, double PMax, double TMin,double TMax) {
  //  init_channels();

  TFeldmanCousinsA fc("mumem::Print",-1);
  fc.SetNExp(100000); // enough for defining 5

  double smin(1.), smax(11.);  // range of tested signal values
  int  npt(20);	       // number of points
  double s[npt], msign[npt];

  double qbgr[100];

  int nbgr = GetNBgrChannels();

  printf("      pmin       pmax       tmin       tmax    %-10s     %-10s %-10s      %-10s Total        Signal         s5          SES         Rmue \n",
	 GetBgrChannel(0)->GetName(),
	 GetBgrChannel(1)->GetName(),
	 GetBgrChannel(2)->GetName(),
	 GetBgrChannel(3)->GetName()
	 );
  
  printf("---------------------------------------------------------------------------------");
  printf("-------------------------------------------------------------------------------\n");

  for (int i=0; i<20; i++) {
    double tmin = TMin-i*10;		// the bin width is 10 ns
    double sig = GetSigIntegral(PMin,PMax,tmin,TMax);
    double ses = GetSES        (PMin,PMax,tmin,TMax);
   
    double bgr_tot = 0;
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      su2020::channel* bgr = GetBgrChannel(ibgr);
      qbgr[ibgr] = bgr->GetIntegral(PMin,PMax,tmin,TMax);
      bgr_tot   += qbgr[ibgr] ;
    }
//-----------------------------------------------------------------------------
// for given bgr_tot, calculate the expected mean 5-sigma sensitivity
// assume bgr_tot < 0.5, [1,11] covers the range in all cases
//-----------------------------------------------------------------------------
    fc.DiscoveryProbMean(bgr_tot, smin, smax, npt, s, msign);

    double s5;                         // signal coresponding to 1-sided "5-sigma"
    fc.SolveFor(5,s,msign,npt,&s5);
    printf("%10.3f %10.3f %10.3f %10.3f",PMin,PMax,tmin,TMax);
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      printf(" %12.5e", qbgr[ibgr]);
    }
    printf(" %12.5e %12.5e %12.5e %12.5e %12.5e\n",bgr_tot,sig,s5,ses,s5*ses);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int mumem::scan_tmax(double PMin, double PMax, double TMin,double TMax) {
  //  init_channels();

  TFeldmanCousinsA fc("mumem::Print",-1);
  fc.SetNExp(100000); // enough for defining 5

  double smin(1.), smax(11.);  // range of tested signal values
  int  npt(20);	       // number of points
  double s[npt], msign[npt];

  double qbgr[100];

  int nbgr = GetNBgrChannels();

  printf("      pmin       pmax       tmin       tmax    %-10s     %-10s %-10s      %-10s Total        Signal         s5          SES         Rmue \n",
	 GetBgrChannel(0)->GetName(),
	 GetBgrChannel(1)->GetName(),
	 GetBgrChannel(2)->GetName(),
	 GetBgrChannel(3)->GetName()
	 );
  
  printf("---------------------------------------------------------------------------------");
  printf("-------------------------------------------------------------------------------\n");

  for (int i=0; i<20; i++) {
    double tmax = TMax-i*10;		// the bin width is 10 ns
    double sig = GetSigIntegral(PMin,PMax,TMin,tmax);
    double ses = GetSES        (PMin,PMax,TMin,tmax);
   
    double bgr_tot = 0;
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      su2020::channel* bgr = GetBgrChannel(ibgr);
      qbgr[ibgr] = bgr->GetIntegral(PMin,PMax,TMin,tmax);
      bgr_tot   += qbgr[ibgr] ;
    }
//-----------------------------------------------------------------------------
// for given bgr_tot, calculate the expected mean 5-sigma sensitivity
// assume bgr_tot < 0.5, [1,11] covers the range in all cases
//-----------------------------------------------------------------------------
    fc.DiscoveryProbMean(bgr_tot, smin, smax, npt, s, msign);

    double s5;                         // signal coresponding to 1-sided "5-sigma"
    fc.SolveFor(5,s,msign,npt,&s5);
    printf("%10.3f %10.3f %10.3f %10.3f",PMin,PMax,TMin,tmax);
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      printf(" %12.5e", qbgr[ibgr]);
    }
    printf(" %12.5e %12.5e %12.5e %12.5e %12.5e\n",bgr_tot,sig,s5,ses,s5*ses);
  }
  return 0;
}

}

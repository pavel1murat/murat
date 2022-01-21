//-----------------------------------------------------------------------------
int test_syst_err() {

  TFeldmanCousinsA fc("mumem::Print",-1);
  fc.SetNExp(10000); // enough for defining 5

  double smin(1.), smax(11.);  // range of tested signal values
  double s[npt], msign[npt];

  double qbgr[100];

  double tmin = TMin-i*10;		// the bin width is 10 ns
  double sig = GetSigIntegral(PMin,PMax,tmin,TMax);
  double ses = GetSES        (PMin,PMax,tmin,TMax);
//-----------------------------------------------------------------------------
// for given bgr_tot, calculate the expected mean 5-sigma sensitivity
// assume bgr_tot < 0.5, [1,11] covers the range in all cases
//-----------------------------------------------------------------------------
  TBgrChannel* ch = new TChannel("emoe",0.1,0.01);
  
  fc.DiscoveryProbMean(ch, smin, smax, npt, s, msign);

  double s5;                         // signal coresponding to 1-sided "5-sigma"
  fc.SolveFor(5,s,msign,npt,&s5);
  printf("%10.3f %10.3f %10.3f %10.3f",PMin,PMax,tmin,TMax);
  for (int ibgr=0; ibgr<nbgr; ibgr++) {
    printf(" %12.5e", qbgr[ibgr]);
  }
  printf(" %12.5e %12.5e %12.5e %12.5e %12.5e\n",bgr_tot,sig,s5,ses,s5*ses);

  return 0;
}

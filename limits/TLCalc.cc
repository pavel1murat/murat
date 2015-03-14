///////////////////////////////////////////////////////////////////////////////
// mode: 1:PYTHIA, 2:HERWIG
///////////////////////////////////////////////////////////////////////////////
#include "limits/TLCalc.hh"
#include "TMath.h"
//-----------------------------------------------------------------------------
TLCalc::TLCalc(int Mode): TNamed() {
  //                M(ZZ)  PT(ZZ)
  double x[8][2] = { {196.,  35.},
		     {190.,  30.},
		     {234.,  10.},
		     {192.,  27.},
					// high mass
		     {321.,  47.},
		     {325., 127.},
		     {329., 111.},
		     {334.,  45.}
  };

  for (int i=0; i<8; i++) {
    fData[i][0] = x[i][0];
    fData[i][1] = x[i][1];
  }

  fMode       = Mode;
  fNMax       = 100000000;
  fMinMass    = 300.;
  fMassWindow = 20.;
  fNMean      = 5.5;    // 8*1.4/2.1*1.03 = 5.5

  InitHistograms();

}

//-----------------------------------------------------------------------------
TLCalc::~TLCalc() {
}


//-----------------------------------------------------------------------------
void TLCalc::InitHistograms() {

  char   name[100];
  double qent;
  TH2F*  h_pt_vs_m;
  TFile* f;
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  fHist.fLhMass     = new TH1F(Form("h_lh_mass_%i",fMode) ,"likelihood MASS",200,-100,0);
  fHist.fLhPt       = new TH1F(Form("h_lh_pt_%i",fMode)   ,"likelihood PT"  ,200,-100,0);
  fHist.fLh2d       = new TH1F(Form("h_lh_2d_%i",fMode)   ,"likelihood 2D"  ,400,-200,0);
  fHist.fNHighClose = new TH2F(Form("nhc_vs_nhm_%i",fMode),"n(high close) vs N(high mass)",
			       50, 0,50, 50,0,50);

  for (int i=0; i<2; i++) {
    sprintf(name,"nev_%i_mode_%i",i,fMode);
    fHist.fNEvents  [i] = new TH1F(name,name,100,0,100);
    sprintf(name,"nhm_%i_mode_%i",i,fMode);
    fHist.fNHighMass[i] = new TH1F(name,name,100,0,100);
  }

  if (fMode == 1) {
//-----------------------------------------------------------------------------
// PYTHIA + TGzzAnaModule...
//-----------------------------------------------------------------------------
    h_pt_vs_m = gh2("hist/zzx/zzx_pythia_4l.hist","ZZ","zz_1/pt_vs_m");
  }
  else if (fMode == 2) {
//-----------------------------------------------------------------------------
// use histogram produced by MC@NLO+HERWIG no acceptance cuts
//-----------------------------------------------------------------------------
    f = TFile::Open("hist/zzx/HERVB_zz.root");
    h_pt_vs_m = (TH2F*) f->Get("h21");
  }

//-----------------------------------------------------------------------------
// prepare probability distributions: 2D
//-----------------------------------------------------------------------------
  fHist.fPtVsMass = (TH2F*) h_pt_vs_m->Clone(Form("zz_1_pt_vs_m_%i",fMode));

  qent = fHist.fPtVsMass->Integral();
  fHist.fPtVsMass->Rebin2D(2,2);
  fHist.fPtVsMass->Scale(1./qent);

  fHist.fMass = h_pt_vs_m->ProjectionX(Form("h1_px_%i",fMode));
  qent = fHist.fMass->GetEntries();
  fHist.fMass->Scale(1./qent);
  // this selects MZZ=325 bin
  fHist.fPt = h_pt_vs_m->ProjectionY(Form("h1_py_%i",fMode),31,32);
  qent = fHist.fPt->GetEntries();
  fHist.fPt->Scale(1./qent);

}

//-----------------------------------------------------------------------------
int TLCalc::calc_lh_2d(int NMax) {
  TObject  *o; 
  double    qent;
  TFile*    f;
  TH2      *h1;

  h1 = fHist.fPtVsMass;

//   TCanvas* c1 = new_slide("c1","c1",2,1,1000,500);
//   TPad*    p1 = (TPad*) c1->GetPrimitive("p1");

//  p1->cd(1);
//   h1->GetXaxis()->SetRangeUser(0,1000);
//   h1->Draw("box");
//-----------------------------------------------------------------------------
// need to generate random numbers according to this histogram
//                    M     PT
//-----------------------------------------------------------------------------
  double  xbin_width, ybin_width, xmin, ymin, m, pt, prob_hm, p;
  double  tot_prob_new;
  double  mass[100], prob[100];
  double  tot_prob;

  int     ix, iy, nn(0), nhm, nlm, nevents, npassed(0);

  xbin_width = h1->GetXaxis()->GetBinWidth(1);  // all bins are the same
  ybin_width = h1->GetYaxis()->GetBinWidth(1);  // all bins are the same

  xmin = h1->GetXaxis()->GetXmin();
  ymin = h1->GetYaxis()->GetXmin();

  fLh2dData = 0;

  for (int i=4; i<8; i++) {
    //  for (int i=1; i<8; i++) {
    ix           = (int) ((fData[i][0]-xmin)/xbin_width + 1);
    iy           = (int) ((fData[i][1]-ymin)/ybin_width + 1);
    tot_prob     = h1->GetBinContent(ix,iy);
    fLh2dData   += 2*TMath::Log(tot_prob);
  }

  printf(" >>> fLh2dData = %10.5f\n",fLh2dData);
//-----------------------------------------------------------------------------
// now generate pseudoexperiments
//-----------------------------------------------------------------------------
  TRandom3  r3;

  double lh;

  for (int ipe=0; ipe<NMax; ipe++) {
    nhm          = 0;
    nlm          = 0;
    tot_prob     = 1;
    prob_hm      = 1;
    tot_prob_new = 1;

    //     nevents   = r3.Poisson(5.8);
    nevents   = r3.Poisson(fNMean);

    fHist.fNEvents[0]->Fill(nevents);

    for (int ievent=0; ievent<nevents; ievent++) {

      h1->GetRandom2(m,pt);

      ix    = (int) ((m -xmin)/xbin_width + 1);
      iy    = (int) ((pt-ymin)/ybin_width + 1);
      p     = h1->GetBinContent(ix,iy);

      if (m > 300) {
	mass[nhm] = m;
	prob[nhm] = p;
	nhm       = nhm+1;
	//	prob_hm   = prob_hm*p ;
	tot_prob_new = tot_prob_new*p;
      }
      else {
	tot_prob  = tot_prob*p;
	nlm       = nlm+1;
      }

      //      tot_prob_new = tot_prob_new*p;
    }

    fHist.fNHighMass[0]->Fill(nhm);

    //    if ((nhm < 4) || (nlm < 4)) {
    //    if ((nevents < 8)) {
    //    if ((nevents < 8) || (nhm != 4)) {
    if (nhm < 4) {
      tot_prob     = 0.5;
      tot_prob_new = 0.5;
    }
    else {
//-----------------------------------------------------------------------------
// check if all masses are within 15 GeV
// first order masses
//-----------------------------------------------------------------------------
      for (int i1=0; i1<nhm-1; i1++) {
	for (int i2=i1+1; i2<nhm; i2++) {
	  if (mass[i1] > mass[i2]) {
	    m        = mass[i1];
	    mass[i1] = mass[i2];
	    mass[i2] = m;

	    p        = prob[i1];
	    prob[i1] = prob[i2];
	    prob[i2] = p;
	  }
	}
      }
//-----------------------------------------------------------------------------
// then count 'nn' - number of events within 20 GeV window
//-----------------------------------------------------------------------------
      nn = 0;
      for (int i1=0; i1<nhm-3; i1++) {
	nn      = 1;
	prob_hm = prob[i1];
	for (int i2=i1+1; i2<nhm; i2++) {
	  if (mass[i2]-mass[i1] < fMassWindow) {
	    nn++;
	    prob_hm = prob_hm*prob[i2];
	  }
	  else {
	    break;
	  }
	}
					// found a solution
 	if (nn >= 4) {
 	  break;
 	}
      }

      //      if (nn >= 4) {
      if (nn >= 4) {
	// tot_prob = tot_prob*prob_hm;
	tot_prob = prob_hm;
      }
      else {
	tot_prob     = 0.05;
	tot_prob_new = 0.05;
      }
    }
    
    //    lh = 2*TMath::Log(tot_prob);
    lh = 2*TMath::Log(tot_prob_new);

    fHist.fLh2d->Fill(lh);

    if (lh < fLh2dData) {
      npassed++;
      printf("--- ipe = %10i ",ipe);
      printf("nevents,nhm,nn = %5i %5i %5i ",nevents,nhm,nn);
      for (int ii=0; ii<nhm; ii++) {
	printf (" %12.5f",mass[ii]) ; 
      }
      printf("\n");

      fHist.fNEvents[1]->Fill(nevents);
      fHist.fNHighMass[1]->Fill(nhm);
      fHist.fNHighClose->Fill(nhm,nn);
    }
  }

  printf(" >>> npassed, ntotal = %10i %10i \n",npassed,NMax);
//   p1->cd(2);
//   gPad->SetLogy(1);
//  fHist.fLh2d->Draw();
  return 0;
}


//-----------------------------------------------------------------------------
int TLCalc::calc_lh_pt(int NMax) {
  TObject  *o; 
  double    qent;
  int       npassed(0);
  TH1*      h1; 
  TFile     *f;

  h1 = fHist.fPt;

  //  TCanvas* c1 = new_slide("calc_lh_pt","calc_lh_pt",2,1,1000,500);
  //  TPad*    p1 = (TPad*) c1->GetPrimitive("p1");

//   p1->cd(1);
//   h1->GetXaxis()->SetRangeUser(0,400);
//   h1->Draw();

  //  return 0;
//-----------------------------------------------------------------------------
// generate random numbers according to this histogram
//-----------------------------------------------------------------------------
  double  xw, xmin, pt, p;

  double  prob;

  int     ix, nevents;

  xw   = h1->GetXaxis()->GetBinWidth(1);  // all bins are the same
  xmin = h1->GetXaxis()->GetXmin();

  fLhPtData = 0;

  for (int i=0; i<8; i++) {
    ix         = (int) ((fData[i][1]-xmin)/xw + 1);
    prob       = h1->GetBinContent(ix);
    fLhPtData += 2*TMath::Log(prob);

    printf("i,data,xmin,xw,ix,prob : %3i %10.3f %10.3f %10.3f %4i  %12.5e\n",
	   i,fData[i][1],xmin,xw,ix,prob);
  }

  printf(" >>> fLhPtData = %10.5f\n",fLhPtData);

  //  return 0;
//-----------------------------------------------------------------------------
// now generate pseudoexperiments
//-----------------------------------------------------------------------------
  TRandom3  r3;

  double lh;

  for (int ipe=0; ipe<NMax; ipe++) {
    //    nevents   = r3.Poisson(5.8);
    nevents   = r3.Poisson(fNMean);
    if (nevents < 8) {
      lh = -0.5;
    }
    else {
      lh        = 0;
      for (int ievent=0; ievent<nevents; ievent++) {
	pt        = h1->GetRandom();
	ix        = (int) ((pt -xmin)/xw + 1);
	p         = h1->GetBinContent(ix);
	lh       += 2*TMath::Log(p);
      }
    }
    
    fHist.fLhPt->Fill(lh);

    if (lh < fLhPtData) {
      npassed += 1;
    }
  }

  printf(" >>> lh_pt: npassed, ntotal = %10i %10i \n",npassed,NMax);
  return 0;
}
//-----------------------------------------------------------------------------
int TLCalc::calc_lh_pt_hm(int NMax) {
  TObject  *o; 
  double    qent;
  int       npassed(0);
  TH1*      h1; 
  TFile     *f;

  h1 = fHist.fPt;

  //  TCanvas* c1 = new_slide("calc_lh_pt","calc_lh_pt",2,1,1000,500);
  //  TPad*    p1 = (TPad*) c1->GetPrimitive("p1");

//   p1->cd(1);
//   h1->GetXaxis()->SetRangeUser(0,400);
//   h1->Draw();

  //  return 0;
//-----------------------------------------------------------------------------
// generate random numbers according to this histogram
//-----------------------------------------------------------------------------
  double  xw, xmin, pt, p;

  double  prob;

  int     ix, nevents;

  xw   = h1->GetXaxis()->GetBinWidth(1);  // all bins are the same
  xmin = h1->GetXaxis()->GetXmin();

  fLhPtData = 0;

  for (int i=4; i<8; i++) {
    ix         = (int) ((fData[i][1]-xmin)/xw + 1);
    prob       = h1->GetBinContent(ix);
    fLhPtData += 2*TMath::Log(prob);

    printf("i,data,xmin,xw,ix,prob : %3i %10.3f %10.3f %10.3f %4i  %12.5e\n",
	   i,fData[i][1],xmin,xw,ix,prob);
  }

  printf(" >>> fLhPtData = %10.5f\n",fLhPtData);

  //  return 0;

//-----------------------------------------------------------------------------
// now generate pseudoexperiments
//-----------------------------------------------------------------------------
  TRandom3  r3;

  double lh;

  for (int ipe=0; ipe<NMax; ipe++) {
    nevents   = 4;
    lh        = 0;
    for (int ievent=0; ievent<nevents; ievent++) {
      pt        = h1->GetRandom();
      ix        = (int) ((pt -xmin)/xw + 1);
      p         = h1->GetBinContent(ix);
      lh       += 2*TMath::Log(p);
    }
    
    fHist.fLhPt->Fill(lh);
    
    if (lh < fLhPtData) {
      npassed += 1;
    }
  }
  
  printf(" >>> lh_pt: npassed, ntotal = %10i %10i \n",npassed,NMax);
  return 0;
}



//-----------------------------------------------------------------------------
int TLCalc::calc_lh_mass(int NMax) {
  TObject  *o; 
  double    qent;
  TH1*      h1;

  h1 = fHist.fMass;

//   TCanvas* c1 = new_slide("calc_lh_mass","calc_lh_mass",2,1,1000,500);
//   TPad*    p1 = (TPad*) c1->GetPrimitive("p1");

//   p1->cd(1);
//   h1->GetXaxis()->SetRangeUser(0,1000);
//   h1->Draw();

  //  return 0;
//-----------------------------------------------------------------------------
// generate random numbers according to this histogram
//-----------------------------------------------------------------------------
  double  xw, xmin, m, mass[100], prob[100], p;

  int     ix, nevents, nhm, nlm, nn, npassed(0);

  xw   = h1->GetXaxis()->GetBinWidth(1);  // all bins are the same
  xmin = h1->GetXaxis()->GetXmin();

  fLhMassData = 0;

  for (int i=4; i<8; i++) {
    m            = fData[i][0];
    ix           = (int) ((m-xmin)/xw + 1);
    p            = h1->GetBinContent(ix);
    fLhMassData += 2*TMath::Log(p);

    printf("i,data,xmin,xw,ix,prob : %3i %10.3f %10.3f %10.3f %4i  %12.5e\n",
	   i,m,xmin,xw,ix,p);
  }

  printf(" >>> fLhMassData = %10.5f\n",fLhMassData);

  //  return 0;
//-----------------------------------------------------------------------------
// now generate pseudoexperiments
//-----------------------------------------------------------------------------
  TRandom3  r3;

  double lh, tot_prob, tot_prob_new, prob_hm;

  for (int ipe=0; ipe<NMax; ipe++) {
    //    nevents   = r3.Poisson(5.8);
    nevents   = r3.Poisson(fNMean);
    lh        = 0;

    nhm          = 0;
    nlm          = 0;
    tot_prob     = 1;
    tot_prob_new = 1;
    prob_hm      = 1;

    fHist.fNEvents[0]->Fill(nevents);

    for (int ievent=0; ievent<nevents; ievent++) {
      m = h1->GetRandom();
      ix    = (int) ((m -xmin)/xw + 1);
      p     = h1->GetBinContent(ix);

      if (m > 300) {
	mass[nhm] = m;
	prob[nhm] = p;
	nhm           = nhm+1;
	tot_prob_new  = tot_prob_new*p;
	//	prob_hm   = prob_hm*p ;
      }
      else {
	tot_prob  = tot_prob*p;
	nlm       = nlm+1;
      }
    }
    
    fHist.fNHighMass[0]->Fill(nhm);

    //    if ((nhm < 4) || (nlm < 4)) {
    //    if ((nevents < 8)) {
    //    if ((nevents < 8) || (nhm != 4)) {
    if (nhm < 4) {
      tot_prob     = 0.5;
      tot_prob_new = 0.5;
    }
    else {
//-----------------------------------------------------------------------------
// check if all masses are within 15 GeV
// first order masses
//-----------------------------------------------------------------------------
      for (int i1=0; i1<nhm-1; i1++) {
	for (int i2=i1+1; i2<nhm; i2++) {
	  if (mass[i1] > mass[i2]) {
	    m        = mass[i1];
	    mass[i1] = mass[i2];
	    mass[i2] = m;

	    p        = prob[i1];
	    prob[i1] = prob[i2];
	    prob[i2] = p;
	  }
	}
      }
//-----------------------------------------------------------------------------
// then count events within 15 GeV window
//-----------------------------------------------------------------------------
      nn = 0;
      for (int i1=0; i1<nhm-3; i1++) {
	nn      = 1;
	prob_hm = prob[i1];
	for (int i2=i1+1; i2<nhm; i2++) {
	  if (mass[i2]-mass[i1] < fMassWindow) {
	    nn++;
	    prob_hm = prob_hm*prob[i2];
	  }
	  else {
	    break;
	  }
	}
	if (nn >= 4) {
	  break;
	}
      }
      
      if (nn >= 4) {
	// tot_prob = tot_prob*prob_hm;
	tot_prob = prob_hm;
      }
      else {
	tot_prob     = 0.05;
	tot_prob_new = 0.05;
      }
    }
    
    //    lh = 2*TMath::Log(tot_prob);
    lh = 2*TMath::Log(tot_prob_new);

    fHist.fLhMass->Fill(lh);

    if (lh < fLhMassData) {
      npassed += 1;
      fHist.fNEvents  [1]->Fill(nevents);
      fHist.fNHighMass[1]->Fill(nhm);
      fHist.fNHighClose->Fill(nhm,nn);
    }
  }

  printf(" >>> npassed, ntotal = %10i %10i \n",npassed,NMax);
  return 0;
}




//-----------------------------------------------------------------------------
int TLCalc::calc_4l_prob(double MinMass, int NMax) {
  TObject  *o; 
  double    qent;
  TH1*      h1;

  h1 = fHist.fMass;
//-----------------------------------------------------------------------------
// generate random numbers according to this histogram
//-----------------------------------------------------------------------------
  double  xw, xmin, m, mass[100], prob[100], p;

  int     ix, nevents, nhm, nlm, nn, npassed(0);

  xw   = h1->GetXaxis()->GetBinWidth(1);  // all bins are the same
  xmin = h1->GetXaxis()->GetXmin();

  fLhMassData = 0;
//-----------------------------------------------------------------------------
// now generate pseudoexperiments
//-----------------------------------------------------------------------------
  TRandom3  r3;

  double lh, tot_prob, tot_prob_new, prob_hm;

  for (int ipe=0; ipe<NMax; ipe++) {
    lh           = 0;
    nhm          = 0;
    nlm          = 0;
    tot_prob     = 1;
    tot_prob_new = 1;
    prob_hm      = 1;

    //    nevents   = r3.Poisson(5.8);
    nevents      = r3.Poisson(fNMean);

    fHist.fNEvents[0]->Fill(nevents);

    for (int ievent=0; ievent<nevents; ievent++) {
      m = h1->GetRandom();
      ix    = (int) ((m -xmin)/xw + 1);
      p     = h1->GetBinContent(ix);

      if (m > MinMass) {
	nhm += 1;
      }
    }

    fHist.fNHighMass[0]->Fill(nhm);
    if (nhm >= 4) {
      npassed += 1;
    }
  }

  printf(" >>> npassed, ntotal = %10i %10i \n",npassed,NMax);
  return 0;
}

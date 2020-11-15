#include "TMath.h"

/*
Float_t protonMass{0.938272};                       // GeV/c2
Float_t protonMass2 = protonMass*protonMass;          
Float_t nucleonMass{0.931494};                      // GeV/c2                
Float_t KineticEnergy{8.0};                         // GeV
*/


//STRIGANOV CODE -------------------------------------------------------STRIGANOV CODE------------------------STRIGANOV CODE-----------------------------------------------------STRIGANOV CODE---------------------



Double_t rTan(Double_t xT, Double_t rootS){
  Double_t B[10] = {-9999999.,.3060,.120,.0552,2.720,0.758,-0.680,1.54,0.594,2.87}; //increase dimension by 1 for fortran->c++
  Double_t sT = 3.752;
  Double_t aH = B[1]*TMath::Exp(-B[2]*xT);
  Double_t bH = B[3]*TMath::Exp(B[4]*xT);
  Double_t cH = B[5] + B[6]*xT + B[7]*xT*xT;
  Double_t dH = B[8]*TMath::Exp(B[9]*xT);
  Double_t q = rootS - sT;
  Double_t r1 = 1. - TMath::Exp(-aH * TMath::Power(q,bH));
  Double_t r2 = TMath::Exp(cH*q - dH);
  Double_t _rTan = 1./(1. - TMath::Exp(-r1*r2));
  // foo << "inside rtan " << aH << " " << bH << " " << cH << " " << dH << " " << q << " " << r1 << " " << r2 << " " << _rTan << std::endl;
  return _rTan;
}
		    
Double_t tan0(Double_t xF, Double_t pT, Double_t rootS){
  // just code A for antiprotons from Striganov

  // added a dummy value in front, making the array one larger than Sergei
  // but then I can copy his indices
  Double_t _tan0{0.};
  // tan ng, J Ph G, doc-db 8403, table 2, pbar line
  Double_t A[13] = {-99999999.,1.05e-04, 10.1, 0.5, 7.90, 0.465, 3.7e-02, 2.31, 1.4e-2, 3.02e-02, 3.19, 0.399, 8.39};
  Double_t mP(0.93827);
  Double_t smx = 3.0*mP;
  Double_t g00 = 3.15;
  Double_t sT = 3.752;

  // foo << "inside tan0 with args " << xF << " " << pT << " " << rootS << std::endl;
  if (rootS <= sT){
    _tan0 = 0;
  }
  if (rootS > sT){
    Double_t s = rootS*rootS;
    Double_t eMax =  (s - smx*smx + mP*mP)/(2.*rootS); // this is exactly the standard eMax for nucleon-nucleon collisions
    // not used in this version, striganov comment Double_t pMax = TMath::Sqrt(eMax*eMax - mP*mP);
    Double_t xM = mP/eMax;
    Double_t pCM = xF*rootS/2.;
    Double_t eCM = TMath::Sqrt(pCM*pCM + mP*mP + pT*pT); // different scaling variable, this is Tan and Ng's x - Xr.
    Double_t x = eCM/eMax;
    //foo << "eCM, emax, x = " << eCM << " " << eMax << " " << x << std::endl;
    if (x < 1){
      //foo << " inside x<1 clause" << std::endl;
      Double_t xT = (x - xM)/(1. - xM);
      Double_t f1 = 0;
      if ( (A[3] - x) > 0 ){
	f1 = A[1]*TMath::Exp(-A[2]*x);
      }
      Double_t f = (g00 - A[1])*(TMath::Power( 1. - x , A[4] )) + f1;
      //foo << "f, x, g00,a1,a4,f1 " << f << " "  << x << " " << g00 << " " << A[1] << " " << A[4] << " " << f1 << std::endl;
      Double_t aH = A[5]*TMath::Exp(-A[6]*x) + A[7]*TMath::Exp(A[8]*x);
      Double_t bH = A[9]*TMath::Exp(-A[10]* (x + A[11])) *TMath::Power((x + A[11]),A[12]);
      if ((rootS - 15) <= 0. ){
	Double_t rLE = rTan(xT,rootS);
	_tan0 = f*TMath::Exp(-(aH+bH*pT)*pT)*rLE;
	//foo << "inside rootS < 15 clause" << std::endl;
	//foo << "f,aH,bH,pT, rle " << f << " " << aH << " " << bH << " " << pT << " " << rLE<< std::endl;
      }
      if ( (rootS - 15) > 0.){
	//foo << "inside rootS>15 clause" << std::endl;
	_tan0 = f*TMath::Exp(-(aH+bH*pT)*pT);
	//foo << "f,aH,bH,pT " << f << " " << aH << " " << bH << " " << pT << std::endl;
      }
    }
  }
  //foo << "leaving tan0 = " << _tan0 << std::endl;
  return _tan0;
}

Double_t tanng1(Double_t xF, Double_t pT, Double_t rootS) {
  //directly from striganov
  Double_t p0 = 0.2;
  //
  // next statement in Sergei's code is an arithmetic if in FORTRAN to allow for more than one particle type, I ignore
  Double_t A;
  Double_t _tanng1;
  //foo << "inside tanng1 with xF,pT, rootS and p0 = " << xF << " " << pT << " " << rootS << " " << p0 << std::endl;
  if ( (pT-p0) <= 0){
    A = tan0(xF,p0,rootS)*TMath::Exp(1.4*p0);
    //foo << "executing first clause " << p0 << " " << xF << " " << rootS << " " << TMath::Exp(1.4*p0) << " " << pT << " " <<  A << std::endl;
    _tanng1 = A*TMath::Exp(-1.4*pT)*1000.*29.; // no idea where the 29 comes from but it seems to get removed later. Go with exact transcription for now
  }
  else{
    _tanng1 = tan0(xF,pT,rootS)*1000.*29.; // again, no idea where the 29 comes from but it seems to get removed later. Go with exact transcription for now
    //foo << " executing second clause " << " " << _tanng1 << std::endl;
  }
  return _tanng1;
}

Double_t xfam(Double_t p0,Double_t p, Double_t tet)
{
  Double_t mP(0.938272); 
  Double_t smx = 3.*mP;
  //foo << "p0,p,tet = " << p0 << " " << p << " " << tet << std::endl;
  // only used in def'n of eMax, which is then not used Double_t smx = 3.*mP;
  Double_t e0 = TMath::Sqrt(p0*p0 + mP*mP);
  Double_t s = 2.*mP*e0 + mP*mP + mP*mP;
  Double_t rootS = TMath::Sqrt(s);
  Double_t eMax = (s + mP*mP - smx*smx)/(2.*rootS);
  Double_t gammaCM = (e0 + mP)/rootS;
  Double_t gammaBCM = p0/rootS;
  Double_t ePi = TMath::Sqrt(p*p + mP*mP);
  Double_t pmp = p*TMath::Cos(tet);
  Double_t pzcm = -gammaBCM*ePi + gammaCM*pmp;
  Double_t _xfam = 2.*pzcm/rootS;
  //foo << "emax, gacm,gabcm, epi,p,pmp,sqs,xfam= " << eMax << " " << gammaCM << " " << gammaBCM << " " << ePi << " " <<
  //  p << " " << pmp << " " << rootS << " " << _xfam << std::endl;
  return _xfam;
}
    
Double_t flrgan1(Double_t pLab, Double_t theta, Double_t xlq){
  Double_t mP = 0.938272;
  Double_t sinTheta = TMath::Sin(theta);
  Double_t _flrgan1(0.);
  Double_t fc(0.);
  //foo << "inside flrgan1 "  << pLab << " " << theta << " " << xlq << " " << sinTheta << " " << xlq-sinTheta 
  //    << std::endl;
  if (xlq >= sinTheta){
    Double_t dsqs = 2.*mP*(pLab - 7.*mP);
    Double_t x01;
    Double_t x02;
    if (pLab <= 12.9){
      x01 = 0.023*(dsqs - 2.93);
      x02 = 0.5598*TMath::Exp(-0.664*dsqs);
    }
    else{
      //
      // this is exactly what he wrote, have to figure out what's going on.  
      x01 = 0.15;
      x02 = 0.0002;
      x01 = 0.2059;
      x02 = 0.00021;
    }
    Double_t c2 =  4.5;
    Double_t c1 = x02;
    Double_t x1 = xlq - x01;
    //      Double_t fc(0.);
    if (x1 >=0){
      fc = TMath::Power(x1,c2)/(TMath::Power(x1,c2) + c1);
    }
    _flrgan1 = 9674.7*TMath::Exp(-10.127*xlq)*TMath::Power((1. - xlq/181.),9)*fc;
  }
  // foo << "pieces = " << fc << " " << TMath::Exp(-10.127*xlq) << " " << TMath::Power((1. - xlq/181.),9) << " " << _flrgan1 << std::endl;
  return _flrgan1;
}

Double_t plmax(Double_t pLab, Double_t theta){
  Double_t _plmax(0.);
  Double_t mP = 0.93827;
  Double_t smx = 3.*mP;
  Double_t e0 = TMath::Sqrt(pLab*pLab + mP*mP);
  Double_t s = 2.*mP*e0 + mP*mP + mP*mP;
  Double_t rootS = TMath::Sqrt(s);
  Double_t eMax = (s + mP*mP - smx*smx)/(2.*rootS);
  //not used Double_t pMax = TMath::Sqrt(eMax*eMax - mP*mP);
  Double_t gacm = (e0 + mP)/rootS;
  Double_t gabcm = pLab/rootS;
  Double_t c2 = -gabcm*gabcm + TMath::Power((gacm/TMath::Cos(theta)),2);
  Double_t c1 = gabcm*eMax;
  Double_t disc = gacm*gacm*(gabcm*gabcm*mP*mP + (eMax*eMax - gacm*gacm*mP*mP)/TMath::Power(TMath::Cos(theta),2)); 
  if (disc >= 0) {_plmax = (c1 + TMath::Sqrt(disc))/c2 ;}
  return _plmax;
}

Double_t getWeight(const Double_t labMomentumInMeV, const Double_t cosTheta, const Double_t pBeamInMeV=8.89*1000, const Int_t returnWeight=1){
  Double_t pBeam = pBeamInMeV/1000.;
  Double_t labMomentum = labMomentumInMeV/1000.;

  //duperray def'n of eInc is kinetic energy of the incoming proton.
  Double_t totalEnergy = TMath::Sqrt(pBeam*pBeam + protonMass2);
  Double_t energy = TMath::Sqrt(labMomentum*labMomentum + protonMass2); // for consistency so I can cut/paste code from strigNorm.cc
  Double_t eInc = totalEnergy - protonMass;
  Double_t mandelS = 2*protonMass*(protonMass + totalEnergy);
  Double_t sqrtS = TMath::Sqrt(mandelS);

  Double_t totalInelastic(1.539e6);                       // xSec totale inelastica
  Double_t scaleFac(1.);
  Double_t sinTheta = TMath::Sqrt(1. - cosTheta*cosTheta);
  Double_t theta = TMath::ACos(cosTheta);
  Double_t xlq =   (TMath::Sqrt(labMomentum*labMomentum + protonMass2) - labMomentum*cosTheta)/protonMass;        
  Double_t pT = labMomentum*sinTheta;               //componente del momento perpendicolare al moto 
  Double_t xfa = xfam(pBeam, labMomentum, theta);
  Double_t firstPiece(0.);
  Double_t secondPiece(0.);

  //costruisco la sezione d'urto invariante secondo il modello di Sergei, primo e secondo pezzo
  if (pBeam <= 12.){
    firstPiece = tanng1(xfa,pT,sqrtS);
  }
  if (pBeam > 12.){
    // default to duperray.  his version has major bug: MeV and GeV units mixed.  Will discuss.
    firstPiece = 0.;
    //foo << " less than 12 " << " xfa = " << xfa << " pT = " << pT << " sqrtS = " << sqrtS << " tanng1 = " << firstPiece << std::endl;

    //  foo << " greater than 12 " << std::endl;
  }
  // now he has a set of choices:
  // 1) p2max = 1.3 with a comment that says used for Mu2e simulation 2011-2012
  // 2) p2max = 2 with a comment that new data (Kiselev?) shows this is better for labMomentum > 1.2 GeV/c only
  // 3) p2max = 2. with no comment
  // the code I have has 1.3
  // but changing from 1.3 to 2 causes discontinuities, other problems...leave at 2 for all data
  Double_t p2max(2.);
  Double_t p2c = 1.1;
  Double_t p0max = plmax(pBeam, theta);
  Double_t p1max = TMath::Max(p2max,p0max);
	
  //foo << " p2max,p0max, p1max = " << p2max << " " << p0max << " " << p1max << std::endl;
  if (labMomentum >= p2c){
    Double_t xlq0 = (TMath::Sqrt(p2c*p2c + protonMass2) - p2c*TMath::Cos(theta))/protonMass;
    //  foo << "in xlq0 clause " << pBeam << 
    //  " " << theta << " " << xlq0 << " " << p2c << std::endl; 
    Double_t f20 = flrgan1(pBeam,theta,xlq0);
    Double_t boundary = (labMomentum - p2c)/(p1max - p2c);
    if (boundary >=1){
      secondPiece = 0.;
    }
    else{
      secondPiece = TMath::Power( (1. - (labMomentum - p2c)/(p1max - p2c)),5)*f20;
      //foo << " ratio = " << (labMomentum - p2c)/(p1max - p2c) << " " << (labMomentum - p2c) << " " << (p1max - p2c) << std::endl;
      //foo << " foo = " << TMath::Power( (1. - (labMomentum - p2c)/(p1max - p2c)),5) << std::endl;
      //foo << " bar = " << labMomentum << " " << p2c << " " << p1max << std::endl;
      //foo << "corresponding in xlq0 clause " << pBeam << " " << theta << " " << xlq0 << " " << f20 << " " << secondPiece << std::endl;
    }
  }
  else{
    //	  foo << "in xlq clause with lab Momentum, p2c " << labMomentum << " " << theta << " " << p2c << " " <<  xlq << std::endl;
    secondPiece = flrgan1(pBeam,theta,xlq);
  }
  //  foo << " first piece = " << firstPiece << " " << " second piece = " << secondPiece << std::endl;
  
  Double_t invariantCrossSection = (firstPiece + secondPiece); // in millibarns, not microbarns
  //         invariantCrossSection = secondPiece;
  //		   invariantCrossSection = firstPiece;
  //
  // 1.539e6 being what Sergei thinks is the inelastic cross-section.  The point of this code is to figure out what needs to be done to normalize properly
  //  Double_t probThisBin = scaleFac*crossSection*labMomentum*labMomentum/totalEnergy/totalInelastic;
  //
  // this gets us d^3 sigma/ dp3 divided by the total cross-section; scaleFac optional, apply here if you want
  // invariant cross = E * d^3 sigma/dp^3; divide by E, multiply by p**2 sin theta for phase space p^2 dp.
  // then since dealing with weights have to multiply by 4 pi and the range in momentum, which I approximate as 0-5 GeV/c
  // integral = <average weight>*interval, average value thm of calculus
  Double_t probThisBin = 4.*TMath::Pi()*5.*invariantCrossSection*labMomentum*labMomentum/energy/totalInelastic;

  //  std::cout << "cos theta, sin theta = " << cosTheta << " " << sinTheta << std::endl;

  //  Double_t rawProb = labMomentum*labMomentum*(crossSection/energy);
  //resultProb  << "pbeam, first piece, second piece, energy = " << pBeam << " "  << firstPiece << " " << secondPiece << " " << energy << std::endl;
  //resultProb << labMomentum << " " << theta << " " << probThisBin*scaleFac/totalInelastic << " " << crossSection << " " << rawProb << std::endl;	  

  //	  	  probThisBin = labMomentum*labMomentum; // pure volume element for integral rho^2 drho
  // resultProb << labMomentum << " " << theta << " " << probThisBin*scaleFac/totalInelastic << " " << crossSection << std::endl;	  
  //  return probThisBin;
  Double_t result = -1.;
  if (returnWeight  == 0) {result = probThisBin;}
  if (returnWeight  == 1) {result = invariantCrossSection;}
  return result;
}




// common source for preparing t-channel likelihood function
// limit inputs -- using SM with SM single top
// as the null hypothesis and SM with SM single top and extra anomalous
// single top as the test hypothesis

const char *ename[15];
const char jesname[]="JES";
const char isrname[]="ISR";
const char fsrname[]="FSR";
const char pdfname[]="PDF";
const char luminame[]="LUMI";
const char btagname[]="BTAG";
const char wbbname[]="WBB";
const char wccname[]="WBB";  // Take these to be all correlated for now -- 38% Method 2
const char mistagname[]="MISTAG";
const char nonwname[]="NONW";
const char nonwflavname[]="NONWFLAV";
const char nonwantiname[]="NONWANTI";
const char bsvname[]="BSV";
const char ttbarname[]="TTBAR";
const char qsquaredname[]="QSQUARED";
const char drjjname[]="DRJJRW";
const char phxename[]="PHXE";
const char etaj2name[]="ETAJ2";

double nps_low[15];
double nps_high[15];
double lowsigma[15];
double highsigma[15];
int i;
double sfact;
int isyst;

TH1 *lowshape[15];
TH1 *highshape[15];

//-----------------------

csm_model* nullhyp = new csm_model();
csm_model* testhyp = new csm_model();
csm_model* nullhyp_pe = new csm_model();
csm_model* testhyp_pe = new csm_model();

/* first index:  
   0:  ttbar 
   1:  wbb
   2:  tchan
   3:  schan
   4:  data
   5:  wc(c)
   6:  mistags
   7:  ww
   8:  wz
   9:  zz
   10: nonw
   11: Z+jets

   Second index:

   0: central value
   1: isr less
   2: isr more
   3: fsr less
   4: fsr more
   5: jes minus
   6: jes plus
   7: bsv opt
   8: bsv pes   (symmetrized)
   9: qsquared +variation
   10: qsquared -variation 
   11: nonW flav+
   12: nonW flav-
   13: nonW jet-electrons
   14: nonW anti-electrons
   15: DRJJ reweight
   16: DRJJ symmetrized
   17: PHXE reweight
   18: PHXE symmetrized
   19: ETAJ2 reweight
   20: ETAJ2 symmetrized
*/


TH1F* histo[12][21];

cout << "reading histos" << endl;

TFile *wglfroot = (TFile*) new TFile("../hbooktmp/ltchan_central.root");
isyst = 0;
histo[0][isyst]  = (TH1F*)  wglfroot->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  wglfroot->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  wglfroot->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  wglfroot->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  wglfroot->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  wglfroot->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  wglfroot->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  wglfroot->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) wglfroot->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) wglfroot->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) wglfroot->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) wglfroot->Get("h120;1")->Clone("zjetscent");


TFile *isrless = (TFile*) new TFile("../hbooktmp/ltchan_isrless.root");
isyst = 1;
histo[0][isyst]  = (TH1F*)  isrless->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  isrless->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  isrless->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  isrless->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  isrless->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  isrless->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  isrless->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  isrless->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) isrless->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) isrless->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) isrless->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) isrless->Get("h120;1")->Clone("zjetscent");


TFile *isrmore = (TFile*) new TFile("../hbooktmp/ltchan_isrmore.root");
isyst = 2;
histo[0][isyst]  = (TH1F*)  isrmore->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  isrmore->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  isrmore->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  isrmore->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  isrmore->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  isrmore->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  isrmore->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  isrmore->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) isrmore->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) isrmore->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) isrmore->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) isrmore->Get("h120;1")->Clone("zjetscent");

TFile *fsrless = (TFile*) new TFile("../hbooktmp/ltchan_fsrless.root");
isyst = 3;
histo[0][isyst]  = (TH1F*)  fsrless->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  fsrless->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  fsrless->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  fsrless->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  fsrless->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  fsrless->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  fsrless->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  fsrless->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) fsrless->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) fsrless->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) fsrless->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) fsrless->Get("h120;1")->Clone("zjetscent");


TFile *fsrmore = (TFile*) new TFile("../hbooktmp/ltchan_fsrmore.root");
isyst = 4;
histo[0][isyst]  = (TH1F*)  fsrmore->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  fsrmore->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  fsrmore->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  fsrmore->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  fsrmore->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  fsrmore->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  fsrmore->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  fsrmore->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) fsrmore->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) fsrmore->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) fsrmore->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) fsrmore->Get("h120;1")->Clone("zjetscent");


TFile *jesm = (TFile*) new TFile("../hbooktmp/ltchan_jesm.root");
isyst = 5;
histo[0][isyst]  = (TH1F*)  jesm->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  jesm->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  jesm->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  jesm->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  jesm->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  jesm->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  jesm->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  jesm->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) jesm->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) jesm->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) jesm->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) jesm->Get("h120;1")->Clone("zjetscent");


TFile *jesp = (TFile*) new TFile("../hbooktmp/ltchan_jesp.root");
isyst = 6;
histo[0][isyst]  = (TH1F*)  jesp->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  jesp->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  jesp->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  jesp->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  jesp->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  jesp->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  jesp->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  jesp->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) jesp->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) jesp->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) jesp->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) jesp->Get("h120;1")->Clone("zjetscent");

// one-sided.  To symmetrize below

TFile *bsvopt = (TFile*) new TFile("../hbooktmp/ltchan_bsvopt.root");
isyst = 7;
histo[0][isyst]  = (TH1F*)  bsvopt->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  bsvopt->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  bsvopt->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  bsvopt->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  bsvopt->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  bsvopt->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  bsvopt->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  bsvopt->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) bsvopt->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) bsvopt->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) bsvopt->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) bsvopt->Get("h120;1")->Clone("zjetscent");

TFile *bsvpes = (TFile*) new TFile("../hbooktmp/ltchan_central.root");
isyst = 8;
histo[0][isyst]  = (TH1F*)  bsvpes->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  bsvpes->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  bsvpes->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  bsvpes->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  bsvpes->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  bsvpes->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  bsvpes->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  bsvpes->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) bsvpes->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) bsvpes->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) bsvpes->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) bsvpes->Get("h120;1")->Clone("zjetscent");

TFile *qsquared1 = (TFile*) new TFile("../hbooktmp/ltchan_qsquared1.root");
isyst = 9;
histo[0][isyst]  = (TH1F*)  qsquared1->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  qsquared1->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  qsquared1->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  qsquared1->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  qsquared1->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  qsquared1->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  qsquared1->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  qsquared1->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) qsquared1->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) qsquared1->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) qsquared1->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) qsquared1->Get("h120;1")->Clone("zjetscent");

// we only have one-sided qsquared systematics -- symmetrize them below.

TFile *qsquared2 = (TFile*) new TFile("../hbooktmp/ltchan_central.root");
isyst = 10;
histo[0][isyst]  = (TH1F*)  qsquared2->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  qsquared2->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  qsquared2->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  qsquared2->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  qsquared2->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  qsquared2->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  qsquared2->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  qsquared2->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) qsquared2->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) qsquared2->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) qsquared2->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) qsquared2->Get("h120;1")->Clone("zjetscent");

TFile *nonwflav1 = (TFile*) new TFile("../hbooktmp/ltchan_nonwflav1.root");
isyst = 11;
histo[0][isyst]  = (TH1F*)  nonwflav1->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  nonwflav1->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  nonwflav1->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  nonwflav1->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  nonwflav1->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  nonwflav1->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  nonwflav1->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  nonwflav1->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) nonwflav1->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) nonwflav1->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) nonwflav1->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) nonwflav1->Get("h120;1")->Clone("zjetscent");

TFile *nonwflav2 = (TFile*) new TFile("../hbooktmp/ltchan_nonwflav2.root");
isyst = 12;
histo[0][isyst]  = (TH1F*)  nonwflav2->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  nonwflav2->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  nonwflav2->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  nonwflav2->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  nonwflav2->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  nonwflav2->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  nonwflav2->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  nonwflav2->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) nonwflav2->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) nonwflav2->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) nonwflav2->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) nonwflav2->Get("h120;1")->Clone("zjetscent");

TFile *nonwanti1 = (TFile*) new TFile("../hbooktmp/ltchan_jetele.root");
isyst = 13;
histo[0][isyst]  = (TH1F*)  nonwanti1->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  nonwanti1->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  nonwanti1->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  nonwanti1->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  nonwanti1->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  nonwanti1->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  nonwanti1->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  nonwanti1->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) nonwanti1->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) nonwanti1->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) nonwanti1->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) nonwanti1->Get("h120;1")->Clone("zjetscent");

TFile *nonwanti2 = (TFile*) new TFile("../hbooktmp/ltchan_antiele.root");
isyst = 14;
histo[0][isyst]  = (TH1F*)  nonwanti2->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  nonwanti2->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  nonwanti2->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  nonwanti2->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  nonwanti2->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  nonwanti2->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  nonwanti2->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  nonwanti2->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) nonwanti2->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) nonwanti2->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) nonwanti2->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) nonwanti2->Get("h120;1")->Clone("zjetscent");

TFile *drjjrw = (TFile*) new TFile("../hbooktmp/ltchan_drjjrw.root");
isyst = 15;
histo[0][isyst]  = (TH1F*)  drjjrw->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  drjjrw->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  drjjrw->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  drjjrw->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  drjjrw->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  drjjrw->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  drjjrw->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  drjjrw->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) drjjrw->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) drjjrw->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) drjjrw->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) drjjrw->Get("h120;1")->Clone("zjetscent");
isyst = 16;
histo[0][isyst]  = (TH1F*)  wglfroot->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  wglfroot->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  wglfroot->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  wglfroot->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  wglfroot->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  wglfroot->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  wglfroot->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  wglfroot->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) wglfroot->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) wglfroot->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) wglfroot->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) wglfroot->Get("h120;1")->Clone("zjetscent");

TFile *phxerw = (TFile*) new TFile("../hbooktmp/ltchan_phxerw.root");
isyst = 17;
histo[0][isyst]  = (TH1F*)  phxerw->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  phxerw->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  phxerw->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  phxerw->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  phxerw->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  phxerw->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  phxerw->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  phxerw->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) phxerw->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) phxerw->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) phxerw->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) phxerw->Get("h120;1")->Clone("zjetscent");
isyst = 18;
histo[0][isyst]  = (TH1F*)  wglfroot->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  wglfroot->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  wglfroot->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  wglfroot->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  wglfroot->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  wglfroot->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  wglfroot->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  wglfroot->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) wglfroot->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) wglfroot->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) wglfroot->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) wglfroot->Get("h120;1")->Clone("zjetscent");

TFile *j2etarw = (TFile*) new TFile("../hbooktmp/ltchan_j2etarw.root");
isyst=19;
histo[0][isyst]  = (TH1F*)  j2etarw->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  j2etarw->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  j2etarw->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  j2etarw->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  j2etarw->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  j2etarw->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  j2etarw->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  j2etarw->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) j2etarw->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) j2etarw->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) j2etarw->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) j2etarw->Get("h120;1")->Clone("zjetscent");
isyst=20;
histo[0][isyst]  = (TH1F*)  wglfroot->Get("h10;1")->Clone("ttcentd");
histo[1][isyst]  = (TH1F*)  wglfroot->Get("h20;1")->Clone("wbbcent");
histo[2][isyst]  = (TH1F*)  wglfroot->Get("h30;1")->Clone("tchcent");
histo[3][isyst]  = (TH1F*)  wglfroot->Get("h40;1")->Clone("schcent");
histo[4][isyst]  = (TH1F*)  wglfroot->Get("h50;1")->Clone("datcent");
histo[5][isyst]  = (TH1F*)  wglfroot->Get("h60;1")->Clone("wcccent");
histo[6][isyst]  = (TH1F*)  wglfroot->Get("h70;1")->Clone("miscent");
histo[7][isyst]  = (TH1F*)  wglfroot->Get("h80;1")->Clone("wwcent ");
histo[8][isyst]  = (TH1F*) wglfroot->Get("h90;1")->Clone("wzcent ");
histo[9][isyst] = (TH1F*) wglfroot->Get("h100;1")->Clone("zzcent ");
histo[10][isyst] = (TH1F*) wglfroot->Get("h110;1")->Clone("nonwcen");
histo[11][isyst] = (TH1F*) wglfroot->Get("h120;1")->Clone("zjetscent");

// normalize syst. variations to the central value and 
// symmetrize errors 

// step 1 -- normalize all syst. variations to unit area
// (isyst=0 = central value)

double histsum;
int nbins;
int ibin;
for (isyst=1;isyst<21;isyst++)
{
  for (i=0;i<12;i++)
    { 
      nbins = histo[i][isyst]->GetNbinsX();
      histsum = 0;
      for (ibin=1;ibin<=nbins;ibin++)
	{
	  histsum += histo[i][isyst]->GetBinContent(ibin);
	}
      if (histsum>0)
	{
	  for (ibin=1;ibin<=nbins;ibin++)
	    {
	      histo[i][isyst]->SetBinContent(ibin,histo[i][isyst]->GetBinContent(ibin)/histsum);
	    }
	}
	       
    }
}

// step 2 : find the fractional change between + and - and apply these to the central value
// do separately for each + and - error

double hc,hcs,hcd;
for (i=0;i<12;i++)
{ 
  nbins = histo[i][0]->GetNbinsX();
  for (ibin=1;ibin<=nbins;ibin++)
    {
      hcs = histo[i][1]->GetBinContent(ibin)+histo[i][2]->GetBinContent(ibin);
      hcd = -histo[i][1]->GetBinContent(ibin)+histo[i][2]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][2]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][1]->SetBinContent(ibin,hc);
	}

      hcs = histo[i][3]->GetBinContent(ibin)+histo[i][4]->GetBinContent(ibin);
      hcd = -histo[i][3]->GetBinContent(ibin)+histo[i][4]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][4]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][3]->SetBinContent(ibin,hc);
	}

      hcs = histo[i][5]->GetBinContent(ibin)+histo[i][6]->GetBinContent(ibin);
      hcd = -histo[i][5]->GetBinContent(ibin)+histo[i][6]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][6]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][5]->SetBinContent(ibin,hc);
	}

      // bsv variation -- symmetrize it and take twice the difference since we
      // only have a 1-sided variation

      hcs = histo[i][7]->GetBinContent(ibin)+histo[i][8]->GetBinContent(ibin);
      hcd = -histo[i][7]->GetBinContent(ibin)+histo[i][8]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+2*hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][8]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-2*hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][7]->SetBinContent(ibin,hc);
	}

      // Q**2 variation -- now have up and down

      hcs = histo[i][9]->GetBinContent(ibin)+histo[i][10]->GetBinContent(ibin);
      hcd = -histo[i][9]->GetBinContent(ibin)+histo[i][10]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][10]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][9]->SetBinContent(ibin,hc);
	}

      // nonW flavor -- have up and down here

      hcs = histo[i][11]->GetBinContent(ibin)+histo[i][12]->GetBinContent(ibin);
      hcd = -histo[i][11]->GetBinContent(ibin)+histo[i][12]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][12]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][11]->SetBinContent(ibin,hc);
	}

      // jet electrons vs anti-electrons.  Symmetrize

      hcs = histo[i][13]->GetBinContent(ibin)+histo[i][14]->GetBinContent(ibin);
      hcd = -histo[i][13]->GetBinContent(ibin)+histo[i][14]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][14]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][13]->SetBinContent(ibin,hc);
	}

      // drjj -- one-sided reweight.  Symmetrize

      hcs = histo[i][15]->GetBinContent(ibin)+histo[i][16]->GetBinContent(ibin);
      hcd = -histo[i][15]->GetBinContent(ibin)+histo[i][16]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+2*hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][16]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-2*hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][15]->SetBinContent(ibin,hc);
	}

      // phx energy -- one-sided reweight.  Symmetrize

      hcs = histo[i][17]->GetBinContent(ibin)+histo[i][18]->GetBinContent(ibin);
      hcd = -histo[i][17]->GetBinContent(ibin)+histo[i][18]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+2*hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][18]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-2*hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][17]->SetBinContent(ibin,hc);
	}

      // Eta Jet 2

      hcs = histo[i][19]->GetBinContent(ibin)+histo[i][20]->GetBinContent(ibin);
      hcd = -histo[i][19]->GetBinContent(ibin)+histo[i][20]->GetBinContent(ibin);
      if (hcs>0)
	{
	  hc = histo[i][0]->GetBinContent(ibin)*(1+2*hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][20]->SetBinContent(ibin,hc);
	  hc = histo[i][0]->GetBinContent(ibin)*(1-2*hcd/hcs);
	  if (hc<0) hc=0;
	  histo[i][19]->SetBinContent(ibin,hc);
	}
    }
}


//-----------------------------------------------------------
// "ttbar background"


cout << "ttbar bg " << endl;

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}

ename[0] = luminame;
nps_low[0] = -0.06;
nps_high[0] = 0.06;

ename[1] = btagname;
nps_low[1] = -0.042;
nps_high[1] = 0.042;

ename[2] = jesname;
nps_low[2] = 0.098;
nps_high[2] = -0.091;
lowsigma[2] = -1;
highsigma[2] = 1;
highshape[2] = histo[0][6];
lowshape[2] = histo[0][5];

ename[3] = pdfname;
nps_low[3] = -0.024;
nps_high[3] = 0.024;

ename[4] = ttbarname;
nps_low[4] = -0.18;
nps_high[4] = 0.18;

ename[5] = drjjname;
lowsigma[5] = -1;
highsigma[5] = 1;
highshape[5] = histo[0][15];
lowshape[5] = histo[0][16];

ename[6] = etaj2name;
lowsigma[6] = -1;
highsigma[6] = 1;
highshape[6] = histo[0][19];
lowshape[6] = histo[0][20];

sfact = 1;

cout << "nullhyp_pe ttbar" <<endl;
nullhyp_pe->add_template(histo[0][0],sfact,7,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
cout << "testhyp_pe ttbar" <<endl;
testhyp_pe->add_template(histo[0][0],sfact,7,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
cout << "nullhyp ttbar" <<endl;
nullhyp->add_template(histo[0][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
cout << "testhyp ttbar" <<endl;
testhyp->add_template(histo[0][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");

/*------------------------------------------------------*/

// "Wbb background"

cout << "wbb bg " << endl;

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}
ename[0] = wbbname;
nps_low[0] = -0.38;
nps_high[0] = 0.38;

ename[1] = qsquaredname;
lowsigma[1] = -1;
highsigma[1] = 1;
lowshape[1] = histo[1][10];
highshape[1] = histo[1][9];

ename[2] = jesname;
nps_low[2] = .070;
nps_high[2] = -.076;

ename[3] = drjjname;
lowsigma[3] = -1;
highsigma[3] = 1;
highshape[3] = histo[1][15];
lowshape[3] = histo[1][16];

ename[4] = etaj2name;
lowsigma[4] = -1;
highsigma[4] = 1;
highshape[4] = histo[1][19];
lowshape[4] = histo[1][20];

sfact = 1;

nullhyp_pe->add_template(histo[1][0],sfact,5,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp_pe->add_template(histo[1][0],sfact,5,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
nullhyp->add_template(histo[1][0],sfact,1,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp->add_template(histo[1][0],sfact,1,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");


//-----------------------------------------------------------
// "Wc(c) background"

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}
ename[0] = wccname;
nps_low[0] = -0.38;
nps_high[0] = 0.38;
sfact = 1;

ename[1] = qsquaredname;
lowsigma[1] = -1;
highsigma[1] = 1;
lowshape[1] = histo[5][10];
highshape[1] = histo[5][9];

ename[2] = jesname;
nps_low[2] = .065;
nps_high[2] = -.058;

ename[3] = drjjname;
lowsigma[3] = -1;
highsigma[3] = 1;
highshape[3] = histo[5][15];
lowshape[3] = histo[5][16];

ename[4] = etaj2name;
lowsigma[4] = -1;
highsigma[4] = 1;
highshape[4] = histo[5][19];
lowshape[4] = histo[5][20];

nullhyp_pe->add_template(histo[5][0],sfact,5,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp_pe->add_template(histo[5][0],sfact,5,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
nullhyp->add_template(histo[5][0],sfact,1,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp->add_template(histo[5][0],sfact,1,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");

/*-----------------------------------------------------*/
// "Mistag background"

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}
ename[0] = mistagname;
nps_low[0] = -0.14;
nps_high[0] = 0.14;
sfact = 1;

ename[1] = jesname;
nps_low[1] = 0;
nps_high[1] = 0;
lowsigma[1] = -1;
highsigma[1] = 1;
lowshape[1] = histo[6][5];
highshape[1] = histo[6][6];

ename[2] = drjjname;
lowsigma[2] = -1;
highsigma[2] = 1;
highshape[2] = histo[6][15];
lowshape[2] = histo[6][16];

ename[3] = etaj2name;
lowsigma[3] = -1;
highsigma[3] = 1;
highshape[3] = histo[6][19];
lowshape[3] = histo[6][20];

// note -- not a Poisson process; events have different weights

nullhyp_pe->add_template(histo[6][0],sfact,4,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp_pe->add_template(histo[6][0],sfact,4,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
nullhyp->add_template(histo[6][0],sfact,1,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp->add_template(histo[6][0],sfact,1,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");


//-----------------------------------------------------------
// "WW background"

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}
ename[0] = luminame;
nps_low[0] = -0.06;
nps_high[0] = 0.06;
ename[1] = btagname;
nps_low[1] = -0.042;
nps_high[1] = 0.042;
sfact = 1;

ename[2] = jesname;
nps_low[2] = -.033;
nps_high[2] = .015;

// drjj doesn't look convincing here

ename[3] = etaj2name;
lowsigma[3] = -1;
highsigma[3] = 1;
highshape[3] = histo[7][19];
lowshape[3] = histo[7][20];

nullhyp_pe->add_template(histo[7][0],sfact,4,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp_pe->add_template(histo[7][0],sfact,4,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
nullhyp->add_template(histo[7][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp->add_template(histo[7][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");

//-----------------------------------------------------------
// "WZ background"

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}
ename[0] = luminame;
nps_low[0] = -0.06;
nps_high[0] = 0.06;
sfact = 1;
ename[1] = btagname;
nps_low[1] = -0.042;
nps_high[1] = 0.042;

// maybe a shape in jes too

ename[2] = jesname;
nps_low[2] = -.031;
nps_high[2] = .015;
lowsigma[2] = -1;
highsigma[2] = 1;
highshape[2] = histo[8][5];
lowshape[2] = histo[8][6];

ename[3] = etaj2name;
lowsigma[3] = -1;
highsigma[3] = 1;
highshape[3] = histo[8][19];
lowshape[3] = histo[8][20];

nullhyp_pe->add_template(histo[8][0],sfact,4,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp_pe->add_template(histo[8][0],sfact,4,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
nullhyp->add_template(histo[8][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp->add_template(histo[8][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");

//-----------------------------------------------------------
// "ZZ background"

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}
ename[0] = luminame;
nps_low[0] = -0.06;
nps_high[0] = 0.06;
ename[1] = btagname;
nps_low[1] = -0.042;
nps_high[1] = 0.042;
sfact = 1;

nullhyp_pe->add_template(histo[9][0],sfact,2,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp_pe->add_template(histo[9][0],sfact,2,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
nullhyp->add_template(histo[9][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp->add_template(histo[9][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");

//-----------------------------------------------------------
// "Z+jets background"

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}
ename[0] = wbbname;
nps_low[0] = -0.38;
nps_high[0] = 0.38;

ename[1] = jesname;
nps_low[1] = -.051;
nps_high[1] = .036;

ename[2] = drjjname;
lowsigma[2] = -1;
highsigma[2] = 1;
highshape[2] = histo[11][15];
lowshape[2] = histo[11][16];

ename[3] = etaj2name;
lowsigma[3] = -1;
highsigma[3] = 1;
highshape[3] = histo[11][19];
lowshape[3] = histo[11][20];

sfact = 1;

nullhyp_pe->add_template(histo[11][0],sfact,4,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp_pe->add_template(histo[11][0],sfact,4,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
nullhyp->add_template(histo[11][0],sfact,1,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp->add_template(histo[11][0],sfact,1,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");


//-----------------------------------------------------------
// "Non-W background "

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}
ename[0] = nonwname;
nps_low[0] = -0.6;
nps_high[0] = 0.6;

ename[1] = nonwflavname;
nps_low[1] = 0;
nps_high[1] = 0;
lowsigma[1] = -1;
highsigma[1] = 1;
lowshape[1] = histo[10][11];
highshape[1] = histo[10][12];

sfact = 1;

nullhyp_pe->add_template(histo[10][0],sfact,2,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp_pe->add_template(histo[10][0],sfact,2,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
nullhyp->add_template(histo[10][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");
testhyp->add_template(histo[10][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,0,"tchan");


//-----------------------------------------------------------
// "t-channel signal"

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}

// skip ISR and FSR shapes -- they are indistinguishable from no shape difference

ename[0] = isrname;
nps_low[0] = -0.01;
nps_high[0] = 0.03;
//lowsigma[0] = -1;
//highsigma[0] = 1;
//lowshape[0] = histo[2][1];
//highshape[0] = histo[2][2];

ename[1] = fsrname;
nps_low[1] = -0.015;
nps_high[1] = 0.0532;
//lowsigma[1] = -1;
//highsigma[1] = 1;
//lowshape[1] = histo[2][3];
//highshape[1] = histo[2][4];

ename[2] = btagname;
nps_low[2] = -0.042;
nps_high[2] = 0.042;

ename[3] = jesname;
nps_low[3] = -0.008;
nps_high[3] = 0.003;
lowsigma[3] = -1;
highsigma[3] = 1;
lowshape[3] = histo[2][5];
highshape[3] = histo[2][6];

ename[4] = luminame;
nps_low[4] = -0.06;
nps_high[4] = 0.06;

ename[5] = pdfname;
nps_low[5] = -0.017;
nps_high[5] = 0.011;

ename[6] = drjjname;
lowsigma[6] = -1;
highsigma[6] = 1;
highshape[6] = histo[2][15];
lowshape[6] = histo[2][16];

ename[7] = etaj2name;
lowsigma[7] = -1;
highsigma[7] = 1;
highshape[7] = histo[2][19];
lowshape[7] = histo[2][20];

sfact = 1;

testhyp_pe->add_template(histo[2][0],sfact,8,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,1,"tchan");
testhyp->add_template(histo[2][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,1,"tchan");

//-----------------------------------------------------------
// "s-channel signal"

for(i=0;i<15;i++)
{
  nps_low[i] = 0;
  nps_high[i] = 0;
  lowsigma[i] = 0;
  highsigma[i] = 0;
  lowshape[i] = 0;
  highshape[i] = 0;
}

// skip the ISR/FSR fits -- their shapes are indistinguishable from flat

ename[0] = isrname;
nps_low[0] = -0.01;
nps_high[0] = 0.0318;
//lowsigma[0] = -1;
//highsigma[0] = 1;
//lowshape[0] = histo[3][1];
//highshape[0] = histo[3][2];

ename[1] = fsrname;
nps_low[1] = -0.0148;
nps_high[1] = 0.0532;
//lowsigma[1] = -1;
//highsigma[1] = 1;
//lowshape[1] = histo[3][3];
//highshape[1] = histo[3][4];

ename[2] = btagname;
nps_low[2] = -0.042;
nps_high[2] = 0.042;

ename[3] = jesname;
nps_low[3] = 0.004;
nps_high[3] = -0.008;
lowsigma[3] = -1;
highsigma[3] = 1;
lowshape[3] = histo[3][5];
highshape[3] = histo[3][6];

ename[4] = luminame;
nps_low[4] = -0.06;
nps_high[4] = 0.06;

ename[5] = pdfname;
nps_low[5] = -0.01;
nps_high[5] = 0.01;

ename[6] = drjjname;
lowsigma[6] = -1;
highsigma[6] = 1;
highshape[6] = histo[3][15];
lowshape[6] = histo[3][16];

ename[7] = etaj2name;
lowsigma[7] = -1;
highsigma[7] = 1;
highshape[7] = histo[3][19];
lowshape[7] = histo[3][20];

sfact = 1;

testhyp_pe->add_template(histo[3][0],sfact,8,ename,nps_low,nps_high,
			 lowshape,lowsigma,highshape,highsigma,0,1,"tchan");
testhyp->add_template(histo[3][0],sfact,0,ename,nps_low,nps_high,
		      lowshape,lowsigma,highshape,highsigma,0,1,"tchan");

nullhyp->set_interpolation_style("tchan",CSM_INTERP_VERTICAL);
nullhyp_pe->set_interpolation_style("tchan",CSM_INTERP_VERTICAL);
testhyp->set_interpolation_style("tchan",CSM_INTERP_VERTICAL);
testhyp_pe->set_interpolation_style("tchan",CSM_INTERP_VERTICAL);

/*
  nullhyp->set_interpolation_style("tchan",CSM_INTERP_HORIZONTAL);
  nullhyp_pe->set_interpolation_style("tchan",CSM_INTERP_HORIZONTAL);
  testhyp->set_interpolation_style("tchan",CSM_INTERP_HORIZONTAL);
  testhyp_pe->set_interpolation_style("tchan",CSM_INTERP_HORIZONTAL);
*/

cout << "end of templates" << endl;


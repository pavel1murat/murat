//

#include "TFile.h"
#include "TTree.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <stdio.h>

void ts1_tooth(const char* Dir="622_0001") {

  char buf[1000];
  TObjArray  files;

  FILE* pipe = gSystem->OpenPipe(Form("ls /mu2e/data/users/murat/datasets/ts1_tooth/%s/root/",Dir),"r");
  while (fgets(buf,10000,pipe)) { 
    files.Add(new TObjString(buf)); 
  }
  gSystem->ClosePipe(pipe);

  TObjString *ostr;
  TString prefix;
  TIter itt(&files);


  int loc = 1;
  while ( (ostr = (TObjString*) itt.Next()) ) {
    char fn[100],full_fn[200];
    sscanf(ostr->String().Data(),"%s",fn);
    sprintf(full_fn,"/mu2e/data/users/murat/datasets/ts1_tooth/%s/root/%s",Dir,fn);
    //    printf("%s\n",full_fn);

    TFile* f = TFile::Open(full_fn);

    TTree* t = (TTree*) f->Get("//stepPointMCDumper/nt");

    TH1F* h = new TH1F("h","PDG code",500,-250,250);

    t->Project("h","hits.pdgId","");
    
    int nep = h->GetBinContent(238);
    int nmp = h->GetBinContent(240);
    int nem = h->GetBinContent(262);
    int nmm = h->GetBinContent(264);

    printf(" loc, ne+, nmu+, ne-, nmu- = %5i %5i %5i %5i %5i \n",loc,nep,nmp,nem,nmm);
    loc += 1;
  }

}

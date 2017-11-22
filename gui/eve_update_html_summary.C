/// \file
/// \ingroup tutorial_eve
/// Html table and event summary for alice_esd.C
///
/// \macro_code
///
/// \author Bertrand Bellenot

#include "TEveEventManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TEveGeoShape.h"
#include "TGHtml.h"
#include "aliesd/AliESDEvent.h"
#include "aliesd/AliESDRun.h"
#include "aliesd/AliESDtrack.h"
#include "aliesd/AliExternalTrackParam.h"
#include "TEveTrackPropagator.h"

#include "murat/gui/eve_HtmlSummary.hh"
#include "murat/gui/eve_HtmlObjTable.hh"

//==============================================================================

TGHtml      *fgHtml        = 0;

//______________________________________________________________________________
void update_html_summary() {
   // Update summary of current event.

   // TEveElement::List_i i;
   // TEveElement::List_i j;
   // Int_t k;
   // TEveElement *el;
   // HtmlObjTable *table;
   // TEveEventManager *mgr = gEve ? gEve->GetCurrentEvent() : 0;
   // if (mgr) {
   //    fgHtmlSummary->Clear("D");
   //    for (i=mgr->BeginChildren(); i!=mgr->EndChildren(); ++i) {
   //       el = ((TEveElement*)(*i));
   //       if (el->IsA() == TEvePointSet::Class()) {
   //          TEvePointSet *ps = (TEvePointSet *)el;
   //          TString ename  = ps->GetElementName();
   //          TString etitle = ps->GetElementTitle();
   //          if (ename.First('\'') != kNPOS)
   //             ename.Remove(ename.First('\''));
   //          etitle.Remove(0, 2);
   //          Int_t nel = atoi(etitle.Data());
   //          table = fgHtmlSummary->AddTable(ename, 0, nel);
   //       }
   //       else if (el->IsA() == TEveTrackList::Class()) {
   //          TEveTrackList *tracks = (TEveTrackList *)el;
   //          TString ename  = tracks->GetElementName();
   //          if (ename.First('\'') != kNPOS)
   //             ename.Remove(ename.First('\''));
   //          table = fgHtmlSummary->AddTable(ename.Data(), 5,
   //                   tracks->NumChildren(), kTRUE, "first");
   //          table->SetLabel(0, "Momentum");
   //          table->SetLabel(1, "P_t");
   //          table->SetLabel(2, "Phi");
   //          table->SetLabel(3, "Theta");
   //          table->SetLabel(4, "Eta");
   //          k=0;
   //          for (j=tracks->BeginChildren(); j!=tracks->EndChildren(); ++j) {
   //             Float_t p     = ((TEveTrack*)(*j))->GetMomentum().Mag();
   //             table->SetValue(0, k, p);
   //             Float_t pt    = ((TEveTrack*)(*j))->GetMomentum().Perp();
   //             table->SetValue(1, k, pt);
   //             Float_t phi   = ((TEveTrack*)(*j))->GetMomentum().Phi();
   //             table->SetValue(2, k, phi);
   //             Float_t theta = ((TEveTrack*)(*j))->GetMomentum().Theta();
   //             table->SetValue(3, k, theta);
   //             Float_t eta   = ((TEveTrack*)(*j))->GetMomentum().Eta();
   //             table->SetValue(4, k, eta);
   //             ++k;
   //          }
   //       }
   //    }
   //    fgHtmlSummary->Build();
   //    fgHtml->Clear();
   //    fgHtml->ParseText((char*)fgHtmlSummary->Html().Data());
   //    fgHtml->Layout();
   // }
}

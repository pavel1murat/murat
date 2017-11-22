#ifndef __murat_gui_eve_HtmlObjTable_hh__
#define __murat_gui_eve_HtmlObjTable_hh__

#include "TObject.h"
#include "TString.h"
#include "TArrayF.h"

////////////////////////////////////////////////////////////////////////////////
class HtmlObjTable : public TObject {
public:                     // make them public for shorter code

   TString   fName;
   Int_t     fNValues;      // number of values
   Int_t     fNFields;      // number of fields
   TArrayF  *fValues;
   TString  *fLabels;
   Bool_t    fExpand;

   TString   fHtml;         // HTML output code

   void Build();
   void BuildTitle();
   void BuildLabels();
   void BuildTable();

public:
   HtmlObjTable(const char *name, Int_t nfields, Int_t nvals, Bool_t exp=kTRUE);
   virtual ~HtmlObjTable();

   void     SetLabel(Int_t col, const char *label) { fLabels[col] = label; }
   void     SetValue(Int_t col, Int_t row, Float_t val) { fValues[col].SetAt(val, row); }
   TString  Html() const { return fHtml; }

   ClassDef(HtmlObjTable, 0);
};

#endif

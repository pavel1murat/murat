///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/eve_HtmlObjTable.hh"

ClassImp(HtmlObjTable)

//______________________________________________________________________________
HtmlObjTable::HtmlObjTable(const char *name, Int_t nfields, Int_t nvals, Bool_t exp) :
   fName(name), fNValues(nvals), fNFields(nfields), fExpand(exp)
{
   // Constructor.

   fValues = new TArrayF[fNFields];
   for (int i=0;i<fNFields;i++)
      fValues[i].Set(nvals);
   fLabels = new TString[fNFields];
}

//______________________________________________________________________________
HtmlObjTable::~HtmlObjTable()
{
   // Destructor.

   delete [] fValues;
   delete [] fLabels;
}

//______________________________________________________________________________
void HtmlObjTable::Build()
{
   // Build HTML code.

   fHtml = "<table width=100% border=1 cellspacing=0 cellpadding=0 bgcolor=f0f0f0> ",

   BuildTitle();
   if (fExpand && (fNFields > 0) && (fNValues > 0)) {
      BuildLabels();
      BuildTable();
   }

   fHtml += "</table>";
}

//______________________________________________________________________________
void HtmlObjTable::BuildTitle()
{
   // Build table title.

   fHtml += "<tr><td colspan=";
   fHtml += Form("%d>", fNFields+1);
   fHtml += "<table width=100% border=0 cellspacing=2 cellpadding=0 bgcolor=6e6ea0>";
   fHtml += "<tr><td align=left>";
   fHtml += "<font face=Verdana size=3 color=ffffff><b><i>";
   fHtml += fName;
   fHtml += "</i></b></font></td>";
   fHtml += "<td>";
   fHtml += "<td align=right> ";
   fHtml += "<font face=Verdana size=3 color=ffffff><b><i>";
   fHtml += Form("Size = %d", fNValues);
   fHtml += "</i></b></font></td></tr>";
   fHtml += "</table>";
   fHtml += "</td></tr>";
}

//______________________________________________________________________________
void HtmlObjTable::BuildLabels()
{
   // Build table labels.

   Int_t i;
   fHtml += "<tr bgcolor=c0c0ff>";
   fHtml += "<th> </th>"; // for the check boxes
   for (i=0;i<fNFields;i++) {
      fHtml += "<th> ";
      fHtml += fLabels[i];
      fHtml += " </th>"; // for the check boxes
   }
   fHtml += "</tr>";
}

//______________________________________________________________________________
void HtmlObjTable::BuildTable()
{
   // Build part of table with values.

   for (int i = 0; i < fNValues; i++) {
      if (i%2)
         fHtml += "<tr bgcolor=e0e0ff>";
      else
         fHtml += "<tr bgcolor=ffffff>";

      TString name = fName;
      name.ReplaceAll(" ", "_");
      // checkboxes
      fHtml += "<td bgcolor=d0d0ff align=\"center\">";
      fHtml += "<input type=\"checkbox\" name=\"";
      fHtml += name;
      fHtml += Form("[%d]\">",i);
      fHtml += "</td>";

      for (int j = 0; j < fNFields; j++) {
         fHtml += "<td width=";
         fHtml += Form("%d%%", 100/fNFields);
         fHtml += " align=\"center\"";
         fHtml += ">";
         fHtml += Form("%1.4f", fValues[j][i]);
         fHtml += "</td>";
      }
      fHtml += "</tr> ";
   }
}


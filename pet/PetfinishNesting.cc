//
// Free function to be used by the nest... functions
//
// $Id: PetfinishNesting.cc,v 1.2 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $
//
// Original author KLG based on nest... functions
//

// C++ includes; if using cout << *info.solid etc...
#include <iostream>
#include <iomanip>

// Mu2e includes
#include "murat/pet/PetfinishNesting.hh"
#include "murat/pet/PetVolumeInfo.hh"
#include "murat/pet/PetG4Helper.hh"

// G4 includes
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

// if using cout << *info.solid etc...
#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

using namespace std;

namespace mu2e {

  void PetfinishNesting(PetVolumeInfo& info,
                     G4Material* material,
                     G4RotationMatrix const* rot,
                     G4ThreeVector const & offset,
                     G4LogicalVolume* parent,
                     int copyNo,
                     bool const isVisible,
                     G4Colour const color,
                     bool const forceSolid,
                     bool const forceAuxEdgeVisible,
                     bool const placePV,
                     bool const doSurfaceCheck,
                     bool const verbose
                     ) {

    PetG4Helper    & _helper = *(art::ServiceHandle<PetG4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    // the code below if activated prints the parameters of the solid
    // being placed
    if (verbose) {
      ios::fmtflags oldfl = cout.flags();
      int const oldpr = cout.precision();
      int const oldwdth = cout.width();
      static int const newpr   = 18;
      static int const newwdth = 28;

      cout.setf(ios::right,ios::adjustfield);
      cout.precision(newpr);
      cout.width(newwdth);

      cout << "Placing " << info.name;
      if (parent!=0) {
        cout << " inside " << parent->GetName()<< endl;
      } else {
        cout << endl;
      }
      cout << scientific << setprecision(newpr) << setw(newwdth) << right << *info.solid << endl;
      cout << "Offset " << offset << endl;
      cout << "Rotation ";
      if (rot != 0) {
        rot->print(cout);
        cout << endl;
        cout << "Tolerance " << rot->getTolerance() << endl;
        cout << "Accessing the rotation matrix " << endl;
        int const vsize = rot->colX().SIZE;
        G4ThreeVector colX(rot->colX());
        G4ThreeVector colY(rot->colY());
        G4ThreeVector colZ(rot->colZ());
        for (int i=0;i!=vsize;++i) {
          cout << " " << setprecision(newpr) << setw(newwdth) << right << colX[i] <<
            " "<< setprecision(newpr) << setw(newwdth) << right << colY[i] <<
            " "<< setprecision(newpr) << setw(newwdth) << right << colZ[i] << endl;
        }
      } else {
        cout << endl << 0 << endl;
      }
      cout << setprecision(oldpr) << setw(oldwdth);
      cout.flags(oldfl);

    }

    info.logical  = new G4LogicalVolume( info.solid, material, info.name);

    // G4 did not get const-ness correctly, thus the const_cast
    info.physical  =  placePV ? new G4PVPlacement( const_cast<G4RotationMatrix*>(rot),
                                                   offset,
                                                   info.logical,
                                                   info.name,
                                                   parent,
                                                   0,
                                                   copyNo,
                                                   doSurfaceCheck)
                              : 0;

    // uncomment for a more thorrow overlap check
    // doSurfaceCheck && info.physical!=0 && info.physical->CheckOverlaps(100000,0.0,false);


    if (!isVisible) {

      info.logical->SetVisAttributes(G4VisAttributes::Invisible);

    } else {

      G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, color));
      visAtt->SetForceSolid(forceSolid);
      // If I do not do this, then the rendering depends on what happens in
      // other parts of the code;  is there a G4 bug that causes something to be
      // unitialized?
      visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      info.logical->SetVisAttributes(visAtt);

    }

    // Save the volume information in case someone else needs to access it by name.
    _helper.addVolInfo(info);

    return;

  }

}

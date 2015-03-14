#ifndef Mu2eG4_PetConstructMaterials_hh
#define Mu2eG4_PetConstructMaterials_hh
//
// Construct materials requested by the run-time configuration system.
//
// $Id: PetConstructMaterials.hh,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//

#include <string>
#include <memory>

// Forward references in global namespace.
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Element;

#include "G4String.hh"

namespace mu2e {

  // Return type of the isNeeded() method.
  class CheckedG4String{
  public:
    CheckedG4String ():
      doit(false),
      name(){
    }
    CheckedG4String ( bool b, G4String const& n):
      doit(b),
      name(n){
    }
    bool doit;
    G4String name;
  };

  // Forward references within mu2e namespace.
  class SimpleConfig;

  class PetConstructMaterials{
  public:

    PetConstructMaterials();
    ~PetConstructMaterials();

    // Construct all of the materials.
    void construct();

  private:

    // Construct materials specific to Mu2e.
    void constructMu2eMaterials( SimpleConfig const& config );

    // Wrapper around FindOrBuildElement.
    G4Element* getElementOrThrow( G4String const& name);

    // Do we need to build this material?
    CheckedG4String isNeeded(  std::vector<std::string> const& V, std::string const& name);

    // Check that a material name is not already in use.
    void uniqueMaterialOrThrow( G4String const& name);

    // Check to see if a string appears within a vector of strings.
    bool isRequested( std::vector<std::string> const& V, std::string const& name);

  };

} // end namespace mu2e
#endif /* Mu2eG4_PetConstructMaterials_hh */


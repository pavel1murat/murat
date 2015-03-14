//
// Geometry and identifier info about the Calorimeter.
//
// Original author B. Echenard
//

// C++ includes
#include <iostream>
#include <algorithm>

// Mu2e includes
#include "murat/pet/PetImager.hh"
// #include "CalorimeterGeom/inc/Disk.hh"

//other includes
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/TwoVector.h"


    int PetImager::caloSectionId(int crystalId) const
    {          
        for (unsigned int i=0;i<_nSections;++i) 
        {
           if (crystalId < section(i).nCrystals()) return i;
           crystalId -= section(i).nCrystals();
        }
        return _nSections;
    }

    int PetImager::localCrystalId(int crystalId) const
    {     
        for (unsigned int i=0;i<_nSections;++i) 
        {
          if (crystalId < section(i).nCrystals()) return crystalId;
          crystalId -= section(i).nCrystals();
        }
        return crystalId;
    }





    unsigned int PetImager::nRO(void) const 
    {
        unsigned total(0);
        for (unsigned int i=0;i<_nSections;++i) total += section(i).nCrystals();
        return total*_nROPerCrystal;
    }

    unsigned int PetImager::nCrystal(void) const 
    {
        unsigned total(0);
        for (unsigned int i=0;i<_nSections;++i) total += section(i).nCrystals();
        return total;
    }







    CLHEP::Hep3Vector PetImager::toCrystalFrame(int CrystalId, CLHEP::Hep3Vector const& pos) const 
    {   

       const mu2e::CaloSection& thisSection = section( caloSectionId(CrystalId) );
       int ic                         = localCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = thisSection.crystal(ic).position();
       crysLocalPos += _crystalShift;

       return (thisSection.rotation())*(pos-thisSection.origin())-crysLocalPos;  
    }


    CLHEP::Hep3Vector PetImager::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
       const mu2e::CaloSection& thisSection = section(sectionId);
       return (thisSection.rotation())*(pos-thisSection.origin());
    }


    CLHEP::Hep3Vector PetImager::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
        const mu2e::CaloSection& thisSection = section(sectionId);
        return thisSection.inverseRotation()*pos + thisSection.origin();
    }


    CLHEP::Hep3Vector PetImager::crystalOrigin(int CrystalId) const 
    {          
       const mu2e::CaloSection& thisSection = section( caloSectionId(CrystalId) );
       int ic                         = localCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = thisSection.crystal(ic).position();
       crysLocalPos += _crystalShift;

       return thisSection.origin() + thisSection.inverseRotation()*crysLocalPos; 
    }


    CLHEP::Hep3Vector PetImager::localCrystalOrigin(int CrystalId) const 
    {          
       const mu2e::CaloSection& thisSection = section( caloSectionId(CrystalId) );
       int ic                         = localCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = thisSection.crystal(ic).position();
       crysLocalPos += _crystalShift;
                   
       return crysLocalPos; 
    }


    CLHEP::Hep3Vector PetImager::crystalAxis(int CrystalId) const 
    {
       const mu2e::CaloSection& thisSection = section( caloSectionId(CrystalId) );
       CLHEP::Hep3Vector vlocal(0,0,1);
       return thisSection.inverseRotation()*vlocal;
    }

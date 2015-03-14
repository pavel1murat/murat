#ifndef CalorimeterGeom_BrainImager_hh
#define CalorimeterGeom_BrainImager_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: BrainImager.hh,v 1.5 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $
//
// Original author R. Bernstein and Rob Kutschke
//

//C++ includes
#include <vector>
#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "murat/pet/PetImager.hh"
#include "CalorimeterGeom/inc/Vane.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


using boost::static_pointer_cast;
using boost::shared_ptr;

class BrainImager: public PetImager {

      
  friend class BrainImagerMaker;


public:
  
  BrainImager(){}
  ~BrainImager(){}


  int                NWedges()   const  { return _nSections; }
  const mu2e::Vane*  Vane(int i) const  { return (const mu2e::Vane*) &section(i); }
  
  int    nCrystalR()            const  { return _nCrystalR; }
  int    nCrystalZ()            const  { return _nCrystalZ; }
  
  double innerRadius ()                 const  {return _rMin;}
  double outherRadius()                 const  {return _rMax;}
  
  bool   isInsideVane(int ivane, CLHEP::Hep3Vector const& pos) const ;
  

        


//keep only for backward compatibility, will disappear in the future.

int crystalByRO(int roid) const              {return (roid/_nROPerCrystal); }
int ROBaseByCrystal(int crystalId) const     {return (crystalId*_nROPerCrystal);}

int crystalVaneByRO(int roid) const { return (roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR);}
int crystalRByRO(int roid) const {return ((roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR))/_nCrystalZ;}
int crystalZByRO(int roid) const {return ((roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR))%_nCrystalZ;}

CLHEP::Hep3Vector crystalOriginByRO(int roid) const { return crystalOrigin(crystalByRO(roid));  }
int vaneByRO(int roid) const {return roid/(_nCrystalZ*_nCrystalR*_nROPerCrystal);}
CLHEP::Hep3Vector crystalAxisByRO(int roid) const { return crystalAxis(crystalByRO(roid));  }

CLHEP::Hep3Vector toVaneFrame(int vaneId, CLHEP::Hep3Vector const& pos) const   {return toSectionFrame(vaneId,pos);}
CLHEP::Hep3Vector fromVaneFrame(int vaneId, CLHEP::Hep3Vector const& pos) const {return fromSectionFrame(vaneId,pos);}

//-----------------------------------------------------------------------------
// overloaded methods
//-----------------------------------------------------------------------------
  virtual double           crystalHalfTrans()     const  {return _crystalHW; }
  virtual double           crystalHalfLength()    const  {return _crystalHL; }
  virtual double           crystalVolume()        const  {return 8*_crystalHW*_crystalHW*_crystalHL;}
  virtual bool             isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const ;        
  virtual int              crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const ;
  virtual double           crystalLongPos(int crystalId, CLHEP::Hep3Vector const& pos) const; 

  virtual std::vector<int> neighbors(int crystalId, int level=1) const;


      private:

          int    _nCrystalZ;
          int    _nCrystalR;
          double _rMin;
          double _rMax;
          double _crystalHL;
          double _crystalHW;        
          double _shieldHalfThickness;
          double _absorberHalfThickness;



   };
#endif /* CalorimeterGeom_BrainImager_hh */

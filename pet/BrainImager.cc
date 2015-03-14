//
// Geometry and identifier info about the BrainImager.
//
//
// $Id: BrainImager.cc,v 1.3 2013/10/15 23:41:14 murat Exp $
// $Author: murat $
// $Date: 2013/10/15 23:41:14 $
//
// Original author R. Bernstein and Rob Kutschke
//
//C++ includes
#include <algorithm>

//mu2e includes
#include "murat/pet/BrainImager.hh"



    bool BrainImager::isInsideVane(int ivane, CLHEP::Hep3Vector const& pos) const 
    {   
	CLHEP::Hep3Vector posInSection = toSectionFrame(ivane, pos);

	double xlim = _crystalHL + _wrapperThickness + _roHalfThickness + 0.5;
	double ylim = _nCrystalR*(_crystalHW+_wrapperThickness+_shellThickness) + 0.5;   
	double zlim = _nCrystalZ*(_crystalHW+_wrapperThickness+_shellThickness) + 0.5;   

	if (posInSection.x() < -xlim || posInSection.x() > xlim ) return false;      
	if (posInSection.y() < -ylim || posInSection.y() > ylim ) return false;      
	if (posInSection.z() < -zlim || posInSection.z() > zlim ) return false;      

	return true;
    }
  
    bool BrainImager::isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const 
    {   
        for (unsigned int ivane=0;ivane<_nSections;++ivane) if (isInsideVane(ivane,pos)) return true;
        return false;    
    }

   
    int BrainImager::crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const 
    {   
        int offset(0);
        for (unsigned int ivane=0;ivane<_nSections;++ivane) {
           if ( isInsideVane(ivane,pos) ) {
                 CLHEP::Hep3Vector posInSection = toSectionFrame(ivane, pos);
                 return offset + Vane(ivane)->idxFromPosition(posInSection.y(),posInSection.z());
  	   }
           offset +=Vane(ivane)->nCrystals();
        }       	  
        return -1;
    }

    std::vector<int> BrainImager::neighbors(int CrystalId, int level) const 
    {

	int iv = caloSectionId(CrystalId);
	int ic = localCrystalId(CrystalId);

        int offset(0);
        for (int i=0;i<iv;++i) offset += Vane(i)->nCrystals();

	std::vector<int> list = Vane(iv)->neighbors(ic,level);
	transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  
	return list;
    }

   double BrainImager::crystalLongPos(int crystalId, CLHEP::Hep3Vector const& pos) const
   {   
	CLHEP::Hep3Vector posInSection = toCrystalFrame(crystalId, pos);
	return posInSection.x();
   }

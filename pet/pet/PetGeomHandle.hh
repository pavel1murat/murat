#ifndef PetGeometryService_GeomHandle_hh
#define PetGeometryService_GeomHandle_hh

//
// A safe pointer to the geometry information for a
// detector component.
//
// $Id: PetGeomHandle.hh,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//

#include "murat/pet/PetGeometryService.hh"

namespace mu2e {
  template <typename DET>
  class PetGeomHandle
  {
  public:
    PetGeomHandle()
    {
      art::ServiceHandle<PetGeometryService> sg;
      _detector = sg->getElement<DET>();
    }
    ~PetGeomHandle() { }

    DET const * operator->() const { return _detector;}
    DET const * get()        const { return _detector;}
    DET const & operator*()  const { return *_detector;}
    DET const * operator->() { return _detector;}
    DET const & operator*()  { return *_detector;}

  private:
    PetGeomHandle(const PetGeomHandle&);
    PetGeomHandle& operator=(const PetGeomHandle&);

    // unnecessary
    DET* operator&();

    DET* _detector;
  };
}

#endif /* PetGeometryService_GeomHandle_hh */

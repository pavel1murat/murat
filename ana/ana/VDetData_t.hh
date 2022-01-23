#ifndef __murat_ana_VDetData_t__
#define __murat_ana_VDetData_t__
//------------------------------------------------------------------------------
//  IXLocal and IZLocal are either (1,3) or (3,1) 
//  vector normal to the VD is pointed along IZLocal
//------------------------------------------------------------------------------
namespace murat {
  struct VDetData_t {
    int    fID;
    int    fIXLocal;			// orientation of the local axis to calculate Pt wrt to the field
    int    fIZLocal;
    double fX      ;		        // for X dist
    double fZ      ;
    double fBField ;                    // magnitude of the B-field, assumed to be orthogonal to the detector plane
    double fPhiXZ  ;                    // angle of the local unit Z vector in the global system
  };
}
#endif

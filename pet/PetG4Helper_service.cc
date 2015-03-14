//
// PetG4Helper plugin.
//
// $Id: PetG4Helper_service.cc,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//

// Framework includes

// Mu2e includes
#include "murat/pet/PetG4Helper.hh"

using namespace std;

namespace mu2e {

  PetG4Helper::PetG4Helper(fhicl::ParameterSet const& iPS,
			   art::ActivityRegistry&iRegistry){
  }

  PetG4Helper::~PetG4Helper(){
  }

  // Return the volume info mapped to the given key, throw if the key does not exist.
  PetVolumeInfo& PetG4Helper::locateVolInfo( const std::string key){
    std::map<std::string,PetVolumeInfo>::iterator i = _volumeInfoList.find(key);
    if ( i == _volumeInfoList.end() ){
      throw cet::exception("GEOM")
        << "PetG4Helper::locateVolInfo cannot find the volume named: "
        << key
        << "\n";
    }
    return i->second;
  } // end of PetG4Helper::locateVolInfo

  // If the key already exists, throw. Otherwise add the (key, value) pair
  // to the map.
  void PetG4Helper::addVolInfo( const PetVolumeInfo& info ){
    std::map<std::string,PetVolumeInfo>::iterator i = _volumeInfoList.find(info.name);
    if ( i != _volumeInfoList.end() ){
      throw cet::exception("GEOM")
        << "locateVolInfo already has the key: "
        << info.name
        << "\n";
    }
    _volumeInfoList[info.name] = info;
  } // end of PetG4Helper::addVolInfo

} // end namespace mu2e

//using mu2e::PetG4Helper;
DEFINE_ART_SERVICE(mu2e::PetG4Helper);

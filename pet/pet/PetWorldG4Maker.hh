#ifndef PETWORLDG4MAKER_HH
#define PETWORLDG4MAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class PetWorldG4; }

namespace mu2e {
  class PetWorldG4Maker {
  public:
    static std::unique_ptr<PetWorldG4> make(const SimpleConfig& config);
  };
}

#endif/*WORLDG4MAKER_HH*/

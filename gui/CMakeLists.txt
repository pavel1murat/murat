#
SET (CMAKE_Fortran_COMPILER  gfortran)
ENABLE_LANGUAGE(Fortran)

set (subpkg "gui")
set (dict "${CMAKE_PROJECT_NAME}_${subpkg}_dict")
set (lib  "${CMAKE_PROJECT_NAME}_${subpkg}")

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

ROOT_GENERATE_DICTIONARY(${dict}
  TEvdCalorimeter.hh TEvdDisk.hh TEvdCrystal.hh
  TEvdCrv.hh TEvdCrvSector.hh
  TEvdTracker.hh TEvdStation.hh TEvdPlane.hh TEvdPanel.hh TEvdStraw.hh
  TEvdManagerA.hh TEvdSubdetector.hh
  LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  TEvdCalorimeter.cc TEvdDisk.cc TEvdCrystal.cc
  TEvdCrv.cc TEvdCrvSector.cc
  TEvdTracker.cc TEvdStation.cc TEvdPlane.cc TEvdPanel.cc TEvdStraw.cc
  TEvdManagerA.cc TEvdSubdetector.cc
  ${dict}.cxx
)

target_link_libraries(${lib}
  Stntuple_alg Stntuple_geom Stntuple_obj Stntuple_base
  ConfigTools
  GeometryService
  fhiclcpp  cetlib
  ROOT::EG ROOT::Core ROOT::Matrix ROOT::Tree ROOT::Graf ROOT::Geom ROOT::Gpad)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/${CMAKE_PROJECT_NAME})

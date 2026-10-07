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
  TEvdCaloNode.hh TEvdCaloView.hh
  TEvdCrv.hh TEvdCrvSector.hh TEvdCrvModule.hh TEvdCrvLayer.hh TEvdCrvCounter.hh
  TEvdDetNode.hh
  TEvdTracker.hh TEvdStation.hh TEvdPlane.hh TEvdPanel.hh TEvdStraw.hh
  TEvdManager.hh TEvdSubdetector.hh
  TEvdView.hh TEvdVisNode.hh
  LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  TEvdCalorimeter.cc TEvdDisk.cc TEvdCrystal.cc
  TEvdCaloNode.cc TEvdCaloView.cc
  TEvdCrv.cc TEvdCrvSector.cc  TEvdCrvModule.cc TEvdCrvLayer.cc TEvdCrvCounter.cc
  TEvdDetNode.cc
  TEvdTracker.cc TEvdStation.cc TEvdPlane.cc TEvdPanel.cc TEvdStraw.cc
  TEvdManager.cc TEvdSubdetector.cc
  TEvdView.cc TEvdVisNode.cc
  ${dict}.cxx
)

target_link_libraries(${lib}
  Stntuple_alg Stntuple_geom Stntuple_obj Stntuple_base
  ConfigTools
  GeometryService
  fhiclcpp  cetlib
  ROOT::EG ROOT::Core ROOT::Matrix ROOT::Tree ROOT::Graf ROOT::Geom ROOT::Gpad ROOT::RGL)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/${CMAKE_PROJECT_NAME})

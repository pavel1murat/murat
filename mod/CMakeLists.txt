#
set (subpkg "mod"                    )
set (dict   "murat_${subpkg}_dict")
set (lib    "murat_${subpkg}"     )

# set (MuHitDisplayModule "Stntuple_MuHitDisplay_module" )
# set (MuHitDisplayDict   "Stntuple_MuHitDisplay_module_dict" )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

# set (include_files
#   THistModule.hh TModule.hh
#   TAnaRint.hh
# )
# 
# ROOT_GENERATE_DICTIONARY(${MuHitDisplayDict} MuHitDisplay_module.hh
#   LINKDEF ${subpkg}/MuHitDisplay_module_Linkdef.h
#   OPTIONS -D__CLING__
# )
# set_property(SOURCE ${MuHitDisplayDict}.cxx PROPERTY generated)
# 
# ROOT_GENERATE_DICTIONARY(${dict} ${include_files}
#   LINKDEF ${subpkg}/Linkdef.h
#   OPTIONS -D__CLING__
# )
# set_property(SOURCE ${dict}.cxx PROPERTY generated)
# 
cet_build_plugin(DumpGeometry art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_gui
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  Offline::CRVReco
  Offline::CosmicRayShieldGeom
  Offline::CRVConditions
  Offline::DataProducts
  Offline::GeometryService
  Offline::RecoDataProducts
  artdaq-core-mu2e::Data
  art_root_io::TFileService_service
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
)


cet_add_compiler_flags(-g -O0)

install(TARGETS ${DumpGeometry})

# install_headers()
# install_source()

message( " --Stntuple/${subpkg}: CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})

install(DIRECTORY ${subpkg} DESTINATION include/murat)

#
SET (CMAKE_Fortran_COMPILER  gfortran)
ENABLE_LANGUAGE(Fortran)

set (subpkg "ana")
set (dict "${CMAKE_PROJECT_NAME}_${subpkg}_dict")
set (lib  "${CMAKE_PROJECT_NAME}_${subpkg}")

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

ROOT_GENERATE_DICTIONARY(${dict}
  CrvHist_t.hh lhr_rejection.hh mva_data.hh poisson_fluct.hh prob_dist.hh rpc_timing.hh
  TAnaModule.hh TBeamFlashAnaModule.hh TCalAnaModule.hh TCaloDose.hh TCaloTimeRes.hh
  TClusterAnaModule.hh TColl1DoseAnaModule.hh TColl3DoseAnaModule.hh TCosmicsAnaModule.hh
  TDegraderRpcAnaModule.hh TDioCalib.hh TDioCalibModule.hh TDoseAnaModule.hh TEmuAnaModule.hh
  TFilterModule.hh TFitDioRcSpectrum.hh TFunAnaModule.hh TG4ValidationModule.hh TGenAnaModule.hh
  THelixAnaModule.hh TMomscaleAnaModule.hh TMuonStopAnaModule.hh TPidAnaModule.hh TPipenuAnaModule.hh
  TPolFitSpectrum.hh TrackHist_t.hh TSimpAnaModule.hh TSpmcAnaModule.hh TStrawHitAnaModule.hh
  TTrackAnaModule.hh TTrackCompModule.hh TTrackerDose.hh TTrackPidAnaModule.hh TTrackRecoEffAnaModule.hh
  TTrackStrawHitAnaModule.hh TTriggerAnaModule.hh TTsMisalignment.hh TValCalPatRecModule.hh
  TVDetAnaModule.hh VDetData_t.hh
  LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  poisson_fluct.cc prob_dist.cc rpc_timing.cc lhr_rejection.cc mva_data.cc TValCalPatRecModule.cc
  TVDetAnaModule.cc TTsMisalignment.cc TTrackerDose.cc TTriggerAnaModule.cc TTrackCompModule.cc
  TTrackPidAnaModule.cc TTrackRecoEffAnaModule.cc TTrackStrawHitAnaModule.cc 
  TTrackAnaModule.cc TSpmcAnaModule.cc TStrawHitAnaModule.cc TPolFitSpectrum.cc TSimpAnaModule.cc
  TMomscaleAnaModule.cc TMuonStopAnaModule.cc TPidAnaModule.cc TPipenuAnaModule.cc TG4ValidationModule.cc
  TGenAnaModule.cc THelixAnaModule.cc TEmuAnaModule.cc TFilterModule.cc TFitDioRcSpectrum.cc
  TFunAnaModule.cc TDegraderRpcAnaModule.cc TDioCalib.cc TDioCalibModule.cc TDoseAnaModule.cc
  TColl1DoseAnaModule.cc TColl3DoseAnaModule.cc TCosmicsAnaModule.cc TCalAnaModule.cc TCaloDose.cc
  TCaloTimeRes.cc TClusterAnaModule.cc TBeamFlashAnaModule.cc SConscript TAnaModule.cc
  InitVirtualDetectors.cc
  ${dict}.cxx
)

target_link_libraries(${lib}
  Stntuple_val Stntuple_loop Stntuple_alg Stntuple_geom Stntuple_obj Stntuple_base
  ROOT::TMVA ROOT::EG ROOT::Core ROOT::Matrix ROOT::Tree ROOT::Graf ROOT::Gpad)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/${CMAKE_PROJECT_NAME})

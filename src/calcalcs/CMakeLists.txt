add_library(pismcalcalcs calcalcs.c utCalendar2_cal.c)
target_link_libraries (pismcalcalcs ${UDUNITS2_LIBRARIES})
set_target_properties(pismcalcalcs PROPERTIES COMPILE_FLAGS "-w")

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismcalcalcs DESTINATION ${Pism_LIB_DIR})
endif()

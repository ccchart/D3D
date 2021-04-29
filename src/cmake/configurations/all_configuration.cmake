# project(all)

# Specify the modules to be included

include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/dflowfm_configuration.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/d_waq_configuration.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/d_waves_configuration.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/flow2d3d_configuration.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/dimr_configuration.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/tools_configuration.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/configurations/tests_configuration.cmake)

project(all)

# Include directories that are required as references
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Set the directory of where the source code is located
set(src_path src) 

# Gather the source files
# set(headers include/swanout1_deltares.inc)

get_fortran_source_files(${src_path} icepack_sourcefiles)
list(FILTER icepack_sourcefiles EXCLUDE REGEX ".*version.txt$")


# Define library
set(library_name icepack)
set(library_source_list     ${headers}
                            ${icepack_sourcefiles})
add_library(${library_name} ${library_source_list})

# Add dependencies
#if (UNIX)
#    # the `pkg_check_modules` function is created with this call
#    find_package(PkgConfig REQUIRED)
#
#    # these calls create special `PkgConfig::<MODULE>` variables
#    pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)
#    pkg_check_modules(NETCDF_FTN REQUIRED IMPORTED_TARGET netcdf-fortran)
#
#    target_link_libraries(${library_name}
#         PkgConfig::NETCDF
#         PkgConfig::NETCDF_FTN)
#
#
#    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
#endif(UNIX)

#if (WIN32)
#    set(oss_dependencies    netcdf4
#                            netcdff)
#    oss_include_libraries(${library_name} ${oss_dependencies})
#    target_link_libraries(${library_name} ${oss_dependencies})
#    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include"
#                        ${mpich2_path})
#endif(WIN32)


# Define additional compilation flags
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_aerosol.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_age.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_algae.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_atmo.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_brine.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_firstyear.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_flux.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_fsd.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_intfc.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_isotope.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_itd.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_kinds.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_mechred.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_meltpond_cesm.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_meltpond_lvl.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_meltpond_topo.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_mushy_physics.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_ocean.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_orbital.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_parameters.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_shortwave.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_snow.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_therm_0layer.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_therm_bl99.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_therm_itd.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_therm_mushy.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_therm_shared.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_therm_vertical.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_tracers.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_warnings.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_wavefracspec.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_zbgc.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_zbgc_shared.F90
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/icepack_zsalinity.F90
                            PROPERTIES COMPILE_OPTIONS "${file_preprocessor_flag}")

# Set additional compilation properties, specific for Debug mode
target_compile_options(${library_name} PRIVATE "$<$<CONFIG:Debug>:${check_bounds_flag}>")

# Define preprocessor definitions 
#if (WIN32)
#    target_compile_definitions(${library_name} PRIVATE "WIN32;USE_MPI")
#endif(WIN32)

#if (UNIX)
#    target_compile_definitions(${library_name} PRIVATE "USE_MPI")
#endif(UNIX)

# Create the folder structure in vfproj
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${library_source_list})
set_target_properties (${library_name} PROPERTIES FOLDER third_party_open/icepack_lib)

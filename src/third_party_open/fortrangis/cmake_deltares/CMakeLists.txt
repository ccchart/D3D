# Set the directory which contains all the source files
set(src_path "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Define the fortranc library
set(fortranc_path ${src_path}/libfortranc)

set(fortranc_sources ${fortranc_path}/fortranc.F90
                        ${fortranc_path}/iso_varying_string.F90)
                  
if (UNIX)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake.in ${fortranc_path}/config.h)  
    set(configFile ${fortranc_path}/config.h) 
endif(UNIX)

set(fortranc_library_name fortranc)
add_library(${fortranc_library_name} ${fortranc_sources})

# Set additional compilation properties, specific for Debug mode
target_compile_options(${fortranc_library_name} PRIVATE "$<$<CONFIG:Debug>:${check_pointers_flag}>")
target_compile_options(${fortranc_library_name} PRIVATE "$<$<CONFIG:Debug>:${floating_point_exception_flag}>")

# Create the folder structure in project
source_group("${visual_studio_source_group_name}" FILES ${fortranc_sources})

# Create the folder structure in vfproj
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/.." FILES ${fortranc_sources})
set_target_properties (${fortranc_library_name} PROPERTIES FOLDER third_party_open/fortrangis)


# Define the fortrangis library
set(fortrangis_path ${src_path}/libfortrangis)
set(fortrangis_sources  ${fortrangis_path}/gdal.F90
                        ${fortrangis_path}/proj.F90
                        ${fortrangis_path}/shapelib.F90)

set(fortrangis_library_name fortrangis)
add_library(${fortrangis_library_name} ${fortrangis_sources})

# Include directories that are necessary
target_include_directories(${fortrangis_library_name} PRIVATE ${src_path}/include)

# Add dependencies
target_link_libraries(${fortrangis_library_name} PRIVATE ${fortranc_library_name})

# Add dependencies of underlying 3rd party libraries
# Use target_link_libraries with INTERFACE to propagate all these dependencies
# to other package that use fortrangis.

# shp intentionally comes first here: otherwise on Windows, the precompiled gdal_i.lib
# will introduce a duplicate definition of Shp* API functions.
target_link_libraries(${fortrangis_library_name} INTERFACE shp)

if (UNIX)
    # the `pkg_check_modules` function is created with this call
    find_package(PkgConfig REQUIRED)

    # these calls create special `PkgConfig::<MODULE>` variables
    pkg_check_modules(PROJ REQUIRED IMPORTED_TARGET proj)
    pkg_check_modules(GDAL REQUIRED IMPORTED_TARGET gdal)

    target_link_libraries(${fortrangis_library_name} INTERFACE PkgConfig::PROJ PkgConfig::GDAL)
endif(UNIX)
if (WIN32)
    target_link_libraries(${fortrangis_library_name} INTERFACE proj)

    # GDAL comes pre-compiled on Windows
    target_link_libraries(${fortrangis_library_name} INTERFACE "gdal_i.lib")
    target_link_directories(${fortrangis_library_name} INTERFACE ${gisinternals_path})
    target_include_directories(${fortrangis_library_name} INTERFACE ${gisinternals_path})
endif(WIN32)

# Propagate preprocessor variable definitions to any other package that uses fortrangis.
target_compile_definitions(${fortrangis_library_name} INTERFACE HAVE_PROJ)
target_compile_definitions(${fortrangis_library_name} INTERFACE HAVE_SHAPELIB)
target_compile_definitions(${fortrangis_library_name} INTERFACE HAVE_GDAL)

# Set additional compilation properties, specific for Debug mode
target_compile_options(${fortrangis_library_name} PRIVATE "$<$<CONFIG:Debug>:${check_pointers_flag}>")
target_compile_options(${fortrangis_library_name} PRIVATE "$<$<CONFIG:Debug>:${floating_point_exception_flag}>")

# Create the folder structure in vfproj
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/.." FILES ${fortran_sources})
set_target_properties (${fortrangis_library_name} PROPERTIES FOLDER third_party_open/fortrangis)


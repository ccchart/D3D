# Set directory which contains all source files
set(src_path src)

if (UNIX)
    #name mangling, these functions are generated by the FortranCInterface cmake package and might be system dependent
    include(FortranCInterface)
    FortranCInterface_HEADER("${CMAKE_CURRENT_SOURCE_DIR}/include/FortranCInterface.h" MACRO_NAMESPACE "FC_")
    #produce a config.h to be included
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/include/config.h)  
    set(configFile "${CMAKE_CURRENT_SOURCE_DIR}/include/config.h") 
    
    include_directories(${NETCDF_INCLUDEDIR}) 
    add_definitions(-DD3D_PLUGIN_EXT=".so")
    add_definitions(-DHAVE_STRUCT_TIMESPEC)
endif(UNIX)

# Gather source files
file(GLOB source_files  ${src_path}/*.c)

# Define executable name
set(executable_name duprol2delwaq)
add_executable(${executable_name} ${source_files})

# Create the folder structure in vfproj
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${source_files})
set_target_properties (${executable_name} PROPERTIES FOLDER tools_gpl/duprol2delwaq)

# post-build
set(install_dir ${CMAKE_BINARY_DIR})
set(build_dir ${CMAKE_BINARY_DIR})

message(STATUS "build_config ${build_config}")

post_build_target( ${executable_name}
                   ${install_dir} 
                   ${build_dir} 
                   ${checkout_src_root} 
                   ${executable_name}  )

install(TARGETS ${executable_name} RUNTIME DESTINATION bin)

# Navigate back to the wave/manager component to retrieve the rc file for setting file properties
set(rc_version_file ${CMAKE_CURRENT_SOURCE_DIR}/../manager/src/wave_version.rc)

# Gather source files
set(library_files ${src_path}/wave_bmi.f90) # Because the .dll and the .exe are defined in the same directory, retrieve the relevant files for the library alone

# Define library
set(library_name wave)
add_library(${library_name} SHARED  ${library_files}
                                    ${rc_version_file})

# Set additional compilation properties
target_compile_options(${library_name} PRIVATE "${extend_source132_flag}")

# Set dependencies
set(library_dependencies    data
                            delftio
                            delftio_shm
                            deltares_common
                            deltares_common_c
                            deltares_common_mpi
                            ec_module
                            gridgeom
                            io
                            io_netcdf
                            kdtree_wrapper
                            kdtree
                            kernel
                            manager
                            nefis
                            netcdf4
                            netcdff
                            triangle_c) 
oss_include_libraries(${library_name} library_dependencies)
target_link_libraries(${library_name} ${library_dependencies})

if (WIN32)
    # Set linker properties
    message(STATUS "Setting linker properties in windows")
    target_link_directories(${library_name}
                            PRIVATE
                            "${checkout_src_root}/third_party_open/netcdf/netCDF 4.6.1/lib"
                            "${checkout_src_root}/third_party_open/pthreads/bin/x64")

    target_link_libraries(${library_name}                                                   
                            "pthreadVC2.lib"
                            "netcdf.lib")

    # Set linker options
    message(STATUS "Setting target_link_options in windows")
    target_link_options(${library_name} PRIVATE ${nologo_flag})
endif(WIN32)

if(UNIX)
    target_link_libraries(${library_name} ${library_dependencies})
endif(UNIX)

# Define how the files should be structured within Visual Studio
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${library_files})
source_group(Resources FILES    ${rc_version_file})
set_target_properties (${library_name} PROPERTIES FOLDER engines_gpl/wave)

# Change the name of the target library to wave.dll
set_target_properties (${library_name} PROPERTIES OUTPUT_NAME wave)

# Set post-build step
set(install_dir ${CMAKE_BINARY_DIR})
set(build_dir ${CMAKE_BINARY_DIR})

post_build_target (${library_name}
                   ${install_dir} 
                   ${build_dir} 
                   ${checkout_src_root} 
                   ${library_name})

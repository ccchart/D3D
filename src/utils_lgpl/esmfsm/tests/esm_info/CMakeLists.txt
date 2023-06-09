set(executable_name esm_info)
set(src_path .)
set(esmfsm_include_path ${CMAKE_CURRENT_SOURCE_DIR}/../../include) # Trace back to the esmfsm path
set(esm_include_path ${CMAKE_CURRENT_SOURCE_DIR}/../../packages/esmfsm_c/include)

if (UNIX)
    #name mangling, these functions are generated by the FortranCInterface cmake package and might be system dependent
    include(FortranCInterface)
    FortranCInterface_HEADER("${CMAKE_CURRENT_SOURCE_DIR}/include/FortranCInterface.h" MACRO_NAMESPACE "FC_")
    #produce a config.h to be included
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/include/config.h)  
    set(configFile "${CMAKE_CURRENT_SOURCE_DIR}/include/config.h") 
    
    add_definitions(-DHAVE_CONFIG_H)
    add_definitions(-DD3D_PLUGIN_EXT=".so")
    add_definitions(-DHAVE_STRUCT_TIMESPEC)
endif(UNIX)

# Gather source files
file(GLOB esm_global_header_files ${esmfsm_include_path}/*.h)
file(GLOB esm_header_files ${esm_include_path}/*.h)
file(GLOB source_files ${src_path}/*.c)

# Define executable
add_executable(${executable_name} ${source_files}
                                  ${esmfsm_include_path}
                                  ${esm_include_path})

# Set dependencies
set(oss_dependencies    deltares_common
                        esmfsm
                        esmfsm_c
                        esmfsm_version_number)
oss_include_libraries(${executable_name} oss_dependencies)
if(UNIX)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(${executable_name} Threads::Threads)
    target_include_directories(${executable_name} 
                               PRIVATE  
                               ${esmfsm_include_path}
                               ${esm_include_path})
endif(UNIX)
if (WIN32)
    target_include_directories(${executable_name} 
                               PRIVATE  
                               ${esmfsm_include_path}
                               ${esm_include_path}
                               "${checkout_src_root}/third_party_open/pthreads/include/x64")
    target_link_directories(${executable_name}
                            PRIVATE
                            "${checkout_src_root}/third_party_open/pthreads/bin/x64")
    target_link_libraries(${executable_name}    ${oss_dependencies}
                                                "pthreadVC2.lib")


    # Set linker options
    message(STATUS "Setting target_link_options in windows")
    target_link_options(${executable_name} PRIVATE ${nologo_flag})

    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        target_link_options(${executable_name} PRIVATE "/LIBPATH: $(IFORT_COMPILER21)/compiler/lib/intel64;")
    endif()
endif(WIN32)
# Set preprocessor definitions
if(WIN32)
    target_compile_definitions(${executable_name} PRIVATE "WIN32")
endif(WIN32)
if(UNIX)
    target_compile_definitions(${executable_name} PRIVATE "HAVE_CONFIG_H")
endif(UNIX)

# Explicitly link the libraries to the target or the compilation will fail as the linker cannot find the Fortran libs
target_link_libraries(${executable_name} ${oss_dependencies})

# Define how the files should be structured within Visual Studio
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${source_files}
                                                    ${esmfsm_include_path}
                                                    ${esm_include_path})
set_target_properties (${executable_name} PROPERTIES FOLDER utils_lgpl/esmfsm/test)

# Not needed on Windows, only on Linux
install(TARGETS ${executable_name} RUNTIME  DESTINATION bin)

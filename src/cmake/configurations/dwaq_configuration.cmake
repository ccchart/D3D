#
# WAQ
#=============
add_subdirectory(${checkout_src_root}/${waq_plugin_wasteload_module} waq_plugin_wasteload)
add_subdirectory(${checkout_src_root}/${waq_utils_c_module} waq_utils_c)
add_subdirectory(${checkout_src_root}/${waq_utils_f_module} waq_utils_f)
add_subdirectory(${checkout_src_root}/${waq_process_module} waq_process)
add_subdirectory(${checkout_src_root}/${waq_kernel_module} waq_kernel)
add_subdirectory(${checkout_src_root}/${waq_io_module} waq_io)
add_subdirectory(${checkout_src_root}/${delwaq_lib_module} delwaq_lib)
add_subdirectory(${checkout_src_root}/${waq_data_module} waq_data)
add_subdirectory(${checkout_src_root}/${delwaq1_module} delwaq1)
add_subdirectory(${checkout_src_root}/${delwaq2_module} delwaq2)
add_subdirectory(${checkout_src_root}/${delwaq_lib_examples_module} delwaq_lib_examples)
add_subdirectory(${checkout_src_root}/${waq_delftio_module} waq_delftio)
add_subdirectory(${checkout_src_root}/${wq_processes_module} wq_processes)

#
# WAQ Tools
#=============
add_subdirectory(${checkout_src_root}/${waqpb_export_module} waqpb_export)
add_subdirectory(${checkout_src_root}/${waqpb_import_module} waqpb_import)
add_subdirectory(${checkout_src_root}/${waqpb_lib_module} waqpb_lib)
add_subdirectory(${checkout_src_root}/${waqmerge_module} waqmerge)
add_subdirectory(${checkout_src_root}/${ddcouple_module} ddcouple)
add_subdirectory(${checkout_src_root}/${agrhyd_module} agrhyd)
add_subdirectory(${checkout_src_root}/${maptonetcdf_module} maptonetcdf)



#
# PART
#=============
add_subdirectory(${checkout_src_root}/${part_data_f_module} part_data_f)
add_subdirectory(${checkout_src_root}/${part_utils_f_module} part_utils_f)
add_subdirectory(${checkout_src_root}/${part_io_f_module} part_io_f)
add_subdirectory(${checkout_src_root}/${part_kernel_f_module} part_kernel_f)
add_subdirectory(${checkout_src_root}/${delpar_module} delpar)


#
# Third party libraries
#=============
add_subdirectory(${checkout_src_root}/${kdtree_module} kdtree2)
add_subdirectory(${checkout_src_root}/${kdtree_wrapper_module} kdtree_wrapper)
add_subdirectory(${checkout_src_root}/${triangle_c_module} triangle_c)
add_subdirectory(${checkout_src_root}/${fortrangis_module} fortrangis)
add_subdirectory(${checkout_src_root}/${shp_module} shp)


if(WIN32)
    add_subdirectory(${checkout_src_root}/${proj_module} proj)
endif(WIN32)

#
# Utils
#=============
# Deltares_common
add_subdirectory(${checkout_src_root}/${deltares_common_module} deltares_common)
add_subdirectory(${checkout_src_root}/${deltares_common_c_module} deltares_common_c)


# netcdf
if(WIN32)
    add_subdirectory(${checkout_src_root}/${netcdf_module} netcdff)
endif()

# io_netcdf
add_subdirectory(${checkout_src_root}/${io_netcdf_module} io_netcdf)

# ec_module
add_subdirectory(${checkout_src_root}/${ec_module} ec_module)

# gridgeom
add_subdirectory(${checkout_src_root}/${gridgeom_module} gridgeom)

# Nefis
add_subdirectory(${checkout_src_root}/${nefis_module} nefis)

# Solvesaphe
add_subdirectory(${checkout_src_root}/${solvesaphe_module} solvesaphe)

# io_hyd
add_subdirectory(${checkout_src_root}/${io_hyd_module} io_hyd)

#
# Linux installation
#=============
if(UNIX)
    add_subdirectory(${checkout_src_root}/${install_waq_module} install_waq)
endif()

# Project name must be at the end of the configuration: it might get a name when including other configurations and needs to overwrite that
project(dwaq)

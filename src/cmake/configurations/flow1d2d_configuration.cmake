project(flow1d2d)

# flow1d2d
# ============
add_subdirectory(${checkout_src_root}/${flow1d2d_module} flow1d2d)
add_subdirectory(${checkout_src_root}/${flow1d2d_api_access_module} flow1d_flowfm_api_access)


# Utils
# =====

# Utils LGPL
# =====
if(NOT TARGET deltares_common)
    add_subdirectory(${checkout_src_root}/${deltares_common_module} deltares_common)
endif()

if(NOT TARGET deltares_common_c)
    add_subdirectory(${checkout_src_root}/${deltares_common_c_module} deltares_common_c)
endif()


# Third party
# ===========

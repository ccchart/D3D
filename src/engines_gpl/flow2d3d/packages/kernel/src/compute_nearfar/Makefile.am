#-------------------------------------------------------------------------------
#   Makefile for delftflow
#
#   Irv.Elshoff@deltares.nl
#   Adri.Mourits@deltares.nl
#   Fedor.Baart@deltares.nl
#   11 jan 2007
#
#   Copyright (C)  Stichting Deltares, 2011-2013.
#-------------------------------------------------------------------------------

include $(top_srcdir)/common.am
noinst_LTLIBRARIES = libKernelComputeNearfar.la

#-------------------------------------------------------------------------------
#Includes
FINCLUDES += \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/include

FMODULES += \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/src/basics \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/src/parallel_mpi \ 
	-I$(top_builddir)/utils_lgpl/inifiles/packages/inifiles/src \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common_mpi/src \
	-I$(top_builddir)/utils_lgpl/constants/packages/mathconsts/src \
	-I$(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src

#-------------------------------------------------------------------------------
#Sources
libKernelComputeNearfar_la_SOURCES = \
         findnmk.f90 \
         desa.f90 \
         dis_qjet.f90 \
         interp_tk.f90 \
         discha_nf.f90 \
         corinp_gen.f90 \
         cormix2flow.f90 \
         jet3d2flow.f90 \
         near_field.f90 \
         wri_cormix.f90 \
         wri_jet3d.f90 \
         corjet2flow.f90 \
         cortim2flow.f90 \
         cortim_no_modules.f90 \
         cortim_xyzs.f90 \
         coupled.f90 \
         ct2jettraj.f90 \
         det_num_dis.f90 \
         init_nfl.f90 \
         wait_a_sec.f90 \
         wri_cortim.f90 \
         wri_tek.f90
         
#-------------------------------------------------------------------------------
#Libraries


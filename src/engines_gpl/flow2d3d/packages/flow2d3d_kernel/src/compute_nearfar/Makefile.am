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
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/flow2d3d_data/include

FMODULES += \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/flow2d3d_data/src/basics \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src \
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
         wri_jet3d.f90

#-------------------------------------------------------------------------------
#Libraries


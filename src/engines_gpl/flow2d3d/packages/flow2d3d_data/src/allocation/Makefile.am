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
noinst_LTLIBRARIES = libDataAllocation.la

#-------------------------------------------------------------------------------
#Includes
FINCLUDES += \
	-I$(top_builddir)/utils_lgpl/esmfsm/include

FMODULES += \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/flow2d3d_data/src/basics \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/flow2d3d_data/src/parallel_mpi \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common_mpi/src \
	-I$(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src

#-------------------------------------------------------------------------------
#Sources
libDataAllocation_la_SOURCES = \
         chnull.f90 \
         decarr.f90 \
         dnull.f90 \
         errpnt.f90 \
         esm_alloc_char.f90 \
         esm_alloc_int.f90 \
         esm_alloc_real.f90 \
         getpointer.f90 \
         gtcpnt.f90 \
         gtdpnt.f90 \
         gtipnt.f90 \
         gtlpnt.f90 \
         gtptrs.f90 \
         gtrpnt.f90 \
         inull.f90 \
         lnull.f90 \
         mkcpnt.f90 \
         mkdpnt.f90 \
         mkfpnt.f90 \
         mkipnt.f90 \
         mklpnt.f90 \
         mkrpnt.f90 \
         rnull.f90

#-------------------------------------------------------------------------------
#Libraries


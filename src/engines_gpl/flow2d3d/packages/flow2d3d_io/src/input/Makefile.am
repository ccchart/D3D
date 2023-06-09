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
noinst_LTLIBRARIES = libIOInput.la

#-------------------------------------------------------------------------------
#Includes
AM_FCFLAGS += $(NETCDF_CFLAGS) $(NETCDF_FORTRAN_CFLAGS)

FINCLUDES += \
	-I$(top_builddir)/utils_lgpl/esmfsm/include \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/flow2d3d_data/include

FMODULES += \
	-I$(top_builddir)/utils_gpl/morphology/packages/morphology_io/src \
	-I$(top_builddir)/utils_gpl/morphology/packages/morphology_kernel/src \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/flow2d3d_data/src/basics \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/flow2d3d_data/src/parallel_mpi \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/include \
	-I$(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src/meteo \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common_mpi/src \
	-I$(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src
#-------------------------------------------------------------------------------
#Sources
libIOInput_la_SOURCES = \
         rdarray.f90 \
         restart_trim_lyrs.f90 \
         split_stations.f90 \
         barfil.f90 \
         bndfil.f90 \
         bubfil.f90 \
         ciofil.f90 \
         ck_dpopt.f90 \
         cmpnum.f90 \
         dimbch.f90 \
         dimbnd.f90 \
         dimbub.f90 \
         dimdis.f90 \
         dimfou.f90 \
         dimiwe.f90 \
         dimmud.f90 \
         dimpro.f90 \
         dimrd.f90 \
         dimsedconst.f90 \
         dimsit.f90 \
         dimstr.f90 \
         dimtrt.f90 \
         drofil.f90 \
         dryfil.f90 \
         edyfil.f90 \
         expblk.f90 \
         fbcorr.f90 \
         restart_trim_flow.f90 \
         flsfil.f90 \
         grdfil.f90 \
         hybfil.f90 \
         icfil.f90 \
         area4.f90 \
         inicut.f90 \
         inimet.f90 \
         inimorlyr.f90 \
         morbndfill.f90 \
         numcmp.f90 \
         rdbedformpar.f90 \
         rdbcg.f90 \
         rdbndd.f90 \
         rdcul.f90 \
         rddept.f90 \
         rdveg3d.f90 \
         rddredge_d3d4.f90 \
         rdfcio.f90 \
         rdfour.f90 \
         rdgrid.f90 \
         rdhyb.f90 \
         rdhyg.f90 \
         rdhyvd.f90 \
         rdibch.f90 \
         rdic.f90 \
         rdirt.f90 \
         rdiwe.f90 \
         rdlaydis.f90 \
         rdmassbal.f90 \
         rdmeteo.f90 \
         rdmud.f90 \
         rdnum.f90 \
         rdprfl.f90 \
         rdproc.f90 \
         rdq2eb.f90 \
         rdrgf.f90 \
         rdroll.f90 \
         rdrund.f90 \
         rdscour.f90 \
         rdsedmortra.f90 \
         rdsite.f90 \
         rdspec.f90 \
         rdstru.f90 \
         rdtimo.f90 \
         rdtrt.f90 \
         rdttar.f90 \
         rdusrc.f90 \
         rdusrf.f90 \
         rdusrp.f90 \
         rdwaqpar.f90 \
         rdxyzo.f90 \
         restart_trim_roller.f90 \
         reacmp.f90 \
         readmd.f90 \
         read_triana.f90 \
         reafou.f90 \
         restart_trim_bdf.f90 \
         restart_trim_fluff.f90 \
         rstfil.f90 \
         srcfil.f90 \
         stafil.f90 \
         strfil.f90 \
         tdfil.f90 \
         trafil.f90 \
         search2.f90

#-------------------------------------------------------------------------------
#Libraries


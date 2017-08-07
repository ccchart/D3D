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
noinst_LTLIBRARIES = libPLICcutCELLS.la

#-------------------------------------------------------------------------------
#Includes
FINCLUDES += \
	-I$(top_builddir)/utils_lgpl/esmfsm/include \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/include
#
FMODULES += \
	-I$(top_builddir)/utils_gpl/morphology/packages/morphology_data/src \
	-I$(top_builddir)/utils_gpl/morphology/packages/morphology_io/src \
	-I$(top_builddir)/utils_gpl/morphology/packages/morphology_kernel/src \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/src/basics \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/src/parallel_mpi \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/include \
	-I$(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src/meteo \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common_mpi/src \
	-I$(top_builddir)/utils_lgpl/delftio/packages/delftio/src/diof90 \
	-I$(top_builddir)/utils_lgpl/delftio/packages/delftio/src/delftio_sync \
	-I$(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src
#-------------------------------------------------------------------------------
#Sources
         libPLICcutCELLS_la_SOURCES = \
         ALPHAvof.f90 \
         AREApoly.f90 \
         A_G_Poly.f90 \
         BC_VOF.f90 \
         BEDslopeCORR.f90 \
         BOUvelFROMdisch.f90 \
         BOUvelFROMdisch2.f90 \
         BOUvelINT.f90 \
         COMPUTEmergingCARATT.f90 \
         COMPUTEtau.f90 \
         COPYFIRSTinSECOND.f90 \
         Compute_inSTENCILuv.f90 \
         DISSvelGHOST.f90 \
         DISSvelGHOST_uzd.f90 \
         DP0isDP.f90 \
         Diff_segm.f90 \
         ERROR_DGETRF.F90 \
         ERROR_DGETRS.F90 \
         EXACTslopeSUB.f90 \
         FORCEdischarge.f90 \
         FRESHqx.f90 \
         FindGhostPoints.f90 \
         Find_BI_PI.f90 \
         INIT_kfuv.f90 \
         INTERPvATuPOINT.f90 \
         IScurvilinear.f90 \
         InterpG_ATs1LOCATION.f90 \
         InterpG_ATu1LOCATION.f90 \
         InterpG_ATv1LOCATION.f90 \
         InterpGhost_s_u_v.f90 \
         IntersSemilineSegm_always.f90 \
         MERGEconc.f90 \
         MERGEneuCARATT.f90 \
         ModCbind.f90 \
         My_intersec.f90 \
         My_intersec_circle.f90 \
         OFFonPERvel.f90 \
         OFFonPERvel_ONLYkf1.f90 \
         PER_GRIDgeometry.f90 \
         PER_NpsiNeta.f90 \
         PER_dp.f90 \
         PER_dpLH.f90 \
         PER_dps.f90 \
         PER_dps_anal.f90 \
         PER_kfscc.f90 \
         PER_poros.f90 \
         PLIC_VOF_INIT.F90 \
         PLIC_VOFdelftSUBR.F90 \
         PNPOLY.f90 \
         POROSupdate.f90 \
         PRINTtestGHOSTaccur.f90 \
         READpoly.f90 \
         RECOMPqkz.f90 \
         REDUCEgsqs.f90 \
         ROTATE.f90 \
         ROTATEback.f90 \
         S1bound_fromCONT_OBSOLETE.f90 \
         SHIFT_xycorU1.f90 \
         SORT4.f90 \
         SUBRsmoothEb.f90 \
         TRANSPOSE_wrapper.F90 \
         U1bound_fromCONT_OBSOLETE.f90 \
         UPDATEbank.f90 \
         UPDATEbankdaCANC.f90 \
         UorVmean_ghosts.f90 \
         VIRTUALlinkAREAS.f90 \
         VIRTUALlinkAVER.f90 \
         VIRTUALlinkAVER_s1.f90 \
         VIRTUALlinkAVER_s1bis.f90 \
         VIRTUALlinkAVER_vert.f90 \
         VIRTUALlinkAVER_volum1.f90 \
         VIRTUALlinkCONC.f90 \
         VIRTUALlinkVERT.f90 \
         VIRTUALlink_gsqsR.f90 \
         W1prescribed.f90 \
         W1qzkZERO.f90 \
         WATERlevelPERIOD.f90 \
         adjacent.f90 \
         allocateWORKARRAYS.f90 \
         bedLOADperiod.f90 \
         caldpu_cc.f90 \
         checkDRY.f90 \
         checkQbnd.f90 \
         computeU1V1accuratezza.f90 \
         compute_kwduv.f90 \
         comvol_onlyVOL.f90 \
         curvat_bis.f90 \
         cutcell_pre_sud_stage1.f90 \
         cutcell_pre_sud_stage2.f90 \
         cutcell_pre_uzd_stage1.f90 \
         cutcell_pre_uzd_stage2.f90 \
         deactivateBANKS.f90 \
         deferred_GRADs1.F90 \
         defineINTEXTper.f90 \
         drychk_cc.f90 \
         extrapPERIODICedgeBANK1.f90 \
         extrapPERIODICedgeBANK2.f90 \
         extrapPERIODICghost.f90 \
         extrapU.F90 \
         extrapUV.F90 \
         extrapV.F90 \
         forceGHOSTsud.f90 \
         forceGHOSTuzd.f90 \
         forcePERvelSUD.f90 \
         forcePERvelUZD.f90 \
         fullyBANKtoFRESH.f90 \
         ghosttype_is_kf.f90 \
         gsqsR_is_gsqs.f90 \
         guu0gvv0.f90 \
         huFRESH.f90 \
         iniFLUXuv.f90 \
         inpolygon.f90 \
         intCELLS.f90 \
         interpG_ATu1LOCATION_Hart.f90 \
         interpG_ATu1LOCATION_exactFREE.f90 \
         interpG_ATv1LOCATION_exactFREE.f90 \
         interpUinV_FROMu1STENCIL_exactFREE.f90 \
         interpVinU_FROMv1STENCIL_exactFREE.f90 \
         kfsuv_ghost.f90 \
         kfuv0_ghost_sud.f90 \
         kfuv1_ghost_sud.f90 \
         kfuvGHOST3.F90 \
         kfuv_ghost_smallcut.f90 \
         locateMN_s1.f90 \
         locateNM_U1.f90 \
         locateNM_V1.f90 \
         mergSMALLbound.f90 \
         min_dist.f90 \
         perCELLvar2D.f90 \
         perVISCdiff.f90 \
         periodSPIRAL.f90 \
         periodicGRADs1.f90 \
         periodic_Eb.f90 \
         periodic_dryWET_poly.f90 \
         periodic_hv.f90 \
         periodic_NxNy.f90 \
         periodic_u0INTv.f90 \
         periodic_uATvPOINTS.f90 \
         periodic_vATuPOINTS.f90 \
         periodic_xGLyGL.f90 \
         periodic_xGyG_U.f90 \
         periodic_xGyG_V.f90 \
         postpr_ghost.f90 \
         recompGHOSTiterUZD.f90 \
         reconVOF.f90 \
         rectPOLinRECTstenc.f90 \
         resetV1toV0subr.f90 \
         s0iss1_ghost.f90 \
         s1iss0_ghost.f90 \
         slopePERIOD.f90 \
         trapz.F90 \
         u0isu1_ghost.f90 \
         u1INTv_is_u0INTv.f90 \
         u1isu0.f90 \
         u1isu0_ghost.f90 \
         u1isuFULL.f90 \
         updateQXK.F90 \
         update_qfilt.f90 \
         uvelocityPERIODsud.f90 \
         v0isv1_ghost.f90 \
         v1isv0_ghost.f90 \
         v1isvFULL.f90 \
         velocFORcurvPERIOD.f90 \
         velocityPERIOD.f90 \
         velocityPERIOD_ghost.f90 \
         virtMERG.f90 \
         virtMERG_neumann.f90 \
         virtMERG_neumann_wrong.f90 \
         virtMERGdisch.f90 \
         volumNL_are_volumes.f90 \
         vvvORD2sub.F90 \
         xGL_dpL.f90 \
         z_aguuCOMP.f90 \
         z_cutgeo.f90 \
         z_defineKGP.f90


#-------------------------------------------------------------------------------
#Libraries


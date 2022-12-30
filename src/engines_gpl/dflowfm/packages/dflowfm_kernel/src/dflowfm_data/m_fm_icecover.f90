!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2021.                                
!                                                                               
!  This file is part of Delft3D (D-Flow Flexible Mesh component).               
!                                                                               
!  Delft3D is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  Delft3D  is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.             
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D",                  
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting 
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!                                                                               
!-------------------------------------------------------------------------------

! $Id$
! $HeadURL$

module m_fm_icecover
use precision
use icecover_module
use icecover_input_module
implicit none

!
! Global data
!
type(icecover_type), target                                :: ice_data                     !< module instance of the icecover data structure
!
real(fp), dimension(:), pointer                            :: ice_af                       !< module pointer to array ice areafrac inside ice_data
real(fp), dimension(:), pointer                            :: ice_h                        !< module pointer to array ice thickness inside ice_data
real(fp), dimension(:), pointer                            :: ice_p                        !< module pointer to array pressure inside ice_data
real(fp), dimension(:), pointer                            :: ice_t                        !< module pointer to array temperature inside ice_data
real(fp), dimension(:), pointer                            :: qh_air2ice                   !< module pointer to array qh_air2ice inside ice_data
real(fp), dimension(:), pointer                            :: qh_ice2wat                   !< module pointer to array qh_ice2wat inside ice_data
real(fp), dimension(:), pointer                            :: snow_h                       !< module pointer to array snow thickness inside ice_data
real(fp), dimension(:), pointer                            :: snow_t                       !< module pointer to array snow temperature inside ice_data

integer, pointer                                           :: ja_aice_read                 !< flag indicating whether ice area fraction is available via EC module
integer, pointer                                           :: ja_hice_read                 !< flag indicating whether ice thickness is available via EC module

logical, pointer                                           :: ice_hisout                   !< module pointer to flag hisout inside ice_data
logical, pointer                                           :: ice_mapout                   !< module pointer to flag mapout inside ice_data

logical, pointer                                           :: ice_apply_pressure           !< module pointer to flag apply_pressure inside ice_data
logical, pointer                                           :: ice_apply_friction           !< module pointer to flag apply_friction inside ice_data
logical, pointer                                           :: ice_reduce_surface_fluxes    !< module pointer to flag reduce_surface_fluxes inside ice_data
logical, pointer                                           :: ice_reduce_waves             !< module pointer to flag reduce_waves inside ice_data
integer, pointer                                           :: ice_modify_winddrag          !< module pointer to flag modify_winddrag inside ice_data

integer, pointer                                           :: ja_icecover                  !< module pointer to modeltype flag inside ice_data that specifies the ice cover model
integer, pointer                                           :: ice_frctp                    !< module pointer to frict_type inside ice_data

real(fp), pointer                                          :: ice_dens                     !< module pointer to ice_dens inside ice_data
real(fp), pointer                                          :: ice_albedo                   !< module pointer to ice_albedo inside ice_data
real(fp), pointer                                          :: ice_conduc                   !< module pointer to ice_conduc inside ice_data
real(fp), pointer                                          :: ice_lh                       !< module pointer to ice_lh inside ice_data
real(fp), pointer                                          :: ice_frcuni                   !< module pointer to frict_val inside ice_data

real(fp), pointer                                          :: snow_albedo                  !< module pointer to snow_albedo inside ice_data
real(fp), pointer                                          :: snow_conduc                  !< module pointer to snow_conduc inside ice_data
real(fp), pointer                                          :: snow_lh                      !< module pointer to snow_lh inside ice_data


contains


!> Nullify/initialize ice data structure.
subroutine fm_ice_null()
!!--declarations----------------------------------------------------------------
    !
    implicit none
    !
    ! Function/routine arguments
    !
    ! NONE
    !
    ! Local variables
    !
    integer                                                    :: istat     !> status flag for allocation
!
!! executable statements -------------------------------------------------------
!
    istat = null_icecover(ice_data)
    call fm_ice_update_all_pointers()
end subroutine fm_ice_null


!> Update all ice data structure.
subroutine fm_ice_update_all_pointers()
!!--declarations----------------------------------------------------------------
    !
    implicit none
    !
    ! Function/routine arguments
    !
    ! NONE
    !
    ! Local variables
    !
    ! NONE
!
!! executable statements -------------------------------------------------------
!
    ja_aice_read => ice_data%areafrac_forcing_available
    ja_hice_read => ice_data%thick_ice_forcing_available

    ja_icecover => ice_data%modeltype
   
    ice_hisout => ice_data%hisout
    ice_mapout => ice_data%mapout
    
    ice_apply_pressure  => ice_data%apply_pressure
    ice_apply_friction  => ice_data%apply_friction
    ice_reduce_waves    => ice_data%reduce_waves
    ice_modify_winddrag => ice_data%modify_winddrag
    
    ice_albedo          => ice_data%ice_albedo
    ice_conduc          => ice_data%ice_conductivity
    ice_lh              => ice_data%ice_latentheat
    ice_dens            => ice_data%ice_dens
    ice_frctp           => ice_data%frict_type
    ice_frcuni          => ice_data%frict_val
    
    snow_albedo         => ice_data%snow_albedo
    snow_conduc         => ice_data%snow_conductivity
    snow_lh             => ice_data%snow_latentheat
    
    call fm_ice_update_spatial_pointers()
end subroutine fm_ice_update_all_pointers


!> Update spatial pointers after (de)allocation
subroutine fm_ice_update_spatial_pointers()
!!--declarations----------------------------------------------------------------
    !
    implicit none
    !
    ! Function/routine arguments
    !
    ! NONE
    !
    ! Local variables
    !
!
!! executable statements -------------------------------------------------------
!
    ice_af => ice_data%areafrac
    ice_h  => ice_data%thick_ice
    ice_t  => ice_data%temp_ice
    ice_p  => ice_data%pressure
    qh_air2ice => ice_data%qh_air2ice
    qh_ice2wat => ice_data%qh_ice2wat
    snow_h => ice_data%thick_snow
    snow_t => ice_data%temp_snow
end subroutine fm_ice_update_spatial_pointers


!> activation of icecover module based on external forcing input
subroutine fm_ice_activate_by_ext_forces(ndx)
!!--declarations----------------------------------------------------------------
    !
    implicit none
    !
    ! Function/routine arguments
    !
    integer                                    , intent(in)    :: ndx       !> number of cells in the D-Flow FM domain
    !
    ! Local variables
    !
    integer                                                    :: istat     !> status flag for allocation
!
!! executable statements -------------------------------------------------------
!
    istat = late_activation_ext_force_icecover(ice_data)
    call fm_ice_alloc(ndx)
end subroutine fm_ice_activate_by_ext_forces


!> Allocate the arrays of ice data structure.
subroutine fm_ice_alloc(ndx)
!!--declarations----------------------------------------------------------------
    !
    implicit none
    !
    ! Function/routine arguments
    !
    integer                                    , intent(in)    :: ndx       !> number of cells in the D-Flow FM domain
    !
    ! Local variables
    !
    integer                                                    :: istat     !> status flag for allocation
!
!! executable statements -------------------------------------------------------
!
    if (associated(ice_af)) return ! don't allocate if already allocated - or should we deallocate and realloc?
    
    istat = alloc_icecover(ice_data, 1, ndx)
    call fm_ice_update_spatial_pointers()
end subroutine fm_ice_alloc


!> Clear the arrays of ice data structure.
subroutine fm_ice_clr()
!!--declarations----------------------------------------------------------------
    !
    implicit none
    !
    ! Function/routine arguments
    !
    ! NONE
    !
    ! Local variables
    !
    integer                                                    :: istat     !> status flag for allocation
!
!! executable statements -------------------------------------------------------
!
    istat = clr_icecover(ice_data)
    call fm_ice_null()
end subroutine fm_ice_clr


!> Read the ice cover module configuration from the mdu file
subroutine fm_ice_read(md_ptr, ierror)
!!--declarations----------------------------------------------------------------
    use dfm_error, only: DFM_WRONGINPUT
    use properties, only: tree_data
    implicit none
    !
    ! Function/routine arguments
    !
    type(tree_data)                            , pointer       :: md_ptr   !> pointer to the input file
    integer                                    , intent(inout) :: ierror   !> D-Flow FM error flag
    !
    ! Local variables
    !
    logical                                                    :: error    !> ice module error flag
!
!! executable statements -------------------------------------------------------
!
    call rd_icecover(ice_data, md_ptr, 'ice',  error)
    call fm_ice_update_spatial_pointers()
    !
    if (error) then
        ierror = DFM_WRONGINPUT
    endif
end subroutine fm_ice_read


!> Report the ice configuration to the diagnostic output.
subroutine fm_ice_echo(mdia)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    integer                                    , intent(in)    :: mdia     !> unit number of diagnostic output
    !
    ! Local variables
    !
    logical                                                    :: error    !> ice module error flag
!
!! executable statements -------------------------------------------------------
!
    error = echo_icecover(ice_data, mdia)
end subroutine fm_ice_echo


!> Update the ice pressure array.
subroutine fm_ice_update_press(ag)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    real(fp)                                   , intent(in)    :: ag       !> gravitational accelaration (m/s2)
    !
    ! Local variables
    !
    ! NONE
!
!! executable statements -------------------------------------------------------
!
    call update_icepress(ice_data, ag)
end subroutine fm_ice_update_press


!> preprocessing for ice cover, because in subroutine HEATUN some ice cover quantities have to be computed
!> this subroutine is comparable with subroutine HEA_ICE.F90 of the Delft3D-FLOW ice module
subroutine preprocess_icecover(n, Qlong_ice, tempwat, wind, timhr)
!!--declarations----------------------------------------------------------------
    use MessageHandling
    use m_flow                         ! test om tair(.) te gebruiken
    use m_flowgeom   , only: ndx
    use m_flowtimes  , only: dts
    use m_physcoef   , only: vonkar
    use m_heatfluxes , only: cpw
    use unstruc_files, only: mdia
    implicit none
    !
    ! Function/routine arguments
    !
    integer                                    , intent (in)   :: n             !> node number
    double precision                           , intent(in)    :: Qlong_ice     !> xx
    double precision                           , intent(in)    :: tempwat       !> xx
    double precision                           , intent(in)    :: wind          !> xx
    double precision                           , intent(in)    :: timhr         !> xx
    !
    ! Local variables
    !
    integer          :: iter
    double precision :: b, p_r, p_rt, kin_vis, t_freeze
    double precision :: b_t, c_tz, tm
    double precision :: conduc, D_t, D_ice, tsi, coef1, coef2, alpha
    logical          :: converged
    double precision :: z00, ustar,hdz, rhow, Qlong  
    
!
!! executable statements -------------------------------------------------------
!
    ! Initialization
    b = 3.0_fp             ! empirical constant in computation of C_tz
    p_r  = 13.0            ! molecular Prandtl number
    p_rt  = 0.85_fp        ! turbulent Prandtl number
    kin_vis = 0.0000018_fp ! kinematic viscosity of sea water
    t_freeze = 0.0_fp      ! freezing temperature of sea water
    rhow     = 1000.0_fp   ! density of water
    z00      = 2e-4_fp     ! Open sea roughness heigth 
    hdz      = 1.0_fp      ! rough estimate
    ustar = 0.025 * wind   ! See Eq. (12.5) ustar = sqrt(C_D) * U_10
    
    select case (ja_icecover)
    case (ICECOVER_KNMI)
        ! follow De Bruin & Wessels (1975)
        ! not implemented yet
        
    case (ICECOVER_SEMTNER)
        ! follow Semtner (1975)
        !
        ! Compute conductivity, depending oo presence of both ice and snow 
        ! 
        if ( snow_h(n) < 0.001_fp ) then
           conduc = ice_conduc 
           D_ice = max (0.01, ice_h(n))
           tsi = ice_t(n)
        else
           conduc = (ice_conduc * snow_conduc)
           D_ice = ( max (0.01_fp, ice_h(n)) * snow_conduc + max (0.01_fp, snow_h(n)) * ice_conduc)
           tsi = snow_t(n)
        endif    
        !
        ! Compute longwave radiation flux from ice surface according to Eq. (7) in (Wang, 2005)
        ! including an iteration proces
        !
        do iter =1,5
           coef1 = Qlong_ice * (tsi + 273.15_fp)**4.0_fp
           coef2 = 4.0_fp * Qlong_ice * (tsi + 273.15_fp)**3.0_fp
           D_t = (qh_air2ice(n) - coef1 - conduc * tsi / D_ice) / (coef2 + conduc / D_ice)
           tsi = tsi + D_t    
           if (abs(D_t) .lt. 1e-2 ) then
              converged = .true.
              Qlong = coef1 + coef2 * D_t
              if (Qlong .lt. 0.0_fp) then
                  Qlong = 0.0_fp
              endif
              if (tsi .gt. 0.0_fp) then
                 !
                 ! in case of melting recompute qbl
                 !
                 Qlong = Qlong_ice * (tsi + 273.15_fp)**4.0_fp
              endif
              !
              ! apply relaxation for stability reasons
              !
              alpha = 0.5_fp
              if (snow_h(n) < 0.001 ) then
                 ice_t(n) = alpha * tsi + (1.0_fp - alpha) * ice_t(n)
              else
                 snow_t(n) = alpha * tsi + (1.0_fp - alpha) * snow_t(n)
              endif
              !
              ! limit ice and snow temperature
              !
              if (ice_h(n) > 0.001) then
                  ice_t(n)  = min (0.0_fp,  ice_t(n))
                  ice_t(n)  = max (-25.0_fp, ice_t(n))
              endif    
              if (snow_h(n) > 0.001) then
                  snow_t(n) = min (0.0_fp,  snow_t(n))
                  snow_t(n) = max (-25.0_fp, snow_t(n))
              endif    
              !
              qh_air2ice(n) = qh_air2ice(n) - Qlong
              !
              ! no freezing in case of air temperatures above zero
              !
              if (tair(n) .gt. 0.0_fp .and. qh_air2ice(n) .lt. 0.0_fp) then
                 qh_air2ice(n) = 0.0_fp
              endif 
              !
              ! no melting in case of air temperatures below zero
              !
              if (tair(n) .lt. 0.0_fp .and. qh_air2ice(n) .gt. 0.0_fp) then
                 qh_air2ice(n) = 0.0_fp
              endif 
              goto 123  ! jump out of the iteration proces         
           endif
        enddo
123     continue
        !
        if (.not. converged) then
            !! write (lundia,*) 'Ice iteration not converged for NM =',nm
        endif    
        !
        ! Compute ice to water flux according to Wang (2015)
        !
        ! Calculate the molecular sublayer correction b_t
        !
        b_t  = b * sqrt(z00 * ustar / kin_vis ) * (p_r)**0.666
        !
        ! Calculate heat transfer coefficient c_tz
        !
        c_tz = ustar / ( b_t + p_rt * log (hdz/z00) / vonkar )
        !
        ! Calculate heat flux out of the ocean
        !
        qh_ice2wat(n) = rhow * cpw * c_tz * ( t_freeze - max(0.01_fp,tempwat) ) 
        !
        ! extra output for ice testbasin
        ! if (n==25) then
            ! write (msgbuf,'(a,i5,10f10.3)') 'ice fluxes:',n,timhr/24,qh_air2ice(n), qh_ice2wat(n), ice_h(n), ice_t(n), tempwat, tair(n), tsi; call msg_flush()
        ! endif
        !
        ! adaptation if QH_ICE2WAT conform KNMI approach (QH_ICE2WAT = 2.4 W/m2) 
        !
        !! qh_ice2wat(n) = -2.4_fp 
        !
        if ( isnan(qh_ice2wat(n)) ) then
           write (msgbuf,'(a,i5,10f10.3)') 'NAN in PREPROCESS_ICECOVER',n,qh_ice2wat(n); call msg_flush()
        endif
        ! 
    case default
        ! no preparation needed
    end select
end subroutine preprocess_icecover



!> update the ice cover -- initial coding here with full access to D-Flow FM arrays via use statements
!> let's see if we can make it gradually more modular and move functionality to the icecover_module.
subroutine update_icecover()
!!--declarations----------------------------------------------------------------
    use m_flowgeom   , only: ndx
    use m_flowtimes  , only: dts
    use m_wind       , only: tair
    use unstruc_files, only: mdia
    implicit none
    !
    ! Function/routine arguments
    !
    ! NONE
    !
    ! Local variables
    !
    integer          :: n
    double precision :: ice_growth, ice_melt, snow_melt
!
!! executable statements -------------------------------------------------------
!
  
    select case (ja_icecover)
    case (ICECOVER_KNMI)
        ! follow De Bruin & Wessels (1975)
        
    case (ICECOVER_SEMTNER)
        ! follow Semtner (1975)
 
        do n = 1, ndx
           if (tair(n) < 0.0_fp .or. ice_h(n) > 0.001 ) then
               if (qh_air2ice(n) > 0.0_fp) then
                   if ( snow_h(n) < 0.001 ) then
                       ! melting of ice
                       !
                       ice_melt = dts * ( 0.0_fp - qh_air2ice(n) ) / ice_lh
                       ice_h(n) = ice_h(n) + dts * ( -qh_air2ice(n) + qh_ice2wat(n) ) / ice_lh
                   else
                       ! melting of snow
                       !
                       snow_melt  = dts / snow_lh * ( 0.0_fp - qh_air2ice(n) )
                       snow_h(n) = snow_h(n) + (dts / snow_lh) * ( 0.0_fp - qh_air2ice(n) )
                   endif    
               else
                   ! freezing of ice
                   !
                   ice_growth = (dts /ice_lh) * ( -qh_air2ice(n) + qh_ice2wat(n) ) 
                   ice_h(n) = ice_h(n) + ice_growth
                   !
                   ! Maximize ice growth to 10 m
                   !
                   ice_h(n) = min(10.0_fp, ice_h(n)) 
               endif
           endif
        enddo
        
    case default
        ! by default no growth
    end select
end subroutine update_icecover


!> determine effective drag coefficient when ice may be present
pure function fm_ice_drag_effect(ice_af, cdw) result (cdeff)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    real(fp)                                   , intent(in)    :: ice_af   !> ice area fraction (-)
    real(fp)                                   , intent(in)    :: cdw      !> wind drag exerted via open water
    real(fp)                                                   :: cdeff    !> effective wind drag coefficient
    !
    ! Local variables
    !
    ! NONE
!
!! executable statements -------------------------------------------------------
!
    cdeff = ice_drag_effect(ice_data, ice_af, cdw)
end function fm_ice_drag_effect

end module m_fm_icecover

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
real(fp), dimension(:), pointer                            :: qh_air2ice                   !< module pointer to array qh_air2ice inside ice_data
real(fp), dimension(:), pointer                            :: qh_ice2wat                   !< module pointer to array qh_ice2wat inside ice_data
real(fp), dimension(:), pointer                            :: snow_h                       !< module pointer to array snow thickness inside ice_data

integer, pointer                                           :: ja_aice_read                 !< flag indicating whether ice area fraction is available via EC module
integer, pointer                                           :: ja_hice_read                 !< flag indicating whether ice thickness is available via EC module

logical, pointer                                           :: ice_hisout                   !< module pointer to flag hisout inside ice_data
logical, pointer                                           :: ice_mapout                   !< module pointer to flag mapout inside ice_data

logical, pointer                                           :: ice_apply_pressure           !< module pointer to flag apply_pressure inside ice_data
logical, pointer                                           :: ice_apply_friction           !< module pointer to flag apply_friction inside ice_data
logical, pointer                                           :: ice_reduce_surface_fluxes    !< module pointer to flag reduce_surface_fluxes inside ice_data
logical, pointer                                           :: ice_reduce_waves             !< module pointer to flag reduce_waves inside ice_data
logical, pointer                                           :: ice_reduce_wind              !< module pointer to flag reduce_wind inside ice_data

integer, pointer                                           :: ja_icecover                  !< module pointer to modeltype flag inside ice_data that specifies the ice cover model
integer, pointer                                           :: ice_frctp                    !< module pointer to frict_type inside ice_data

real(fp), pointer                                          :: ice_dens                     !< module pointer to ice_dens inside ice_data
real(fp), pointer                                          :: ice_albedo                   !< module pointer to ice_albedo inside ice_data
real(fp), pointer                                          :: ice_frcuni                   !< module pointer to frict_val inside ice_data

real(fp), pointer                                          :: snow_albedo                  !< module pointer to snow_albedo inside ice_data


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
    
    ice_apply_pressure => ice_data%apply_pressure
    ice_apply_friction => ice_data%apply_friction
    ice_reduce_waves   => ice_data%reduce_waves
    ice_reduce_wind    => ice_data%reduce_wind
    
    ice_albedo => ice_data%ice_albedo
    ice_dens   => ice_data%ice_dens
    ice_frctp  => ice_data%frict_type
    ice_frcuni => ice_data%frict_val
    
    snow_albedo => ice_data%snow_albedo
    
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
    ice_p  => ice_data%pressure
    qh_air2ice => ice_data%qh_air2ice
    qh_ice2wat => ice_data%qh_ice2wat
    snow_h => ice_data%thick_snow
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


!> update the ice cover -- initial coding here with full access to D-Flow FM arrays via use statements
!> let's see if we can make it gradually more modular and move functionality to the icecover_module.
subroutine update_icecover()
!!--declarations----------------------------------------------------------------
    use m_flowgeom , only: ndx
    use m_flowtimes, only: dts
    implicit none
    !
    ! Function/routine arguments
    !
    ! NONE
    !
    ! Local variables
    !
    integer :: k
!
!! executable statements -------------------------------------------------------
!
    select case (ja_icecover)
    case (ICECOVER_KNMI)
        ! follow De Bruin & Wessels (1975)
        
    case (ICECOVER_SEMTNER)
        ! follow Semtner (1975)
        do k = 1, ndx
            if (qh_air2ice(k) > 0.0_fp) then
                ! melting
            else
                ! freezing
            endif
            ice_h(k) = ice_h(k) + 0.1_fp * dts ! not sure that dts is actually known before the time step ...
        enddo
        
    case default
        ! by default no growth
    end select
end subroutine update_icecover

end module m_fm_icecover

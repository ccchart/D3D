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

module icecover_module
use precision
use MessageHandling, only: mess, LEVEL_ALL, LEVEL_FATAL
private

!
! public data types
!
public icecover_type

! public parameters
!
integer, parameter, public :: ICECOVER_NONE    = 0 !> no ice cover
integer, parameter, public :: ICECOVER_EXT     = 1 !> externally forced ice cover --> EC module, or BMI?
integer, parameter, public :: ICECOVER_KNMI    = 2 !> ice thickness computed based on De Bruin & Wessels
integer, parameter, public :: ICECOVER_SEMTNER = 3 !> ice thickness computed based on Semtner (1975)
! ... add IcePack?

integer, parameter, public :: FRICT_AS_DRAG_COEFF = 11 ! should be extension of D-Flow FM friction numbers

!
! public routines
!
public null_icecover
public select_icecover_model
public late_activation_ext_force_icecover
public alloc_icecover
public clr_icecover
!public update_icecover
public update_icepress

! ice cover type
type icecover_type
    !
    ! input
    !
    logical  :: hisout                            !> flag indicating whether ice cover should be written to history-file
    logical  :: mapout                            !> flag indicating whether ice cover should be written to map-file
    !
    logical  :: apply_pressure                    !> flag indicating whether pressure of ice cover should be applied
    logical  :: apply_friction                    !> flag indicating whether ice cover friction should be applied
    logical  :: reduce_surface_exchange           !> flag indicating whether precipitation, evaporation and heat exchange should be reduced
    logical  :: reduce_waves                      !> flag indicating whether waves should be reduced
    logical  :: reduce_wind                       !> flag indicating whether wind should be reduced
    !
    integer  :: modeltype                         !> type of the ice cover (one of ICECOVER_...)
    integer  :: frict_type                        !> friction type excerted by the ice cover
    !
    integer  :: areafrac_forcing_available        !> flag indicating whether ice area fraction is available via external forcing
    integer  :: thick_ice_forcing_available       !> flag indicating whether ice thickness is available via external forcing
    !
    real(fp) :: ice_albedo                        !> albedo of ice
    real(fp) :: snow_albedo                       !> albedo of snow
    real(fp) :: ice_dens                          !> ice density
    real(fp) :: frict_val                         !> friction coefficient of ice cover (unit depends on frict_type)
    !
    ! state
    !
    real(fp), dimension(:), pointer :: areafrac   => null() !> area fraction covered by ice (-)
    real(fp), dimension(:), pointer :: thick_ice  => null() !> ice cover thickness (m)
    real(fp), dimension(:), pointer :: thick_snow => null() !> snow cover thickness (m)
    !
    ! extra
    !
    real(fp), dimension(:), pointer :: qh_air2ice => null() !> heat flux from air to ice (?)
    real(fp), dimension(:), pointer :: qh_ice2wat => null() !> heat flux from ice to water (?)
    real(fp), dimension(:), pointer :: pressure   => null() !> pressure exerted by the ice cover (Pa)
end type icecover_type

contains

!> Nullify/initialize an icecover data structure.
function null_icecover(icecover) result(istat)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover  !> data structure containing ice cover data
    integer                                                    :: istat     !> status flag for allocation
    !
    ! Local variables
    !
    ! None
!
!! executable statements -------------------------------------------------------
!
    istat = select_icecover_model(icecover, ICECOVER_NONE)
    !
    ! state
    !
    nullify(icecover%areafrac)
    nullify(icecover%thick_ice)
    nullify(icecover%thick_snow)
    !
    ! extra
    !
    nullify(icecover%qh_air2ice)
    nullify(icecover%qh_ice2wat)
    nullify(icecover%pressure)
end function null_icecover


!> activation of icecover module based on external forcing input
function late_activation_ext_force_icecover(icecover) result(istat)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover  !> data structure containing ice cover data
    integer                                                    :: istat     !> status flag for allocation
    !
    ! Local variables
    !
    ! None
!
!! executable statements -------------------------------------------------------
!
    istat = 0
    if (icecover%modeltype == ICECOVER_EXT) then
       ! icecover already set to externally forced
    elseif (icecover%modeltype == ICECOVER_NONE) then
       ! activate icecover and switch on the pressure effect
       icecover%modeltype = ICECOVER_EXT
       icecover%apply_pressure = .true.
       call mess(LEVEL_ALL, 'Activating ice cover module based on external forcing.')
       ! note: spatial arrays haven't been allocated yet!
    else
       ! don't overrule previously selected icecover ...
       call mess(LEVEL_FATAL, 'Ice cover forcing data conflicts with selected ice cover model.')
    endif
end function late_activation_ext_force_icecover


!> set default values for selected ice cover model and allocate
function select_icecover_model(icecover, modeltype) result(istat)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover  !> data structure containing ice cover data
    integer                                    , intent(in)    :: modeltype !> desired ice cover type
    integer                                                    :: istat     !> status flag for allocation
    !
    ! Local variables
    !
    ! None
!
!! executable statements -------------------------------------------------------
!
    icecover%modeltype                 = modeltype

    icecover%hisout                    = .false.
    icecover%mapout                    = .false.
    
    icecover%areafrac_forcing_available   = 0
    icecover%thick_ice_forcing_available  = 0
    
    if (modeltype == ICECOVER_NONE) then
       icecover%apply_pressure            = .false.
       icecover%apply_friction            = .false.
       icecover%reduce_surface_exchange   = .false.
       icecover%reduce_waves              = .false.
       icecover%reduce_wind               = .false.
    else
       icecover%apply_pressure            = .true.
       icecover%apply_friction            = .false.
       icecover%reduce_surface_exchange   = .false.
       icecover%reduce_waves              = .false.
       icecover%reduce_wind               = .false.
    endif

    icecover%ice_albedo                = 0.75_fp
    icecover%snow_albedo               = 0.9_fp
    icecover%ice_dens                  = 917.0_fp
    icecover%frict_type                = FRICT_AS_DRAG_COEFF
    icecover%frict_val                 = 0.005_fp
    
    if (modeltype == ICECOVER_NONE) then
        istat = clr_icecover(icecover)
    else
        istat = 0
    endif
end function select_icecover_model


!> Allocate the arrays of an icecover data structure.
function alloc_icecover(icecover, nmlb, nmub) result(istat)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover  !> data structure containing ice cover data
    integer                                    , intent(in)    :: nmlb      !> lower bound index for spatial data arrays
    integer                                    , intent(in)    :: nmub      !> upper bound index for spatial data arrays
    integer                                                    :: istat     !> status flag for allocation
    !
    ! Local variables
    !
    ! NONE
!
!! executable statements -------------------------------------------------------
!
    istat = 0
    !
    ! state
    !
    if (icecover%modeltype /= ICECOVER_NONE) then
       if (istat==0) allocate(icecover%areafrac  (nmlb:nmub), STAT = istat)
       if (istat==0) allocate(icecover%thick_ice (nmlb:nmub), STAT = istat)
       if (icecover%modeltype /= ICECOVER_EXT) then
          if (istat==0) allocate(icecover%thick_snow(nmlb:nmub), STAT = istat)
       endif
       if (istat==0) then
          icecover%areafrac  = 0.0_fp
          icecover%thick_ice = 0.0_fp
          if (icecover%modeltype /= ICECOVER_EXT) then
             icecover%thick_snow = 0.0_fp
          endif
       endif
    endif
    !
    ! extra
    !
    if (icecover%modeltype /= ICECOVER_NONE) then
       if (istat==0) allocate(icecover%qh_air2ice(nmlb:nmub), STAT = istat)
       if (istat==0) allocate(icecover%qh_ice2wat(nmlb:nmub), STAT = istat)
       if (istat==0) allocate(icecover%pressure  (nmlb:nmub), STAT = istat)
       if (istat==0) then
          icecover%qh_air2ice = 0.0_fp
          icecover%qh_ice2wat = 0.0_fp
          icecover%pressure   = 0.0_fp
       endif
    endif
end function alloc_icecover


!> Clear the arrays of sedtra_type data structure.
function clr_icecover(icecover) result (istat)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover  !> data structure containing ice cover data
    integer                                                    :: istat     !> status flag for deallocation
    !
    ! Local variables
    !
    ! NONE
!
!! executable statements -------------------------------------------------------
!
    istat = 0
    !
    ! state
    !
    if (associated(icecover%areafrac  )) deallocate(icecover%areafrac  , STAT = istat)
    if (associated(icecover%thick_ice )) deallocate(icecover%thick_ice , STAT = istat)
    if (associated(icecover%thick_snow)) deallocate(icecover%thick_snow, STAT = istat)
    !
    ! extra
    !
    if (associated(icecover%qh_air2ice)) deallocate(icecover%qh_air2ice, STAT = istat)
    if (associated(icecover%qh_ice2wat)) deallocate(icecover%qh_ice2wat, STAT = istat)
    if (associated(icecover%pressure  )) deallocate(icecover%pressure  , STAT = istat)
end function clr_icecover

!--------------- following routines should move to ice kernel ---------------

!> Update the ice pressure array. I hope that we can extract the initial update_icecover from m_fm_icecover to here ...
!subroutine update_icecover(icecover, nm)
!!!--declarations----------------------------------------------------------------
!    !
!    ! Function/routine arguments
!    !
!    type (icecover_type)                       , intent(inout) :: icecover  !> data structure containing ice cover data
!    integer                                    , intent(in)    :: nm        !> Spatial index
!    !
!    ! Local variables
!    !
!!
!!! executable statements -------------------------------------------------------
!!
!    select case (icecover%modeltype)
!    case (ICECOVER_KNMI)
!        ! follow De Bruin & Wessels (1975)
!    case (ICECOVER_SEMTNER)
!        ! follow Semtner (1975)
!    case default
!        ! by default no growth
!    end select
!end subroutine update_icecover


!> Update the ice pressure array.
subroutine update_icepress(icecover, ag)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover  !> data structure containing ice cover data
    real(fp)                                   , intent(in)    :: ag        !> gravitational accelaration (m/s2)
    !
    ! Local variables
    !
    integer                         :: nm        !> Spatial loop index
    real(fp)                        :: ice_dens  !> Local variable for ice density
    real(fp), dimension(:), pointer :: areafrac  !> Pointer to ice area fraction array
    real(fp), dimension(:), pointer :: pressure  !> Pointer to ice pressure array
    real(fp), dimension(:), pointer :: thick_ice !> Pointer to ice thickness array
!
!! executable statements -------------------------------------------------------
!
    areafrac  => icecover%areafrac
    pressure  => icecover%pressure
    thick_ice => icecover%thick_ice
    ice_dens  =  icecover%ice_dens
    do nm = lbound(pressure,1),ubound(pressure,1)
        pressure(nm) = areafrac(nm) * thick_ice(nm) * ice_dens * ag
        ! + optionally snow or is that weight always negligible?
    enddo
end subroutine update_icepress

end module icecover_module

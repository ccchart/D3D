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
!use precision
private

!
! public data types
!
public icecover_type

!
! public routines
!
public null_icecover
public alloc_icecover
public clr_icecover
public update_icepress
!
! puclic constants
!
public ICECOVER_NONE
public ICECOVER_EXT
public fp

integer, parameter :: fp = 8

integer, parameter :: ICECOVER_NONE = 0
integer, parameter :: ICECOVER_EXT  = 1

! ice cover type
type icecover_type
   logical  :: mapout                           !> flag indicating whether ice cover should be written to map-file
   logical  :: clip_waves                       !> flag indicating whether waves need to be clipped
   !
   integer  :: modeltype                        !> type of the ice cover (one of ICECOVER_...)
   !
   real(fp) :: dens                             !> ice density
   !
   real(fp), dimension(:), pointer :: areafrac  !> area fraction covered by ice (-)
   real(fp), dimension(:), pointer :: pressure  !> pressure exerted by the ice cover (Pa)
   real(fp), dimension(:), pointer :: thickness !> ice cover thickness (m)
end type icecover_type

contains

!> Nullify/initialize an icecover data structure.
subroutine null_icecover(icecover)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover
    !
    ! Local variables
    !
    ! None
!
!! executable statements -------------------------------------------------------
!
    icecover%modeltype  = ICECOVER_NONE
    icecover%dens       = 910.0_fp
    icecover%mapout     = .false.
    icecover%clip_waves = .false.
    nullify(icecover%areafrac)
    nullify(icecover%pressure)
    nullify(icecover%thickness)
end subroutine null_icecover

!> Allocate the arrays of an icecover data structure.
function alloc_icecover(icecover, nmlb, nmub) result(istat)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover
    integer                                    , intent(in)    :: nmlb
    integer                                    , intent(in)    :: nmub
    integer                                                    :: istat
    !
    ! Local variables
    !
    ! NONE
!
!! executable statements -------------------------------------------------------
!
                  allocate(icecover%areafrac (nmlb:nmub), STAT = istat)
    if (istat==0) allocate(icecover%pressure (nmlb:nmub), STAT = istat)
    if (istat==0) allocate(icecover%thickness(nmlb:nmub), STAT = istat)
end function alloc_icecover


!> Clear the arrays of sedtra_type data structure.
function clr_icecover(icecover) result (istat)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover
    integer                                                    :: istat
    !
    ! Local variables
    !
    ! NONE
!
!! executable statements -------------------------------------------------------
!
    istat = 0
    if (associated(icecover%areafrac )) deallocate(icecover%areafrac , STAT = istat)
    if (associated(icecover%pressure )) deallocate(icecover%pressure , STAT = istat)
    if (associated(icecover%thickness)) deallocate(icecover%thickness, STAT = istat)
end function clr_icecover


!> Update the ice pressure array.
subroutine update_icepress(icecover, ag)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    type (icecover_type)                       , intent(inout) :: icecover
    real(fp)                                   , intent(in)    :: ag       !> gravitational accelaration (m/s2)
    !
    ! Local variables
    !
    integer                         :: nm        !> Spatial loop index
    real(fp)                        :: density   !> Local variable for ice density
    real(fp), dimension(:), pointer :: areafrac  !> Pointer to ice area fraction array
    real(fp), dimension(:), pointer :: pressure  !> Pointer to ice pressure array
    real(fp), dimension(:), pointer :: thickness !> Pointer to ice thickness array
!
!! executable statements -------------------------------------------------------
!
    areafrac  => icecover%areafrac
    pressure  => icecover%pressure
    thickness => icecover%thickness
    density   =  icecover%dens
    do nm = lbound(pressure,1),ubound(pressure,1)
        pressure(nm) = areafrac(nm) * thickness(nm) * ag * density
    enddo
end subroutine update_icepress

end module icecover_module


module m_fm_icecover
use icecover_module
implicit none

!
! Global data
!
type(icecover_type), target                                :: ice_data
!
real(fp), dimension(:), pointer                            :: ice_af
real(fp), dimension(:), pointer                            :: ice_h
real(fp), dimension(:), pointer                            :: ice_p
integer                                                    :: ja_icecover = ICECOVER_NONE
integer                                                    :: ja_aice_read = 0
integer                                                    :: ja_hice_read = 0

contains


!> Nullify/initialize ice data structure.
subroutine fm_null_icecover()
!!--declarations----------------------------------------------------------------
    !
    implicit none
    !
    ! Function/routine arguments
    !
    ! None
    !
    ! Local variables
    !
    ! None
!
!! executable statements -------------------------------------------------------
!
    call null_icecover(ice_data)
    nullify(ice_af)
    nullify(ice_h)
    nullify(ice_p)
    ja_icecover = ICECOVER_NONE
end subroutine fm_null_icecover


!> Allocate the arrays of ice data structure.
subroutine fm_alloc_icecover(ndx)
!!--declarations----------------------------------------------------------------
    !
    implicit none
    !
    ! Function/routine arguments
    !
    integer                                    , intent(in)    :: ndx
    !
    ! Local variables
    !
    integer                                                    :: istat
!
!! executable statements -------------------------------------------------------
!
    if (associated(ice_af)) return ! don't allocate if already allocated
    istat = alloc_icecover(ice_data, 1, ndx)
    ice_af => ice_data%areafrac
    ice_h  => ice_data%thickness
    ice_p  => ice_data%pressure
end subroutine fm_alloc_icecover


!> Clear the arrays of ice data structure.
subroutine fm_clr_icecover()
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
    integer                                                    :: istat
!
!! executable statements -------------------------------------------------------
!
    istat = clr_icecover(ice_data)
    call fm_null_icecover()
    ja_icecover = ICECOVER_NONE
    ja_aice_read = 0
    ja_hice_read = 0
end subroutine fm_clr_icecover

!> Update the ice pressure array.
subroutine fm_update_icepress(ag)
!!--declarations----------------------------------------------------------------
    implicit none
    !
    ! Function/routine arguments
    !
    real(fp)                                   , intent(in)    :: ag       !> gravitational accelaration (m/s2)
    !
    ! Local variables
    !
    ! None
!
!! executable statements -------------------------------------------------------
!
    call update_icepress(ice_data, ag)
end subroutine fm_update_icepress

end module m_fm_icecover

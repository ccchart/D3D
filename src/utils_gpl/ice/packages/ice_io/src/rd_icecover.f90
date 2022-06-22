module icecover_input_module
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2021.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id$
!  $HeadURL$
!!--description-----------------------------------------------------------------

private

!
! functions and subroutines
!
public rd_icecover
public echo_icecover

contains

!> Read ice cover parameters (Note: the meteo module should already have been initialized)
subroutine rd_icecover(icecover, md_ptr, chapter, error)
!!--declarations----------------------------------------------------------------
    use precision
    use icecover_module, only: icecover_type, alloc_icecover, select_icecover_model, ICECOVER_NONE, ICECOVER_EXT, ICECOVER_KNMI, ICECOVER_SEMTNER, FRICT_AS_DRAG_COEFF
    use MessageHandling, only: mess, LEVEL_ERROR
    use properties
    !
    implicit none
    !
!
! Arguments
!
    type (icecover_type)         , intent(inout) :: icecover !> ice cover data structure containing the data read
    type(tree_data)              , pointer       :: md_ptr   !> pointer to the input file
    character(len=*)             , intent(in)    :: chapter  !> chapter name of the ice section
    logical                      , intent(out)   :: error    !> flag indicating an execution error
!
! Local variables
!
    logical                                      :: success
    integer                                      :: istat
    integer                                      :: model
    character(256)                               :: filename
    character(256)                               :: tmp
!
!! executable statements -------------------------------------------------------
!
    error = .false.
    model =  ICECOVER_NONE
    
    tmp = ' '
    call prop_get(md_ptr,chapter,'IceCoverModel',tmp)
    call str_lower(tmp,len(tmp))
    select case (tmp)
    case ('none', ' ')
       ! default selected
    case ('external')
       model =  ICECOVER_EXT
    case ('knmi')
       model =  ICECOVER_KNMI
    case ('semtner','deltares')
       model =  ICECOVER_SEMTNER
    case default
       call mess(LEVEL_ERROR, 'invalid ice cover model "'//trim(tmp)//'"')
       error = .true.
       ! still want to properly initialize the icecover module, so don't return immediately
    end select
    istat = select_icecover_model(icecover, model)
    if (error) return
    !
    ! Process flags
    !
    call prop_get(md_ptr,chapter,'ApplyPressure',icecover%apply_pressure)
    call prop_get(md_ptr,chapter,'ApplyFriction',icecover%apply_friction)
    call prop_get(md_ptr,chapter,'ReduceSurfExch',icecover%reduce_surface_exchange)
    call prop_get(md_ptr,chapter,'ReduceWaves',icecover%reduce_waves)
    call prop_get(md_ptr,chapter,'ReduceWind',icecover%reduce_wind)
    !
    ! Parameters
    !
    call prop_get(md_ptr,chapter,'IceDensity',icecover%dens_ice)
    !
    tmp = ' '
    call prop_get(md_ptr,chapter,'IceFricType',tmp)
    call str_lower(tmp,len(tmp))
    select case (tmp)
    case ('cdrag', ' ')
       ! default selected
       icecover%frict_type = FRICT_AS_DRAG_COEFF
    case default
       call mess(LEVEL_ERROR, 'invalid ice cover friction type "'//trim(tmp)//'", only "cdrag" supported.')
       error = .true.
       return
    end select
    call prop_get(md_ptr,chapter,'IceFricValue',icecover%frict_val)
    !
    select case (icecover%modeltype)
    case (ICECOVER_EXT)
       ! No extra parameters for external forcing of ice cover
    case default
    
    end select
    !
    ! Output flags
    !
    call prop_get(md_ptr,chapter,'AddIceToHis',icecover%hisout)
    call prop_get(md_ptr,chapter,'AddIceToMap',icecover%mapout)
    !
    !istat = alloc_icecover(icecover, nmlb, nmub)
    !if (istat /= 0) then
    !   call mess(LEVEL_ERROR, 'error while allocation memory for ice cover')
    !   error = .true.
    !   return
    !endif
end subroutine rd_icecover


!> Write ice cover parametersto diagnostic
function echo_icecover(icecover, lundia) result (error)
!!--declarations----------------------------------------------------------------
    use precision
    use icecover_module, only: icecover_type, alloc_icecover, ICECOVER_NONE, ICECOVER_EXT, ICECOVER_KNMI, ICECOVER_SEMTNER, FRICT_AS_DRAG_COEFF
    use MessageHandling, only: mess, LEVEL_ERROR
    !
    implicit none
!
! Arguments
!
    type (icecover_type)         , intent(inout) :: icecover !> ice cover data structure containing the data read
    integer                      , intent(in)    :: lundia   !> unit number of diagnostics file
    logical                                      :: error    !> flag indicating an execution error
!
! Local variables
!
    character(45)                                                     :: txtput1
    character(120)                                                    :: txtput2
!
!! executable statements -------------------------------------------------------
!
    !
    ! don't print any ice messages if there is no ice cover
    !
    error = .false.
    if (icecover%modeltype == ICECOVER_NONE) return
    !
    write (lundia, '(a)' ) '*** Start  of ice cover input'
    
    txtput1 = '  Ice Cover Model'
    select case (icecover%modeltype)
    case (ICECOVER_EXT)
        txtput2 = 'external'
    case (ICECOVER_KNMI)
        txtput2 = 'KNMI: De Bruin & Wessels (1975)'
    case (ICECOVER_SEMTNER)
        txtput2 = 'Semtner (1975)'
    case default
        txtput2 = 'unknown'
    end select
    write (lundia, '(3a)') txtput1, ': ', txtput2
    
    if (icecover%modeltype == ICECOVER_EXT) then
       txtput1 = '  Area Fraction Forcing'
       call write_logical(lundia, txtput1, icecover%areafrac_forcing_available /= 0)

       txtput1 = '  Ice Thickness Forcing'
       call write_logical(lundia, txtput1, icecover%thick_ice_forcing_available /= 0)
       
       if (icecover%areafrac_forcing_available == 0 .and. icecover%thick_ice_forcing_available == 0) then
          call mess(LEVEL_ERROR, 'icecover set to external forcing but no forcing data found')
          error = .true.
       elseif (icecover%areafrac_forcing_available == 0) then
          call mess(LEVEL_ERROR, 'icecover set to external forcing but no area fraction forcing found')
          error = .true.
       elseif (icecover%thick_ice_forcing_available == 0) then
          call mess(LEVEL_ERROR, 'icecover set to external forcing but no ice thickness forcing found')
          error = .true.
       endif
    endif
    !
    ! process flags
    !
    txtput1 = '  Apply Pressure'
    call write_logical(lundia, txtput1, icecover%apply_pressure)

    txtput1 = '  Apply Friction'
    call write_logical(lundia, txtput1, icecover%apply_friction)

    txtput1 = '  Reduce Surface Exchange'
    call write_logical(lundia, txtput1, icecover%reduce_surface_exchange)

    txtput1 = '  Reduce Waves'
    call write_logical(lundia, txtput1, icecover%reduce_waves)

    txtput1 = '  Reduce Wind'
    call write_logical(lundia, txtput1, icecover%reduce_wind)

    !
    ! parameters
    !
    txtput1 = '  Ice Density'
    write (lundia, '(2a,e20.4)') txtput1, ': ', icecover%dens_ice

    txtput1 = '  Ice Cover Friction Type'
    select case (icecover%frict_type)
    case (FRICT_AS_DRAG_COEFF)
        txtput2 = 'drag coefficient'
    case default
        txtput2 = 'unknown'
    end select
    write (lundia, '(3a)') txtput1, ': ', txtput2
    txtput1 = '  Ice Cover Friction Value'
    write (lundia, '(2a,e20.4)') txtput1, ': ', icecover%frict_val

    !
    ! output flags
    !
    txtput1 = '  History File Output'
    call write_logical(lundia, txtput1, icecover%hisout)
    
    txtput1 = '  Map File Output'
    call write_logical(lundia, txtput1, icecover%mapout)

    write (lundia, '(a)' ) '*** End    of ice cover input'
    write (lundia, *)
end function echo_icecover


!> Support routine to write logical flags to diagnostic file
subroutine write_logical(lundia, txtput1, option)
!
! Arguments
!
    integer                      , intent(in)    :: lundia   !> unit number of diagnostics file
    character(*)                 , intent(in)    :: txtput1  !> base string
    logical                      , intent(in)    :: option   !> logical option to report
!
! Local variables
!
!    NONE
!
!! executable statements -------------------------------------------------------
!
    if (option) then
        write (lundia, '(3a)') txtput1, ': Yes'
    else
        write (lundia, '(3a)') txtput1, ': No'
    endif
end subroutine write_logical

end module icecover_input_module
!!  Copyright (C)  Stichting Deltares, 2012-2023.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.

!> This module implements the statistical output in D-Flow FM.
module m_statistical_output
   use MessageHandling
   
private

   public realloc
   public dealloc

   !> Realloc memory cross-section definition or cross-sections
   interface realloc
      module procedure realloc_stat_output
   end interface

   !> Free the memory of cross-section definition or cross-sections
   interface dealloc
      module procedure dealloc_stat_output
   end interface dealloc

   enum, bind(c)
      enumerator :: SO_CURRENT = 1
      enumerator :: SO_AVERAGE = 2
      enumerator :: SO_MAX     = 3
      enumerator :: SO_MIN     = 4
   end enum

   !> Derived type for the statistical output items. 
   type, public :: t_statistical_output_item
      character(len=Idlen)             :: name                 !< Name of the output item on the NETCDF file.      
      character(len=Idlen)             :: long_name            !< Long name of the output item on the NETCDF file.      
      character(len=Idlen)             :: standard_name        !< Standard name of the output item on the NETCDF file.                     
      integer                          :: location_specifier   !< Specifies the location, where the variable is specified (e.g. UNC_LOC_S or UNC_LOC_U)
      integer                          :: operation_id         !< Specifies the kind of operation to perform on the output variable.
      integer                          :: total_steps_count    !< Number of elements in moving average.
      integer                          :: current_step         !< Latest entry in the work array. MOD(current_step+1,total_steps_count) is the next 
                                                               !< item to remove.   
      double precision, pointer, dimension(:)   :: output      !< Array that is to be written to the Netcdf file. In case the current values are
                                                               !< required this variable points to the basic variable (e.g. s1).
                                                               !< Otherwise during the simulation the intermediate results are stored.
      double precision, pointer, dimension(:)   :: input       !< In case a statistical operation is requested. This pointer points to the
                                                               !< basic variable.
      double precision, pointer, dimension(:)   :: moving_average  !< In case a moving average is requested. This pointer points to the
                                                               !< actual average values.
      double precision, pointer, dimension(:,:) :: samples     !< In case a moving average is requested. This pointer points to the
                                                               !< work array, where the different samples are stored.

   end type t_statistical_output_item

   !> Derived type to store the cross-section set
   type, public :: t_statistical_output_set
      integer                                                :: size = 0                  !< Actual size of cross-section set
      integer                                                :: growsby = 200             !< Increment for cross-section set
      integer                                                :: count= 0                  !< Actual number of cross-section sets
      type(t_statistical_output_item), pointer, dimension(:) :: statout                   !< Current cross-section
   end type t_statistical_output_set

contains

   subroutine realloc_stat_output(statoutput)
      ! Modules
      use m_alloc

      implicit none
      ! Input/output parameters
      type(t_statistical_output_set), intent(inout)   :: statoutput !< Current cross-section definition
      
   
      ! Local variables
      integer                   :: ierr
      type(t_statistical_output_item), pointer, dimension(:)    :: oldstats
   
      ! Program code
   
      if (statoutput%size > 0) then
         oldstats=>statoutput%statout
      endif
   
      if (statoutput%growsBy <=0) then
         statoutput%growsBy = 200
      endif
      allocate(statoutput%statout(statoutput%size+statoutput%growsBy),stat=ierr)
      call aerr('statoutput%statout(statoutput%size+statoutput%growsBy)',ierr,statoutput%size+statoutput%growsBy)
   
      if (statoutput%size > 0) then
         statoutput%statout(1:statoutput%size) = oldstats(1:statoutput%size)
         deallocate(oldstats)
      endif
      statoutput%size = statoutput%size+statoutput%growsBy
   end subroutine realloc_stat_output

   subroutine dealloc_stat_output(statoutput)
      implicit none
      ! Input/output parameters
      type(t_statistical_output_set), intent(inout)   :: statoutput !< Current cross-section definition

      if (statoutput%size> 0) then
         deallocate(statoutput%statout)
      endif
   end subroutine dealloc_stat_output

end module m_statistical_output

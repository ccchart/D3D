!> This module implements the statistical output in D-Flow FM.
module m_input_items
   use MessageHandling

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

private
   
   !> Derived type for the statistical output items. 
   type t_statistical_output_item
      character(len=Idlen)             :: key             !< Key of the input item in the MDU file (e.g. wrimap_s1).                       
      character(len=Idlen)             :: name            !< Name of the output item on the NETCDF file.      
      character(len=Idlen)             :: long_name       !< Long name of the output item on the NETCDF file.      
      character(len=Idlen)             :: standard_name   !< Standard name of the output item on the NETCDF file.                     
      character(len=Idlen)             :: input_value     !< Retrieved value from the MDU file for the given key.             
      integer                          :: location specifier !< Specifies the location, where the variable is specified (e.g. UNC_LOC_S or UNC_LOC_U)
   end type t_statistical_output_item

   !> Derived type to store the cross-section set
   type, public :: t_statistical_output_set
      integer                                                :: size = 0                  !< Actual size of cross-section set
      integer                                                :: growsBy = 200            !< Increment for cross-section set
      integer                                                :: count= 0                  !< Actual number of cross-section sets
      type(t_statistical_output_item), pointer, dimension(:) :: statout                     !< Current cross-section
   end type t_statistical_output_set

contains

   subroutine realloc_stat_output(statoutput)
      ! Modules

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
         statoutput%statout(1:statoutput%size) = oldDefs(1:statoutput%size)
         deallocate(oldstats)
      endif
      statoutput%size = statoutput%size+statoutput%growsBy
   end subroutine realloc

   subroutine dealloc_stat_output(statoutput)
      implicit none
      ! Input/output parameters
      type(t_statistical_output_set), intent(inout)   :: statoutput !< Current cross-section definition

      if (statoutput%size> 0) then
         deallocate(statoutput)
      endif
   end subroutine dealloc_stat_output

end module m_input_items
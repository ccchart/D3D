!> This module contains the definition of the input items for the [output] section of the 
!! MDU file.
module m_input_items
   use MessageHandling

private
   
   !> Derived type for the input items, defining one entry [output] section of the MDU file. 
   type t_input_items
      character(len=Idlen)             :: key             !< Key of the input item in the MDU file (e.g. wrimap_s1).                       
      character(len=Idlen)             :: name            !< Name of the output item on the NETCDF file.      
      character(len=Idlen)             :: long_name       !< Long name of the output item on the NETCDF file.      
      character(len=Idlen)             :: standard_name   !< Standard name of the output item on the NETCDF file.                     
      character(len=Idlen)             :: input_value     !< Retrieved value from the MDU file for the given key.             
      integer                          :: location_specifier !< Specifies the location, where the variable is specified (e.g. UNC_LOC_S or UNC_LOC_U)
   end type t_input_items

contains

end module m_input_items
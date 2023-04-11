!> This module contains the definition of the input items for the [output] section of the 
!! MDU file.
module m_output_config
   use MessageHandling

   public t_output_quantity_config
private
   
   !> Derived type for the input items, defining one entry [output] section of the MDU file. 
   type t_output_quantity_config
      character(len=Idlen)             :: key             !< Key of the input item in the MDU file (e.g. wrimap_s1).                       
      character(len=Idlen)             :: name            !< Name of the output item on the NETCDF file.      
      character(len=Idlen)             :: long_name       !< Long name of the output item on the NETCDF file.      
      character(len=Idlen)             :: unit            !< unit of the output item on the NETCDF file.      
      character(len=Idlen)             :: standard_name   !< Standard name of the output item on the NETCDF file.                     
      character(len=Idlen)             :: input_value     !< Original user-provided input valuestring (unparsed) (<<key>> = <<input value>>.         
      integer                          :: location_specifier !< Specifies the location, where the variable is specified (One of UNC_LOC_CN, UNC_LOC_S, 
                                                             !< UNC_LOC_U, UNC_LOC_L, UNC_LOC_S3D, UNC_LOC_U3D, UNC_LOC_W, UNC_LOC_WU)
   end type t_output_quantity_config

contains

end module m_output_config
!> This module contains the definition of the input items for the [output] section of the 
!! MDU file.
module m_output_config
   use MessageHandling
   use coordinate_reference_system

   integer, parameter, public :: UNC_LOC_CN  = 1  !< Data location: corner point.
   integer, parameter, public :: UNC_LOC_S   = 2  !< Data location: pressure point.
   integer, parameter, public :: UNC_LOC_U   = 3  !< Data location: horizontal velocity point.
   integer, parameter, public :: UNC_LOC_L   = 13 !< Data location: horizontal net link.
   integer, parameter, public :: UNC_LOC_S3D = 4  !< Data location: pressure point in all layers.
   integer, parameter, public :: UNC_LOC_U3D = 5  !< Data location: horizontal velocity point in all layers.
   integer, parameter, public :: UNC_LOC_W   = 6  !< Data location: vertical velocity point on all layer interfaces.
   integer, parameter, public :: UNC_LOC_WU  = 16 !< Data location: vertical viscosity point on all layer interfaces.
   integer, parameter, public :: UNC_LOC_WB        = 21 !< Data location: his file water balance
   integer, parameter, public :: UNC_LOC_SOSI      = 22 !< Data location: his file sources and sinks
   integer, parameter, public :: UNC_LOC_GENSTRU   = 23 !< Data location: his file general structure data
   integer, parameter, public :: UNC_LOC_DAMBREAK  = 24 !< Data location: his file dambreak data
   integer, parameter, public :: UNC_LOC_PUMP      = 24 !< Data location: his file pump data
   
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
      integer                          :: num_additional_attributes  !< number of additional attributes
      type(nc_attribute), pointer      :: additional_attributes(:)   !< optional additional attributes for this entity
   end type t_output_quantity_config

contains

end module m_output_config
module fm_statistical_output
   use m_output_config
   use m_statistical_output
   implicit none
   
   type(t_output_quantity_config), allocatable, dimension(:) :: out_quan_conf_his
   type(t_output_quantity_config), allocatable, dimension(:) :: out_quan_conf_map
   type(t_output_quantity_config), allocatable, dimension(:) :: out_quan_conf_classmap

   type(t_output_variable_set), allocatable :: out_variable_set_his
   type(t_output_variable_set), allocatable :: out_variable_set_map
   type(t_output_variable_set), allocatable :: out_variable_set_classmap
   
contains
   
end module fm_statistical_output
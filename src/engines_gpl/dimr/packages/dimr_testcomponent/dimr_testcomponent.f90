!  dimr_testcomponent.f90 
!
!  FUNCTIONS/SUBROUTINES exported from dimr_testcomponent.dll:
!  dimr_testcomponent - subroutine 
!
module bmi
  use iso_c_binding
  implicit none
  integer(c_int), BIND(C, name="MAXSTRLEN") :: MAXSTRLEN = 1024
  
contains
  
  subroutine main() bind(C, name="main")
     !DEC$ ATTRIBUTES DLLEXPORT :: main
     ! Somehow intel fortran compiler expects a main routine in the dll, it is required since interactor is used (and win calls)
  end subroutine main

  !> The initialize() function accepts a string argument that
  !! gives the name (and path) of its "main input file", called
  !! a configuration file. This function should perform all tasks
  !! that are to take place before entering the model's time loop.
  integer(c_int) function initialize(c_config_file) result(c_iresult) bind(C, name="initialize")
    !DEC$ ATTRIBUTES DLLEXPORT :: initialize
    use iso_c_binding, only: c_char
    character(kind=c_char),intent(in)    :: c_config_file(MAXSTRLEN)
    character(len=strlen(c_config_file)) :: config_file
    c_iresult = 0 ! TODO: is this return value BMI-compliant?
  end function initialize


  integer function update(dt) bind(C, name="update")
    !DEC$ ATTRIBUTES DLLEXPORT :: update
    use iso_c_binding, only: c_double
    real(c_double), value, intent(in) :: dt

    update = 0
  end function update


  integer function finalize() bind(C, name="finalize")
    !DEC$ ATTRIBUTES DLLEXPORT :: finalize

    finalize = 0
  end function finalize

   subroutine get_var(c_var_name, x) bind(C, name="get_var")
     !DEC$ ATTRIBUTES DLLEXPORT :: get_var
     ! Return a pointer to the variable
     use iso_c_binding, only: c_double, c_char, c_loc
     use iso_c_utils
     character(kind=c_char), intent(in) :: c_var_name(*) !< Variable name. May be slash separated string "name/item/field": then get_compound_field is called.
     type(c_ptr), intent(inout) :: x
     integer(c_int), target, allocatable, save :: xi(:,:)
   
     integer :: i, j, k
   
   
     ! The fortran name of the attribute name
     character(len=strlen(c_var_name)) :: var_name
     character(len=strlen(c_var_name)) :: tmp_var_name
     character(len=strlen(c_var_name)) :: varset_name, item_name, field_name !< For parsing compound variable names.
     double precision, dimension(:), allocatable :: uabs
   
     ! Store the name
     var_name = char_array_to_string(c_var_name)
   
     ! Please be conservative in adding variables here. Most variables
     ! can be computed outside.
     ! You can generate extra variables here.
     select case(var_name)
     case("uabs")
!         if (.not.allocated(uabs))  allocate(uabs(size(ucx)))
!         x = c_loc(uabs)
        continue
     end select
   end subroutine get_var
   
   subroutine set_var(c_var_name, xptr) bind(C, name="set_var")
     !DEC$ ATTRIBUTES DLLEXPORT :: set_var
     ! Return a pointer to the variable
     use iso_c_binding, only: c_double, c_char, c_loc, c_f_pointer
   
     character(kind=c_char), intent(in) :: c_var_name(*)
     type(c_ptr), value, intent(in) :: xptr
     character(len=strlen(c_var_name)) :: var_name
     character(kind=c_char), pointer :: valuestr
   
     ! Store the name and pointer to the value
     var_name = char_array_to_string(c_var_name, strlen(c_var_name))
     call c_f_pointer(xptr, valuestr)
     
   end subroutine set_var
   
   
   
   subroutine dimr_testcomponent
   
     ! Expose subroutine dimr_testcomponent to users of this DLL
     !
     !DEC$ ATTRIBUTES DLLEXPORT::dimr_testcomponent
   
     ! Variables
   
    ! Body of dimr_testcomponent
   
   end subroutine dimr_testcomponent

   
! Utility functions, move these to interop module
! Make functions pure so they can be used as input arguments.
   integer(c_int) pure function strlen(char_array)
    character(c_char), intent(in) :: char_array(MAXSTRLEN)
    integer :: inull, i
    strlen = 0
    do i = 1, size(char_array)
       if (char_array(i) .eq. C_NULL_CHAR) then
          strlen = i-1
          exit
       end if
    end do
   end function strlen

   pure function char_array_to_string(char_array, length)
      integer(c_int), intent(in) :: length
      character(c_char),intent(in) :: char_array(length)
      character(len=length) :: char_array_to_string
      integer :: i
      do i = 1, length
         char_array_to_string(i:i) = char_array(i)
      enddo
   end function char_array_to_string
   
end module bmi


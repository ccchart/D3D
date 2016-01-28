! write_to_matlab.f90 --
!     Module to facilitate communication with MATLAB
!
module write_to_matlab
    implicit none

    interface to_matlab
        module procedure to_matlab_int
        module procedure to_matlab_real
        module procedure to_matlab_double
        module procedure to_matlab_string
    end interface to_matlab

    character(len=1), parameter, private :: tab = achar(9)

contains
subroutine to_matlab_int( lun, key, value )
    integer, intent(in) :: lun
    character(len=*), intent(in) :: key
    integer, intent(in) :: value

    write( lun, '(2a,i0)' ) trim(key), tab, value
end subroutine to_matlab_int

subroutine to_matlab_real( lun, key, value )
    integer, intent(in) :: lun
    character(len=*), intent(in) :: key
    real, intent(in) :: value

    character(len=20) :: string

    write( string, '(e15.7)' ) value
    write( lun, '(3a)' ) trim(key), tab, adjustl(string)
end subroutine to_matlab_real

subroutine to_matlab_double( lun, key, value )
    integer, intent(in) :: lun
    character(len=*), intent(in) :: key
    real(kind=kind(1.0d0)), intent(in) :: value

    character(len=30) :: string

    write( string, '(e23.17)' ) value
    write( lun, '(3a)' ) trim(key), tab, adjustl(string)
end subroutine to_matlab_double

subroutine to_matlab_string( lun, key, value )
    integer, intent(in) :: lun
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value

    write( lun, '(3a)' ) trim(key), tab, adjustl(value)
end subroutine to_matlab_string

end module write_to_matlab

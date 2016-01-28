! Demo program for xml_process
!
! The program reads an XML file and extracts a few pieces of information
!
! Note:
! In this example we do not check the structure of the
! XML file.
!
module getdata
    use tables_module
    use read_xml_primitives
    implicit none

    type(tables)                      :: list_of_tables
    type(table)                       :: current_table
    logical                           :: found_modulename, found_item, found_data
    character(len=40)                 :: table_name
    character(len=40), dimension(100) :: column_name
    integer :: number_columns

contains

! startfunc --
!     Called at the start of an element
!
recursive subroutine startfunc( tag, attribs, error )
    character(len=*)                 :: tag
    character(len=*), dimension(:,:) :: attribs
    logical                          :: error

    error = .false.

    select case( tag )
        case( 'modulename' )
            !
            ! We need to store the module name
            !
            found_modulename = .true.

        case( 'variables' )
            !
            ! Reset the list
            !
            number_columns = 0

        case( 'item' )
            found_item = .true. ! The data that follow: column name

        case( 'datafield' )
            found_data = .true. ! Read the data that follow as a list of values

        case default
            ! Ignore anything else
    end select
end subroutine

recursive subroutine datafunc( tag, data, error )
    character(len=*)                 :: tag
    character(len=*), dimension(:)   :: data
    logical                          :: error
    real(dp), dimension(:), pointer  :: table_values

    type(xml_parse)                  :: dummy_info         ! Not actually used
    character(len=1), dimension(2,1) :: attribs            ! Not used either
    logical                          :: has_values, endtag ! Ignored

    error = .false.

    !
    ! Store the data in the right variable
    !
    if ( found_modulename ) then
        table_name = data(1)
    endif
    if ( found_item ) then
        number_columns = number_columns + 1
        column_name(number_columns) = data(1)
    endif
    if ( found_data ) then
        !
        ! This is a bit tricky - rely on the read_xml_prims module
        !
        call read_xml_double_array( dummy_info, "Dummy", endtag, attribs, 0, data, size(data), table_values, has_values )

        call add( current_table, table_values, error )
        deallocate( table_values )
        if ( error ) then
            write(*,*) 'Wrong number of values? Reading data for ', table_name
        else
            call add( list_of_tables, table_name, current_table )
        endif
    endif

end subroutine

recursive subroutine endfunc( tag, error )
    character(len=*)               :: tag
    logical                        :: error

    integer                        :: i

    select case( tag )
        case( 'modulename' )
            found_modulename = .false.

        case( 'variables' )
            !
            ! Create the table
            !
            call define_table( current_table, column_name(1:number_columns) )

        case( 'item' )
            found_item = .false.

        case( 'datafield' )
            found_data = .false.

        case default
            ! Ignore anything else
    end select

end subroutine

end module getdata
program waqbininputdump
    implicit none
    !
    character(len=100)                             :: filename
    character(len=100)                             :: dummy
    character(len=40), dimension(4)                :: title
    character(len=20), dimension(:)  , allocatable :: parameter_name
    character(len=20), dimension(:)  , allocatable :: location_name
    integer          , dimension(:)  , allocatable :: location_index
    real             , dimension(:), allocatable :: data
    integer                                        :: time
    integer                                        :: nosys
    integer                                        :: noseg
    integer                                        :: lu
    integer                                        :: i
    integer                                        :: ierr
    !
    ! Body
    call getarg(1,filename)
    write(*,*) "Filename: ", trim(filename)
    call getarg(2,dummy)
    read(dummy,*) noseg
    write(*,*) "noseg: ", noseg
    open(newunit=lu, file = filename, access = 'stream', status = 'old' )
    !read( lu ) title
    !write(*,*) title
    !read( lu ) nosys, noseg
    !write(*,*) nosys, noseg
    !allocate( parameter_name(nosys) )
    !allocate( location_name(noseg), location_index(noseg) )
    allocate( data(noseg) )
    !read( lu ) parameter_name
    !write(*,*) parameter_name
    !read( lu ) ( location_index(i), location_name(i), i = 1,noseg )
    !write(*,*) ( location_index(i), location_name(i), i = 1,noseg )
    do
        read( lu, iostat = ierr ) time, data
        write(*,*)      "time: ", time
        write(*,*)                      data
        ! Or more explicitly:
        ! read( lu, iostat = ierr ) time, (data(isys,iseg) ,isys = 1,nosys), iseg = 1,noseg)
        !
        if ( ierr /= 0 ) then
            exit
        endif
    enddo
endprogram waqbininputdump

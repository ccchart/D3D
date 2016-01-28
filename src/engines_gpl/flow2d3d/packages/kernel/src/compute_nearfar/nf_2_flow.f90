subroutine nf_2_flow( filename,x_jet,y_jet,z_jet,s_jet,h_jet,b_jet, no_val)
    use precision
    integer, intent(inout)                   :: no_val
    real(fp), dimension(no_val), intent(out) :: x_jet, y_jet, z_jet, s_jet, h_jet, b_jet
    character(len=*), intent(in)             :: filename
    
    integer                 :: i, lun, ierr
    
    open( newunit=lun, file = filename )
    
    do i = 1,no_val
        read( lun, *, iostat = ierr ) x_jet(i), y_jet(i), z_jet(i), s_jet(i), h_jet(i), b_jet(i)
        if ( ierr /= 0 ) then
            no_val = i - 1
            exit
        endif
    enddo
    
    close( lun )
end subroutine nf_2_flow
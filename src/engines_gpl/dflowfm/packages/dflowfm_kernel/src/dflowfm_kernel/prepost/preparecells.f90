subroutine preparecells(md_netfile, jaidomain, jaiglobal_s, ierr)

 use network_data
 use unstruc_netcdf
 use dfm_error
 use unstruc_display, only: jareinitialize
 use gridoperations
 use unstruc_messages

 implicit none
 character(len=*), intent(in) :: md_netfile  !< net filename
 integer,          intent(in) :: jaidomain   !< read subdomain numbers (1) or not (0)
 integer,          intent(in) :: jaiglobal_s !< read global cell numbers (1) or not (0)
 integer,          intent(out):: ierr
 integer                      ::k, L, nv, flag, c, i, k1, nc1, j

 ierr = DFM_NOERR

 call readcells(md_netfile, ierr, jaidomain, jaiglobal_s, jareinitialize)
 if (ierr /= DFM_NOERR) then
    ! Cells could not be read, so make sure to call findcells on call site.
    return
 end if

 ! generate lne and lnn
 if(allocated(lne)) deallocate(lne)
 allocate(lne(2, numl))
 if(allocated(lnn)) deallocate(lnn)
 allocate(lnn(numl))
 lne = 0
 lnn = 0
 ! 2D
 do c = 1, nump
    nv = netcell(c)%n
    do k = 1, nv
       L = netcell(c)%lin(k)
       if (lne(1, L) .ne. c .and. lne(2, L) .ne. c) then
         j = lnn(L) + 1
         lnn(L) = j
         lne(j, L) = c
       endif
    enddo
 enddo
 ! 1D
 do c = nump+1, nump1d2d
    k = netcell(c)%nod(1)
    nv = nmk(k)
    do i = 1, nv
       L = nod(k)%lin(i)
       if (lne(1, L) .ne. -c .and. lne(2, L) .ne. -c ) then
          if (lnn(L) >= 2) then
             write (msgbuf, '(a,i0,a,i0,a)') 'Illegal: net link ', L, ' already has 2 2D cells, and is also connected to 1D netnode ', k, '.'
             call qnerror(trim(msgbuf), '', '')
          else
             ! the orientation of lne(1:2,L) should be the same as net link orientation kn(1:2,L). See subroutine is_1d_boundary_candidate().
             if (k == kn(1,L)) then
                j = 1
             else
                j = 2
             end if

             lnn(L) = lnn(L)+1 ! Only #j is filled for now, but it will always end up to be 2.
             lne(j, L) = -c
          end if
       endif
       if (kn(3, L) == 3  .or. kn(3, L) == 4  .or. kn(3, L) == 5 .or. kn(3, L) == 7) then   ! If 1d link L enters a 2d cell
          k1 = kn(1, L)
          if (k1 == k) then
             k1 = kn(2, L)
          endif
          nc1 = 0
          call incells(xk(k1), yk(k1), nc1)
          if (nc1 .ne. 0) then
             if (lnn(L) >= 2) then
                write (msgbuf, '(a,i0,a,i0,a)') 'Illegal: net link ', L, ' already has 2 2D cells, and is also connected to 1D netnode ', k, '.'
                call qnerror(trim(msgbuf), '', '')
             else
                ! the orientation of lne(1:2,L) should be the same as net link orientation kn(1:2,L).
                if (k == kn(1,L)) then
                   j = 2
                else
                   j = 1
                end if

                lne(j, L) = nc1 ! j here is the position of the 2D netcell.
                lnn(L) = 2 ! Must always end up as 2 here.
             end if
          endif
       endif
    enddo
 enddo

 call update_cell_circumcenters()

 end subroutine preparecells

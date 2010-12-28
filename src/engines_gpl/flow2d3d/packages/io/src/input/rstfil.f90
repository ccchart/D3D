subroutine rstfil(lundia    ,error     ,restid    ,lturi     ,mmax      , &
                & nmaxus    ,kmax      ,lstsci    ,ltur      , &
                & s1        ,u1        ,v1        ,r1        ,rtur1     , &
                & umnldf    ,vmnldf    ,kfu       ,kfv       , &
                & dp        ,gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011.                                     
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!!--description-----------------------------------------------------------------
!
!    Function: Reads initial field condition records from an
!              unformatted (single precision) restart file
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use dfparall
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !
    real(fp)              , pointer :: tstart
    integer               , pointer :: julday
!
! Global variables
!
    integer                                                                    , pointer     :: mfg
    integer                                                                    , pointer     :: mlg
    integer                                                                    , pointer     :: nfg
    integer                                                                    , pointer     :: nlg
    integer                                                                    , pointer     :: mmaxgl
    integer                                                                    , pointer     :: nmaxgl
    integer                                                                    , intent(in)  :: kmax   !  Description and declaration in iidim.f90
    integer                                                                    , intent(in)  :: lstsci !  Description and declaration in iidim.f90
    integer                                                                    , intent(in)  :: ltur   !  Description and declaration in iidim.f90
    integer                                                                    , intent(out) :: lturi  !  Description and declaration in tricom.igs
    integer                                                                                  :: lundia !  Description and declaration in inout.igs
    integer                                                                    , intent(in)  :: mmax   !  Description and declaration in iidim.f90
    integer                                                                                  :: nmax   !  Description and declaration in iidim.f90
    integer                                                                    , intent(in)  :: nmaxus !  Description and declaration in iidim.f90
    integer , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                            :: kfu    !  Description and declaration in iidim.f90
    integer , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                            :: kfv    !  Description and declaration in iidim.f90
    logical                                                                                  :: error  !!  Flag=TRUE if an error is encountered
    real(fp), dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(out) :: dp     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(out) :: s1     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(out) :: umnldf !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)              , intent(out) :: vmnldf !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub, 0:kmax, ltur), intent(out) :: rtur1  !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub, kmax)        , intent(out) :: u1     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub, kmax)        , intent(out) :: v1     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub, kmax, lstsci), intent(out) :: r1     !  Description and declaration in rjdim.f90
    character(*)                                                                             :: restid !!  Run identification of the restart file. If RESTID = non-blank then current simulation will use this file for setting the initial conditions
!
! Local variables
!
    integer                                              :: idate
    integer                                              :: iocond  ! IO status for reading 
    integer                                              :: itime
    integer                                              :: k       ! Help var. 
    integer                                              :: l       ! Help var. 
    integer                                              :: lengl
    integer                                              :: lenlo
    integer                                              :: lfile   ! Help var. specifying the length of character variables for files 
    integer                                              :: lrid    ! Help var., length of restid 
    integer                                              :: luntmp  ! Unit number file 
    integer                                              :: m       ! Help var. 
    integer                                              :: n       ! Help var. 
    integer                                              :: newlun
    integer       , dimension(4,0:nproc-1)               :: iarrc   ! array containing collected grid indices
    integer       , dimension(0:nproc-1)                 :: mf      ! first index w.r.t. global grid in x-direction
    integer       , dimension(0:nproc-1)                 :: ml      ! last index w.r.t. global grid in x-direction
    integer       , dimension(0:nproc-1)                 :: nf      ! first index w.r.t. global grid in y-direction
    integer       , dimension(0:nproc-1)                 :: nl      ! last index w.r.t. global grid in y-direction
    logical                                              :: ex
    logical                                              :: ex_nfs
    real(sp)      , dimension(:,:,:,:)     , allocatable :: sbuff   !!  Single precision buffer to read from file
    real(sp)      , dimension(:,:)         , allocatable :: sbuffg
    character(256)                                       :: message
    character(300)                                       :: filtmp  ! File name restart file 300 = 256 + a bit
!
!! executable statements -------------------------------------------------------
!
    julday      => gdp%gdinttim%julday
    tstart      => gdp%gdexttim%tstart
    !
    mfg         => gdp%gdparall%mfg
    mlg         => gdp%gdparall%mlg
    nfg         => gdp%gdparall%nfg
    nlg         => gdp%gdparall%nlg
    mmaxgl      => gdp%gdparall%mmaxgl
    nmaxgl      => gdp%gdparall%nmaxgl
    !
    error = .false.
    !
    ! test file existence, first 'tri-rst.<restid>.idate.itime'
    !
    write(lundia, '(a)') '*** Start of restart messages'
    !
    call timdat(julday    ,tstart*60 ,idate     ,itime     )
    call noextspaces(restid    ,lrid      )
    filtmp = 'tri-rst.' // restid(1:lrid)
    write (filtmp(8 + lrid + 1:8 + lrid + 16), '(a1,i8.8,a1,i6.6)') &
        & '.', idate, '.', itime
    lfile = 8 + lrid + 16
    inquire (file = filtmp(1:lfile), exist = ex)
    ex_nfs = .false.
    if (.not.ex) then
       !
       ! test file existance, second try 'tri-rst.<restid>'
       !
       inquire (file = filtmp(1:8 + lrid), exist = ex)
       if (.not.ex) then
          !
          ! none of these two files exists
          ! check new option: it may be a reference to a TRIM file.
          !
          call flow_nefis_restart(lundia    ,error     ,restid    ,lturi     ,mmax      , &
                                & nmaxus    ,kmax      ,lstsci    ,ltur      , &
                                & s1        ,u1        ,v1        ,r1        ,rtur1     , &
                                & umnldf    ,vmnldf    ,kfu       ,kfv       , &
                                & dp        ,ex_nfs    ,gdp       )
          if (error .and. .not.ex_nfs) then
             call prterr(lundia    ,'G004'    , &
             & filtmp(:lfile) // ', ' // filtmp(:8 + lrid) // ' and ' // filtmp(9:8+lrid) // '.dat/.def', &
             & gdp       )
          endif
       else
          lfile = 8 + lrid
       endif
    endif
    if (.not.error .and. .not.ex_nfs) then
       !
       ! restart file found
       !
       write(lundia, '(a)') 'Restarting from ' // filtmp(:lfile)
       !
       ! Allocate temporary single precision array for the ENTIRE domain
       !
       !      allocate(sbuff(nmaxus, mmax, 0:kmax, max(1, lstsci, ltur)))
       allocate(sbuff(nmaxgl, mmaxgl, 0:kmax, max(1, lstsci, ltur)))
       !
       ! the restart file is opened and read by the master
       !
       if (inode == master) then
          luntmp = newlun(gdp)
          open (luntmp, file = filtmp(1:lfile), form = 'unformatted',              &
               & status = 'old')
       endif
       !
       ! gather grid indices of all subdomains
       !
       call dfgather_grddim(lundia, nfg, nlg, mfg, mlg, nmaxgl, mmaxgl, &
          &                 nf, nl, mf, ml, iarrc, lengl, lenlo, gdp )
       call dfbroadc ( iarrc, 4*nproc, dfint, gdp )
       call dfbroadc ( nf, nproc, dfint, gdp )
       call dfbroadc ( nl, nproc, dfint, gdp )
       call dfbroadc ( mf, nproc, dfint, gdp )
       call dfbroadc ( ml, nproc, dfint, gdp )
       ! read restart values and distributed via a broadcast to the slaves
       ! per nmaxus mmax values in s1 array
       ! NOTE: nmaxus and mmax equal nmaxgl and mmaxgl, respectively (for entire domain)
       !       in case of parallel runs. Moreover, buffer is associated with entire domain
       !
       if (inode == master) then
          read (luntmp, iostat = iocond) ((sbuff(n, m, 1, 1), m = 1, mmaxgl), n = 1, nmaxgl)
       endif
       call dfbroadc(iocond, 1, dfint, gdp)
       if (iocond /= 0) then
          if (iocond < 0) then
             call prterr(lundia    ,'G006'    ,filtmp(:lfile)       )
          else
             call prterr(lundia    ,'G007'    ,filtmp(:lfile)       )
          endif
          error = .true.
          goto 200
       endif
       !
       ! send buffer to other nodes
       !
       call dfbroadc(sbuff, mmaxgl*nmaxgl*(kmax+1)*max(1,lstsci,ltur), dfreal, gdp)
       do m = 1, mmaxgl
         do n = 1, nmaxgl
            if (isnan(sbuff(n,m,1,1))) then
               write(message,*)'NaN found in restart file for s1 at (n,m)= (',n,',',m,')'
               call prterr(lundia, 'U021', trim(message))
               call d3stop(1, gdp)
            endif
         enddo
       enddo
       !
       ! put copies of parts of s1 for each subdomain
       !
       if (parll) then
          do m = mf(inode-1), ml(inode-1)
             do n = nf(inode-1), nl(inode-1)
                s1(n-nfg+1,m-mfg+1) = sbuff(n,m,1,1)
             enddo
          enddo
       else
          s1(1:nmaxus,1:mmax) = sbuff(1:nmaxus,1:mmax,1,1)
       endif
       call dfexchg ( s1, 1, 1, dfloat, gdp )
       !
       ! per layer k: nmaxus mmax values in u1 array
       !
       do k = 1, kmax
          if (inode == master) then
             read (luntmp, iostat = iocond) ((sbuff(n, m, k, 1), m = 1, mmaxgl), n = 1,nmaxgl)
          endif
          call dfbroadc(iocond, 1, dfint, gdp)
          if (iocond /= 0) then
             if (iocond < 0) then
                call prterr(lundia    ,'G006'    ,filtmp(:lfile)       )
             else
                call prterr(lundia    ,'G007'    ,filtmp(:lfile)       )
             endif
             error = .true.
             goto 200
          endif
       enddo
       !
       ! send buffer to other nodes
       !
       call dfbroadc(sbuff, mmaxgl*nmaxgl*(kmax+1)*max(1,lstsci,ltur), dfreal, gdp)
       do k = 1, kmax
          do m = 1, mmaxgl
            do n = 1, nmaxgl
               if (isnan(sbuff(n,m,k,1))) then
                  write(message,*)'NaN found in restart file for u1 at (n,m,k)= (',n,',',m,',',k,')'
                  call prterr(lundia, 'U021', trim(message))
                  call d3stop(1, gdp)
               endif
            enddo
          enddo
       enddo
       !
       ! put copies of parts of u1 for each subdomain
       !
       if (parll) then
          do m = mf(inode-1), ml(inode-1)
             do n = nf(inode-1), nl(inode-1)
                u1(n-nfg+1,m-mfg+1,1:kmax) = sbuff(n,m,1:kmax,1)
             enddo
          enddo
       else
          u1(1:nmaxus,1:mmax,1:kmax) = sbuff(1:nmaxus,1:mmax,1:kmax,1)
       endif
       call dfexchg ( u1, 1, kmax, dfloat, gdp )
       !
       ! per layer k: nmaxus mmax values in v1 array
       !
       do k = 1, kmax
          if (inode == master) then
             read (luntmp, iostat = iocond) ((sbuff(n, m, k, 1), m = 1, mmaxgl), n = 1, nmaxgl)
          endif
          call dfbroadc(iocond, 1, dfint, gdp)
          if (iocond /= 0) then
             if (iocond < 0) then
                call prterr(lundia    ,'G006'    ,filtmp(:lfile)       )
             else
                call prterr(lundia    ,'G007'    ,filtmp(:lfile)       )
             endif
             error = .true.
             goto 200
          endif
       enddo
       !
       ! send buffer to other nodes
       !
       call dfbroadc(sbuff, mmaxgl*nmaxgl*(kmax+1)*max(1,lstsci,ltur), dfreal, gdp)
       do k = 1, kmax
          do m = 1, mmaxgl
            do n = 1, nmaxgl
               if (isnan(sbuff(n,m,k,1))) then
                  write(message,*)'NaN found in restart file for v1 at (n,m,k)= (',n,',',m,',',k,')'
                  call prterr(lundia, 'U021', trim(message))
                  call d3stop(1, gdp)
               endif
            enddo
          enddo
       enddo
       !
       ! put copies of parts of v1 for each subdomain
       !
       if (parll) then
          do m = mf(inode-1), ml(inode-1)
             do n = nf(inode-1), nl(inode-1)
                v1(n-nfg+1,m-mfg+1,1:kmax) = sbuff(n,m,1:kmax,1)
             enddo
          enddo
       else
          v1(1:nmaxus,1:mmax,1:kmax) = sbuff(1:nmaxus,1:mmax,1:kmax,1)
       endif
       call dfexchg ( v1, 1, kmax, dfloat, gdp )
       !
       ! per constituent l: kmax nmaxus mmax values in r1 array
       ! only Salinity, Temperature, real constituents and secondary
       ! flow; no turbulence
       !
       if (lstsci > 0) then
          do l = 1, lstsci
             do k = 1, kmax
                if (inode==master) then
                   read (luntmp, iostat = iocond) ((sbuff(n, m, k, l), m = 1, mmaxgl), n = 1, nmaxgl)
                endif
                call dfbroadc(iocond, 1, dfint, gdp)
                if (iocond /= 0) then
                   if (iocond < 0) then
                      call prterr(lundia    ,'G006'    ,filtmp(:lfile)       )
                   else
                      call prterr(lundia    ,'G007'    ,filtmp(:lfile)       )
                   endif
                   error = .true.
                   goto 200
                endif
             enddo
          enddo
          !
          ! send buffer to other nodes
          !
          call dfbroadc(sbuff, mmaxgl*nmaxgl*(kmax+1)*max(1,lstsci,ltur), dfreal, gdp)
          do l = 1, lstsci
             do k = 1, kmax
                do m = 1, mmaxgl
                  do n = 1, nmaxgl
                     if (isnan(sbuff(n,m,k,l))) then
                        write(message,*)'NaN found in restart file for r1 at (n,m,k,l)= (',n,',',m,',',k,',',l,')'
                        call prterr(lundia, 'U021', trim(message))
                        call d3stop(1, gdp)
                     endif
                  enddo
                enddo
             enddo
          enddo
          !
          ! put copies of parts of r1 for each subdomain
          !
          if (parll) then
             do m = mf(inode-1), ml(inode-1)
                do n = nf(inode-1), nl(inode-1)
                   r1(n-nfg+1,m-mfg+1,1:kmax,1:lstsci) = sbuff(n,m,1:kmax,1:lstsci)
                enddo
             enddo
             do l = 1, lstsci
                call dfexchg ( r1(:,:,:,l), 1, kmax, dfloat, gdp )
             enddo
          else
             r1(1:nmaxus,1:mmax,1:kmax,1:lstsci) = sbuff(1:nmaxus,1:mmax,1:kmax,1:lstsci)
          endif
       endif
       !
       ! Per turbulence l: 0:kmax nmaxus mmax values in rtur1 array
       ! If no turbulence arrays on restart file then rtur1 will be
       ! initialized in INITUR
       ! If only K on restart file EPS will be calculated in INITUR
       !
       if (ltur>0) then
          lturi = 0
          do l = 1, ltur
             do k = 0, kmax
                if (inode==master) then
                   read (luntmp, iostat = iocond) ((sbuff(n, m, k, l), m = 1, mmaxgl) , n = 1, nmaxgl)
                endif
                call dfbroadc(iocond, 1, dfint, gdp)
                if (iocond /= 0) then
                   if (iocond < 0) then
                      lturi = ltur
                      if (l==2) lturi = -ltur
                   else
                      call prterr(lundia    ,'G007'    ,filtmp(:lfile)       )
                      error = .true.
                   endif
                   goto 200
                endif
             enddo
          enddo
          !
          ! send buffer to other nodes
          !
          call dfbroadc(sbuff, mmaxgl*nmaxgl*(kmax+1)*max(1,lstsci,ltur), dfreal, gdp)
          do l = 1, ltur
             do k = 1, kmax
                do m = 1, mmaxgl
                  do n = 1, nmaxgl
                     if (isnan(sbuff(n,m,k,l))) then
                        write(message,*)'NaN found in restart file for rtur1 at (n,m,k,l)= (',n,',',m,',',k,',',l,')'
                        call prterr(lundia, 'U021', trim(message))
                        call d3stop(1, gdp)
                     endif
                  enddo
                enddo
             enddo
          enddo
          !
          ! put copies of parts of rtur1 for each subdomain
          !
          if (parll) then
             do m = mf(inode-1), ml(inode-1)
                do n = nf(inode-1), nl(inode-1)
                   rtur1(n-nfg+1,m-mfg+1,0:kmax,1:ltur) = sbuff(n,m,0:kmax,1:ltur)
                enddo
             enddo
             do l = 1, ltur
                call dfexchg ( rtur1(:,:,:,l), 0, kmax, dfloat, gdp )
             enddo
          else
             rtur1(1:nmaxus,1:mmax,0:kmax,1:ltur) = sbuff(1:nmaxus,1:mmax,0:kmax,1:ltur)
          endif
       endif
       !
       ! read filtered velocity components to allow restarts
       ! using subgrid viscosity model
       !
       if (inode == master) then
          read (luntmp, iostat = iocond, end = 200) ((sbuff(n, m, 1, 1), m = 1, mmaxgl), n = 1, nmaxgl)
       endif
       !
       ! send buffer to other nodes
       !
       call dfbroadc(sbuff, mmaxgl*nmaxgl*(kmax+1)*max(1,lstsci,ltur), dfreal, gdp)
       do m = 1, mmaxgl
         do n = 1, nmaxgl
            if (isnan(sbuff(n,m,1,1))) then
               write(message,*)'NaN found in restart file for umnldf at (n,m)= (',n,',',m,')'
               call prterr(lundia, 'U021', trim(message))
               call d3stop(1, gdp)
            endif
         enddo
       enddo
       !
       ! put copies of parts of umnldf for each subdomain
       !
       if (parll) then
          do m = mf(inode-1), ml(inode-1)
             do n = nf(inode-1), nl(inode-1)
                umnldf(n-nfg+1,m-mfg+1) = sbuff(n,m,1,1)
             enddo
          enddo
       else
          umnldf(1:nmaxus,1:mmax) = sbuff(1:nmaxus,1:mmax,1,1)
       endif
       call dfexchg ( umnldf, 1, 1, dfloat, gdp )
       !
       if (inode == master) then
          read (luntmp, iostat = iocond) ((sbuff(n, m, 1, 1), m = 1, mmaxgl), n = 1, nmaxgl)
       endif
       !
       ! send buffer to other nodes
       !
       call dfbroadc(sbuff, mmaxgl*nmaxgl*(kmax+1)*max(1,lstsci,ltur), dfreal, gdp)
       do m = 1, mmaxgl
         do n = 1, nmaxgl
            if (isnan(sbuff(n,m,1,1))) then
               write(message,*)'NaN found in restart file for vmnldf at (n,m)= (',n,',',m,')'
               call prterr(lundia, 'U021', trim(message))
               call d3stop(1, gdp)
            endif
         enddo
       enddo
       !
       ! put copies of parts of vmnldf for each subdomain
       !
       if (parll) then
          do m = mf(inode-1), ml(inode-1)
             do n = nf(inode-1), nl(inode-1)
                vmnldf(n-nfg+1,m-mfg+1) = sbuff(n,m,1,1)
             enddo
          enddo
       else
          vmnldf(1:nmaxus,1:mmax) = sbuff(1:nmaxus,1:mmax,1,1)
       endif
       call dfexchg ( vmnldf, 1, 1, dfloat, gdp )
       !
       ! stop reading file
       !
       !
       ! close file
       !
  200  continue
       if (inode == master) close (luntmp)
       deallocate(sbuff)
    endif
    write (lundia, '(a)') '*** End   of restart messages'
    write (lundia, *)
end subroutine rstfil

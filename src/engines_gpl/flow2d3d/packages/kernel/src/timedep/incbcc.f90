subroutine incbcc(lundia    ,timnow    ,zmodel    ,nmax      ,mmax      , &
                & kmax      ,nto       ,nrob      ,lstsc     ,noroco    , &
                & tprofc    ,itbcc     ,mnbnd     ,nob       ,kstp      , &
                & kfsmin    ,kfsmax    ,rob       ,rbnd      ,guu       , &
                & gvv       ,dps       ,s0        ,sig       ,procbc    , &
                & zstep     ,dzs1      ,zk        ,gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2012.                                
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
!  $Id$
!  $HeadURL$
!!--description-----------------------------------------------------------------
!
!    Function: Carry out interpolation in space, determine time
!              increments and updates the time dependent const.
!              BC (if NTO > 0 and LSTSC > 0)
! Method used: - At each time step the increment values are added
!                to update the values at the vertices of the
!                boundaries
!              - Hereafter space interpolation is applied to cal-
!                culate the boundary values at each point.
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use dfparall
    use globaldata
    use m_openda_exchange_items, only : get_openda_buffer
    !
    implicit none
    !
    ! Enumeration
    !
    integer, parameter :: start_pivot = 1
    integer, parameter ::   end_pivot = 2
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    integer                 , pointer :: ltem
    integer                 , pointer :: lsal
    real(fp)                , pointer :: tstop
    integer                 , pointer :: itfinish
    integer                 , pointer :: lunbcc
    logical                 , pointer :: newstp
    real(fp), dimension(:,:), pointer :: dist_pivot_part
!
! Global variables
!
    integer                                                                                :: kmax   !  Description and declaration in esm_alloc_int.f90
    integer                                                                                :: lstsc  !  Description and declaration in dimens.igs
    integer                                                                                :: lundia !  Description and declaration in inout.igs
    integer                                                                                :: mmax   !  Description and declaration in esm_alloc_int.f90
    integer                                                                                :: nmax   !  Description and declaration in esm_alloc_int.f90
    integer                                                                  , intent(in)  :: noroco !  Description and declaration in esm_alloc_int.f90
    integer                                                                  , intent(in)  :: nrob   !  Description and declaration in esm_alloc_int.f90
    integer                                                                                :: nto    !  Description and declaration in esm_alloc_int.f90
    integer      , dimension(7, nto)                                         , intent(in)  :: mnbnd  !  Description and declaration in esm_alloc_int.f90
    integer      , dimension(5, nto, lstsc)                                                :: itbcc  !  Description and declaration in esm_alloc_int.f90
    integer      , dimension(8, nrob)                                        , intent(in)  :: nob    !  Description and declaration in esm_alloc_int.f90
    integer      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: kfsmax !  Description and declaration in esm_alloc_int.f90
    integer      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)       , intent(in)  :: kfsmin !  Description and declaration in esm_alloc_int.f90
    integer      , dimension(nrob, lstsc)                                                  :: kstp   !  Description and declaration in esm_alloc_int.f90
    logical                                                                  , intent(in)  :: zmodel !  Description and declaration in procs.igs
    real(fp)                                                                               :: timnow !!  Current timestep (multiples of dt)
    real(fp)     , dimension(2, nto, lstsc)                                                :: zstep  !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(4, nto, kmax, lstsc)                                          :: procbc !  Description and declaration in esm_alloc_real.f90
    real(prec)   , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)       , intent(in)  :: dps    !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)       , intent(in)  :: s0     !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)       , intent(in)  :: guu    !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)       , intent(in)  :: gvv    !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub, kmax) , intent(in)  :: dzs1   !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(kmax)                                           , intent(in)  :: sig    !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(0:kmax)                                         , intent(in)  :: zk
    real(fp)     , dimension(kmax, max(lstsc, 1), 2, noroco)                               :: rbnd   !  Description and declaration in esm_alloc_real.f90
    real(fp)     , dimension(kmax, nrob)                                                   :: rob    !  Description and declaration in esm_alloc_real.f90
    character(10), dimension(nto, lstsc)                                                   :: tprofc !  Description and declaration in esm_alloc_char.f90
!
! Local variables
!
    integer                                      :: i
    integer                                      :: incx    ! Nr. of grid points-1 (in the x-dir.) between the begin and the end point of an opening section 
    integer                                      :: incy    ! Nr. of grid points-1 (in the y-dir.) between the begin and the end point of an opening section 
    integer                                      :: istsc   ! Index number of constituent 
    integer                                      :: ito     ! Index number of open boundary loc. 
    integer                                      :: j       ! Loop variable 
    integer                                      :: k       ! Loop counter over KMAX 
    integer                                      :: kp      ! First array index of array CIRC pointing to the nr. of row/column in array IROCOL
    integer                                      :: kpc     ! First array index of array CIRC pointing to the column number in arra IROCOL 
    integer                                      :: kpp     ! Hulp varible 
    integer                                      :: kpr     ! First array index of array CIRC pointing to the row number in array IROCOL 
    integer                                      :: kq      ! Second array index of array CIRC pointing to the nr. of row/column in array IROCOL 
    integer                                      :: kqc     ! Second array index of array CIRC pointing to the column number in arra IROCOL 
    integer                                      :: kqq     ! Hulp varible 
    integer                                      :: kqr     ! Second array index of array CIRC pointing to the row number in array IROCOL 
    integer                                      :: l       ! Loop variable for constituents 
    integer                                      :: maxinc  ! Max. of (INCX,INCY,1) 
    integer                                      :: mend    ! End coord. (in the x-dir.) of an open bound. section 
    integer                                      :: mp
    integer                                      :: msta    ! Starting coord. (in the x-dir.) of an open bound. section 
    integer                                      :: n       ! Loop variable 
    integer                                      :: n1      ! Pointer var. relating NOB to MNBND 
    integer                                      :: nend    ! End coord. (in the y-dir.) of an open bound. section 
    integer                                      :: np
    integer                                      :: nsta    ! Starting coord. (in the y-dir.) of an open bound. section
    integer                                      :: lb      ! lowerboundary of loopcounter
    integer                                      :: ub      ! upperboundary of loopcounter
    integer                                      :: mgg     ! M-coord. of the actual open boundary point, which may differ from the ori- ginal position due to grid staggering
    integer                                      :: ngg     ! N-coord. of the actual open boundary point, which may differ from the ori- ginal position due to grid staggering
    integer                                      :: posrel  ! code denoting the position of the open boundary, related to the complete grid
    logical                                      :: first   ! Flag = TRUE in case a time-dependent file is read for the 1-st time 
    real(fp)                                     :: dist    ! Real distance between an open bounda- ry point to the begin point of the related opening section 
    real(fp)                                     :: distx   ! Incremental distance (in the x-dir.) between two consecutive open boundary points belonging to the same section 
    real(fp)                                     :: disty   ! Incremental distance (in the x-dir.) between two consecutive open boundary points belonging to the same section 
    real(fp)                                     :: frac    ! Fraction between DIST and the total length of an opening section 
    real(fp)                                     :: h0
    real(fp)                                     :: reldep
    real(fp)                                     :: sigjmp
    real(fp)                                     :: timscl  ! Multiple factor to create minutes from read times 
    real(fp)                                     :: totl    ! Actual length of an openbnd. section 
    real(fp)                                     :: zcord
    real(fp)                                     :: zjmp
    real(fp), dimension(:, :, :, :), allocatable :: procbc_old
!
!! executable statements -------------------------------------------------------
!
    ltem            => gdp%d%ltem
    lsal            => gdp%d%lsal
    newstp          => gdp%gdincbcc%newstp
    lunbcc          => gdp%gdluntmp%lunbcc
    itfinish        => gdp%gdinttim%itfinish
    tstop           => gdp%gdexttim%tstop
    dist_pivot_part => gdp%gdbcdat%dist_pivot_part
    if (lsal .ne. 0 .or. ltem .ne. 0) then
       allocate(procbc_old(2, nto, kmax, max(lsal,ltem)))
    endif
    !
    ! Initialization local parameters
    ! TIMSCL will not been used in UPDBCC
    !
    first  = .false.
    timscl = 1.0
    !
    ! Update values if necessary
    !
    do ito = 1, nto
       do istsc = 1, lstsc
          !
          ! When using several instances, it is necessary to do a partial rewind of the file since
          ! a newly selected instance has an older time stamp. 
          ! Luckily, this file (see updbcc) can be read by record.
          ! So we reset the record pointer to the start record if timnow < itbcc(1,..)
          !
          if (timnow < real(itbcc(1,ito,istsc),fp)) then
              itbcc(5, ito, istsc) = itbcc(3, ito, istsc) - 1

             first = .true.
             call updbcc(lunbcc    ,lundia    ,first     ,itbcc     ,ito       , &
                       & istsc     ,timnow    ,itfinish  ,timscl    , &
                       & nto       ,kmax      ,lstsc     ,procbc    ,tprofc    , &
                       & zstep     ,gdp       )
             !          
             ! Now, we have to update procbc(1,:) . Its current value is still the value in the
             ! file at itbcc(1,:)
             !
             do k = 1, kmax
                procbc(1,ito,k,istsc) = procbc(1,ito,k,istsc) +                                                 &
                                      & 2.0_fp*(timnow-0.5_fp - real(itbcc(1,ito,istsc)))*procbc(3,ito,k,istsc) 
                procbc(2,ito,k,istsc) = procbc(2,ito,k,istsc) +                                                 &
                                      & 2.0_fp*(timnow-0.5_fp - real(itbcc(1,ito,istsc)))*procbc(4,ito,k,istsc)               
             enddo 
          endif
          !
          if (timnow > real(itbcc(2, ito, istsc),fp)) then
             zstep(1, ito, istsc) = zstep(2, ito, istsc)
             itbcc(1, ito, istsc) = itbcc(2, ito, istsc)
             !
             ! Read new time step and concentrations data
             !
             call updbcc(lunbcc    ,lundia    ,first     ,itbcc     ,ito       , &
                       & istsc     ,timnow    ,itfinish  ,timscl    , &
                       & nto       ,kmax      ,lstsc     ,procbc    ,tprofc    , &
                       & zstep     ,gdp       )
             newstp = .true.
          endif
          !
          ! interpolation in time
          !
          do k = 1, kmax
             procbc(1, ito, k, istsc) = procbc(1, ito, k, istsc)                &
                                      & + procbc(3, ito, k, istsc)
             procbc(2, ito, k, istsc) = procbc(2, ito, k, istsc)                &
                                      & + procbc(4, ito, k, istsc)
          enddo
       enddo
       !
       ! Overrule procbc(1:2,:) with OpenDA input, if present
       ! this is done for salinity and temperature, if present
       ! note: store this value and reset procbc at the end of this routine.
       ! This is to prevent noise-on-noise since procbc(1:2) is never read from file again
       !
       if (lsal .ne. 0) then
          procbc_old(1:2,ito,1:kmax,lsal) = procbc(1:2,ito,1:kmax,lsal)
          call get_openda_buffer('bound_salt', ito, 2, kmax, procbc(1:2,ito,1:kmax,lsal))
       endif
       if (ltem .ne. 0) then
          procbc_old(1:2,ito,1:kmax,ltem) = procbc(1:2,ito,1:kmax,ltem)
          call get_openda_buffer('bound_temp', ito, 2, kmax, procbc(1:2,ito,1:kmax,ltem))    
       endif
    enddo
    !
    ! for all "conservative" constituents
    !
    do l = 1, lstsc
       !
       ! calculate space fraction of opening for function values in hob
       !
       do n = 1, nrob
          n1     = nob(8, n)
          msta   = mnbnd(1, n1)
          nsta   = mnbnd(2, n1)
          mend   = mnbnd(3, n1)
          nend   = mnbnd(4, n1)
          posrel = mnbnd(7, n1)
          incx   = mend - msta
          incy   = nend - nsta
          maxinc = max(1, abs(incx), abs(incy))
          incx   = incx/maxinc
          incy   = incy/maxinc
          !
          ! use coord. at zeta points to calculate the interpolated
          ! values from the two utmost points of the opening
          !
          if (parll) then
             !
             ! If the start point/pivot is outside this partition,
             ! add the distance from that point/pivot to the first point inside this partition
             ! to totl and distl
             !
             totl = dist_pivot_part(start_pivot,n1)
          else
             totl = 0.0_fp
          endif
          dist = totl
          frac  = 0.0_fp
          do j = 1, maxinc
             msta  = msta + incx
             nsta  = nsta + incy
             !
             ! Pythagoras is used to calculate the distance from xz,yz(nsta-incy,msta-incx) to xz,yz(nsta,msta),
             ! using distx = |xz,yz(nsta     ,msta-incx) - xz,yz(nsta,msta     )|
             !       disty = |xz,yz(nsta-incy,msta-incx) - xz,yz(nsta,msta-incx)|
             ! Assumption: the grid is more or less rectangular locally
             !
             if (incx == 0) then
                distx = 0.0_fp
             else
                if (posrel == 4) then
                   !
                   ! Boundary on top side of the domain
                   ! The correct gvv is one index down (staggered grid)
                   !
                   ngg = nsta - 1
                else
                   ngg = nsta
                endif
                !
                ! General case: incx > 1
                ! The distance in the x-direction between the two zeta points is:
                ! 0.5*gvv(lowest_point) + sum(gvv(intermediate_points)) + 0.5*gvv(highest_point)
                !
                distx = 0.5_fp * (real(gvv(ngg,msta-incx),fp)+real(gvv(ngg,msta),fp))
                lb = min(msta-incx,msta)
                ub = max(msta-incx,msta)
                do i=lb+1,ub-1
                   distx = distx + real(gvv(ngg,i),fp)
                enddo
             endif
             distx = distx*distx
             if (incy == 0) then
                disty = 0.0_fp
             else
                if (posrel == 3) then
                   !
                   ! Boundary on right side of the domain
                   ! The correct guu is one index down (staggered grid)
                   !
                   mgg = msta - incx - 1
                else
                   mgg = msta - incx
                endif
                !
                ! General case: incy > 1
                ! The distance in the y-direction between the two zeta points is:
                ! 0.5*guu(lowest_point) + sum(guu(intermediate_points)) + 0.5*guu(highest_point)
                !
                disty = 0.5_fp * (real(guu(nsta-incy,mgg),fp)+real(guu(nsta,mgg),fp))
                lb = min(nsta-incy,nsta)
                ub = max(nsta-incy,nsta)
                do i=lb+1,ub-1
                   disty = disty + real(guu(i,mgg),fp)
                enddo
             endif
             disty = disty*disty
             totl  = totl + sqrt(distx + disty)
             if (msta==nob(1, n) .and. nsta==nob(2, n)) then
                dist = totl
             endif
          enddo
          if (parll) then
             !
             ! If the end point/pivot is outside this partition,
             ! add the distance from that point/pivot to the last point inside this partition
             ! to totl
             !
             totl = totl + dist_pivot_part(end_pivot,n1)
          endif
          if (totl > 0.) then
             frac = dist/totl
          endif
          !
          ! Interpolate in space constituents along surface and bottom
          !
          do k = 1, kmax
             ! rob   (k,n) = VALPNT(n1,k,l,frac)
             rob(k, n) = procbc(1, n1, k, l)                                    &
                       & + frac*(procbc(2, n1, k, l) - procbc(1, n1, k, l))
          enddo
          !
          ! calculate rbnd array
          ! where nob (4,n) := sort opening and
          !       nob (5,n) := row (conform irocol table)
          !       nob (6,n) := sort opening and
          !       nob (7,n) := column (conform irocol table)
          !
          kpr = nob(4, n)
          kqr = nob(5, n)
          kpc = nob(6, n)
          kqc = nob(7, n)
          kp  = kpr
          kq  = kqr
          kpp = kpc
          kqq = kqc
          !
          ! If row boundary not available then use col boundary
          !
          if (kp*kq == 0) then
             kp  = kpc
             kq  = kqc
             kpp = kpr
             kqq = kqr
          endif
          !
          ! Interpolation in vertical direction
          !
          if (tprofc(n1, l)(1:7) == 'uniform') then
             !
             ! uniform
             ! For all KMAX layers ROB is equal (stored in ROB(1,N))
             !
             do k = 1, kmax
                rbnd(k, l, kp, kq) = rob(1, n)
             enddo
          elseif (tprofc(n1, l)(1:6) == 'linear') then
             !
             ! Linear
             ! Values of bottom and surface will be used to interpolated
             ! for all KMAX layers
             !
             if (zmodel) then
                zcord = 0.
                do k = kfsmin(nsta, msta), kmax
                   if (k==kfsmin(nsta, msta)) then
                      zcord = zcord + .5*dzs1(nsta, msta, k)
                   else
                      zcord = zcord +                                           &
                            & .5*(dzs1(nsta, msta, k - 1) + dzs1(nsta, msta, k))
                   endif
                   reldep = zcord/(real(dps(nsta, msta),fp) + s0(nsta, msta))
                   rbnd(k, l, kp, kq) = rob(kmax, n)                            &
                                      & + reldep*(rob(1, n) - rob(kmax, n))
                enddo
             else
                do k = 1, kmax
                   rbnd(k, l, kp, kq) = rob(kmax, n) + (1. + sig(k))            &
                                      & *(rob(1, n) - rob(kmax, n))
                enddo
             endif
          elseif (tprofc(n1, l)(1:4) == 'step') then
             !
             ! step function
             ! The two values of the Step function are written in
             ! ROB(1,N) and ROB(KMAX,N)
             !
             if (newstp) then
                mp     = nob(1, n)
                np     = nob(2, n)
                h0     = real(dps(np, mp),fp) + s0(np, mp)
                zjmp   = -zstep(1, n1, l)
                sigjmp = zjmp/max(h0, 0.01_fp)
                do k = 1, kmax
                   if (zmodel) then
                      if (zjmp > zk(k)) then
                         kstp(n, l) = k
                         goto 610
                      endif
                   elseif (sigjmp > sig(k)) then
                      kstp(n, l) = max(k - 1, 1)
                      goto 610
                   else
                   endif
                enddo
                kstp(n, l) = kmax
             endif
  610        continue
             do k = 1, kmax
                if (k <= kstp(n, l)) then
                   rbnd(k, l, kp, kq) = rob(1, n)
                else
                   rbnd(k, l, kp, kq) = rob(kmax, n)
                endif
             enddo
          elseif (tprofc(n1, l) == '3d-profile') then
             !
             ! 3D-profile function
             ! The user defined input is stored in ROB for all KMAX
             ! layers
             !
             do k = 1, kmax
                rbnd(k, l, kp, kq) = rob(k, n)
             enddo
          else
          endif
          !
          ! If oblique boundary is defined then fill the
          ! computed RBND also in other direction
          !
          do k = 1, kmax
             if (kpp*kqq /= 0) then
                rbnd(k, l, kpp, kqq) = rbnd(k, l, kp, kq)
             endif
          enddo
       enddo
    enddo
    newstp = .false.
    !
    ! Now reset procbc(1:2)
    !
    if (lsal .ne. 0) then
       procbc(1:2,1:nto,1:kmax,lsal) = procbc_old(1:2,1:nto,1:kmax,lsal)
    endif
    if (ltem .ne. 0) then
       procbc(1:2,1:nto,1:kmax,ltem) = procbc_old(1:2,1:nto,1:kmax,ltem)
    endif    
    if (lsal .ne. 0 .or. ltem .ne. 0) then
       deallocate(procbc_old)
    endif
end subroutine incbcc

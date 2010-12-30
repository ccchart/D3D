subroutine z_forfil(nmmax     ,kmax      ,lstsci    , &
                  & lsecfl    ,lsal      ,ltem      ,icx       ,icy       , &
                  & nst       ,forfuv    ,forfww    , &
                  & kfs       ,kcs       ,kcu       ,kcv       ,kfsmin    , &
                  & kfsmax    ,kfuz1     ,kfvz1     ,kfsz1     ,idifu     , &
                  & dzs1      ,r0        , &
                  & r1        ,rmneg     ,volum1    ,dicww     ,w1        , &
                  & sigdif    ,sigmol    ,gdp       )
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
!    Function: Filters concentrations. Numerical diffusion is
!              added in points where the concentration is nega-
!              tive (FORESTER-filter). The diffusion is two
!              dimensional along sigma-planes.
!              For the spiral motion intensity there is no
!              filtering applied
!              In the vertical direction the solution is made
!              monotone for salinity and temperature. Wiggles
!              generated by the central differences in the
!              vertical are removed.
!              First horizontal filter to remove negative values
!              then vertical to make monotone solution.
!              The end solution is physically realistic.
!              Vertical filter only for Peclet number > 2.0
! Method used: Forester-filter is used. (C.K. Forester, Higher
!              Order Monotonic Convective Difference Schemes,
!              J. Computational Phys. 23, 1977)
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    integer                , pointer :: lundia
    real(fp)               , pointer :: vicmol
!
! Global variables
!
    integer                                                   , intent(in)  :: icx    !!  Increment in the X-dir., if ICX= NMAX
                                                                                      !!  then computation proceeds in the X-
                                                                                      !!  dir. If icx=1 then computation pro-
                                                                                      !!  ceeds in the Y-dir.
    integer                                                   , intent(in)  :: icy    !!  Increment in the Y-dir. (see ICX)
    integer                                                   , intent(in)  :: kmax   !  Description and declaration in iidim.f90
    integer                                                   , intent(in)  :: lsal   !  Description and declaration in dimens.igs
    integer                                                   , intent(in)  :: lsecfl !  Description and declaration in dimens.igs
    integer                                                   , intent(in)  :: lstsci !  Description and declaration in iidim.f90
    integer                                                   , intent(in)  :: ltem   !  Description and declaration in dimens.igs
    integer                                                   , intent(in)  :: nmmax  !  Description and declaration in dimens.igs
    integer                                                   , intent(in)  :: nst    !!  Current time step number
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)                               :: idifu  !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)                 , intent(in)  :: kcs    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)                 , intent(in)  :: kcu    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)                 , intent(in)  :: kcv    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)                 , intent(in)  :: kfs    !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)                 , intent(in)  :: kfsmax !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub)                 , intent(in)  :: kfsmin !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)           , intent(in)  :: kfsz1  !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)           , intent(in)  :: kfuz1  !  Description and declaration in iidim.f90
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax)           , intent(in)  :: kfvz1  !  Description and declaration in iidim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)        , intent(in)  :: dzs1   !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)                      :: volum1 !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, 0:kmax)      , intent(in)  :: dicww  !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, 0:kmax)      , intent(in)  :: w1     !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax, lstsci)              :: r0     !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax, lstsci)              :: r1     !  Description and declaration in rjdim.f90
    real(fp)  , dimension(lstsci)                                           :: sigdif !  Description and declaration in rjdim.f90
    real(fp)  , dimension(lstsci)                             , intent(in)  :: sigmol !  Description and declaration in rjdim.f90
    real(fp)  , dimension(lstsci)                                           :: rmneg  !  Description and declaration in rjdim.f90
    character(1)                                              , intent(in)  :: forfuv !  Description and declaration in tricom.igs
    character(1)                                              , intent(in)  :: forfww !  Description and declaration in tricom.igs
!
! Local variables
!
    integer       :: ifil
    integer       :: itfil
    integer       :: m
    integer       :: n
    integer       :: k
    integer       :: kd
    integer       :: ku
    integer       :: l
    integer       :: maxfil ! Upper limit for the number of times that Forester filter for negative concentrations is applied 
    integer       :: ndm
    integer       :: nm
    integer       :: nmd
    integer       :: nmu
    integer       :: num
    real(fp)      :: coef
    real(fp)      :: cofndm ! Coefficient for NM and NDM
    real(fp)      :: cofnmd ! Coefficient for NM and NMD
    real(fp)      :: cofnmu ! Coefficient for NM and NMU
    real(fp)      :: cofnum ! Coefficient for NM and NUM
    real(fp)      :: dr
    real(fp)      :: dr1
    real(fp)      :: dr2
    real(fp)      :: dz1
    real(fp)      :: dz2
    real(fp)      :: volndm
    real(fp)      :: volnmd
    real(fp)      :: volnmu
    real(fp)      :: volnum
    real(fp)      :: peclz
    character(80) :: errtxt
!
!! executable statements -------------------------------------------------------
!
    lundia   => gdp%gdinout%lundia
    vicmol   => gdp%gdphysco%vicmol
    !
    ! Forester filter for UV only if FORFUV = 'Y'
    !
    if (forfuv=='Y') then
       !
       ! define value for 'smoothing'.
       ! for most constituents -1.e-2 is small enough
       !
       do l = 1, lstsci
          rmneg(l) = -1.0e-2_fp
       enddo
       maxfil = 100
       if (lsecfl>0) rmneg(lsecfl) = -1.0e-8_fp
       !
       ! loop over the constituents
       !
       do l = 1, lstsci
          !
          ! iteration loop over computational grid
          !
          do k = 1, kmax
             !
             ! Forester filter initialisation
             !
             do nm = 1, nmmax
                idifu(nm) = 0
             enddo
             !
             ! Forester filter number of iterations
             !
             do itfil = 1, maxfil
                ifil = 0
                !
                do nm = 1, nmmax
                   r0(nm, k, l) = r1(nm, k, l)
                   if (r0(nm, k, l)<rmneg(l) .and. kcs(nm)*kfsz1(nm, k)==1) then
                      idifu(nm) = 1
                      ifil = 1
                   endif
                enddo
                if (ifil==0) then
                   exit
                endif
                nmd = -icx
                ndm = -icy
                nmu = icx
                num = icy
                do nm = 1, nmmax
                   nmd = nmd + 1
                   ndm = ndm + 1
                   nmu = nmu + 1
                   num = num + 1
                   if (kcs(nm)*kfsz1(nm, k)==1) then
                      volnmu = min(1.0_fp, volum1(nmu, k)/volum1(nm, k))
                      cofnmu = 0.125*(idifu(nmu) + idifu(nm))*kfuz1(nm,k)*volnmu
                      volnmd = min(1.0_fp, volum1(nmd, k)/volum1(nm, k))
                      cofnmd = 0.125*(idifu(nmd) + idifu(nm))*kfuz1(nmd,k)*volnmd
                      volnum = min(1.0_fp, volum1(num, k)/volum1(nm, k))
                      cofnum = 0.125*(idifu(num) + idifu(nm))*kfvz1(nm,k)*volnum
                      volndm = min(1.0_fp, volum1(ndm, k)/volum1(nm, k))
                      cofndm = 0.125*(idifu(ndm) + idifu(nm))*kfvz1(ndm,k)*volndm
                      !
                      ! corrections for subdomain interfaces:
                      !
                      if (kcu(nm) ==3) cofnmu = 0.0
                      if (kcu(nmd)==3) cofnmd = 0.0
                      if (kcv(nm) ==3) cofnum = 0.0
                      if (kcv(ndm)==3) cofndm = 0.0
                      !
                      r1(nm, k, l) = r0(nm, k, l)                               &
                                   & *(1 - cofnmu - cofnmd - cofndm - cofnum)   &
                                   & + r0(nmu, k, l)*cofnmu + r0(nmd, k, l)     &
                                   & *cofnmd + r0(num, k, l)                    &
                                   & *cofnum + r0(ndm, k, l)*cofndm
                   endif ! computational wet points only
                enddo    ! nm-loop
                !
                ! test if number of iteration steps for filtering is
                ! exceeded
                !
                if (itfil==maxfil) then
                   write (errtxt, '(a,i2,i3,i9)') &
                        & 'Negative concentrations for  l, k, nst:',l, k, nst
                   call prterr(lundia    ,'U190'    ,errtxt    )
                endif
             enddo  ! itfil-loop negative concentrations
          enddo     ! k-loop
       enddo        ! l-loop
    endif           ! forfuv
    !
    ! Forester filter for Z-direction only if FORFWW = 'Y'
    !
    if (forfww=='Y') then
       do l = 1, lstsci
          !
          ! Vertical Forester filter only for salt and temperature
          !
          if (l/=lsal .and. l/=ltem) then
             cycle
          endif
          !
          ! iteration loop over computational grid
          !
          do nm = 1, nmmax
             if (kfs(nm)*kcs(nm)==1) then
                !
                ! Filter for vertical wiggles
                !
                rmneg(l) = 1.0e-6_fp
                maxfil   = 1000
                do itfil = 1, maxfil
                   ifil = 0
                   do k = kfsmin(nm) + 1, kfsmax(nm) - 1
                      kd = max(kfsmin(nm), k - 1)
                      ku = min(kfsmax(nm), k + 1)
                      if (r1(nm, k, l)>max(r1(nm, kd, l), r1(nm, ku, l)) + rmneg(l)) then
                         !
                         ! local maximum; filter applied for layers > 1 mm
                         !
                         dr1 = r1(nm, k, l) - r1(nm, kd, l)
                         dr2 = r1(nm, k, l) - r1(nm, ku, l)
                         if (dr2>dr1) then
                            dr    = min(dr1, 0.5*dr2)
                            dz1   = dzs1(nm, k)
                            dz2   = dzs1(nm, ku)
                            coef  = min(dz1, dz2, 1.0_fp)
                            peclz = abs(w1(nm,k)) * coef                         &
                                  & / (dicww(nm,k)/sigdif(l) + vicmol/sigmol(l))
                            !
                            !  only filter if Peclet number > 2 (numerical generated wiggle)
                            !
                            if (coef > 0.001_fp .and. peclz > 2.0_fp) then
                               ifil = 1
                               r1(nm, k , l) = r1(nm, k , l) - coef*dr/dz1
                               r1(nm, ku, l) = r1(nm, ku, l) + coef*dr/dz2
                            endif
                         else
                            dr = min(0.5*dr1, dr2)
                            dz1 = dzs1(nm, k)
                            dz2 = dzs1(nm, kd)
                            coef = min(dz1, dz2, 1.0_fp)
                            peclz = abs(w1(nm,k)) * coef                         &
                                  & / (dicww(nm,k)/sigdif(l) + vicmol/sigmol(l))
                            !
                            !  only filter if Peclet number > 2 (numerical generated wiggle)
                            !
                            if (coef > 0.001_fp .and. peclz > 2.0_fp  ) then
                               ifil = 1
                               r1(nm, k , l) = r1(nm, k , l) - coef*dr/dz1
                               r1(nm, kd, l) = r1(nm, kd, l) + coef*dr/dz2
                            endif
                         endif
                      elseif (r1(nm, k, l)<min(r1(nm, kd, l), r1(nm, ku, l)) - rmneg(l)) then
                         !
                         ! local minimum; filter applied for layers > 1 mm
                         !
                         dr1 = r1(nm, k, l) - r1(nm, kd, l)
                         dr2 = r1(nm, k, l) - r1(nm, ku, l)
                         if (dr2<dr1) then
                            dr = max(dr1, 0.5*dr2)
                            dz1 = dzs1(nm, k)
                            dz2 = dzs1(nm, ku)
                            coef = min(dz1, dz2, 1.0_fp)
                            peclz = abs(w1(nm,k)) * coef                         &
                                  & / (dicww(nm,k)/sigdif(l) + vicmol/sigmol(l))
                            !
                            !  only filter if Peclet number > 2 (numerical generated wiggle)
                            !
                            if (coef > 0.001_fp.and. peclz > 2.0_fp) then
                               ifil = 1
                               r1(nm, k , l) = r1(nm, k , l) - coef*dr/dz1
                               r1(nm, ku, l) = r1(nm, ku, l) + coef*dr/dz2
                            endif
                         else
                            dr = max(0.5*dr1, dr2)
                            dz1 = dzs1(nm, k)
                            dz2 = dzs1(nm, kd)
                            coef = min(dz1, dz2, 1.0_fp)
                            peclz = abs(w1(nm,k)) * coef                         &
                                  & / (dicww(nm,k)/sigdif(l) + vicmol/sigmol(l))
                            !
                            !  only filter if Peclet number > 2 (numerical generated wiggle)
                            !
                            if (coef > 0.001_fp.and. peclz > 2.0_fp) then
                               ifil = 1
                               r1(nm, k , l) = r1(nm, k , l) - coef*dr/dz1
                               r1(nm, kd, l) = r1(nm, kd, l) + coef*dr/dz2
                            endif
                         endif
                      else
                      endif
                   enddo ! k-loop
                   if (ifil==0) then
                      exit
                   endif
                   !
                   ! test if number of iteration steps for vertical filtering is
                   ! exceeded
                   !
                   if (itfil==maxfil) then
                      m = nm/icx + 1
                      n = (nm - (m-1)*icx)/icy
                      write (errtxt, '(a,i2,i4,i4,i9)') &
                           & 'Vertical wiggle for  l, m ,n, nst:', l, m, n, nst
                      call prterr(lundia    ,'U190'    ,errtxt    )
                   endif
                enddo    ! itfil-loop wiggle in vertical for constituents
             endif    ! computational wet points
          enddo       ! nm-loop
       enddo          ! l-loop
    endif             ! forfww
end subroutine z_forfil

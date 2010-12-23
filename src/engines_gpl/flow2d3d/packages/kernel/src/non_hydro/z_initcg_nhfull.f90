subroutine z_initcg_nhfull(aak       ,bbk       ,cck       ,aak2      ,cck2      , &
                  & bbka      ,bbkc      ,ddk       ,kmax      ,icx       , &
                  & icy       ,nmmax     ,nst       ,kfsz1     ,dinv      , &
                  & pbbk      ,pbbkc     ,pnhcor    ,rj        ,p0        , &
                  & p1        ,s1        ,gdp       )
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
! Initialisation for CG iterative solver.
! CG-method
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !
    integer     , pointer :: m1_nhy
    integer     , pointer :: m2_nhy
    integer     , pointer :: n1_nhy
    integer     , pointer :: n2_nhy
    integer     , pointer :: ifirst
    character(8), pointer :: precon
!
! Global variables
!
    integer                                                       :: icx
    integer                                                       :: icy
    integer                                                       :: kmax   !  Description and declaration in iidim.f90
    integer                                                       :: nmmax  !  Description and declaration in dimens.igs
    integer                                                       :: nst
    integer , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: kfsz1  !  Description and declaration in iidim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: aak
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: aak2
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: bbk
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: bbka
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: bbkc
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: cck
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: cck2
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: ddk
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: dinv
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: p1     !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)               :: p0    !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: pbbk
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: pbbkc
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: pnhcor !  Description and declaration in rjdim.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)              :: rj
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub)                    :: s1     !  Description and declaration in rjdim.f90
!
! Local variables
!
    integer :: ddb
    integer :: icxy
    integer :: k
    integer :: m
    integer :: ndelta
    integer :: ndm
    integer :: nm
    integer :: nmd
    integer :: nmst
    integer :: nmstart
    integer :: nmu
    integer :: num
!
!! executable statements -------------------------------------------------------
!
    m1_nhy  => gdp%gdnonhyd%m1_nhy
    m2_nhy  => gdp%gdnonhyd%m2_nhy
    n1_nhy  => gdp%gdnonhyd%n1_nhy
    n2_nhy  => gdp%gdnonhyd%n2_nhy
    precon  => gdp%gdnonhyd%precon
    ifirst  => gdp%gdz_initcg%ifirst
    !
    ddb  = gdp%d%ddbound
    icxy = max(icx, icy)
    !
    ndelta  = n2_nhy - n1_nhy
    nmstart = (n1_nhy+ddb) + (m1_nhy-1+ddb)*icxy
    if (precon == 'none') then
       !
       ! The preconditioner is the unity matrix.
       !
       do m = m1_nhy, m2_nhy
          nmst = nmstart + (m-m1_nhy)*icxy
          do nm = nmst, nmst+ndelta
             do k = 1, kmax
                if (kfsz1(nm,k) /= 0) then
                   dinv(nm,k) = 1.0_fp
                endif
             enddo
          enddo
       enddo
    endif
    if ((precon == 'diag') .or. (precon == 'tridiag')) then
       !
       ! If D is the diagonal part of A, the vector dinv is equal to the 
       ! diagonal of D^(-1/2). Even for 'tridiag' a diagonal scaling is applied.
       !
       do m = m1_nhy, m2_nhy
          nmst = nmstart + (m-m1_nhy)*icxy
          do nm = nmst, nmst+ndelta
             do k = 1, kmax
                if (kfsz1(nm,k) /= 0) then
                   dinv(nm,k) = 1.0_fp / sqrt(bbk(nm,k))
                endif
             enddo
          enddo
       enddo
    endif
    if (precon == 'ilu') then
       !
       ! All information of ILU is stored in dinv, which is D^(-1/2).
       ! We define A = D^(-1/2) * A * D^(-1/2).
       ! The preconditioner is then given by M = (I + L_A) * A * (I + U_A)
       ! with L_A and U_A the strictly lower and upper triangular parts of A.
       !
       call z_ilu_nhfull(aak       ,bbk       ,cck       ,aak2      ,cck2      , &
                       & bbka      ,bbkc      ,kmax      ,icx       , &
                       & icy       ,nmmax     ,kfsz1     , &
                       & dinv      ,gdp)
    endif
    !
    ! column scaling
    !
    do m = m1_nhy, m2_nhy
       nmst = nmstart + (m-m1_nhy)*icxy
       do nm = nmst, nmst+ndelta
          do k = 1, kmax
             if (kfsz1(nm,k) /= 0) then
                ndm        = nm - icy
                nmd        = nm - icx
                nmu        = nm + icx
                num        = nm + icy
                aak (nm,k) = aak (nm,k) * dinv(nmd,k)
                bbk (nm,k) = bbk (nm,k) * dinv(nm ,k)
                cck (nm,k) = cck (nm,k) * dinv(nmu,k)
                ddk (nm,k) = ddk (nm,k) * dinv(nm ,k)
                aak2(nm,k) = aak2(nm,k) * dinv(ndm,k)
                cck2(nm,k) = cck2(nm,k) * dinv(num,k)
             endif
          enddo
          do k = 2, kmax
             if (kfsz1(nm,k) /= 0) then
                bbka(nm,k) = bbka(nm,k) * dinv(nm,k-1)
             endif
          enddo
          do k = 1, kmax-1
             if (kfsz1(nm,k) /= 0) then
                bbkc(nm,k) = bbkc(nm,k) * dinv(nm,k+1)
             endif
          enddo
       enddo
    enddo
    !
    ! row scaling
    !
    do m = m1_nhy, m2_nhy
       nmst = nmstart + (m-m1_nhy)*icxy
       do nm = nmst, nmst+ndelta
          do k = 1, kmax
             if (kfsz1(nm,k) /= 0) then
                aak (nm,k) = aak (nm,k) * dinv(nm,k)
                bbk (nm,k) = bbk (nm,k) * dinv(nm,k)
                cck (nm,k) = cck (nm,k) * dinv(nm,k)
                aak2(nm,k) = aak2(nm,k) * dinv(nm,k)
                cck2(nm,k) = cck2(nm,k) * dinv(nm,k)
                bbka(nm,k) = bbka(nm,k) * dinv(nm,k)
                bbkc(nm,k) = bbkc(nm,k) * dinv(nm,k)
             endif
          enddo
       enddo
    enddo
    if (precon == 'tridiag') then
        !
        ! compute the LU decomposition (and store in pbbk/pbbkc)
        !
        call z_lu(bbka      ,bbk       ,bbkc      ,kmax      ,icx       , &
                & icy       ,nmmax     ,kfsz1     ,pbbk      ,pbbkc     , &
                & gdp       )
    endif
    !
    ! initialize solution pnhcor at the first time step, otherwise phcor
    ! of the last computation is still available as a first guess
    !
    if (ifirst == 1) then
       do nm = 1, nmmax
          do k = 1, kmax
             pnhcor(nm,k) = 0.0_fp
          enddo
       enddo
       ifirst = 0
    endif
end subroutine z_initcg_nhfull

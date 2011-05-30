subroutine z_precon(bbka      ,bbkc      ,pbbk      ,pbbkc     ,kmax      , &
                  & icx       ,icy       ,nmmax     ,kfsz1     ,rj        , &
                  & p1        ,gdp       )
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
!  $Id$
!  $HeadURL$
!!--description-----------------------------------------------------------------
!
! Computes the preconditioning
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
    !
    integer , pointer :: m1_nhy
    integer , pointer :: m2_nhy
    integer , pointer :: n1_nhy
    integer , pointer :: n2_nhy
!
! Global variables
!
    integer                                         , intent(in) :: icx
    integer                                         , intent(in) :: icy
    integer                                         , intent(in) :: kmax  !  Description and declaration in esm_alloc_int.f90
    integer                                                      :: nmmax !  Description and declaration in dimens.igs
    integer, dimension(gdp%d%nmlb:gdp%d%nmub, kmax) , intent(in) :: kfsz1 !  Description and declaration in esm_alloc_int.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax), intent(in) :: bbka
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax), intent(in) :: bbkc
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)             :: p1    !  Description and declaration in esm_alloc_real.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax), intent(in) :: pbbk
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)             :: pbbkc
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax), intent(in) :: rj
!
! Local variables
!
    integer :: ddb
    integer :: icxy
    integer :: k
    integer :: m
    integer :: ndelta
    integer :: nm
    integer :: nmst
    integer :: nmstart
!
!! executable statements -------------------------------------------------------
!
    m1_nhy  => gdp%gdnonhyd%m1_nhy
    m2_nhy  => gdp%gdnonhyd%m2_nhy
    n1_nhy  => gdp%gdnonhyd%n1_nhy
    n2_nhy  => gdp%gdnonhyd%n2_nhy
    !
    ddb  = gdp%d%ddbound
    icxy = max(icx, icy)
    !
    ndelta  = n2_nhy - n1_nhy
    nmstart = (n1_nhy + ddb) + (m1_nhy - 1 + ddb)*icxy
    !
    do m = m1_nhy, m2_nhy
       nmst = nmstart + (m - m1_nhy)*icxy
       do nm = nmst, nmst + ndelta
          !
          ! in case of no preconditioning:
          !  do k = 1, kmax
          !     if (kfsz1(nm, k)/=0) then
          !        p1(nm, k) = rj(nm, k)
          !     endif
          !  enddo
          !
          !
          ! in case of ADI preconditioning in vertical direction:
          !
          if (kfsz1(nm,1) /= 0) then
             p1(nm,1) = rj(nm,1)
          endif
          do k = 2, kmax
             if (kfsz1(nm,k) /= 0) then
                p1(nm,k) = ( rj(nm,k) - bbka(nm,k) * p1(nm,k-1) ) * pbbk(nm,k) 
             endif
          enddo
          do k = kmax-1,1,-1
             if (kfsz1(nm,k) /= 0) then
                p1(nm,k) = p1(nm,k) - pbbkc(nm,k) * p1(nm,k+1)
             endif
          enddo
       enddo
    enddo
end subroutine z_precon

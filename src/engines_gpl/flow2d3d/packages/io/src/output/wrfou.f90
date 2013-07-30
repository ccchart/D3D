subroutine wrfou(nmax      ,mmax      ,nmaxus    ,kmax      ,lmax      , &
               & nofou     ,runid     ,dtsec     ,versio    ,namcon    , &
               & kcs       ,xz        ,yz        ,alfas     ,xcor      , &
               & ycor      ,kfu       ,kfv       ,itdate    ,gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2013.                                
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
!    Function: - open fourier analysis output file
!              - writes results of fourier analysis to output
!                file
!              - closes fourier analysis output file
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use globaldata
    use dfparall
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    character(1)   , dimension(:)        , pointer :: foutyp
!
! Global variables
!
    integer                                                                                 :: itdate !  Reference time in YYYYMMDD
    integer                                                                                 :: kmax   !  Description and declaration in esm_alloc_int.f90
    integer                                                                                 :: lmax   !  Description and declaration in dimens.igs
    integer                                                                                 :: mmax   !  Description and declaration in esm_alloc_int.f90
    integer                                                                                 :: nmax   !  Description and declaration in esm_alloc_int.f90
    integer                                                                                 :: nmaxus !  Description and declaration in esm_alloc_int.f90
    integer                                                                                 :: nofou  !  Description and declaration in dimens.igs
    integer       , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: kcs    !  Description and declaration in esm_alloc_int.f90
    integer       , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: kfu    !  Description and declaration in esm_alloc_int.f90
    integer       , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: kfv    !  Description and declaration in esm_alloc_int.f90
    real(fp)                                                                                :: dtsec  !!  Integration time step [in seconds]
    real(fp)      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: alfas  !  Description and declaration in esm_alloc_real.f90
    real(fp)      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: xcor   !  Description and declaration in esm_alloc_real.f90
    real(fp)      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: xz     !  Description and declaration in esm_alloc_real.f90
    real(fp)      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: ycor   !  Description and declaration in esm_alloc_real.f90
    real(fp)      , dimension(gdp%d%nlb:gdp%d%nub, gdp%d%mlb:gdp%d%mub)                     :: yz     !  Description and declaration in esm_alloc_real.f90
    character(*)                                                                            :: runid  !!  Run identification of this simulation
    character(20) , dimension(lmax)                                                         :: namcon !  Description and declaration in esm_alloc_char.f90
    character(5)                                                               , intent(in) :: versio !!  Version nr. of the current package
!
!
! Local variables
!
    integer                        :: ifou         ! Local teller for fourier functions 
    integer                        :: lrid         ! Length of RUNID character string 
    integer                        :: lunfou
    integer        , external      :: newlun
    character(256)                 :: fixid        ! fixed size version of runid, needed for character concatenation 
    character(256)                 :: filename     ! fixed size version of runid, needed for character concatenation 
!
!
!! executable statements -------------------------------------------------------
!
!
    foutyp => gdp%gdfourier%foutyp
    !
    ! define length of runid and put in fixed size array
    ! size is tested in iniid
    !
    filename = "fourier." // trim(runid)
    if ( parll ) then
       write(filename,'(2a,i3.3)') trim(filename), '-', inode
    endif
    !
    ! Open output file 'fourier.'runid
    !
    lunfou = newlun(gdp)
    open (lunfou, file = filename, status = 'unknown')
    !
    ! Write all requested fourier function output until IFOU > NOFOU
    !
    write (lunfou, '(a,a,a)') '*** Delft3D-FLOW utility FOUMOD *** version ',   &
                            & versio, ' ***'
    ifou = 1
    !
    ! Write requested fourier function output for scalar quantity
    !
    do 
       if (ifou > nofou) exit
       !
       if (foutyp(ifou)=='s') then
          call wrfous(nmax      ,mmax      ,nmaxus    ,kmax      ,lmax      , &
                    & nofou     ,ifou      ,lunfou    ,dtsec     ,namcon    , &
                    & kcs       ,xz        ,yz        ,xcor      ,ycor      , &
                    & kfu       ,kfv       ,itdate    ,gdp       )
          ifou = ifou + 1
       else
          !
          ! Write requested fourier function output for vector quantity
          ! and add 2 to ifou afterwards
          !
          call wrfouv(nmax      ,mmax      ,nmaxus    ,kmax      ,nofou     , &
                    & ifou      ,lunfou    ,dtsec     ,kcs       ,xz        , &
                    & yz        ,alfas     ,xcor      ,ycor      ,kfu       , &
                    & kfv       ,itdate    ,gdp       )
          !
          ifou = ifou + 2
       endif
    enddo
    !
    !
    ! Close fourier output file
    !
    close (lunfou)
end subroutine wrfou

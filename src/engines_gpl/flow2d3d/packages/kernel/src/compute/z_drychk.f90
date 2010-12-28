subroutine z_drychk(idry      ,j         ,nmmaxj    ,nmmax     ,kmax      , &
                  & nfltyp    ,icx       ,icy       ,kfu       ,kfv       , &
                  & kfs       ,kcs       ,kfuz1     ,kfvz1     ,kfsz1     , &
                  & kfumin    ,kfumax    ,kfvmin    ,kfvmax    ,kfsmin    , &
                  & kfsmax    ,kfsmx0    ,s1        ,s0        ,dps       , &
                  & qxk       ,qyk       ,w1        ,dzu1      ,dzv1      , &
                  & dzs1      ,zk        ,gdp       )
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
!    Function: This subroutine checks for drying in water level
!              points. In case the point is dry, all surrounding
!              mask arrays (KFU and KFV) are set to zero and sub-
!              sequently SUD computation will be repeated
! Method used:
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
    ! They replace the  include igd / include igp lines
    !
    integer  , pointer :: lundia
    real(fp) , pointer :: dzmin
!
! Global variables
!
    integer                                             , intent(in)  :: icx    !!  Increment in the X-dir., if ICX= NMAX then computation proceeds in the X-dir. If icx=1 then computation proceeds in the Y-dir.
    integer                                             , intent(in)  :: icy    !!  Increment in the Y-dir. (see ICX)
    integer                                             , intent(out) :: idry   !!  Flag set to 1 if a dry point is detected in routine DRYCHK after SUD is completed
    integer                                                           :: j      !!  Begin pointer for arrays which have been transformed into 1D arrays. Due to the shift in the 2nd (M-) index, J = -2*NMAX + 1
    integer                                             , intent(in)  :: kmax   !  Description and declaration in iidim.f90
    integer                                                           :: nfltyp !  Description and declaration in iidim.f90
    integer                                             , intent(in)  :: nmmax  !  Description and declaration in dimens.igs
    integer                                                           :: nmmaxj !  Description and declaration in dimens.igs
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)        , intent(in)  :: kcs    !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfs    !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfsmax !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)        , intent(in)  :: kfsmin !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfsmx0 !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfu    !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfumax !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfumin !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfv    !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfvmax !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub)                      :: kfvmin !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)  , intent(out) :: kfsz1  !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)  , intent(out) :: kfuz1  !  Description and declaration in iidim.f90
    integer   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)  , intent(out) :: kfvz1  !  Description and declaration in iidim.f90
    real(prec), dimension(gdp%d%nmlb:gdp%d%nmub)        , intent(in)  :: dps    !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub)        , intent(in)  :: s1     !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub)        , intent(in)  :: s0     !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, 0:kmax)              :: w1     !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)  , intent(out) :: dzs1   !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)                :: dzu1   !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)                :: dzv1   !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)  , intent(out) :: qxk    !  Description and declaration in rjdim.f90
    real(fp)  , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)  , intent(out) :: qyk    !  Description and declaration in rjdim.f90
    real(fp)  , dimension(0:kmax), intent(in) :: zk
!
! Local variables
!
    integer                :: k
    integer                :: ndm
    integer                :: nm
    integer                :: nmd
    integer, dimension(1)  :: nm_s1max1
    integer, dimension(1)  :: nm_s1max2
    logical                :: flood
    logical                :: zmodel
    real(fp)               :: s1max1
    real(fp)               :: s1max2
    real(fp)               :: zdiff
    character(300)         :: errmsg
!
!! executable statements -------------------------------------------------------
!
    lundia  => gdp%gdinout%lundia
    dzmin   => gdp%gdzmodel%dzmin
    !
    ! DZMIN has been set in Z-iniflm
    !
    idry = 0
    !
    do nm = 1, nmmax
       nmd = nm - icx
       ndm = nm - icy
       if (kfu(nm)==1 .or. kfu(nmd)==1 .or. kfv(nm)==1 .or. kfv(ndm)==1) then
          if (s1(nm) <= - real(dps(nm),fp)) then
             kfu(nm) = 0
             kfu(nmd) = 0
             kfv(nm) = 0
             kfv(ndm) = 0
             do k = 1, kmax  ! evt. kfsmin, kmax
                kfuz1(nm, k) = 0
                kfuz1(nmd, k) = 0
                kfvz1(nm, k) = 0
                kfvz1(ndm, k) = 0
                qxk(nm, k) = 0.0
                qxk(nmd, k) = 0.0
                qyk(nm, k) = 0.0
                qyk(ndm, k) = 0.0
                dzs1(nm, k) = 0.0
             enddo
             idry = 1
          endif
       endif
    enddo
    !
    ! CHECK FOR FOUR DRY VELOCITY POINTS
    !
    s1max1 = -999.0_fp
    s1max2 = -999.0_fp
    do nm = 1, nmmax
       kfsmax(nm) = -1
       !
       ! determination number of layers at old time level in case of flooding.
       !
       if (kcs(nm)>0) then
          nmd = nm - icx
          ndm = nm - icy
          kfs(nm) = max(kfu(nm), kfu(nmd), kfv(nm), kfv(ndm))
       endif
       !
       if (kfs(nm)==1) then
          !
          ! 15-3-2007 change to allow S1 > ZK(KMAX), needed for NH-models
          !
          if (kfsmx0(nm)== - 1) then
             if (s0(nm)>zk(kmax) ) then
                kfsmx0(nm) = kmax
             else
                do k = kfsmin(nm), kmax
                   if ( zk(k) + dzmin>=s0(nm) ) then
                      kfsmx0(nm) = k
                      exit
                   endif
                enddo
             endif
          endif
          !
          ! 15-3-2007 change to allow S1 > ZK(KMAX), needed for NH-models
          !
          if ( s1(nm)>zk(kmax) ) then
             kfsmax(nm)=kmax
          else
             do k = kfsmin(nm), kmax
                if ( zk(k) + dzmin>=s1(nm) ) then
                   kfsmax(nm)=k
                   exit
                endif
             enddo
          endif
       endif ! determination number of layers at old time level in case of flooding.
       !
       if (kfs(nm)==1) then
          do k = kfsmin(nm), kmax
             if (k<=kfsmax(nm)) then
                kfsz1(nm, k) = 1
             else
                kfsz1(nm, k) = 0
             endif
          enddo
       else
          kfsmax(nm) = -1
          do k = max(kfsmin(nm),1), kmax
             kfsz1(nm, k) = 0
          enddo
       endif
    enddo
    !
    ! issue warning if maximum water level is above zk(kmax) (ZTOP)
    !
    !
    s1max1 = maxval(s1)
    s1max2 = maxval(s1(1:nmmax))
    nm_s1max1 = maxloc(s1)
    nm_s1max2 = maxloc(s1(1:nmmax))
    if (s1max1 > zk(kmax)) then
       write (errmsg, '(a,g10.3,2a,i0,a,i0,2a)') '1: Maximum water level ', s1max1, &
                    & ', (m) exceeds ZTOP specified in input', &
                    & ', for nm = ', nm_s1max1, ', and icx = ', icx, &
                    & '; changing ZTOP to above this value', &
                    & ' is strongly advised'
       call prterr(lundia, 'U190', trim(errmsg))
    endif
    if (s1max2 > zk(kmax)) then
       write (errmsg, '(a,g10.3,2a,i0,a,i0,2a)') '2: Maximum water level ', s1max2, &
                    & ', (m) exceeds ZTOP specified in input', &
                    & ', for nm = ', nm_s1max2, ', and icx = ', icx, &
                    & '; changing ZTOP to above this value', &
                    & ' is strongly advised'
       call prterr(lundia, 'U190', trim(errmsg))
    endif
    !
    ! Recalculate DZS1
    ! Reset all DZS1 to for all inactive points above KFSMAX
    !
    do nm = 1, nmmax
       if (kfs(nm)==1) then
          do k = kfsmin(nm), kfsmax(nm)
             if (kfsmin(nm)==kfsmax(nm)) then
                dzs1(nm, k) = real(dps(nm),fp) + s1(nm)
             elseif (k==kfsmin(nm)) then
                dzs1(nm, k) = zk(k) + real(dps(nm),fp)
             elseif (k==kfsmax(nm)) then
                dzs1(nm, k) = s1(nm) - zk(k - 1)
             else
                dzs1(nm, k) = zk(k) - zk(k - 1)
             endif
          enddo
          do k = kfsmax(nm) + 1, kmax
             dzs1(nm, k) = 0.
          enddo
       endif
    enddo
    !
    ! A "trick" to ensure that "wet" points that were dry
    ! obtains a velocity (see also Z_CHECKU)
    !
    do nm = 1, nmmax
       if (kfsmax(nm)>kfsmx0(nm)) then
          do k = kfsmx0(nm), kfsmax(nm)
             w1(nm, k) = w1(nm, kfsmx0(nm))
          enddo
       endif
    enddo
end subroutine z_drychk

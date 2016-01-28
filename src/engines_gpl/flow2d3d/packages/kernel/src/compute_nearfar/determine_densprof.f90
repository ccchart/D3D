subroutine determine_densprof(kmax      ,thick     ,s0        ,dps       ,rho       ,ha        ,hd        ,&
                             &stype1    ,stype2    ,rhoam     ,rhoas     ,rhoab     ,hint      ,drohj     )

!!--copyright-------------------------------------------------------------------
! Copyright (c) 2009, WL | Delft Hydraulics. All rights reserved.
!!--disclaimer------------------------------------------------------------------
! This code is part of the Delft3D software system. WL|Delft Hydraulics has
! developed c.q. manufactured this code to its best ability and according to the
! state of the art. Nevertheless, there is no express or implied warranty as to
! this software whether tangible or intangible. In particular, there is no
! express or implied warranty as to the fitness for a particular purpose of this
! software, whether tangible or intangible. The intellectual property rights
! related to this software code remain with WL|Delft Hydraulics at all times.
! For details on the licensing agreement, we refer to the Delft3D software
! license and any modifications to this license, if applicable. These documents
! are available upon request.
!!--version information---------------------------------------------------------
! $Author$
! $Date$
! $Revision$
!!--description-----------------------------------------------------------------
!
!    Function: Determines density profile following CORMIX definition
!
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
!
! Global variables
!
    integer                                                     , intent(in) :: kmax
    real(fp)                                                    , intent(in) :: s0
    real(fp)                                                    , intent(in) :: dps
    real(fp)                                                    , intent(in) :: ha
    real(fp)                                                    , intent(in) :: hd
    real(fp)   , dimension(kmax)                                , intent(in) :: rho
    real(fp)   , dimension(kmax)                                , intent(in) :: thick
    real(fp)                                                    , intent(out):: rhoam
    real(fp)                                                    , intent(out):: rhoas
    real(fp)                                                    , intent(out):: rhoab
    real(fp)                                                    , intent(out):: hint
    real(fp)                                                    , intent(out):: drohj
    
    character*1                                                 , intent(out):: stype1
    character*1                                                 , intent(out):: stype2
!
! Local variables
!
    integer                                :: k
    integer                                :: kgrad
    real(fp)                               :: dengra
    real(fp)                               :: d_diff
    real(fp)                               :: maxgrad
    real(fp)                               :: thck
    real(fp) , dimension(:) , allocatable  :: h1
    real(fp) , dimension(:) , allocatable  :: rhoa
!
!
!! executable statements -------------------------------------------------------
!
    !
    allocate (h1   (kmax) )
    allocate (rhoa (kmax) )
    !
    ! Compute parameters density profile C,
    ! start with computing heigths and densities at the cell centers
    ! of the position were the ambient conditions are taken from
    ! (switch positive direction from positive downward to positive upward)
    !

    h1(1)   = 0.5_fp * thick(kmax) * (s0+dps)
    rhoa(1) = rho(kmax)
    !
    do k = 2, kmax
       thck     = 0.5_fp * (thick(kmax-k+2) + thick(kmax-k+1))
       h1 (k)   = h1(k-1) + thck*(s0 + dps)
       rhoa (k) = rho(kmax-k+1)
    enddo

    !
    ! Determine the density gradients,
    ! determine maxmimum density gradient and its location, but first,
    ! ensure stable density profile
    !

    do k = 1, kmax - 1
       if (rhoa(k) < rhoa(k+1)) then
          rhoa(k+1) = rhoa(k)
        endif
    enddo

    maxgrad = -1.0e36_fp
    do k = 1, kmax - 1
       dengra = (rhoa(k) - rhoa(k+1))
       if (dengra > maxgrad) then
          maxgrad = dengra
          kgrad   = k
       endif
    enddo

    !
    ! Determine Profile type C parameters
    !

    rhoab = rhoa(1)
    rhoas = 0.0_fp

    do k = kgrad + 1, kmax
       rhoas = rhoas + rhoa(k)/(kmax - kgrad)
    enddo

    hint  = 0.5_fp*(h1(kgrad + 1) + h1 (kgrad))
    drohj = min(rhoa(kgrad) - rhoa(kgrad + 1), (rhoab - rhoas)-0.01_fp)

    !
    ! Adjust hint such that it is accepted by cormix
    !

    if (hint > 0.89*ha .or. hint > 0.89_fp*hd) then
       hint = min(0.85_fp*ha,0.85_fp*hd)
    endif

    if (hint < 0.41*ha .or. hint < 0.41_fp*hd) then
       hint = max(0.45_fp*ha,0.45_fp*hd)
    endif


    !
    ! Determine profile type
    !

    d_diff = rhoa(1) - rhoa(kmax)
    if (d_diff < 0.2_fp) then
       stype1 = 'U'
       rhoam = 0.0_fp
       do k = 1, kmax
          rhoam = rhoam + rho(k)*thick(k)
       enddo
    else
       stype1 = 'S'
       if (maxgrad < 0.5_fp*d_diff) then
          rhoas  = rho(1)
          rhoab  = rho(kmax)
          stype2 = 'A'
       else
          stype2 = 'C'
       endif
    endif
    
    deallocate (h1   )
    deallocate (rhoa )

end subroutine determine_densprof

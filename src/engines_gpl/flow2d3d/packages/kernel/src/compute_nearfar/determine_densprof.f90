subroutine determine_densprof(kmax       ,thick      ,s0       ,dps    ,rho       , &
                            & ha         ,hd         ,stype1   ,stype2 ,rhoam     , &
                            & rhoas      ,rhoab      ,hint     ,drohj  , &
                            & kfsmin_amb ,kfsmax_amb ,dzs0_amb ,zmodel )

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
    integer                                                     , intent(in)  :: kmax
    integer                                                     , intent(in)  :: kfsmin_amb
    integer                                                     , intent(in)  :: kfsmax_amb 
    real(fp)   , dimension(kmax)                                , intent(in)  :: dzs0_amb
    real(fp)                                                    , intent(in)  :: s0
    real(fp)                                                    , intent(in)  :: dps
    real(fp)                                                    , intent(in)  :: ha
    real(fp)                                                    , intent(in)  :: hd
    real(fp)   , dimension(kmax)                                , intent(in)  :: rho
    real(fp)   , dimension(kmax)                                , intent(in)  :: thick
    real(fp)                                                    , intent(out) :: rhoam
    real(fp)                                                    , intent(out) :: rhoas
    real(fp)                                                    , intent(out) :: rhoab
    real(fp)                                                    , intent(out) :: hint
    real(fp)                                                    , intent(out) :: drohj
    logical                                                     , intent(in)  :: zmodel
    
    character*1                                                 , intent(out) :: stype1
    character*1                                                 , intent(out) :: stype2
!
! Local variables
!
    integer                                :: k
    integer                                :: k0
    integer                                :: k1
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
    h1   = 0.0_fp
    rhoa = 0.0_fp
    !
    ! Compute parameters density profile C,
    ! start with computing heigths and densities at the cell centers
    ! of the position where the ambient conditions are taken from
    ! (switch positive direction from positive downward to positive upward)
    !
    if (.not. zmodel) then
        k0 = 1
        k1 = kmax
        !
        h1(k0)   = 0.5_fp * thick(kmax) * (s0+dps)
        rhoa(k0) = rho(kmax)
        !
        do k = 2, kmax
           thck     = 0.5_fp * (thick(kmax-k+2) + thick(kmax-k+1))
           h1 (k)   = h1(k-1) + thck*(s0 + dps)
           rhoa (k) = rho(kmax-k+1)
        enddo
    else
        k0 = kfsmin_amb
        k1 = kfsmax_amb
        !
        h1(k0)   = 0.5_fp * dzs0_amb(k0)
        rhoa(k0) = rho(k0)
        !
        do k = k0+1, k1
           thck     = 0.5_fp * (dzs0_amb(k) + dzs0_amb(k-1))
           h1 (k)   = h1(k-1) + thck*(s0 + dps)
           rhoa (k) = rho(k)
        enddo
    endif
    !
    ! Determine the density gradients,
    ! determine maxmimum density gradient and its location, but first,
    ! ensure stable density profile
    !
    do k = k0, k1 - 1
       if (rhoa(k) < rhoa(k+1)) then
          rhoa(k+1) = rhoa(k)
        endif
    enddo
    !
    maxgrad = -1.0e36_fp
    do k = k0, k1 - 1
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
    !
    do k = kgrad + 1, k1
       rhoas = rhoas + rhoa(k)/(1 - kgrad)
    enddo
    !
    hint  = 0.5_fp*(h1(kgrad + 1) + h1 (kgrad))
    drohj = min(rhoa(kgrad) - rhoa(kgrad + 1), (rhoab - rhoas)-0.01_fp)
    !
    ! Adjust hint such that it is accepted by cormix
    !
    if (hint > 0.89_fp*ha .or. hint > 0.89_fp*hd) then
       hint = min(0.85_fp*ha,0.85_fp*hd)
    endif
    !
    if (hint < 0.41_fp*ha .or. hint < 0.41_fp*hd) then
       hint = max(0.45_fp*ha,0.45_fp*hd)
    endif
    !
    ! Determine profile type
    !
    d_diff = rhoa(1) - rhoa(kmax)
    if (d_diff < 0.2_fp) then
       stype1 = 'U'
       rhoam = 0.0_fp
       if (.not. zmodel) then
          do k = k0, k1
             rhoam = rhoam + rho(k)*thick(k)
          enddo
       else
          do k = k0, k1
             rhoam = rhoam + rho(k)*dzs0_amb(k)/(s0 + dps)
          enddo
       endif
    else
       stype1 = 'S'
       if (maxgrad < 0.5_fp*d_diff) then
          rhoas  = rho(k0)
          rhoab  = rho(k1)
          stype2 = 'A'
       else
          stype2 = 'C'
       endif
    endif
    
    deallocate (h1   )
    deallocate (rhoa )

end subroutine determine_densprof

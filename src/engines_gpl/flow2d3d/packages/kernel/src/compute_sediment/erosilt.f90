subroutine erosilt(thick    ,kmax     ,ws       , &
                 & wstau    ,entr     ,dicww    ,seddif   ,lundia   , &
                 & h0       ,h1       ,um       ,vm       ,uuu      ,vvv      , &
                 & taub     ,error    ,fixfac   , &
                 & frac     ,sinkse   ,sourse   ,oldmudfrac, flmd2l , tcrdep  , &
                 & tcrero   ,eropar   ,iform    , &
                 & numintpar,numrealpar,numstrpar,dllfunc ,dllhandle, &
                 & intpar   ,realpar  ,strpar   )
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
!    Function: Computes sediment fluxes at the bed using
!              the Partheniades-Krone formulations.
!              Arrays SOURSE and SINKSE are filled
!              Array seddif id filled with dicww for mud
! Method used:
!
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
    !
    include 'sedparams.inc'
    !
    integer                                                   , intent(in)  :: kmax
    integer                                                                 :: lundia   !  Description and declaration in inout.igs
    real(fp)                                                  , intent(in)  :: entr
    real(fp)  , dimension(0:kmax)                             , intent(in)  :: ws
    real(fp)                                                  , intent(out) :: wstau
    real(fp)  , dimension(0:kmax)                             , intent(in)  :: dicww
    real(fp)  , dimension(0:kmax)                             , intent(out) :: seddif
    real(fp)  , dimension(kmax)                               , intent(in)  :: thick
    real(fp)                                                  , intent(in)  :: h0
    real(fp)                                                  , intent(in)  :: h1
    real(fp)                                                  , intent(in)  :: um
    real(fp)                                                  , intent(in)  :: uuu
    real(fp)                                                  , intent(in)  :: vm
    real(fp)                                                  , intent(in)  :: vvv
    real(fp)                                                  , intent(in)  :: taub
    logical                                                   , intent(out) :: error
    real(fp)                                                  , intent(in)  :: fixfac
    real(fp)                                                  , intent(in)  :: frac
    real(fp)                                                  , intent(out) :: sinkse
    real(fp)                                                  , intent(out) :: sourse
    logical                                                   , intent(in)  :: oldmudfrac
    logical                                                   , intent(in)  :: flmd2l
    real(fp)                                                  , intent(in)  :: tcrdep
    real(fp)                                                  , intent(in)  :: tcrero
    real(fp)                                                  , intent(in)  :: eropar
    integer                                                   , intent(in)  :: iform
    integer, dimension(numintpar)   , intent(inout):: intpar
    integer                         , intent(in)   :: numintpar
    integer                         , intent(in)   :: numrealpar
    integer                         , intent(in)   :: numstrpar
    real(hp), dimension(numrealpar) , intent(inout):: realpar
    character(256), dimension(numstrpar), intent(inout):: strpar
    character(256)                  , intent(in)   :: dllfunc
    integer                         , intent(in)   :: dllhandle
!
! Local variables
!
    integer  :: k
    real(fp) :: sour
    real(fp) :: sink
    real(fp) :: taum
    real(fp) :: thick0
    real(fp) :: thick1

    ! Interface to dll is in High precision!
    !
    real(fp)          :: ee
    real(hp)          :: sink_dll
    real(hp)          :: sour_dll
    integer           :: ierror
    integer, external :: perf_function_erosilt
    character(256)    :: errmsg
    character(256)    :: message     ! Contains message from
!
!! executable statements ------------------
!
    ee     = exp(1.0_fp)
    error  = .false.
    !
    ! Calculate total (possibly wave enhanced) roughness
    !
    thick0 = thick(kmax) * h0
    thick1 = thick(kmax) * h1
    !
    ! Bed transport following Partheniades and Krone
    ! but in case of fluid mud, source term is determined by
    ! fluid mud part (sourmu). Information is passed via entr()
    ! maximum erosion is sediment available at bed (ignores sediment
    ! settling during the current morphological timestep)
    ! In case of fluid mud the maximum erosion is determined in sourmu
    ! of the fluid mud module. So ignore this check when fluid mud.
    ! Also, taum is not required in the formulation since whether or not
    ! and how much entrainment occurs is entirely handled by the sourmu
    ! routine.
    !
    ! For 3D model set sediment diffusion coefficient
    ! NOTE THAT IF ALGEBRAIC OR K-L TURBULENCE MODEL IS USED THEN WAVES
    ! ONLY AFFECT THE VERTICAL TURBULENT MIXING VIA THE ENHANCED BED
    ! ROUGHNESS
    !
    if (kmax > 1) then
       do k = 1, kmax
          seddif(k) = dicww(k)
       enddo
    endif
    !
    ! calculation both for mud and floc
    !
    if (flmd2l) then
       !
       ! maximum erosion is sediment available at bed
       ! (ignores sediment settling during the current morphological timestep)
       !
       sour = entr
       if (tcrdep > 0.0) then
          sink = max(0.0_fp , 1.0-taub/tcrdep)
       else
          sink = 0.0
       endif
    else
       if (iform == -1) then
          !
          ! Default Partheniades-Krone formula
          !
          taum = max(0.0_fp, taub/tcrero - 1.0)
          sour = eropar * taum
          if (tcrdep > 0.0) then
             sink = max(0.0_fp , 1.0-taub/tcrdep)
          else
             sink = 0.0
          endif
       elseif (iform == 15) then
          !
          ! User defined formula in DLL
          ! Input parameters are passed via realpar/intpar/strpar-arrays
          !
          realpar( 2) = real(um     ,hp)
          realpar( 3) = real(vm     ,hp)
          realpar( 4) = real(sqrt(um*um + vm*vm),hp)
          realpar( 5) = real(uuu    ,hp)
          realpar( 6) = real(vvv    ,hp)
          realpar( 7) = real(sqrt(uuu*uuu + vvv*vvv),hp)
          if (kmax>1) then
             realpar( 8) = real(h1*thick(kmax)/2.0_fp,hp)
          else
             realpar( 8) = real(h1/ee,hp)
          endif
          !
          ! Initialisation of output variables of user defined transport formulae
          !
          sink_dll    = 0.0_hp
          sour_dll    = 0.0_hp
          message     = ' '
          !
          ! psem/vsem is used to be sure this works fine in DD calculations
          !
          call psemlun
          ierror = perf_function_erosilt(dllhandle       , dllfunc           , &
                                         intpar          , numintpar         , &
                                         realpar         , numrealpar        , &
                                         strpar          , numstrpar         , &
                                         sink_dll        , sour_dll          , &
                                         message)
          call vsemlun
          if (ierror /= 0) then
             write(errmsg,'(a,a,a)') 'Cannot find function "',trim(dllfunc),'" in dynamic library.'
             call prterr (lundia,'U021', trim(errmsg))
             error = .true.
             return
          endif
          if (message /= ' ') then
             write (lundia,'(a,a,a)') '*** ERROR Message from user defined erosion/deposition formulae ',trim(dllfunc),' :'
             write (lundia,'(a,a  )') '          ', trim(message)
             write (lundia,'(a    )') ' '
             error = .true.
             return
          endif
          !
          ! Output parameters
          !
          sour    = real(sour_dll,fp)
          sink    = real(sink_dll,fp)
       endif
    endif
    !
    wstau         = ws(kmax) * sink
    if (.not.flmd2l) then
       if (oldmudfrac) then
          sour = fixfac * sour
       else
          sour = fixfac * frac * sour
       endif
    endif
    sourse = sour / thick0
    sinkse = wstau / thick1
end subroutine erosilt

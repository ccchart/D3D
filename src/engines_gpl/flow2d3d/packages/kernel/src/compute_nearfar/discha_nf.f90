subroutine discha_nf(kmax      ,lstsci    ,nmmax     ,kfs       ,sour      ,sink    , &
                   & volum1    ,volum0    ,r0        ,thick     ,kfsmn0   ,kfsmx0   , gdp   )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2016.                                
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
!   Function: The discharges resulting from near field simulation (jet3d/nrfield etc) are
!             added to the sink and source terms of the continuity equation.
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

    integer                          ,pointer :: no_dis
    integer , dimension(:)           ,pointer :: m_intake
    integer , dimension(:)           ,pointer :: n_intake
    integer , dimension(:)           ,pointer :: k_intake
    real(fp), dimension(:)           ,pointer :: q_diff
    real(fp), dimension(:,:,:)       ,pointer :: disnf
    real(fp), dimension(:,:,:,:)     ,pointer :: sournf
    logical                          ,pointer :: zmodel


!
! Global variables
!
    integer                                                 , intent(in)  :: kmax
    integer,  dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kfs    ! Description and declaration in esm_alloc_int.f90
    integer                                                 , intent(in)  :: lstsci ! Description and declaration in dimens.igs
    integer                                                               :: nmmax  ! Description and declaration in dimens.igs
    integer , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kfsmx0 ! Description and declaration in esm_alloc_int.f90
    integer , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kfsmn0 ! Description and declaration in esm_alloc_int.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax, lstsci), intent(in)  :: r0     ! Description and declaration in esm_alloc_real.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax, lstsci)              :: sink   ! Description and declaration in esm_alloc_real.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax, lstsci)              :: sour   ! Description and declaration in esm_alloc_real.f90
    real(fp), dimension(kmax)                               , intent(in)  :: thick  ! Description and declaration in esm_alloc_real.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)        , intent(in)  :: volum0 ! Description and declaration in esm_alloc_real.f90
    real(fp), dimension(gdp%d%nmlb:gdp%d%nmub, kmax)        , intent(in)  :: volum1 ! Description and declaration in esm_alloc_real.f90
!
! Local variables
!
    integer                               :: idis
    integer                               :: nm
    integer                               :: lcon
    integer                               :: k
    integer                               :: nm_intake
    real(fp)                              :: q_tot
    real(fp) , dimension(:) , allocatable :: total_mass
!
!! executable statements -------------------------------------------------------
!
    no_dis    => gdp%gdnfl%no_dis
    disnf     => gdp%gdnfl%disnf
    sournf    => gdp%gdnfl%sournf
    m_intake  => gdp%gdnfl%m_intake
    n_intake  => gdp%gdnfl%n_intake
    k_intake  => gdp%gdnfl%k_intake
    q_diff    => gdp%gdnfl%q_diff
    zmodel    => gdp%gdprocs%zmodel
    !
    ! Fill sinks for difu
    ! Determine total subtracted mass for negative discharges
    !
    if (.not. zmodel) then
       !
       ! Sigma-model
       !
       do idis = 1, no_dis
          allocate (total_mass(lstsci))
          !
          total_mass  = 0.0_fp
          q_tot       = 0.0_fp
          !
          do k = 1, kmax
             do nm = 1, nmmax
                if (disnf(nm,k,idis) > 0.0_fp) then
                   q_tot = q_tot + disnf(nm,k,idis)
                endif
                if (disnf(nm,k,idis) < 0.0_fp) then
                   do lcon = 1,lstsci
                      sink(nm, k, lcon) = sink(nm, k, lcon) -             &
                                        & disnf(nm,k,idis)/volum1(nm, k)
                      total_mass(lcon)  = total_mass(lcon)  - disnf(nm,k,idis)*r0(nm,k,lcon)
                   enddo
                endif
             enddo
          enddo
       
          !
          ! Fill sinks for coupled intake
          !
          ! (Op verzoek van Robin uitgezet. Onttrekkening dus als normaal, ontkoppeld, onttrekkingspunt
          !
          !         if (m_intake(idis) > 0) then
          !            call n_and_m_to_nm(n_intake(idis), m_intake(idis), nm_intake, gdp)
          !            do lcon = 1, lstsci
          !               if (k_intake(idis) /=0) then
          !                  sink(nm_intake, k_intake(idis), lcon) = sink(nm_intake, k_intake(idis), lcon) +             &
          !                                                        & q_diff(idis)/volum1(nm_intake, k_intake(idis))
          !               else
          !                  do k = 1, kmax
          !                     sink(nm_intake, k, lcon) = sink(nm_intake, k, lcon) +             &
          !                                              & thick(k)*q_diff(idis)/volum1(nm_intake, k)
          !                  enddo
          !               endif
          !            enddo
          !         endif
          !
          ! Fill sour array for difu
          ! Add total subtracted mass and initial discharge amount (sournf)
          !
          do lcon = 1,lstsci
             do k = 1, kmax
                do nm = 1, nmmax
                   if (sournf(nm,k,lcon,idis) > 0.0_fp) then
                      ! if (disnf(nm,k,idis) < 0.0_fp) then
                      !    sour(nm,k,lcon) = sour(nm,k,lcon) + sournf(nm,k,lcon,idis)/volum0(nm,k)
                      ! else
                      sour(nm, k, lcon) = sour(nm, k, lcon)  +                        &
                                        & (sournf(nm,k,lcon,idis) + disnf(nm,k,idis)*total_mass(lcon)/q_tot)/&
                                        & volum0(nm, k)
                      ! endif
                   endif
                enddo
             enddo
          enddo
          !
          deallocate (total_mass)
          !
       enddo
    else
       !
       ! Z-model
       !
       do idis = 1, no_dis
          allocate (total_mass(lstsci))
          !
          total_mass  = 0.0_fp
          q_tot       = 0.0_fp
          !
          do nm = 1, nmmax
             do k = kfsmn0(nm), kfsmx0(nm)
                if (disnf(nm,k,idis) > 0.0_fp) then
                   q_tot = q_tot + disnf(nm,k,idis)
                endif
                if (disnf(nm,k,idis) < 0.0_fp) then
                   do lcon = 1,lstsci
                      sink(nm, k, lcon) = sink(nm, k, lcon) -             &
                                        & disnf(nm,k,idis) !/volum1(nm, k)
                      total_mass(lcon)  = total_mass(lcon)  - disnf(nm,k,idis)*r0(nm,k,lcon)
                   enddo
                endif
             enddo
          enddo
          !
          ! Fill sinks for coupled intake
          !
          ! (Op verzoek van Robin uitgezet. Onttrekkening dus als normaal, ontkoppeld, onttrekkingspunt
          !
          !         if (m_intake(idis) > 0) then
          !            call n_and_m_to_nm(n_intake(idis), m_intake(idis), nm_intake, gdp)
          !            do lcon = 1, lstsci
          !               if (k_intake(idis) /=0) then
          !                  sink(nm_intake, k_intake(idis), lcon) = sink(nm_intake, k_intake(idis), lcon) +             &
          !                                                        & q_diff(idis)/volum1(nm_intake, k_intake(idis))
          !               else
          !                  do k = 1, kmax
          !                     sink(nm_intake, k, lcon) = sink(nm_intake, k, lcon) +             &
          !                                              & thick(k)*q_diff(idis)/volum1(nm_intake, k)
          !                  enddo
          !               endif
          !            enddo
          !         endif
          !
          ! Fill sour array for difu
          ! Add total subtracted mass and initial discharge amount (sournf)
          !
          do lcon = 1,lstsci
             do nm = 1, nmmax
                 do k = kfsmn0(nm), kfsmx0(nm)
                   if (sournf(nm,k,lcon,idis) > 0.0_fp) then
                      ! if (disnf(nm,k,idis) < 0.0_fp) then
                      !    sour(nm,k,lcon) = sour(nm,k,lcon) + sournf(nm,k,lcon,idis)/volum0(nm,k)
                      ! else
                      sour(nm, k, lcon) = sour(nm, k, lcon)  +                        &
                                        & (sournf(nm,k,lcon,idis) + disnf(nm,k,idis)*total_mass(lcon)/q_tot) !/&
                                        !& volum0(nm, k)
                      ! endif
                   endif
                enddo
             enddo
          enddo
          !
          deallocate (total_mass)
          !
       enddo
    endif
end subroutine discha_nf

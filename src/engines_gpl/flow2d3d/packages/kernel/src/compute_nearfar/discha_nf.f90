subroutine discha_nf(kmax      ,lstsci    ,nmmax     ,kfs       ,sour      ,sink    , &
                   & volum1    ,volum0    ,r0        ,thick     ,kfsmn0   ,kfsmx0   , &
                   & kcs       ,gdp   )
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
    use mathconsts
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
    real(fp), dimension(:,:,:)       ,pointer :: disnf_intake
    real(fp), dimension(:,:,:)       ,pointer :: disnf_entr
    real(fp), dimension(:,:,:,:)     ,pointer :: sournf
    logical                          ,pointer :: zmodel
!
! Constants
!
    integer, parameter :: Q_TOT = 0

!
! Global variables
!
    integer                                                 , intent(in)  :: kmax
    integer,  dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kfs    ! Description and declaration in esm_alloc_int.f90
    integer                                                 , intent(in)  :: lstsci ! Description and declaration in dimens.igs
    integer                                                               :: nmmax  ! Description and declaration in dimens.igs
    integer , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kcs    !  Description and declaration in esm_alloc_int.f90
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
    integer                                 :: ierror
    integer                                 :: idis
    integer                                 :: nm, n, m
    integer                                 :: lcon
    integer                                 :: k
    integer                                 :: nm_intake
    real(fp) , dimension(:,:) , allocatable :: tot_mass
!
!! executable statements -------------------------------------------------------
!
    no_dis       => gdp%gdnfl%no_dis
    disnf        => gdp%gdnfl%disnf
    disnf_intake => gdp%gdnfl%disnf_intake
    disnf_entr   => gdp%gdnfl%disnf_entr
    sournf       => gdp%gdnfl%sournf
    m_intake     => gdp%gdnfl%m_intake
    n_intake     => gdp%gdnfl%n_intake
    k_intake     => gdp%gdnfl%k_intake
    q_diff       => gdp%gdnfl%q_diff
    zmodel       => gdp%gdprocs%zmodel
    !
    allocate (tot_mass(0:lstsci,no_dis), stat=ierror)
    tot_mass  = 0.0_fp
    !
    ! Fill tot_mass and sink for difu
    !
    if (.not. zmodel) then
       !
       ! Sigma-model
       !
       do idis = 1, no_dis
          ! Determine total subtracted mass from negative discharges in disnf_entr (or disnf)
          ! tot_mass(lcon, idis): mass of constituent lcon, removed from the system by negative discharges, summed over the full domain, for this discharge
          ! q_tot           : discharge volume        , removed from the system by negative discharges, summed over the full domain, for this discharge
          ! Parallel: The loop below performs a summation of the total nearfield-farfield discharges (per diffusor idis)
          ! and of the total mass per constituent and per diffusor idis for all points.
          ! This requires an allreduce for parallel runs or the summation must be done previously and stored
          !
          do nm = 1, nmmax
             if (kcs(nm) <= 0) cycle
             do k = 1, kmax
                call nm_to_n_and_m(nm, n, m, gdp)
                !if (abs(disnf_entr(nm,k,idis)) > eps_fp .or. abs(disnf(nm,k,idis)) > eps_fp) write(gdp%gdinout%lundia,'(4i4,2f20.10)') nm, n, m, k, disnf_entr(nm,k,idis), disnf(nm,k,idis)
                if (disnf_entr(nm,k,idis) < 0.0_fp) then
                   tot_mass(Q_TOT, idis) = tot_mass(Q_TOT, idis) - disnf_entr(nm,k,idis)
                   do lcon = 1,lstsci
                      sink(nm, k, lcon) = sink(nm, k, lcon) - disnf_entr(nm,k,idis)/volum1(nm, k)
                      tot_mass(lcon, idis)  = tot_mass(lcon, idis)  - disnf_entr(nm,k,idis)*r0(nm,k,lcon)
                   enddo
                endif
                if (disnf_intake(nm,k,idis) < 0.0_fp) then
                   do lcon = 1,lstsci
                      sink(nm, k, lcon) = sink(nm, k, lcon) - disnf_intake(nm,k,idis)/volum1(nm, k)
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
          ! Add amount discharged via diffuser (sournf),
          ! and the amount subtracted via entrainment (disnf_entr * r0),  weighted with the fraction for this cell (disnf_entr/tot_mass(Q_TOT, idis))
          !
       enddo
    else
       !
       ! Z-model
       !
       do idis = 1, no_dis
          ! Determine total subtracted mass from negative discharges in disnf_entr (or disnf)
          ! tot_mass(lcon, idis): mass of constituent lcon, removed from the system by negative discharges, summed over the full domain, for this discharge
          ! tot_mass(Q_TOT, idis)           : discharge volume        , removed from the system by negative discharges, summed over the full domain, for this discharge
          ! Parallel: The loop below performs a summation of the total nearfield-farfield discharges (per diffusor idis)
          ! and of the total mass per constituent and per diffusor idis for all points.
          ! This requires an allreduce for parallel runs or the summation must be done previously and stored
          !
          do nm = 1, nmmax
             if (kcs(nm) <= 0) cycle
             do k = kfsmn0(nm), kfsmx0(nm)
                if (disnf_entr(nm,k,idis) < 0.0_fp) then
                   tot_mass(Q_TOT, idis) = tot_mass(Q_TOT, idis) - disnf_entr(nm,k,idis)
                   do lcon = 1,lstsci
                      ! Note that for the z-model, the sinks are not divided by the volume! (compare with sigma)
                      sink(nm, k, lcon) = sink(nm,k,lcon)  - disnf_entr(nm,k,idis)
                      tot_mass(lcon, idis)  = tot_mass(lcon, idis) - disnf_entr(nm,k,idis)*r0(nm,k,lcon)
                   enddo
                endif
                if (disnf_intake(nm,k,idis) < 0.0_fp) then
                   do lcon = 1,lstsci
                      ! Note that for the z-model, the sinks are not divided by the volume! (compare with sigma)
                      sink(nm, k, lcon) = sink(nm,k,lcon) - disnf_intake(nm,k,idis)
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
          ! Add amount discharged via diffuser (sournf),
          ! and the amount subtracted via entrainment (disnf_entr * r0),  weighted with the fraction for this cell (disnf_entr/tot_mass(Q_TOT, idis))
          !
       enddo
    endif
    !
    ! TODO: exchange tot_mass with other partitions
    do idis = 1, no_dis
          write(gdp%gdinout%lundia,'(a,3f20.10)') "tot_mass(Q_TOT, idis), tot_mass(1,2):", (tot_mass(lcon, idis), lcon=0,lstsci)
    enddo
    !
    ! Fill sour (using tot_mass) for difu
    !
    if (.not. zmodel) then
       !
       ! Sigma-model
       !
       do idis = 1, no_dis
          do lcon = 1,lstsci
             do k = 1, kmax
                do nm = 1, nmmax
                   if (disnf(nm,k,idis) > 0.0_fp) then
                      sour(nm, k, lcon) = sour(nm, k, lcon)                                                         &
                                        & + (sournf(nm,k,lcon,idis) + disnf_entr(nm,k,idis)*tot_mass(lcon, idis)/tot_mass(Q_TOT, idis)) &
                                        &   / volum0(nm, k)
                   endif
                enddo
             enddo
          enddo
       enddo
    else
       !
       ! Z-model
       !
       do idis = 1, no_dis
          do lcon = 1,lstsci
             do nm = 1, nmmax
                 do k = kfsmn0(nm), kfsmx0(nm)
                   if (disnf(nm,k,idis) > 0.0_fp) then
                      ! Note that for the z-model, the sources are not divided by the volume!
                      sour(nm, k, lcon) = sour(nm, k, lcon)                          &
                                        & + (sournf(nm,k,lcon,idis) + disnf_entr(nm,k,idis)*tot_mass(lcon, idis)/tot_mass(Q_TOT, idis))
                   endif
                enddo
             enddo
          enddo
       enddo
    endif
    !
    deallocate (tot_mass, stat=ierror)
end subroutine discha_nf

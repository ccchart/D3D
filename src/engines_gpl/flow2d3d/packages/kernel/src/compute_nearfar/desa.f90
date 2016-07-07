subroutine desa(nlb     ,nub     ,mlb     ,mub        ,kmax       , &
              & lstsci  ,no_dis  ,lsal    ,ltem       , &
              & idis    ,thick   , &
              & kcs     ,xz      ,yz      ,alfas      , &
              & dps     ,s0      ,r0      ,kfsmn0     ,kfsmx0     , &
              & dzs0    ,disnf   ,sournf  ,nf_src_momu,nf_src_momv, &
              & linkinf ,gdp     )
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
!    Function: Converts Cosumo output to delft3d sources
!              following the DESA methodology of Joseph Lee
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
    integer ,dimension(:)          , pointer :: m_intake
    integer ,dimension(:)          , pointer :: n_intake
    integer ,dimension(:)          , pointer :: k_intake
    integer                        , pointer :: lstsc
    integer                        , pointer :: lunscr
    integer                        , pointer :: lundia
	integer                        , pointer :: nf_const_operator
    real(fp)                       , pointer :: grdang
    real(fp)                       , pointer :: nf_q_source
    real(fp)                       , pointer :: nf_q_intake
    real(fp), dimension(:)         , pointer :: nf_const 
    real(fp), dimension(:,:)       , pointer :: nf_intake
    real(fp), dimension(:,:)       , pointer :: nf_sink  
    real(fp), dimension(:,:)       , pointer :: nf_sour  
    real(fp),dimension(:)          , pointer :: q_diff
    real(fp),dimension(:,:)        , pointer :: const_diff
	logical                        , pointer :: nf_src_mom
    logical , dimension(:)         , pointer :: flbcktemp
    logical                        , pointer :: zmodel
!
! Parameters
!
    integer, parameter :: IX    = 1
    integer, parameter :: IY    = 2
    integer, parameter :: IZ    = 3
    integer, parameter :: IS    = 4
    integer, parameter :: IH    = 5
    integer, parameter :: IW    = 6
    integer, parameter :: IUMAG = 7
    integer, parameter :: IUDIR = 8
!
! Global variables
!
    integer                                                    , intent(in)    :: nlb
    integer                                                    , intent(in)    :: nub
    integer                                                    , intent(in)    :: mlb
    integer                                                    , intent(in)    :: mub
    integer                                                    , intent(in)    :: no_dis
    integer                                                    , intent(in)    :: idis     !  Description and declaration in tricom.igs
    integer                                                    , intent(in)    :: kmax     !  Description and declaration in tricom.igs
    integer                                                    , intent(in)    :: lstsci   !  Description and declaration in tricom.igs
    integer                                                    , intent(in)    :: lsal     !  Description and declaration in tricom.igs
    integer                                                    , intent(in)    :: ltem     !  Description and declaration in tricom.igs
    integer    , dimension(nlb:nub,mlb:mub)                    , intent(in)    :: kcs      !  Description and declaration in esm_alloc_real.f90
    integer    , dimension(nlb:nub,mlb:mub)                    , intent(in)    :: kfsmn0   !  Description and declaration in esm_alloc_real.f90
    integer    , dimension(nlb:nub,mlb:mub)                    , intent(in)    :: kfsmx0   !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(10)                                 , intent(inout) :: linkinf
    real(fp)   , dimension(nlb:nub,mlb:mub)                    , intent(in)    :: s0       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(nlb:nub,mlb:mub,kmax)               , intent(in)    :: dzs0     !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(nlb:nub,mlb:mub,kmax,lstsci)        , intent(in)    :: r0       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(nlb:nub,mlb:mub)                    , intent(in)    :: xz       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(nlb:nub,mlb:mub)                    , intent(in)    :: yz       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(nlb:nub,mlb:mub)                    , intent(in)    :: alfas    !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(kmax)                               , intent(in)    :: thick    !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(nlb:nub,mlb:mub,kmax,no_dis)                        :: disnf
    real(fp)   , dimension(nlb:nub,mlb:mub,kmax,lstsci,no_dis)                 :: sournf
    real(fp)   , dimension(nlb:nub,mlb:mub,kmax,no_dis)                        :: nf_src_momu
    real(fp)   , dimension(nlb:nub,mlb:mub,kmax,no_dis)                        :: nf_src_momv
    real(prec) , dimension(nlb:nub,mlb:mub)                    , intent(in)    :: dps      !  Description and declaration in esm_alloc_real.f90
!
! Local variables
!
    integer                              :: ierror
    integer                              :: irow
    integer                              :: idum
    integer                              :: iidis
    integer                              :: k
    integer                              :: k_irow
    integer                              :: k_last
    integer                              :: k_start
    integer                              :: k_top_tmp
    integer                              :: k_down_tmp
    integer                              :: k_end_top
    integer                              :: k_end_down
    integer                              :: lcon
    integer                              :: n
    integer                              :: m
    integer                              :: ndis_track
    integer                              :: sink_cnt
    integer                              :: sour_cnt
    integer                              :: intake_cnt
    integer                              :: n_irow
    integer                              :: m_irow
    integer                              :: n_start
    integer                              :: m_start
    integer                              :: n_end
    integer                              :: m_end
    integer                              :: n_last
    integer                              :: m_last
    integer                              :: n_tmp
    integer                              :: m_tmp
    integer                              :: k_tmp
    real(fp)                             :: conc
    real(fp)                             :: conc_intake
    real(fp)                             :: dis_dil
    real(fp)                             :: dis_tot
    real(fp)                             :: dis_per_intake
    real(fp)                             :: q1
    real(fp)                             :: q2
    real(fp)                             :: pi
    real(fp)                             :: thick_tot
    real(fp)                             :: ang_end
    real(fp)                             :: dx
    real(fp)                             :: dy
    real(fp)                             :: hhi
    real(fp)                             :: wght
    real(fp)                             :: xstart
    real(fp)                             :: xend
    real(fp)                             :: ystart
    real(fp)                             :: yend
    real(fp)                             :: momu_tmp
    real(fp)                             :: momv_tmp
    real(fp),dimension(:), allocatable   :: momu
    real(fp),dimension(:), allocatable   :: momv
    real(fp),dimension(:), allocatable   :: weight
    integer, dimension(:), allocatable   :: n_dis
    integer, dimension(:), allocatable   :: m_dis
    integer, dimension(:), allocatable   :: k_dis
!
!! executable statements -------------------------------------------------------
!
    lstsc             => gdp%d%lstsc
    lunscr            => gdp%gdinout%lunscr
    lundia            => gdp%gdinout%lundia
    m_intake          => gdp%gdnfl%m_intake
    n_intake          => gdp%gdnfl%n_intake
    k_intake          => gdp%gdnfl%k_intake
    grdang            => gdp%gdtricom%grdang
    q_diff            => gdp%gdnfl%q_diff
    const_diff        => gdp%gdnfl%const_diff
    flbcktemp         => gdp%gdheat%flbcktemp
    zmodel            => gdp%gdprocs%zmodel
    nf_intake         => gdp%gdnfl%nf_intake
    nf_q_source       => gdp%gdnfl%nf_q_source
    nf_q_intake       => gdp%gdnfl%nf_q_intake
    nf_const_operator => gdp%gdnfl%nf_const_operator
    nf_const          => gdp%gdnfl%nf_const 
    nf_sink           => gdp%gdnfl%nf_sink  
    nf_sour           => gdp%gdnfl%nf_sour  
    nf_src_mom        => gdp%gdnfl%nf_src_mom
    !
    dis_dil   = 0.0_fp
    dis_tot   = 0.0_fp
    pi        = acos(-1.0_fp)
    !
    disnf      (nlb:nub,mlb:mub, 1:kmax, idis)          = 0.0_fp
    sournf     (nlb:nub,mlb:mub, 1:kmax, 1:lstsci,idis) = 0.0_fp
    nf_src_momu(nlb:nub,mlb:mub, 1:kmax, idis)          = 0.0_fp
    nf_src_momv(nlb:nub,mlb:mub, 1:kmax, idis)          = 0.0_fp
    !
    ! Handle sinks and sources
    !
    sink_cnt = size(nf_sink,1)
    sour_cnt = size(nf_sour,1)
    !
    if (sink_cnt > 0) then
       !
       ! Get characteristics starting point
       !
       call findnmk(nlb          ,nub          ,mlb          ,mub    ,xz     ,yz     , &
                  & dps          ,s0           ,kcs          ,thick  ,kmax   , &
                  & nf_sink(1,IX),nf_sink(1,IY),nf_sink(1,IZ),n_start,m_start,k_start, &
                  & kfsmn0       ,kfsmx0       ,dzs0         ,zmodel ,gdp    )
       n_last = n_start
       m_last = m_start
       k_last = k_start
       !
       ! Get characteristics end      point
       !
       call findnmk(nlb          ,nub          ,mlb          ,mub    ,xz     , yz      , &
                  & dps          ,s0           ,kcs          ,thick  ,kmax   , &
                  & nf_sour(1,IX),nf_sour(1,IY),nf_sour(1,IZ),n_end  ,m_end  ,k_end_top, &
                  & kfsmn0       ,kfsmx0       ,dzs0         ,zmodel ,gdp    )
       !
       ! For postprocessing store begin and end coordinates of the plume trajectory
       !
       linkinf( 7) = n_start
       linkinf( 8) = m_start
       linkinf( 9) = n_end
       linkinf(10) = m_end
       !
       ! Cycle over sink points in Cosumo output file
       !
       do irow = 2, sink_cnt
          !
          ! Get position of point
          !
          call findnmk(nlb             ,nub             ,mlb             ,mub   ,xz    ,yz    , &
                     & dps             ,s0              ,kcs             ,thick ,kmax  ,  &
                     & nf_sink(irow,IX),nf_sink(irow,IY),nf_sink(irow,IZ),n_irow,m_irow,k_irow, &
                     & kfsmn0          ,kfsmx0          ,dzs0            ,zmodel,gdp   )
          if (n_irow==0 .or. m_irow==0 .or. k_irow==0) then
             n_irow  = n_last
             m_irow  = m_last
             k_irow  = k_last
          endif
          n_last   = n_irow
          m_last   = m_irow
          k_last   = k_irow
          !
          ! Fill dis_nf array, Desa Method:
          ! For all non-end-points:
          ! Substract the amount of water corresponding with the dilution
          ! Keep track of total amounts of water, salt in order to discharge the correct
          ! amounts at the end of the near field
          !
          if (n_last/=n_end .or. m_last/=m_end .or. k_last/=k_end_top) then
             dis_dil = 1.0_fp * (nf_sink(irow,IS)-nf_sink(irow-1,IS)) * nf_q_source
             dis_tot = dis_tot + dis_dil
             disnf(n_last,m_last,k_last,idis) = disnf(n_last,m_last,k_last,idis) - dis_dil
          endif
       enddo
    endif
    !
    if (sour_cnt > 0) then
       if (sour_cnt==1 .and. sink_cnt>0) then
          ! (Single) source point:
          ! Determine the relative thickness over which to distribute the diluted discharge:
          ! Define the line through the source point, perpendicular to the line connecting the
          ! last sink point with the source point. Define the line piece on this line, using
          ! the specified source-width. Walk with 1000 steps over this line piece and check in what
          ! nm cell you are. With this information, you can define the nm-points to distribute the
          ! sources over and their relative weights.
          !
          call findnmk(nlb          ,nub          ,mlb                        ,mub   ,xz   , yz      , &
                     & dps          ,s0           ,kcs                        ,thick ,kmax , &
                     & nf_sour(1,IX),nf_sour(1,IY),nf_sour(1,IZ)-nf_sour(1,IH),n_end ,m_end,k_end_top, &
                     & kfsmn0       ,kfsmx0       ,dzs0                       ,zmodel,gdp  )
          call findnmk(nlb          ,nub          ,mlb                        ,mub   ,xz   ,yz        , &
                     & dps          ,s0           ,kcs                        ,thick ,kmax , &
                     & nf_sour(1,IX),nf_sour(1,IY),nf_sour(1,IZ)+nf_sour(1,IH),n_end ,m_end,k_end_down, &
                     & kfsmn0       ,kfsmx0       ,dzs0                       ,zmodel,gdp  )
          !
          ! Determine grid cells over which to distribute the diluted discharge, begin and and of horizontal distribution area
          !
          ang_end = atan2( (nf_sour(1,IY)-nf_sink(sink_cnt,IY)) , (nf_sour(1,IX)-nf_sink(sink_cnt,IX)) )
          dx      = -1.0_fp*nf_sour(1,IW)*cos(pi/2.0_fp - ang_end)
          dy      =  1.0_fp*nf_sour(1,IW)*sin(pi/2.0_fp - ang_end)
          !
          !      dx = 0.0_fp
          !      dy = 0.0_fp
          !
          xstart   = nf_sour(1,IX) + dx
          ystart   = nf_sour(1,IY) + dy
          xend     = nf_sour(1,IX) - dx
          yend     = nf_sour(1,IY) - dy
          !
          ! Determine grid cell numbers over which to distribute the diluted discharge
          !
          allocate (n_dis (1000), stat=ierror)
          allocate (m_dis (1000), stat=ierror)
          if (allocated(k_dis)) deallocate(k_dis, stat=ierror)
          allocate (weight(1000), stat=ierror)
          allocate (momu  (1000), stat=ierror)
          allocate (momv  (1000), stat=ierror)
          n_dis      = 0
          m_dis      = 0
          weight     = 0.0_fp
          ndis_track = 1
          call findnmk(nlb   ,nub   ,mlb   ,mub   ,xz   ,yz  , &
                     & dps   ,s0    ,kcs   ,thick ,kmax , &
                     & xstart,ystart,0.0_fp,n_tmp ,m_tmp,idum, &
                     & kfsmn0,kfsmx0,dzs0  ,zmodel,gdp  )
          n_dis (1) = n_tmp
          m_dis (1) = m_tmp
          weight(1) = 0.001_fp
          !
          ! Momentum
          !
          if (nf_src_mom) then
             call magdir_to_uv(alfas(n_dis(1),m_dis(1)), grdang          , &
                             & nf_sour(1,IUMAG)        , nf_sour(1,IUDIR), momu_tmp, momv_tmp)
             momu(1) = momu_tmp * 0.001_fp
             momv(1) = momv_tmp * 0.001_fp
          endif
          !
          dx = (xend - xstart)/999.0_fp
          dy = (yend - ystart)/999.0_fp
          !
          do iidis = 1, 999
             call findnmk(nlb              ,nub              ,mlb   ,mub   ,xz   ,yz  , &
                        & dps              ,s0               ,kcs   ,thick ,kmax , &
                        & xstart + iidis*dx,ystart + iidis*dy,0.0_fp,n_tmp ,m_tmp,idum, &
                        & kfsmn0           ,kfsmx0           ,dzs0  ,zmodel,gdp  )
             if (n_tmp/=n_dis(ndis_track) .or. m_tmp/=m_dis(ndis_track)) then
                ndis_track             = ndis_track + 1
                n_dis(ndis_track)      = n_tmp
                m_dis(ndis_track)      = m_tmp
             endif
             weight(ndis_track) = weight(ndis_track) + 0.001_fp
             !
             ! Momentum
             !
             if (nf_src_mom) then
                call magdir_to_uv(alfas(n_dis(ndis_track),m_dis(ndis_track)), grdang          , &
                                & nf_sour(1,IUMAG)                          , nf_sour(1,IUDIR), momu_tmp, momv_tmp)
                momu(ndis_track) = momu(ndis_track) + momu_tmp * 0.001_fp
                momv(ndis_track) = momv(ndis_track) + momv_tmp * 0.001_fp
             endif
          enddo
       else
          ! Multiple source points defined or no sink points defined
          !
          wght = 1.0_fp / real(sour_cnt,fp)
          allocate (n_dis (sour_cnt), stat=ierror)
          allocate (m_dis (sour_cnt), stat=ierror)
          allocate (k_dis (sour_cnt), stat=ierror)
          allocate (weight(sour_cnt), stat=ierror)
          allocate (momu  (sour_cnt), stat=ierror)
          allocate (momv  (sour_cnt), stat=ierror)
          n_dis      = 0
          m_dis      = 0
          weight     = 0.0_fp
          ndis_track = 1
          call findnmk(nlb          ,nub          ,mlb          ,mub   ,xz   ,yz   , &
                     & dps          ,s0           ,kcs          ,thick ,kmax , &
                     & nf_sour(1,IX),nf_sour(1,IY),nf_sour(1,IZ),n_tmp ,m_tmp,k_tmp, &
                     & kfsmn0       ,kfsmx0       ,dzs0         ,zmodel,gdp  )
          n_dis (1) = n_tmp
          m_dis (1) = m_tmp
          k_dis (1) = k_tmp
          weight(1) = wght
          !
          ! Momentum
          !
          if (nf_src_mom) then
             !
             ! Watch out: the discharge volume is distributed (multiplied with wght) over all source points
             ! The momentum is NOT distributed (multiplied with wght), but specified for each individual source point
             !
             call magdir_to_uv(alfas(n_dis(1),m_dis(1)), grdang          , &
                             & nf_sour(1,IUMAG)        , nf_sour(1,IUDIR), momu_tmp, momv_tmp)
             momu(1) = momu_tmp
             momv(1) = momv_tmp
          endif
          do iidis = 2, sour_cnt
             call findnmk(nlb              ,nub              ,mlb              ,mub   ,xz   ,yz   , &
                        & dps              ,s0               ,kcs              ,thick ,kmax , &
                        & nf_sour(iidis,IX),nf_sour(iidis,IY),nf_sour(iidis,IZ),n_tmp ,m_tmp,k_tmp, &
                        & kfsmn0           ,kfsmx0           ,dzs0             ,zmodel,gdp  )
             if (n_tmp/=n_dis(ndis_track) .or. m_tmp/=m_dis(ndis_track)) then
                ndis_track             = ndis_track + 1
                n_dis(ndis_track)      = n_tmp
                m_dis(ndis_track)      = m_tmp
                k_dis(ndis_track)      = k_tmp
             endif
             weight(ndis_track) = weight(ndis_track) + wght
             !
             ! Momentum
             !
             if (nf_src_mom) then
                call magdir_to_uv(alfas(n_dis(ndis_track),m_dis(ndis_track)), grdang              , &
                                & nf_sour(iidis,IUMAG)                      , nf_sour(iidis,IUDIR), momu_tmp, momv_tmp)
                momu(ndis_track) = momu(ndis_track) + momu_tmp
                momv(ndis_track) = momv(ndis_track) + momv_tmp
             endif
          enddo
       endif
       !
       ! Distribute sources discharges horizontal and vertical
       !
       thick_tot = 0.0_fp
       !
       if (.not. zmodel) then
          do iidis = 1, ndis_track
             n    = n_dis(iidis)
             m    = m_dis(iidis)
             wght = weight(iidis)
             if (allocated(k_dis)) then
                k_end_top  = k_dis(iidis)
                k_end_down = k_dis(iidis)
             endif
             do k = k_end_top, k_end_down
                if (disnf(n,m,k,idis) == 0.0_fp) then
                   thick_tot = thick_tot + wght*thick(k)
                endif
             enddo
          enddo
          do iidis = 1, ndis_track
             n    = n_dis(iidis)
             m    = m_dis(iidis)
             wght = weight(iidis)
             if (allocated(k_dis)) then
                k_end_top  = k_dis(iidis)
                k_end_down = k_dis(iidis)
             endif
             do k = k_end_top, k_end_down
                if (disnf(n,m,k,idis) == 0.0_fp) then
                   disnf(n,m,k,idis) = disnf(n,m,k,idis) + (nf_q_source+dis_tot)/(thick_tot/(wght*thick(k)))
                   do lcon = 1, lstsc
                      if ( flbcktemp(lcon) ) then
                         !
                         ! Background temperature: discharge with the temeprature last time step in discharge point
                         !
                         sournf(n,m,k,lcon,idis) = nf_q_source * r0(n,m,k,lcon) &
                                                 & / (thick_tot/(wght*thick(k)))
                      else
                         if (nf_const_operator == NFLCONSTOPERATOR_ABS) then
                            conc = nf_const(lcon)
                         else
                            !
                            ! nf_const_operator = NFLCONSTOPERATOR_EXC:
                            ! Add it to the concentration at the intake location
                            !
                            call coupled(nlb           ,nub           ,mlb           ,mub   ,conc_intake, &
                                       & r0            ,kmax          ,lstsci        ,lcon  ,thick      , &
                                       & m_intake(idis),n_intake(idis),k_intake(idis),s0    ,dps        , &
                                       & dzs0          ,kfsmn0        ,kfsmx0        ,zmodel,gdp        )
                            conc = nf_const(lcon) + conc_intake
                         endif
                         sournf(n,m,k,lcon,idis) = nf_q_source * conc &
                                                 & / (thick_tot/(wght*thick(k)))
                      endif
                   enddo
                   if (nf_src_mom) then
                      nf_src_momu(n,m,k,idis) = nf_src_momu(n,m,k,idis) + momu(iidis)
                      nf_src_momv(n,m,k,idis) = nf_src_momv(n,m,k,idis) + momv(iidis)
                   endif
                endif
             enddo
          enddo
       else
          !
          ! Z-model
          !
          do iidis = 1, ndis_track
             n    = n_dis(iidis)
             m    = m_dis(iidis)
             wght = weight(iidis)
             if (allocated(k_dis)) then
                k_end_top  = k_dis(iidis)
                k_end_down = k_dis(iidis)
             endif
             hhi  = 1.0_fp / max( s0(n,m)+real(dps(n,m),fp) , 0.01_fp )
             do k = k_end_top, k_end_down, -1
                if (k < kfsmn0(n,m)) cycle
                if (k > kfsmx0(n,m)) cycle
                if (disnf(n,m,k,idis) == 0.0_fp) then
                   thick_tot = thick_tot + wght*dzs0(n,m,k)*hhi
                endif
             enddo
          enddo
          do iidis = 1, ndis_track
             n    = n_dis(iidis)
             m    = m_dis(iidis)
             wght = weight(iidis)
             if (allocated(k_dis)) then
                k_end_top  = k_dis(iidis)
                k_end_down = k_dis(iidis)
             endif
             hhi  = 1.0_fp / max( s0(n,m)+real(dps(n,m),fp) , 0.01_fp )
             do k = k_end_top, k_end_down, -1
                if (k < kfsmn0(n,m)) cycle
                if (k > kfsmx0(n,m)) cycle
                if (disnf(n,m,k,idis) == 0.0_fp) then
                   disnf(n,m,k,idis) = disnf(n,m,k,idis) + (nf_q_source+dis_tot)/(thick_tot/(wght*dzs0(n,m,k)*hhi))
                   do lcon = 1, lstsci
                      if ( flbcktemp(lcon) ) then
                         !
                         ! Background temperature: discharge with the temeprature last time step in discharge point
                         !
                         sournf(n,m,k,lcon,idis) = nf_q_source * r0(n,m,k,lcon) &
                                                 & / (thick_tot/(wght*dzs0(n,m,k)*hhi))
                      else
                         if (nf_const_operator == NFLCONSTOPERATOR_ABS) then
                            conc = nf_const(lcon)
                         else
                            !
                            ! nf_const_operator = NFLCONSTOPERATOR_EXC:
                            ! Add it to the concentration at the intake location
                            !
                            call coupled(nlb           ,nub           ,mlb           ,mub   ,conc_intake, &
                                       & r0            ,kmax          ,lstsci        ,lcon  ,thick      , &
                                       & m_intake(idis),n_intake(idis),k_intake(idis),s0    ,dps        , &
                                       & dzs0          ,kfsmn0        ,kfsmx0        ,zmodel,gdp        )
                            conc = nf_const(lcon) + conc_intake
                         endif
                         sournf(n,m,k,lcon,idis) = nf_q_source * conc &
                                                 & / (thick_tot/(wght*dzs0(n,m,k)*hhi))
                      endif
                   enddo
                   if (nf_src_mom) then
                      nf_src_momu(n,m,k,idis) = nf_src_momu(n,m,k,idis) + momu(iidis)
                      nf_src_momv(n,m,k,idis) = nf_src_momv(n,m,k,idis) + momv(iidis)
                   endif
                endif
             enddo
          enddo
       endif
       !
       deallocate(n_dis , stat=ierror)
       deallocate(m_dis , stat=ierror)
       if (allocated(k_dis)) deallocate(k_dis, stat=ierror)
       deallocate(weight, stat=ierror)
       deallocate(momu  , stat=ierror)
       deallocate(momv  , stat=ierror)
    endif
    !
    ! Handle intakes
    !
    intake_cnt = size(nf_intake,1)
    if (intake_cnt > 0) then
       dis_per_intake = nf_q_intake / real(intake_cnt,fp)
       do irow = 1, intake_cnt
          call findnmk(nlb               ,nub               ,mlb               ,mub   ,xz    ,yz    , &
                     & dps               ,s0                ,kcs               ,thick ,kmax  ,  &
                     & nf_intake(irow,IX),nf_intake(irow,IY),nf_intake(irow,IZ),n_irow,m_irow,k_irow, &
                     & kfsmn0            ,kfsmx0            ,dzs0              ,zmodel,gdp   )
          disnf(n_irow,m_irow,k_irow,idis) = disnf(n_irow,m_irow,k_irow,idis) - dis_per_intake
          do lcon = 1, lstsci
             sournf(n_irow,m_irow,k_irow,lcon,idis) = sournf(n_irow,m_irow,k_irow,lcon,idis) - dis_per_intake * r0(n_irow,m_irow,k_irow,lcon)
          enddo
       enddo
    endif
end subroutine desa

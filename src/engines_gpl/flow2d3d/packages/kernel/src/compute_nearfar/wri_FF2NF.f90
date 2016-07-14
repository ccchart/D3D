subroutine wri_FF2NF(nlb     ,nub      ,mlb      ,mub       ,kmax   , &
                   & lstsci  ,lsal     ,ltem     ,idensform ,idis   , &
                   & time    ,taua     ,saleqs   ,temeqs    ,thick  , &
                   & sig     ,zk       ,kfu      ,kfv       , &
                   & alfas   ,s0       ,s1       ,u0        ,v0     , &
                   & r0      ,rho      ,dps      ,xz        ,yz     , &
                   & kfsmn0  ,kfsmx0   ,dzs0     ,filename  ,namcon , gdp    )
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
! NONE
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use properties
    use dfparall
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
    integer                        , pointer :: lundia
    integer                        , pointer :: lstsc
    integer ,dimension(:)          , pointer :: m_diff
    integer ,dimension(:)          , pointer :: n_diff
    integer ,dimension(:)          , pointer :: no_amb
    integer ,dimension(:,:)        , pointer :: m_amb
    integer ,dimension(:,:)        , pointer :: n_amb
    integer ,dimension(:)          , pointer :: m_intake
    integer ,dimension(:)          , pointer :: n_intake
    integer ,dimension(:)          , pointer :: k_intake
    real(fp),dimension(:,:)        , pointer :: x_amb
    real(fp),dimension(:,:)        , pointer :: y_amb
    real(fp),dimension(:)          , pointer :: x_intake
    real(fp),dimension(:)          , pointer :: y_intake
    real(fp),dimension(:)          , pointer :: q_diff
    real(fp),dimension(:,:)        , pointer :: const_diff
    real(fp),dimension(:)          , pointer :: d0
    real(fp),dimension(:)          , pointer :: h0
    real(fp),dimension(:)          , pointer :: sigma0
    integer                        , pointer :: lunsrc
    logical                        , pointer :: zmodel
!
! Global variables
!
    integer                                               , intent(in)  :: nlb
    integer                                               , intent(in)  :: nub
    integer                                               , intent(in)  :: mlb
    integer                                               , intent(in)  :: mub
    integer                                               , intent(in)  :: idis
    integer                                               , intent(in)  :: kmax
    integer                                               , intent(in)  :: lstsci
    integer                                               , intent(in)  :: lsal
    integer                                               , intent(in)  :: ltem
    integer                                               , intent(in)  :: idensform
    integer    , dimension(nlb:nub,mlb:mub)               , intent(in)  :: kfu        ! Description and declaration in esm_alloc_int.f90
    integer    , dimension(nlb:nub,mlb:mub)               , intent(in)  :: kfv        ! Description and declaration in esm_alloc_int.f90
    integer    , dimension(nlb:nub,mlb:mub)               , intent(in)  :: kfsmx0     ! Description and declaration in esm_alloc_int.f90
    integer    , dimension(nlb:nub,mlb:mub)               , intent(in)  :: kfsmn0     ! Description and declaration in esm_alloc_int.f90
    real(fp)                                              , intent(out) :: taua
    real(fp)                                              , intent(in)  :: saleqs
    real(fp)                                              , intent(in)  :: temeqs
    real(fp)   , dimension(nlb:nub,mlb:mub)               , intent(in)  :: alfas
    real(fp)   , dimension(nlb:nub,mlb:mub, kmax)         , intent(in)  :: dzs0       ! Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(nlb:nub,mlb:mub)               , intent(in)  :: s0
    real(fp)   , dimension(nlb:nub,mlb:mub, kmax)         , intent(in)  :: rho
    real(fp)   , dimension(nlb:nub,mlb:mub, kmax)         , intent(in)  :: u0
    real(fp)   , dimension(nlb:nub,mlb:mub, kmax)         , intent(in)  :: v0
    real(fp)   , dimension(nlb:nub,mlb:mub, kmax,lstsci)  , intent(in)  :: r0
    real(fp)   , dimension(kmax)                          , intent(in)  :: thick
    real(fp)   , dimension(kmax)                          , intent(in)  :: sig      !  Vertical coordinates of cell interfaces (SIGMA-MODEL)
    real(fp)   , dimension(0:kmax)                        , intent(in)  :: zk       !  Vertical coordinates of cell interfaces (Z-MODEL)
    real(prec) , dimension(nlb:nub,mlb:mub)               , intent(in)  :: dps
    real(fp)   , dimension(nlb:nub,mlb:mub)               , intent(in)  :: s1  
    real(fp)   , dimension(nlb:nub,mlb:mub)               , intent(in)  :: xz
    real(fp)   , dimension(nlb:nub,mlb:mub)               , intent(in)  :: yz
    character(256), dimension(3)                          , intent(in)  :: filename
    character(20), dimension(lstsci)                      , intent(in)  :: namcon   !  Description and declaration in esm_alloc_char.f90
!
! Local variables
!
    integer                                :: i
    integer                                :: ierror
    integer                                :: ilen
    integer                                :: ipnt
    integer                                :: istat
    integer                                :: k
    integer                                :: kenmu
    integer                                :: kenmv
    integer                                :: m
    integer                                :: md
    integer                                :: n
    integer                                :: nd
    integer                 , external     :: newlun
    integer                                :: npnt
    integer                                :: luntmp
    integer , dimension(:,:), allocatable  :: nm
    real(fp)                               :: deg2rad
    real(fp)                               :: csalfa
    real(fp)                               :: snalfa
    real(fp)                               :: drohj
    real(fp), dimension(:), allocatable    :: ha
    real(fp)                               :: hd
    real(fp)                               :: hint
    real(fp)                               :: pi
    real(fp)                               :: rad2deg
    real(fp)                               :: rhoam
    real(fp)                               :: rhoas
    real(fp)                               :: rhoab
    real(fp)                               :: taurel
    real(fp)                               :: time
    real(fp)                               :: uuu
    real(fp)                               :: ua
    real(fp)                               :: umag
    real(fp)                               :: vvv
    real(fp)                               :: sal
    real(fp)                               :: temp
    real(fp)                               :: dummy
    real(fp)                               :: add
    real(fp)                               :: rho0
    real(fp), dimension(:,:), allocatable  :: dzs0_amb
    real(fp), dimension(:,:), allocatable  :: xy
    character(1)                           :: tab
    character(1)                           :: stype1
    character(1)                           :: stype2
    character(3)                           :: c_inode
    character(12)                          :: crhoam
    character(12)                          :: crhoas
    character(12)                          :: crhoab
    character(12)                          :: chint
    character(12)                          :: cdrohj
    character(12)                          :: ctaua
    character(12)                          :: inttostring
    character(20)                          :: rundat       ! Current date and time containing a combination of DATE and TIME
    character(1000)                        :: string
    character(256), external               :: windows_path
    type(tree_data)              , pointer :: outfile_ptr
    type(tree_data)              , pointer :: cosumo_ptr
    type(tree_data)              , pointer :: subgrid_ptr
    type(tree_data)              , pointer :: node_ptr
    type(tree_data)              , pointer :: subnode_ptr
    type(tree_data)              , pointer :: data_ptr
!
!! executable statements -------------------------------------------------------
!
    lundia         => gdp%gdinout%lundia
    lstsc          => gdp%d%lstsc
    m_diff         => gdp%gdnfl%m_diff
    n_diff         => gdp%gdnfl%n_diff
    no_amb         => gdp%gdnfl%no_amb
    m_amb          => gdp%gdnfl%m_amb
    n_amb          => gdp%gdnfl%n_amb
    m_intake       => gdp%gdnfl%m_intake
    n_intake       => gdp%gdnfl%n_intake
    k_intake       => gdp%gdnfl%k_intake
    x_amb          => gdp%gdnfl%x_amb
    y_amb          => gdp%gdnfl%y_amb
    x_intake       => gdp%gdnfl%x_intake
    y_intake       => gdp%gdnfl%y_intake
    q_diff         => gdp%gdnfl%q_diff
    const_diff     => gdp%gdnfl%const_diff
    d0             => gdp%gdnfl%d0
    h0             => gdp%gdnfl%h0
    sigma0         => gdp%gdnfl%sigma0
    zmodel         => gdp%gdprocs%zmodel
    !
    write(c_inode(1:3),'(i3.3)') inode
    !
    pi      = acos(-1.0_fp)
    rad2deg = 180.0_fp / pi
    deg2rad = pi / 180.0_fp
    tab     = char(9)
    !
    ! Read the general diffusor characteritics from cormix input file
    ! Parallel: Only the master partition executes this
    ! the n_diff(idis), m_diff(idis), etc. will be global indices
    !
    allocate(dzs0_amb(kmax,no_amb(idis)), stat=ierror)
    allocate(ha      (     no_amb(idis)), stat=ierror)
    if (zmodel) then
       do i = 1, no_amb(idis)
          dzs0_amb(:,i) = dzs0(n_amb(idis,i), m_amb(idis,i), :)
       enddo
    else
       dzs0_amb = -999.0_fp
    endif
    !
    ! Compute the depths
    !
    do i = 1, no_amb(idis)
       ha(i) = s0(n_amb(idis,i), m_amb(idis,i))+real(dps(n_amb(idis,i), m_amb(idis,i)),fp)
    enddo
    hd = s0(n_diff(idis), m_diff(idis))+real(dps(n_diff(idis), m_diff(idis)),fp)
    !
    ! Compute depth averaged velocity magnitude and direction
    !
    uuu = 0.0_fp  
    vvv = 0.0_fp
    if (.not. zmodel) then
       !
       ! Sigma-model
       !
       do k = 1, kmax
          uuu      = uuu + 0.5_fp * (u0(n_amb(idis,1), m_amb(idis,1) ,k) + u0(n_amb(idis,1), m_amb(idis,1)-1 ,k))*thick(k)
          vvv      = vvv + 0.5_fp * (v0(n_amb(idis,1), m_amb(idis,1) ,k) + v0(n_amb(idis,1)-1, m_amb(idis,1) ,k))*thick(k)
       enddo
    else
       !
       ! Z-model
       ! We now take the velocity at the k-level of the corresponding cell centre.
       ! If we loop over the kfumn0 to kfumx0 of the velocity points (corresponding to n_amb(idis), m_amb(idis) and n_amb(idis), m_amb(idis)-1),
       ! and divide by dzu0/hu and dzv0/hv, would it then be more accurate?
       !
       do k = kfsmn0(n_amb(idis,1), m_amb(idis,1)), kfsmx0(n_amb(idis,1), m_amb(idis,1))
          uuu      = uuu + 0.5_fp * (u0(n_amb(idis,1), m_amb(idis,1) ,k) + u0(n_amb(idis,1), m_amb(idis,1)-1 ,k))*dzs0(n_amb(idis,1), m_amb(idis,1),k)/max(ha(1), 0.01_fp)
          vvv      = vvv + 0.5_fp * (v0(n_amb(idis,1), m_amb(idis,1) ,k) + v0(n_amb(idis,1)-1, m_amb(idis,1) ,k))*dzs0(n_amb(idis,1), m_amb(idis,1),k)/max(ha(1), 0.01_fp)
       enddo
    endif
    umag = sqrt (uuu*uuu + vvv*vvv)
    taua = atan2(vvv,uuu)*rad2deg + alfas(n_amb(idis,1), m_amb(idis,1))
    taua = mod(taua + 360.0_fp,360.0_fp)
    ua   = umag
    !
    ! Density profile classification (Cormixtype)
    !
    call determine_densprof(kmax           ,thick          ,s0(n_amb(idis,1),m_amb(idis,1)), &
                          & real(dps(n_amb(idis,1),m_amb(idis,1)),fp),  rho(n_amb(idis,1),m_amb(idis,1),:) , &
                          & ha(1)          ,hd             ,stype1         ,stype2               ,rhoam         , &
                          & rhoas          ,rhoab          ,hint           ,drohj                , &
                          & kfsmn0(n_amb(idis,1), m_amb(idis,1)),  kfsmx0(n_amb(idis,1),m_amb(idis,1)) ,dzs0_amb    ,zmodel         )
    !
    ! Compute the density of the discharged water
    !
    if (lsal /= 0) then
       sal  = const_diff(idis,2)
       call coupled(nlb           ,nub           ,mlb           ,mub   ,add   , &
                  & r0            ,kmax          ,lstsci        ,lsal  ,thick , &
                  & m_intake(idis),n_intake(idis),k_intake(idis),s0    ,dps   , &
                  & dzs0          ,kfsmn0        ,kfsmx0        ,zmodel,gdp   )
       sal = sal + add
    else
       sal = saleqs
    endif
    if (ltem /= 0) then
       temp = const_diff(idis,1)
       call coupled(nlb           ,nub           ,mlb           ,mub   ,add   , &
                  & r0            ,kmax          ,lstsci        ,ltem  ,thick , &
                  & m_intake(idis),n_intake(idis),k_intake(idis),s0    ,dps   , &
                  & dzs0          ,kfsmn0        ,kfsmx0        ,zmodel,gdp   )
       temp = temp + add
    else
       temp = temeqs
    endif
    select case (idensform)
       case( dens_Eckart )
          call dens_eck    (temp, sal ,rho0, dummy, dummy)
       case( dens_Unesco)
          call dens_unes   (temp, sal ,rho0, dummy, dummy)
       case( dens_NaClSol)
          call dens_nacl   (temp, sal ,rho0, dummy, dummy)
    end select
    !
    ! Make character strings from all requested input
    !
    if (stype1 == 'U') then
       write(crhoam(1:12),'(f12.3)') rhoam
       crhoas = '-'
       crhoab = '-'
       stype2 = '-'
       chint  = '-'
       cdrohj = '-'
    else
       crhoam ='-'
       write (crhoas(1:12),'(f12.3)') rhoas
       write (crhoab(1:12),'(f12.3)') rhoab
       if (stype2 == 'C') then
          write (chint (1:12),'(f12.3)') hint
          write (cdrohj(1:12),'(f12.3)') drohj
       else
          chint  = '-'
          cdrohj = '-'
       endif
    endif
    !
    ! sigma0 given as direction relative to north in stead of main flow direction; 0, pointing to east, 90 pointing to north etc.
    ! ctaua is port direction relative to main flow direction
    !
    taua = mod(taua,360.0_fp)
    taurel = mod(sigma0(idis) - taua + 360.0_fp,360.0_fp)
    if (taurel > 179.0_fp .and. taurel < 181.0_fp) then
       taurel = 179.0_fp
    endif
    write (ctaua (1:12),'(f12.3)') taurel
    !
    ! Fill new tree with data to be written
    !
    call tree_create('COSUMO Input, created by Delft3D-FLOW', outfile_ptr)
    call tree_put_data(outfile_ptr, transfer(trim(adjustl(filename(1))),node_value), 'STRING:XML')
    call tree_create_node(outfile_ptr, '?xml version="1.0" encoding="utf-8"?', node_ptr)
    call tree_create_node(outfile_ptr, 'COSUMO', cosumo_ptr)
    call tree_create_node(cosumo_ptr, 'fileVersion', node_ptr)
    call tree_put_data(node_ptr, transfer("0.3",node_value), 'STRING:XMLDATA')
    call tree_create_node(cosumo_ptr, 'comm', node_ptr)
    !
    ! Filenames should always be written in Windows style, even on Linux,
    ! Because Cosumo is reading/using it and runs on Windows
    !
    call tree_create_node(node_ptr, 'Filename', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(windows_path(filename(1)))),node_value), 'STRING:XMLDATA')
    call tree_create_node(node_ptr, 'waitForFile', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(windows_path(filename(2)))),node_value), 'STRING:XMLDATA')
    call tree_create_node(node_ptr, 'FFrundir', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(windows_path(filename(3)))),node_value), 'STRING:XMLDATA')
    string = trim(gdp%runid) // '.mdf'
    call tree_create_node(node_ptr, 'FFinputFile', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    call tree_create_node(node_ptr, 'FFuniqueID', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(gdp%uniqueid)),node_value), 'STRING:XMLDATA')
    !
    call tree_create_node(cosumo_ptr, 'SubgridModel', subgrid_ptr)
    write(string,'(i0)') idis
    call tree_create_node(subgrid_ptr, 'SubgridModelNr', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') time/60.0_fp
    call tree_create_node(subgrid_ptr, 'TIME', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    call tree_create_node(subgrid_ptr, 'constituentsNames', subnode_ptr)
    write(string,'(i0)') lstsc
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLNUMDATALINES')
    do i=1,lstsc
       write(inttostring,'(i0)') i
       call tree_create_node(subnode_ptr, trim(inttostring), data_ptr)
       call tree_put_data(data_ptr, transfer(trim(namcon(i)),node_value), "STRING:XMLDATALINE")
    enddo
    !
    call tree_create_node(subgrid_ptr, 'cormix', node_ptr)
    write(string,'(e24.17)') ha(1)
    call tree_create_node(node_ptr, 'HA', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') hd
    call tree_create_node(node_ptr, 'HD', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') max(ua,0.001_fp)
    call tree_create_node(node_ptr, 'UA', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(i0)') idis
    call tree_create_node(node_ptr, 'UorS', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(stype1)),node_value), 'STRING:XMLDATA')
    write(string,'(i0)') idis
    call tree_create_node(node_ptr, 'RHOAM', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(crhoam)),node_value), 'STRING:XMLDATA')
    !
    call tree_create_node(node_ptr, 'STYPE', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(stype2)),node_value), 'STRING:XMLDATA')
    !
    call tree_create_node(node_ptr, 'RHOAS', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(crhoas)),node_value), 'STRING:XMLDATA')
    !
    call tree_create_node(node_ptr, 'RHOAB', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(crhoab)),node_value), 'STRING:XMLDATA')
    !
    call tree_create_node(node_ptr, 'HINT', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(chint)),node_value), 'STRING:XMLDATA')
    !
    call tree_create_node(node_ptr, 'DROHJ', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(cdrohj)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') q_diff(idis)
    call tree_create_node(node_ptr, 'Q0', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') rho0
    call tree_create_node(node_ptr, 'RHO0', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') d0(idis)
    call tree_create_node(node_ptr, 'D0', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    !
    call tree_create_node(node_ptr, 'PHI', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(ctaua)),node_value), 'STRING:XMLDATA')
    call tree_create_node(node_ptr, '<!-- S1: zero=reference level, down is positive -->', subnode_ptr)
    write(string,'(e24.17)') - s1(n_diff(idis), m_diff(idis))
    call tree_create_node(node_ptr, 'S1', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') dps(n_diff(idis), m_diff(idis))
    call tree_create_node(node_ptr, 'h_dps', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') xz(n_diff(idis), m_diff(idis))
    call tree_create_node(node_ptr, 'x_diff', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') yz(n_diff(idis), m_diff(idis))
    call tree_create_node(node_ptr, 'y_diff', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    write(string,'(e24.17)') taua
    call tree_create_node(node_ptr, 'taua', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLDATA')
    !
    ! Intake
    !
    call tree_create_node(subgrid_ptr, 'FFIntake', node_ptr)
    call tree_create_node(node_ptr, '<!-- Z: zero=reference level, down is positive -->', subnode_ptr)
    allocate(nm(1,2), stat=ierror)
    allocate(xy(1,2), stat=ierror)
    nm(1,1) = n_intake(idis)
    nm(1,2) = m_intake(idis)
    xy(1,1) = x_intake(idis)
    xy(1,2) = y_intake(idis)
    !
    call writePointInfoToFF2NF()
    !
    deallocate(nm, stat=ierror)
    deallocate(xy, stat=ierror)
    !
    ! Ambients
    !
    call tree_create_node(subgrid_ptr, 'FFAmbient', node_ptr)
    call tree_create_node(node_ptr, '<!-- Z: zero=reference level, down is positive -->', subnode_ptr)
    allocate(nm(no_amb(idis),2), stat=ierror)
    allocate(xy(no_amb(idis),2), stat=ierror)
    do i=1,no_amb(idis)
       nm(i,1) = n_amb(idis,i)
       nm(i,2) = m_amb(idis,i)
       xy(i,1) = x_amb(idis,i)
       xy(i,2) = y_amb(idis,i)
    enddo
    !
    call writePointInfoToFF2NF()
    !
    deallocate(nm, stat=ierror)
    deallocate(xy, stat=ierror)
    !
    write(*,'(3a)') "Writing file '", trim(filename(1)), "' ..."
    luntmp = newlun(gdp)
    open (luntmp, file=trim(filename(1)), status='new', iostat=istat)
    if (istat /= 0) then
       write(lundia,'(3a)') "ERRROR: file '", trim(filename(1)), "' already exists."
       call d3stop(1,gdp)
    endif
    call prop_write_xmlfile(luntmp, outfile_ptr, 0, istat)
    close (luntmp)
    !
    deallocate(dzs0_amb, stat=ierror)         
    deallocate(ha      , stat=ierror)         
    call tree_destroy(outfile_ptr)



contains


subroutine writePointInfoToFF2NF
    integer :: k_top
    integer :: k_down
    integer :: k_incr
    !
    npnt = size(nm,1)
    write(string,'(i0)') kmax * npnt
    call tree_create_node(node_ptr, 'XYZ', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLNUMDATALINES')
    do ipnt=1,npnt
       if (zmodel) then
          k_top  = kfsmx0(nm(ipnt,1),nm(ipnt,2))
          k_down = kfsmn0(nm(ipnt,1),nm(ipnt,2))
          k_incr = -1
       else
          k_top  = 1
          k_down = kmax
          k_incr = 1
       endif
       hd = s0(nm(ipnt,1),nm(ipnt,2))+real(dps(nm(ipnt,1),nm(ipnt,2)),fp)
       do k = k_top, k_down, k_incr
          if (zmodel) then
             write(string,'(3(e24.17,x))') xy(ipnt,1), xy(ipnt,2), -(zk(k-1)+zk(k))/2.0_fp
          else
             write(string,'(3(e24.17,x))') xy(ipnt,1), xy(ipnt,2), -hd*sig(k)
          endif
          write(inttostring,'(i0)') (ipnt-1)*k + k
          call tree_create_node(subnode_ptr, trim(inttostring), data_ptr)
          call tree_put_data(data_ptr, transfer(trim(string),node_value), "STRING:XMLDATALINE")
       enddo
    enddo
    !
    write(string,'(i0)') npnt
    call tree_create_node(node_ptr, 'waterLevel', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLNUMDATALINES')
    do ipnt=1,npnt
       write(string,'(e24.17)') - s0(nm(ipnt,1), nm(ipnt,2))
       write(inttostring,'(i0)') ipnt
       call tree_create_node(subnode_ptr, trim(inttostring), data_ptr)
       call tree_put_data(data_ptr, transfer(trim(string),node_value), "STRING:XMLDATALINE")
    enddo
    !
    write(string,'(i0)') kmax * npnt
    call tree_create_node(node_ptr, 'XYvelocity', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLNUMDATALINES')
    do ipnt=1,npnt
       !
       ! Copied from wave::flow2wav
       !
       n      = nm(ipnt,1)
       m      = nm(ipnt,2)
       nd     = max(1, n - 1)
       md     = max(1, m - 1)
       kenmu  = max(1, kfu(n,m)+kfu(n ,md))
       kenmv  = max(1, kfv(n,m)+kfv(nd,m ))
       csalfa = cos(alfas(n,m)*deg2rad)
       snalfa = sin(alfas(n,m)*deg2rad)
       if (zmodel) then
          k_top  = kfsmx0(nm(ipnt,1),nm(ipnt,2))
          k_down = kfsmn0(nm(ipnt,1),nm(ipnt,2))
          k_incr = -1
       else
          k_top  = 1
          k_down = kmax
          k_incr = 1
       endif
       do k = k_top, k_down, k_incr
          uuu = (u0(n,m,k)*real(kfu(n,m),fp) + u0(n ,md,k)*real(kfu(n ,md),fp)) / real(kenmu,fp)
          vvv = (v0(n,m,k)*real(kfv(n,m),fp) + v0(nd,m ,k)*real(kfv(nd,m ),fp)) / real(kenmv,fp)
          write(string,'(2(e24.17,x))') uuu*csalfa - vvv*snalfa, uuu*snalfa + vvv*csalfa
          write(inttostring,'(i0)') (ipnt-1)*k + k
          call tree_create_node(subnode_ptr, trim(inttostring), data_ptr)
          call tree_put_data(data_ptr, transfer(trim(adjustl(string)),node_value), "STRING:XMLDATALINE")
       enddo
    enddo
    !
    write(string,'(i0)') kmax * npnt
    call tree_create_node(node_ptr, 'rho', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLNUMDATALINES')
    do ipnt=1,npnt
       if (zmodel) then
          k_top  = kfsmx0(nm(ipnt,1),nm(ipnt,2))
          k_down = kfsmn0(nm(ipnt,1),nm(ipnt,2))
          k_incr = -1
       else
          k_top  = 1
          k_down = kmax
          k_incr = 1
       endif
       do k = k_top, k_down, k_incr
          write(string,'(e24.17)') rho(nm(ipnt,1), nm(ipnt,2), k)
          write(inttostring,'(i0)') (ipnt-1)*k + k
          call tree_create_node(subnode_ptr, trim(inttostring), data_ptr)
          call tree_put_data(data_ptr, transfer(trim(adjustl(string)),node_value), "STRING:XMLDATALINE")
       enddo
    enddo
    !
    write(string,'(i0)') kmax * npnt
    call tree_create_node(node_ptr, 'constituents', subnode_ptr)
    call tree_put_data(subnode_ptr, transfer(trim(adjustl(string)),node_value), 'STRING:XMLNUMDATALINES')
    do ipnt=1,npnt
       if (zmodel) then
          k_top  = kfsmx0(nm(ipnt,1),nm(ipnt,2))
          k_down = kfsmn0(nm(ipnt,1),nm(ipnt,2))
          k_incr = -1
       else
          k_top  = 1
          k_down = kmax
          k_incr = 1
       endif
       do k = k_top, k_down, k_incr
          string = ' '
          do i=1,lstsc
             write(string,'(a,e24.17,x)') trim(string), r0(nm(ipnt,1), nm(ipnt,2), k, i)
          enddo
          write(inttostring,'(i0)') (ipnt-1)*k + k
          call tree_create_node(subnode_ptr, trim(inttostring), data_ptr)
          call tree_put_data(data_ptr, transfer(trim(adjustl(string)),node_value), "STRING:XMLDATALINE")
       enddo
    enddo
end subroutine writePointInfoToFF2NF
end subroutine wri_FF2NF

function windows_path(inpath) result(outpath)
!
! return value
!
character(256) :: outpath
!
! arguments
!
character(*), intent(in)  :: inpath
!
! locals
!
integer      :: i
character(1) :: bslash = '\'
character(1) :: fslash = '/'
!
! body
if (inpath(1:1)==fslash .and. inpath(3:3)== fslash) then
   ! Replace /p by p:
   outpath = inpath
   outpath(1:1) = outpath(2:2)
   outpath(2:2) = ':'
else if (inpath(1:4)=='/opt') then
   ! Replace /opt by p:\h6\opt
   outpath = "p:\h6\opt"
   outpath(10:) = inpath(5:)
else if (inpath(1:4)=='/mnt') then
   ! Replace /mnt by d:
   outpath = "d:"
   outpath(3:) = inpath(5:)
else
   outpath = inpath
endif
!
! Replace / by \
do i=1,len_trim(outpath)
   if (outpath(i:i) == fslash) then
      outpath(i:i) = bslash
   endif
enddo
end function windows_path

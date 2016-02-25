subroutine near_field(u0     ,v0     ,rho    ,thick  , &
                    & kmax   ,alfas  ,dps    ,s0     , &
                    & lstsci ,lsal   ,ltem   ,xz     , &
                    & yz     ,nmmax  ,kcs    ,kcs_nf , &
                    & r0     ,time   ,saleqs ,temeqs , &
                    & s1     ,kfsmn0 ,kfsmx0 ,dzs0   , &
                    & gdp   )
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
!    Function: Converts flow results to cormix input
!              nog wel doen:
!              1) interpoleren naar equidistante laagverdeling jet3d
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    !
    use globaldata
    use dfparall
    use xml_data_discharge_def
    use xmlparse
    use getdata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    ! They replace the  include igd / include igp lines
    !
    integer                            , pointer :: idensform
    integer                            , pointer :: ifirst
    integer                            , pointer :: no_dis
    integer       , dimension(:)       , pointer :: m_diff
    integer       , dimension(:)       , pointer :: n_diff
    integer       , dimension(:)       , pointer :: m_amb
    integer       , dimension(:)       , pointer :: n_amb
    real(fp)      , dimension(:)       , pointer :: q_diff
    real(fp)      , dimension(:)       , pointer :: t0_diff
    real(fp)      , dimension(:)       , pointer :: s0_diff
    real(fp)      , dimension(:)       , pointer :: rho0_diff
    real(fp)      , dimension(:)       , pointer :: d0
    real(fp)      , dimension(:)       , pointer :: h0
    real(fp)      , dimension(:)       , pointer :: sigma0
    real(fp)      , dimension(:)       , pointer :: theta0
    real(fp)      , dimension(:,:,:)   , pointer :: disnf
    real(fp)      , dimension(:,:,:,:) , pointer :: sournf

    character(256), dimension(:,:)     , pointer :: basecase
    character(256)                     , pointer :: nflmod
    integer                            , pointer :: lundia
!
! Global variables
!
    integer                                                    , intent(in)  :: kmax     !  Description and declaration in
    integer                                                    , intent(in)  :: lstsci
    integer                                                    , intent(in)  :: lsal
    integer                                                    , intent(in)  :: ltem
    integer                                                    , intent(in)  :: nmmax
    real(fp)                                                   , intent(in)  :: time
    real(fp)                                                   , intent(in)  :: saleqs
    real(fp)                                                   , intent(in)  :: temeqs
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kcs      !
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)                            :: kcs_nf   !
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kfsmn0   !  Description and declaration in esm_alloc_int.f90
    integer    , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: kfsmx0   !  Description and declaration in esm_alloc_int.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: alfas    !
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: s0       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: s1       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: xz       !  Description and declaration in
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: yz       !  Description and declaration in
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)        , intent(in)  :: dzs0     !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)        , intent(in)  :: rho      !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)        , intent(in)  :: u0       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)        , intent(in)  :: v0       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax,lstsci) , intent(in)  :: r0       !  Description and declaration in esm_alloc_real.f90
    real(fp)   , dimension(kmax)                               , intent(in)  :: thick    !  Description and declaration in esm_alloc_real.f90
    real(prec) , dimension(gdp%d%nmlb:gdp%d%nmub)              , intent(in)  :: dps      !  Description and declaration in esm_alloc_real.f90
!
! Local variables
!
    integer                         :: i
    integer                         :: idis
    integer                         :: istat
    real(fp)                        :: flwang
    real(fp)                        :: signx
    real(fp)                        :: taua
    real(fp),dimension(8)           :: linkinf
    logical                         :: corend
    logical                         :: first_time
    logical                         :: error_reading
    character*3                     :: c_inode
    character*3                     :: c_idis
    character*256, dimension(3)     :: filename
    character*14                    :: cctime
    real(fp)                        :: xstart
    real(fp)                        :: xend
    real(fp)                        :: ystart
    real(fp)                        :: yend
    
    integer, parameter              :: no_jet_max = 10000
    integer                         :: no_val
    real(fp), dimension(no_jet_max) :: x_jet, y_jet, z_jet, s_jet, h_jet, b_jet
    
    character(len=40), dimension(2,10) :: attribs
    character(len=80), dimension(200)  :: data
    logical                            :: error
    integer                            :: ii
    integer                            :: jj
    type(table)                        :: mytable
    
!
!! executable statements -------------------------------------------------------
!

    idensform      => gdp%gdphysco%idensform
    nflmod         => gdp%gdnfl%nflmod
    no_dis         => gdp%gdnfl%no_dis
    !
    ! Update local pointers
    !
    m_diff         => gdp%gdnfl%m_diff
    n_diff         => gdp%gdnfl%n_diff
    m_amb          => gdp%gdnfl%m_amb
    n_amb          => gdp%gdnfl%n_amb
    q_diff         => gdp%gdnfl%q_diff
    t0_diff        => gdp%gdnfl%t0_diff
    s0_diff        => gdp%gdnfl%s0_diff
    rho0_diff      => gdp%gdnfl%rho0_diff
    d0             => gdp%gdnfl%d0
    h0             => gdp%gdnfl%h0
    sigma0         => gdp%gdnfl%sigma0
    theta0         => gdp%gdnfl%theta0
    basecase       => gdp%gdnfl%basecase
    disnf          => gdp%gdnfl%disnf
    sournf         => gdp%gdnfl%sournf
    lundia         => gdp%gdinout%lundia
    !    
    write(c_inode,'(i3.3)') inode
    !
    ! Convert flow results (velocities and densities) at (mdiff,ndiff) to nearfield input
    ! and write cormix/nrfield input file
    !
    select case (nflmod)
       case('corjet')
          !
          ! Convert flow results to input for cormix and write to input file
          !
          call wri_cormix(u0    ,v0    ,rho   ,thick ,kmax  ,dps   , &
                        & s0    ,alfas ,gdp                        )
          !
          ! Do the Cormix simulation
          !
          corend = .false.
          call util_system('corjet.bat')
          do while (.not. corend)
             inquire (file='corjetrun.end',exist=corend)
          enddo
          !
          ! Finally convert cormix results to flow input
          !
          call corjet2flow(thick  ,kmax   ,dps   ,s0   ,disnf    ,sournf , &
                         & lstsci ,lsal   ,ltem  ,xz   ,yz       ,nmmax  , &
                         & kcs    ,time   ,gdp   )
       case('cortime')
          !
          ! Read the general information from the corinp.dat file every time a cortime simulation is requested.
          ! This allows for restarting of cormix on a different pc (request Robin Morelissen)
          !
          call corinp_gen(idensform,gdp)
          !
          ! Convert flow results to input for cormix and write to input file
          !
          write(cctime,'(f14.3)') time/60.0_fp
          !
          do idis = 1, no_dis
             filename(1) = trim(basecase(idis,1))//'.cmx'
             filename(2) = trim(basecase(idis,2))//'_time-step-'//trim(adjustl(cctime))//'.prd'
             filename(3) = trim(basecase(idis,2))//'_trajectory_time-step-'//trim(adjustl(cctime))//'_'//c_inode//'.tek'
             !
             call wri_cortim(u0    ,v0    ,rho   ,thick ,kmax  ,dps    , &
                           & s0    ,alfas ,time  ,taua         ,r0     , &
                           & lstsci,lsal  ,ltem  ,idensform    ,saleqs , &
                           & temeqs,idis  ,filename(1)         ,linkinf, &
                           & kfsmn0,kfsmx0,dzs0  ,gdp    )
             !
             ! Wait for the Cortime simulation to end (use existance of output file as indicator)
             !
             !            corend   = .false.
             !            do while (.not. corend)
             !               inquire (file=filename(2),exist=corend)
             !            enddo
             !
             ! Finally convert cortime results to flow input
             !
             call wait_until_finished(filename(2),gdp)
             !
             call cortim2flow(thick  ,kmax   ,dps   ,s0   ,r0       ,         &
                            & lstsci ,lsal   ,ltem  ,xz   ,yz       ,nmmax  , &
                            & kcs    ,filename      ,taua           ,idis   , &
                            & linkinf,gdp           )
          enddo
          !
       case('generic')
          !
          ! Read the general information from the nff2ff.xml file every time a cortime simulation is requested.
          ! This allows for restarting of cormix on a different pc (request Robin Morelissen)
          !    
                call corinp_gen2(idensform,gdp)
                !
                ! Convert flow results to input for cormix and write to input file
                !
                write(cctime,'(f14.3)') time/60.0_fp
                !
                do idis = 1, no_dis
                   !
                   ! Add SubGridModel number to filename to prevent overwriting
                   !
                   write(c_idis,'(i3.3)') idis
                   !
                   filename(1) =trim(gdp%gdnfl%base_path)//'FF2NF_'//trim(gdp%runid)//'_'//c_inode//'_SubMod'//c_idis//'_'//trim(adjustl(cctime))//'.txt'
                   filename(2) =trim(basecase(idis,1))//'COSUMO\NF2FF\NF2FF_'//trim(gdp%runid)//'_'//c_inode//'_SubMod'//c_idis//'_'//trim(adjustl(cctime))//'.txt'
                   filename(3) =trim(basecase(idis,1))
                   !
                   ! You should get the filenames (the dirs) from the COSUMOsettings.xml
                   !
                   call wri_FF2NF(u0       ,v0      ,rho       ,thick  ,kmax   ,dps    , &
                                & s0       ,alfas   ,time      ,taua   ,r0     ,lstsci , &
                                & lsal     ,ltem    ,idensform ,saleqs ,temeqs ,idis   , &
                                & filename ,s1        ,xz     ,yz     , &
                                & kfsmn0   ,kfsmx0  ,dzs0      ,gdp    )
                   !
                   call wait_until_finished(filename(2),gdp)
                   !
                   no_val=size(x_jet)
                   call nf_2_flow(filename(2),x_jet,y_jet,z_jet,s_jet,h_jet,b_jet, no_val)               
                   !
                   ! Fill sources and sinks following the Desa Method of Prof. Lee
                   !
                   call desa(x_jet   ,y_jet   ,z_jet   ,s_jet   ,no_val  , &
                           & kcs     ,xz      ,yz      ,dps     ,s0      , &
                           & nmmax   ,thick   ,kmax    ,lstsci  ,lsal    , &
                           & ltem    ,h_jet   ,b_jet   ,idis    , &
                           & xstart  ,xend    ,ystart  ,yend    ,r0      , &
                           & linkinf ,kfsmn0  ,kfsmx0  ,dzs0    ,gdp     )
               enddo
               !
               !          call read_xml_file_discharge_def( "nf2ff.xml", lurep = lundia, errout = error_reading )
               !          do idis = 1,size(discharges)
               !              write(*,*) discharges(idis)%name
               !          enddo 
               !          
               !          
               !          write(cctime,'(f14.3)') time/60.0_fp
               !          file_nf_inp = discharges(idis)%inputdir//discharges(idis)%name//trim(gdp%runid)//'_'//c_inode//'_'//trim(adjustl(cctime))//'.coupledinp'
               !          file_nf_out = discharges(idis)%outputdir//discharges(idis)%name//trim(gdp%runid)//'_'//c_inode//'_'//trim(adjustl(cctime))//'.coupledout'
               !
               !
               !
               !
               !          do idis = 1, no_dis
               !
               !             !
               !             ! Convert flow results to input for nearfield model andwait for near field simulation
               !             ! to be finished
               !             !
               !
               !             call wri_nf_inp(u0    ,v0    ,rho   ,thick ,kmax  ,dps    , &
               !                           & s0    ,alfas ,time  ,taua         ,r0     , &
               !                           & lstsci,lsal  ,ltem  ,idensform    ,saleqs , &
               !                           & temeqs,idis  ,file_nf_inp         ,linkinf, &
               !                           & discharges, gdp   )
               !                           !(file_nf_inp, ?????????????)
               !                           
               !             call wait_until_finished(file_nf_out)
               !             call nf_2_flow          (file_nf_out,x_jet,y_jet,z_jet,s_jet,bh_jet,bv_jet)
               !
               !             call desa(x_jet   ,y_jet    ,z_jet   ,s_jet   ,no_val  , &
               !            &               kcs     ,xz       ,yz      ,dps     ,s0      , &
               !            &               nmmax   ,thick    ,kmax    ,lstsci  ,lsal    , &
               !            &               ltem    ,bv_jet  ,bh_jet   ,idis    ,r0      , &
               !            &               gdp     )
       case ('jet3d')
          !
          ! Convert flow results to input for jet3d and write to input file
          !
          call wri_jet3d(u0    ,v0    ,rho    ,thick ,kmax      ,dps   , &
                       & s0    ,alfas ,flwang ,signx ,idensform ,time  ,gdp   )
          !
          ! Do the Jet3d simulation
          !
          corend = .false.
          call util_system('jet3d.bat')
          do while (.not. corend)
             inquire (file='jet3drun.end',exist=corend)
          enddo
          !
          ! Finally convert Jet3D results to flow input
          !    The appoach followed is the DESA approach as suggested by Prof. Lee
          !    Clean up the mess afterwards
          !
          call jet3d2flow(thick  ,kmax   ,dps  ,s0   ,r0       ,        &
                        & lstsci ,lsal   ,ltem ,xz   ,yz       ,nmmax  ,&
                        & kcs    ,flwang ,signx,time ,linkinf  ,gdp  )
          call util_system('del str3dinp.xxx')
          call util_system('del str3dprt.xxx')
          call util_system('del str3dtek.xxx')
          call util_system('del jet3drun.end')
       case ('nrfield')
          !
          ! Nothing (yet)
          !
    end select
end subroutine near_field

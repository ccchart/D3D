subroutine wri_cortim(u0    ,v0    ,rho    ,thick  ,kmax   ,dps    ,&
                    & s0    ,alfas ,time   ,taua   ,r0     ,lstsci ,&
                    & lsal  ,ltem  ,idensform      ,saleqs ,temeqs ,&
                    & idis  ,filename      ,linkinf,gdp   )

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
!    Function: Writes input for cortim
!
! Method used:
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
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
    integer ,dimension(:)          , pointer :: m_diff
    integer ,dimension(:)          , pointer :: n_diff
    integer ,dimension(:)          , pointer :: m_amb
    integer ,dimension(:)          , pointer :: n_amb
    integer ,dimension(:)          , pointer :: m_intake
    integer ,dimension(:)          , pointer :: n_intake
    integer ,dimension(:)          , pointer :: k_intake

    real(fp),dimension(:)          , pointer :: q_diff
    real(fp),dimension(:)          , pointer :: t0_diff
    real(fp),dimension(:)          , pointer :: s0_diff
    real(fp),dimension(:)          , pointer :: d0
    real(fp),dimension(:)          , pointer :: h0
    real(fp),dimension(:)          , pointer :: sigma0
    integer                        , pointer :: lunsrc
!
! Global variables
!
    integer                                                     , intent(in) :: idis
    integer                                                     , intent(in) :: kmax
    integer                                                     , intent(in) :: lstsci
    integer                                                     , intent(in) :: lsal
    integer                                                     , intent(in) :: ltem
    integer                                                     , intent(in) :: idensform
    real(fp)                                                    , intent(out):: taua
    real(fp)   , dimension(8)                                   , intent(out):: linkinf
    real(fp)                                                    , intent(in) :: saleqs
    real(fp)                                                    , intent(in) :: temeqs
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)               , intent(in) :: alfas
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub)               , intent(in) :: s0
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)         , intent(in) :: rho
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)         , intent(in) :: u0
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax)         , intent(in) :: v0
    real(fp)   , dimension(gdp%d%nmlb:gdp%d%nmub, kmax,lstsci)  , intent(in) :: r0
    real(fp)   , dimension(kmax)                                , intent(in) :: thick
    real(prec) , dimension(gdp%d%nmlb:gdp%d%nmub)               , intent(in) :: dps
    character*256                                                            :: filename
!
! Local variables
!
    integer                                :: ilen
    integer                                :: k
    integer                                :: kgrad
    integer                                :: nm_diff
    integer                                :: nmd_diff
    integer                                :: ndm_diff
    integer                                :: nm_amb
    integer                                :: nmd_amb
    integer                                :: ndm_amb
    integer                 , external     :: newlun
    integer                                :: luntmp
    real(fp)                               :: deg2rad
    real(fp)                               :: dengra
    real(fp)                               :: drohj
    real(fp)                               :: d_diff
    real(fp)                               :: ha
    real(fp)                               :: hd
    real(fp)                               :: hint
    real(fp)                               :: maxgrad
    real(fp)                               :: pi
    real(fp)                               :: rad2deg
    real(fp)                               :: rhoam
    real(fp)                               :: rhoas
    real(fp)                               :: rhoab
    real(fp)                               :: taurel
    real(fp)                               :: thck
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
    real(fp) , dimension(:) , allocatable  :: h1
    real(fp) , dimension(:) , allocatable  :: rhoa
    character*1                            :: tab
    character*1                            :: stype1
    character*1                            :: stype2
    character*3                            :: c_inode
    character*12                           :: ctime
    character*12                           :: ctime1
    character*12                           :: cha
    character*12                           :: chd
    character*12                           :: ch0
    character*12                           :: cua
    character*12                           :: crhoam
    character*12                           :: crhoas
    character*12                           :: crhoab
    character*12                           :: chint
    character*12                           :: cdrohj
    character*12                           :: cq0
    character*12                           :: crho0
    character*12                           :: cd0
    character*12                           :: ctaua
    logical                                :: linkinp
!
!
!! executable statements -------------------------------------------------------
!
    m_diff         => gdp%gdnfl%m_diff
    n_diff         => gdp%gdnfl%n_diff
    m_amb          => gdp%gdnfl%m_amb
    n_amb          => gdp%gdnfl%n_amb
    m_intake       => gdp%gdnfl%m_intake
    n_intake       => gdp%gdnfl%n_intake
    k_intake       => gdp%gdnfl%k_intake
    q_diff         => gdp%gdnfl%q_diff
    s0_diff        => gdp%gdnfl%s0_diff
    t0_diff        => gdp%gdnfl%t0_diff
    d0             => gdp%gdnfl%d0
    h0             => gdp%gdnfl%h0
    sigma0         => gdp%gdnfl%sigma0
    !
    write(c_inode(1:3),'(i3.3)') inode
    !
    pi      = acos(-1.0_fp)
    rad2deg = 180.0_fp / pi
    deg2rad = pi / 180.0_fp
    tab     = char(9)
    !
    allocate (h1   (kmax) )
    allocate (rhoa (kmax) )
    !
    ! Read the general diffusor characteritics from cormix input file
    !
    call n_and_m_to_nm(n_diff(idis)    , m_diff(idis)     , nm_diff  , gdp)
    call n_and_m_to_nm(n_diff(idis) - 1, m_diff(idis)     , ndm_diff , gdp)
    call n_and_m_to_nm(n_diff(idis)    , m_diff(idis) - 1 , nmd_diff , gdp)
    call n_and_m_to_nm(n_amb(idis)     , m_amb(idis)      , nm_amb   , gdp)
    call n_and_m_to_nm(n_amb(idis)  - 1, m_amb(idis)      , ndm_amb  , gdp)
    call n_and_m_to_nm(n_amb(idis)     , m_amb(idis)  - 1 , nmd_amb  , gdp)

    !
    ! Compute the depths
    !

    ha = s0(nm_amb)+real(dps(nm_amb),fp)
    hd = s0(nm_diff)+real(dps(nm_diff),fp)

    !
    ! Compute depth averaged velocity magnitude and direction
    !

    uuu = 0.0_fp
    vvv = 0.0_fp

    do k = 1, kmax
       uuu      = uuu + 0.5_fp * (u0(nm_amb ,k) + u0(nmd_amb ,k))*thick(k)
       vvv      = vvv + 0.5_fp * (v0(nm_amb ,k) + v0(ndm_amb ,k))*thick(k)
    enddo

    umag = sqrt (uuu*uuu + vvv*vvv)
    taua = atan2(vvv,uuu)*rad2deg + alfas(nm_amb)
    taua = mod(taua + 360.0_fp,360.0_fp)
    ua   = umag


    !
    ! Compute parameters density profile C,
    ! start with computing heigths and densities at the cell centers
    ! of the position were the ambient conditions are taken from
    ! (switch positive direction from positive downward to positive upward)
    !

    h1(1)   = 0.5_fp * thick(kmax) * (s0(nm_amb)+real(dps(nm_amb),fp))
    rhoa(1) = rho(nm_amb,kmax)
    !
    do k = 2, kmax
       thck     = 0.5_fp * (thick(kmax-k+2) + thick(kmax-k+1))
       h1 (k)   = h1(k-1) + thck*(s0(nm_amb) + real(dps(nm_amb),fp))
       rhoa (k) = rho(nm_amb, kmax-k+1)
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
          rhoam = rhoam + rho(nm_amb,k)*thick(k)
       enddo
    else
       stype1 = 'S'
       if (maxgrad < 0.5_fp*d_diff) then
          rhoas  = rho(nm_amb,1)
          rhoab  = rho(nm_amb,kmax)
          stype2 = 'A'
       else
          stype2 = 'C'
       endif
    endif

    !
    ! Compute the density of the discharged water
    !

    sal  = s0_diff(idis)
    temp = t0_diff(idis)
    if (lsal /= 0) then
       call coupled (add,r0,kmax,lstsci,lsal,thick,m_intake(idis),n_intake(idis),k_intake(idis),gdp)
       sal = sal + add
    else
       sal = saleqs
    endif

    if (ltem /= 0) then
       call coupled (add,r0,kmax,lstsci,ltem,thick,m_intake(idis),n_intake(idis),k_intake(idis),gdp)
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
    ! Write Cortime input file
    !

    luntmp = newlun(gdp)
    linkinp   = .true.

    do while (linkinp)
       inquire (file=trim(gdp%gdnfl%base_path)//'cortime_'//trim(gdp%runid)//'_'//c_inode//'.linkinp',exist=linkinp)
    enddo

    open (luntmp,file=trim(gdp%gdnfl%base_path)//'cortime_'//trim(gdp%runid)//'_'//c_inode//'.linkinp',status='new')
    write (luntmp,'(''CorTime v7.0'')')
    write (luntmp,'()')
    write (luntmp,'(''File name='',a1,a )') tab,'cortime_'//trim(gdp%runid)//'_'//c_inode//'.linkinp'
    write (luntmp,'(''Base case='',a1,a )') tab,trim(filename)
    write (luntmp,'(''Node='',a1,a )') tab,trim(gdp%runid)
    write (luntmp,'(''TOTSTEP='' ,a1,i1)') tab,1
    write (luntmp,'()')
    write (luntmp,'(27(a,a1),a)') 'TIME' , tab, 'HA'   , tab, 'HD'   , tab, 'UA'     , tab, &
   &                              'UorS' , tab, 'RHOAM', tab, 'STYPE', tab, 'RHOAS'  , tab, &
   &                              'RHOAB', tab, 'HINT' , tab, 'DROHJ', tab, 'Q0'     , tab, &
   &                              'C0'   , tab, 'RHO0' , tab, 'Gamma', tab, 'Sigma'  , tab, &
   &                              'D0'   , tab, 'B0'   , tab, 'H0'   , tab, 'PollTyp', tab, &
   &                              'L1Sub', tab, 'L1Den', tab, 'L2Sub', tab, 'L2Den'  , tab, &
   &                              'L3Sub', tab, 'L3Den', tab, 'Distb', tab, 'PHI'

    !
    ! Make character strings from all requested input
    !

    write (ctime (1:12),'(f12.3)') time/60.0_fp
    write (cha   (1:12),'(f12.3)') ha
    write (chd   (1:12),'(f12.3)') hd
    write (cua   (1:12),'(f12.3)') max(ua,0.001_fp)
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

    write (cq0   (1:12),'(f12.8)') q_diff(idis)
    write (crho0 (1:12),'(f12.3)') rho0
    write (ch0   (1:12),'(f12.3)') h0(idis)
    write (cd0   (1:12),'(f12.3)') d0(idis)

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

    ctime = adjustl(ctime )
    cha   = adjustl(cha   )
    chd   = adjustl(chd   )
    cua   = adjustl(cua   )
    crhoam= adjustl(crhoam)
    crhoas= adjustl(crhoas)
    crhoab= adjustl(crhoab)
    chint = adjustl(chint )
    cdrohj= adjustl(cdrohj)
    cq0   = adjustl(cq0   )
    crho0 = adjustl(crho0 )
    ch0   = adjustl(ch0   )
    cd0   = adjustl(cd0   )
    ctaua = adjustl(ctaua )

    write (luntmp,'(27(a,a1),a)') trim(ctime)  , tab, trim(cha)    , tab, trim(chd)    , tab, trim(cua)      , tab, &
   &                              trim(stype1) , tab, trim(crhoam) , tab, trim(stype2) , tab, trim(crhoas)   , tab, &
   &                              trim(crhoab) , tab, trim(chint)  , tab, trim(cdrohj) , tab, trim(cq0)      , tab, &
   &                              '1.0'        , tab, trim(crho0)  , tab, '-'          , tab, '-'            , tab, &
   &                              trim(cd0)    , tab, '-'          , tab, '-'          , tab, '-'            , tab, &
   &                              '-'          , tab, '-'          , tab, '-'          , tab, '-'            , tab, &
   &                              '-'          , tab, '-'          , tab, '-'          , tab, trim(ctaua)

    close (luntmp)
    !
    deallocate (h1)
    deallocate (rhoa)
    !
    ! Store some general information to write to the plume trajectory file
    !
    linkinf(1) = sal
    linkinf(2) = temp
    linkinf(3) = rho0
    linkinf(4) = ua
    linkinf(5) = taua
    linkinf(6) = taurel
    !
end subroutine wri_cortim

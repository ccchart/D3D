!!  Copyright (C)  Stichting Deltares, 2012-2019.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.

module part10_mod
!
!  module declarations
!
contains
      subroutine part10 ( lgrid  , volume , flow   , dx     , dy     ,      &
                          area   , angle  , nmax   , mnmaxk , idelt  ,      &
                          nopart , npart  , mpart  , xpart  , ypart  ,      &
                          zpart  , iptime , rough  , drand  , lgrid2 ,      &
                          lgrid3 ,                                          &
                          wvelo  , wdir   , decays , wpart  , pblay  ,      &
                          npwndw , vdiff  , nosubs , dfact  , modtyp ,      &
                          t0buoy , abuoy  , kpart  , mmax   , layt   ,      &
                          wsettl , depth  , ldiffz , ldiffh , lcorr  ,      &
                          acomp  , ipc    , accur  , xcor   , ycor   ,      &
                          tcktot , lun2   , alpha  , mapsub , nfract ,      &
                          taucs  , tauce  , chezy  , rhow   , lsettl ,      &
                          mstick , nstick , ioptdv , cdisp  , dminim ,      &
                          fstick , defang , floil  , xpart0 , ypart0 ,      &
                          xa0    , ya0    , xa     , ya     , npart0 ,      &
                          mpart0 , za     , locdep , dps    , nolay  ,      &
                          vrtdsp , stickdf, subst  , nbmax  , nconn  ,      &
                          conn   , tau    , caltau , nboomint , nboomtint,  &
                          iboomset , iboomtset, boomdepth,                  &
                          tyboom , efboom , tyboomt , efboomt, xpolboom ,   &
                          ypolboom , xipolboom, yipolboom, nrowsboom ,      &
                          itime  , v_swim , d_swim )
!
!    CALCULATES PARTICLE MOTION FROM ADVECTION, DISPERSION AND SETTLING
!                 (per time step)
!

!     system administration : frank kleissen


!     created               : january 1990, by l. postma


!     modified              : August 1990 by r.j.diependaal
!                             uses a vertical dispersion coefficient
!                             that is averaged over the vertical
!                             August 1990 by l.postma
!                             particle displacement is calculated in a
!                             more sophisticated manner (cross-overs to other
!                             segments induce a new calculation!)
!                           : February 1991 by a. hendriks and a. markus
!                             account for oil particles ( z coordinate
!                             larger than 1, and a two-layer system,
!                             where the two layers are recognised by
!                             positive and negative z coordinates and
!                             a separate parameter in drand)
!                           : March 1991 by a. markus (no oil anymore)
!                           : July 1992 by r.j. vos, for 8 substances
!                             with identical transport but different
!                             decay
!               delparv2.20:  April  1994 by r.j. vos, for buoyancy
!                             with more than 1 discharge and only in
!                             upper layer !!!
!               delparv3.00:  May    1994 by l.postma, smooth loading
!                             for continuous releases
!                             June   1994 by r.j. vos, uniform sampling
!                             in x and y direction for horizontal
!                             random walk
!               delparv3.10:  May    1996 by r.j. vos  3d version
!                             for the time being no pred-corr and uniform diff.
!                             including predictor-corrector of dunsbergen
!               delparv3.21:  Okt    1996 by r.j. vos  3d version
!               delparv3.22:  Nov    1996 by r.j. vos: correction pred-corr scheme 3d version
!               delparv3.40:  Nov    1997 by r.j. vos: 3d version with oil
!                             oil is transported as floating when any of the fractions
!                             of oil within the particle has a mass larger zero.
!               delparv3.43:  July   1998 by r.j. vos: settling/erosion at the bed
!                             a particle at the bed gets kp = 0
!               delparv3.50:  Sept   1998 by r.j. vos: mass may stick to land
!                                                      by horizontal diffusion or settling
!               delparv3.60:  April  1999 by r.j. vos: constant vertical diffusivity added
!                                                      sticking oil: definitions changed
!                                                      sticking parameter introduced
!               delparv3.60:  July     2000 by r.j. vos: when floating oil
!                                                         sticks, then prevent dispersed to
!                                                         overwrite the stick array
!                             July     2007 by Leo Postma: redesign of this routine, thin dams
!                                                          for advection and diffusion
!                             ????     ???? by Jan van Beek: unofficial DD support
!                             June     2011 by Leo Postma: support OMP parallelism
!                             October  2011 by Leo Postma: official support Domain Decomposition
!                             August   2015 by F.M. Kleissen: including moving oil booms

!     logical unit numbers  : standard output: error messages
!                             lun2           : error and debug messages


!     subroutines called    : stop_exit  - stops with return code (typical 1 for stop on error)
!                             p10cor  - function unknow
!                             part07  - not documented here
!                             part11  - not documented here


!     functions   called    : rnd     - random number generator
!                             vdisp   - not documented here
!
!
!**   modules and externals
!     note:
!       random function rnd() must be declared external, as it is an
!       intrinsic function for the lahey fortran95 compiler under linux
!
!  data definition module(s)
!
      use precision_part        ! single/double precision
      use timers           ! performance timer
      use typos
      use normal_mod
      use pinpok_mod
      use pmov_boom_mod  
      use p10cor_mod
      use grid_search_mod
      use spec_feat_par    ! special feature parameters
!
      implicit none             ! force explicit typing

      real   (sp), external      :: rnd                 ! random number function

!**   parameters used for dimensioning

      integer(ip), intent(in)    :: layt                ! number of layers of hydr. database
      integer(ip), intent(in)    :: mmax                ! second grid dimension
      integer(ip), intent(in)    :: mnmaxk              ! total number of active grid cells
      integer(ip), intent(in)    :: nmax                ! first grid dimension
      integer(ip), intent(in)    :: nopart              ! total number of particles
      integer(ip), intent(in)    :: nosubs              ! number of substances per particle

!**   other parameters

      integer(ip), intent(in)    :: idelt               ! time step size in seconds
      integer(ip), intent(in)    :: ioptdv              ! if 0 constant vertical diffusion &
                                                        ! if 1 depth averaged algebraic model
                                                        ! if 2 from file with constant added plus scaling plus minimum
      integer(ip), intent(in)    :: ipc                 ! if > 1 predictor corrector method used &
                                                        ! if   5 something special happens
      integer(ip), pointer    :: lgrid ( : , : )     ! grid with active grid numbers, negatives for open boundaries
      integer(ip), pointer    :: lgrid2( :, : )      ! total grid with grid numbers
      integer(ip), pointer       :: lgrid3( :, : )         ! total active grid with grid numbers
      integer(ip), intent(in)    :: lun2                ! unit number debug in formation file
      integer(ip), pointer    :: mapsub( : )         ! index for substances, used for oil
      integer(ip), intent(in)    :: modtyp              ! 1 = tracer model              &
                                                        ! 2 = 2-layer temperature model &
                                                        ! 3 = obsolete                  &
                                                        ! 4 = oil model                 &
                                                        ! 5 = 1-layer temperature model
      integer(ip), intent(in)    :: nfract              ! nr of oil fractions, each fraction 3 substances &
                                                        ! floating, dispersed and sticked
      integer(ip), intent(in)    :: npwndw              ! first active particle
      integer(ip), intent(in)    :: nstick              ! number of sticked particles !!!!! in !!!!!
      integer(ip), pointer :: iptime( : )         ! particle age in seconds
      integer(ip), pointer :: kpart ( : )         ! third grid index of the particles
      integer(ip), pointer :: mpart ( : )         ! second grid index of the particles
      integer(ip), pointer :: mpart0( : )         ! second grid index particles for previous time step
      integer(ip), pointer :: mstick( : )         ! which active substances can sticking (>0) and what is inactive partner?
                                                        ! j = mstick(i), j = inactive, i = active ; if j = 0 no sticking
                                                        ! if j is negative then i itself is sticking
      integer(ip), intent(inout) :: nolay               ! number of layers == layt
      integer(ip), pointer :: npart ( : )         ! first  grid index of the particles
      integer(ip), pointer :: npart0( : )         ! first  grid index particles for previous time step
      integer(ip), pointer :: floil ( : )         ! contains values 1 or 0
!
      logical    , intent(in)    :: acomp               ! use an analytical function for umagi
      logical    , intent(in)    :: lcorr               ! if true, apply the corrector step
      logical    , intent(in)    :: ldiffh              ! horizontal diffusion is on/off
      logical    , intent(in)    :: ldiffz              ! vertical diffusion is on/off
      logical    , intent(in)    :: lsettl              ! if on deposition/erosion may occur at the bed
!
      real   (sp), pointer    :: abuoy ( : )         ! dispersion coefficient buoyancy
      real   (sp), intent(in)    :: accur               ! accuracy limit taylor expansion
      real   (sp), intent(in)    :: alpha               ! scale factor vertical diffusion
      real   (sp), pointer    :: angle ( : )         ! angle with horizontal
      real   (sp), pointer    :: area  ( :  )        ! horizontal surface areas
      real   (sp), intent(in)    :: cdisp               ! constant vertical diffusivity in m2/s
      real   (sp), intent(in)    :: dminim              ! minimum value vertical diffusivity in m2/s
      real   (sp), pointer    :: decays( :  )        ! decay coefficients of substances in 1/day
      real   (sp), pointer    :: depth ( :  )        ! depth of segment  (in m)
      real   (sp), intent(in)    :: drand ( 3 )         ! drand(1) is 2.0*sqrt(dt*a) in D = a*t^b &
                                                        ! drand(2) is 0.5*b          in D = a*t^b &
                                                        ! drand(3) is winddrag coefficient in %
      real   (sp), pointer    :: flow  ( :    )      ! all flows
      real   (sp), pointer    :: fstick( : )         ! part of mass that sticks to land
      real   (sp), intent(in)    :: pblay               ! relative thickness of bottom layer ( 0.0 - 1.0 )
      real   (sp), intent(in)    :: rhow                ! density of water in g/l (= kg/m3)
      real   (sp), intent(in)    :: rough               ! bottom roughness length (m)
      real   (sp), pointer    :: t0buoy( : )         ! t0-coeff buoyancy
      real   (sp), intent(in)    :: tauce               ! critical shear stress for erosion (pa)
      real   (sp), intent(in)    :: taucs               ! critical shear stress for sedimentation (pa)
      real   (sp), pointer    :: volume( : )         ! volumes of all computational elements
      real   (dp), intent(in)    :: wdir (:)            ! wind direction from north
      real   (dp), intent(in)    :: wvelo(:)            ! wind velocity
      real   (sp), pointer    :: xcor  ( : )         ! bottom coordinate x
      real   (sp), pointer    :: ycor  ( : )         ! bottom coordinate y
      real   (sp), intent(inout) :: chezy               ! chezy coefficient (is set to 50 if .le. 1.0)
      real   (sp), intent(inout) :: defang              ! deflection angle for 3d oil simulations
                                                        ! enters in degrees, becomes radians
      real   (sp), pointer :: dfact ( : )         ! decay factor ( is exp(-decays*t) )
      real   (sp), pointer :: dps   ( : )         ! depths
      real   (sp), pointer :: dx    ( : )         ! dx of the computational elements
      real   (sp), pointer :: dy    ( : )         ! dy of the computational elements
      real   (sp), pointer :: locdep( :, : )      ! depth per layer
      real   (sp), pointer :: tcktot( : )         ! relative thickness of the layers
      real   (sp), pointer :: wpart ( :, :)       ! weight factors of the subs per particle
      real   (sp), pointer :: wsettl( : )         ! settling per particle
      real   (sp), pointer :: xa    ( : )         !
      real   (sp), pointer :: xa0   ( : )         !
      real   (sp), pointer :: xpart ( : )         ! x-value (0.0-1.0) first  direction within grid cell
      real   (sp), pointer :: xpart0( : )         ! x of particles for previous time step
      real   (sp), pointer :: ya    ( : )         !
      real   (sp), pointer :: ya0   ( : )         !
      real   (sp), pointer :: ypart ( : )         ! y-value (0.0-1.0) second direction within grid cell
      real   (sp), pointer :: ypart0( : )         ! y of particles for previous time step
      real   (sp), pointer :: za    ( : )         !
      real   (sp), pointer :: zpart ( : )         ! z-value (0.0-1.0) third  direction within grid cell
      real   (sp), pointer :: vdiff ( : )         ! vertical diffusion - work array
      real   (sp), pointer :: vrtdsp( :,: )       ! storage of vert disp info for debugging
      character(len=*), pointer  :: subst ( : )         ! substance names per substance &
                                                        ! and per layer for the plo file
      integer(ip), intent(in   ) :: stickdf             ! if 1 oil sticks at drying flats
      integer(ip), intent(in   ) :: nbmax         ! highest regular open boundary number
      integer(ip), intent(in   ) :: nconn         ! number of interdomain connections
      type( pnt ), intent(in   ) :: conn (nconn)  ! array with interdomain connections
      real   (sp), pointer       :: tau   ( : )   ! tau
      logical    , intent(in   ) :: caltau        ! if true, tau must be calculated

      integer(ip), intent(in   ) :: nboomint      ! number of boom introductions
      integer(ip), intent(in   ) :: nboomtint     ! number of moving boom introductions
      integer(ip), pointer :: iboomset(:)         ! timing of boom introduction
      integer(ip), pointer       :: iboomtset(:)  ! timing of moving boom introduction
      integer(ip), intent(in   ) :: tyboom        ! type of boom effectiveness parameter
      integer(ip), intent(in   ) :: tyboomt       ! type of moving boom effectiveness parameter
      real     ( sp), pointer  :: boomdepth (:)   ! depth of screen of boom
      real     ( sp), pointer  :: efboom (:,:)  ! effectiveness parameter of boom per oil type
      real     ( sp), pointer  :: efboomt (:,:)   ! effectiveness parameter of moving boom per oil type
      real     ( sp), pointer  :: xpolboom (:,:)  ! x-coordinates of boom polygon
      real     ( sp), pointer  :: ypolboom (:,:)  ! y-coordinates of boom polygon
      real     ( sp), pointer  :: xipolboom (:,:)  ! interpolated x-coordinates of moving boom polygon
      real     ( sp), pointer  :: yipolboom (:,:)  ! interpolated y-coordinates of moving boom polygon
      integer  ( ip), pointer  :: nrowsboom (:)  ! length of boom polygon

      integer(ip), intent(in   ) :: itime

      real   (sp), pointer :: v_swim( : )         ! horizontal swimming velocity m/s
      real   (sp), pointer :: d_swim( : )         ! horizontal swimming direction (degree)

!**   local parameters


      integer(ip)   :: icounz                  ! count the number of vertical bounces
      integer(ip)   :: ic                      ! segment number of active grid lgrid3
      integer(ip)   :: icvis                   ! strange construct to avoid endless operation, cyclic counter
      integer(ip)   :: icvist                  ! strange construct to avoid endless operation, total counter
      integer(ip)   :: idep                    ! layer offset in the hydrodynamic arrays
      integer(ip)   :: idepm1                  ! layer offset - one layer in the hydrodynamic arrays
      integer(ip)   :: idx                     ! boolean x direction of transport -1, 0, +1
      integer(ip)   :: idy                     ! boolean y direction of transport -1, 0, +1
      integer(ip)   :: idz                     ! boolean z direction of transport -1, 0, +1
      integer(ip)   :: ierror                  ! needed to call part07
      integer(ip)   :: ifract                  ! loop counter for nfract
      integer(ip)   :: inside                  ! to indicate whether a particle is within a polygon
      integer(ip)   :: ipart                   ! loop counter particle loop
      integer(ip)   :: ipcgo                   ! unclear helpvariable
      integer(ip)   :: isub                    ! loop counter nosubs
      integer(ip)   :: itdelt                  ! delta-t of the particle for smooth loading
      integer(ip)   :: ivisit( 4 )             ! part of the strange construct of icvis
      integer(ip)   :: jsub                    ! pointer containing mstick(isub)
      integer(ip)   :: kd                      ! loop counter of vertical layers
      integer(ip)   :: kp                      ! k of the particle
      integer(ip)   :: kpp                     ! local k of the particle
      integer(ip)   :: mp                      ! m of the particle
      integer(ip)   :: n0                      ! segment number 2d
      integer(ip)   :: n03d                    ! segment number 3d
      integer(ip)   :: n0new                   ! new value of segment number 2d
      integer(ip)   :: n1                      ! one back from n0 in first index
      integer(ip)   :: n2                      ! one back from n0 in second index
      integer(ip)   :: ninact                  ! summation counter inactive particles
      integer(ip)   :: np                      ! nr of this particle
      integer(ip)   :: npolbounce              ! poygon row number for bouncing particles
      integer(ip)   :: npolinside              ! poygon row number for particles moved due to moving boom
      integer(ip)   :: nstpar                  ! summation counter sticking particles
      integer(ip)   :: nopart_sed              ! no. of particles settled into bed layer (per time step)
      integer(ip)   :: nopart_ero              ! no. of particles eroded from bed layer (per time step)
      logical       :: coriol                  ! = abs( defang ) .ge. 1.0e-6
      logical       :: dstick                  ! logical that determines sticking
      logical       :: ldispo                  ! vertical diffusion is on (true) or off (false) for oil particle
      logical       :: lstick                  ! keeps the possibility to stick
      logical       :: oilmod                  ! = modtyp .eq. 4
      logical       :: threed                  ! = layt .gt. 1
      logical       :: twolay                  ! = modtyp .eq. 2
      logical       :: extract                 ! extract particle or not
      real(dp)      :: a                       ! tsja
      real(sp)      :: codef                   ! cosine of deflection angle
      real(sp)      :: dxx                     ! delta x with deflected oil
      real(sp)      :: dyy                     ! delta y with deflected oil
      real(sp)      :: nrms                    ! number of vertical diffusion values in dsprms
      real(sp)      :: sidef                   ! sine   of deflection angle
      real(sp)      :: xx                      ! new x with deflected oil
      real(sp)      :: yy                      ! new y with deflected oil
      real(sp)      :: abuac                   ! actual value of abuoy(ipart) ( * sqrt(ddfac)
      real(sp)      :: c1                      ! depth curve of the wind component
      real(sp)      :: c2                      ! depth curve of the velocities
      real(sp)      :: c2g                     ! = grav / chezy / chezy
      real(sp)      :: c3                      ! something undeterminded yet, fixed at 1.0
      real(sp)      :: cdrag                   ! wind drag in fraction = drand(3)
      real(sp)      :: cf                      ! help variable
      real(sp)      :: chi0                    ! variable in correction process
      real(sp)      :: chi1                    ! variable in correction process
      real(sp)      :: dax                     ! delta of diffusive spreading x
      real(sp)      :: day                     ! delta of diffusive spreading y
      real(sp)      :: ddfac                   ! control variable for smooth loading
      real(sp)      :: deltt                   ! real value of itdelt, remainder of time step
      real(sp)      :: deppar                  ! depth of particle from the surface
      real(sp)      :: depth1                  ! total depth of the old segment
      real(sp)      :: depth2                  ! total depth of the new segment
      real(sp)      :: depthl                  ! depth of the volume(segment)
      real(sp)      :: disp                    ! vertical diffusion
      real(sp)      :: dran1                   ! actual value of drand(1) or drand(1)*sqrt(ddfac)
      real(sp)      :: dred                    ! depth reduction factor for 2 layer modelling
      real(sp)      :: dspmax                  ! maximum vertical diffusion value
      real(sp)      :: dspmin                  ! minimum vertical diffusion value
      real(dp)      :: dsprms                  ! sum of squares of vertical diffusion value
      real(sp)      :: dvz                     ! vertical displacement due to vertical diffusion
      real(sp)      :: dvzs                    ! vertical displacement due to settling
      real(sp)      :: dvzt                    ! total vertical displacement
      real(sp)      :: dxp                     ! dx at location of the particle
      real(sp)      :: dyp                     ! dy at location of the particle
      real(sp)      :: f1                      ! helpvariable rough / depth
      real(sp)      :: h0wav                   ! wave height, as calculated
      real(sp)      :: rnorm                   ! normal distributed random number
      real(sp)      :: rnorm1                   ! normal distributed random number
      real(sp)      :: pstick                  ! used in fraction computation with oil
      real(sp)      :: ptlay                   ! relative fraction of top layer in twolayer system is 1.0 - pblay
      real(sp)      :: rtim                    ! timestep size that can be set minimum of x, y and z
      real(sp)      :: rtim1                   ! unclear help variable associated with rtim
      real(sp)      :: rtimx                   ! timestep size that can be set in x direction
      real(sp)      :: rtimy                   ! timestep size that can be set in y direction
      real(sp)      :: rtimz                   ! timestep size that can be set in z direction
      real(sp)      :: sangl                   ! angle of the computational element
      real(sp)      :: sq6                     ! = sqrt(6.0)
      real(sp)      :: t0                      ! helpvariable with initial time step buoyant spreading
      real(sp)      :: tdis1                   !helpvariable for lognormal distribution
      real(sp)      :: tdis2                   !helpvariable for lognormal distribution
      real(sp)      :: tp                      ! real value of iptime(ipart)
      real(sp)      :: trp                     ! horizontal random walk
      real(sp)      :: twopi                   ! 2 * pi
      real(sp)      :: ubstar                  ! bottom  friction velocity
      real(sp)      :: ubstar_b                ! bottom  friction velocity at sediment bed
      real(sp)      :: uecrit                  ! critical shearstress velocity for erosion
      real(sp)      :: umagi                   ! variable for strange correction process
      real(sp)      :: umagp                   ! variable for strange correction process
      real(sp)      :: uscrit                  ! critical shearstress velocity for sedimentation
      real(sp)      :: uwstar                  ! surface friction velocity
      real(sp)      :: vol                     ! value of volume (speed)
      real(sp)      :: vvx                     ! particle velocity in x dir.
      real(sp)      :: vvy                     ! particle velocity in y dir.
      real(sp)      :: vvz                     ! particle velocity in z dir.
      real(sp)      :: vx                      ! particle velocity in x dir. normalized
      real(sp)      :: vx0                     ! particle velocity in x dir. 0
      real(sp)      :: vx1                     ! particle velocity in x dir. 1
      real(sp)      :: vxnew                   ! particle velocity in x dir. corrected
      real(sp)      :: vxr                     ! particle velocity in x dir. real
      real(sp)      :: vxw                     ! particle velocity in x dir. winddriven
      real(sp)      :: vy                      ! particle velocity in y dir. normalized
      real(sp)      :: vy0                     ! particle velocity in y dir. 0
      real(sp)      :: vy1                     ! particle velocity in y dir. 1
      real(sp)      :: vynew                   ! particle velocity in y dir. corrected
      real(sp)      :: vyr                     ! particle velocity in y dir. real
      real(sp)      :: vyw                     ! particle velocity in y dir. winddriven
      real(sp)      :: vz                      ! particle velocity in z dir. normalized
      real(sp)      :: vz0                     ! particle velocity in z dir. 0
      real(sp)      :: vz1                     ! particle velocity in z dir. 1
      real(sp)      :: vznew                   ! particle velocity in z dir. corrected
      real(sp)      :: vzs                     ! vertical velocity due to settling
      real(sp)      :: wdirr                   ! is wind direction in radians
      real(sp)      :: wstick                  ! used in fraction computation with oil
      real(sp)      :: wsum                    ! used to sum the weights of a particle
      real(sp)      :: wrem (10)               ! used removed mass of a particle (maximum of 10 fractions)
      real(sp)      :: xnew                    ! new x value
      real(sp)      :: xp                      ! x of the particle
      real(sp)      :: xpold                   ! old x of the particle needed for the correction routine
      real(sp)      :: ynew                    ! new y value
      real(sp)      :: yp                      ! y of the particle
      real(sp)      :: ypold                   ! old y of the particle needed for the correction routine
      real(sp)      :: znew                    ! new z value
      real(sp)      :: zp                      ! z of the particle
      real(sp)      :: zp2                     ! help variable to correct z of the particle for twolay
      real(sp)      :: zpold                   ! old z of the particle needed for the correction routine
      integer(ip)   :: n03d2                   ! help variable sequence number of cell in next layer
      real   (rp)   :: disp2                   ! help variable diffusion in next layer
      real   (rp)   :: pbounce                 ! probability that a particle will NOT bounce
      integer(ip)   :: ddshift                 ! is it n (1) or m (2) direction
      real(sp)      :: xpnew                   ! potential new x of the particle
      real(sp)      :: ypnew                   ! potential new y of the particle
      real(sp)      :: xrest                   ! residual x offset
      real(sp)      :: yrest                   ! residual y offset
      real(sp)      :: dxpnew                  ! dxp new potential cell
      real(sp)      :: dypnew                  ! dyp new potential cell
      integer(ip)   :: mpnew                   ! potential new m of the particle
      integer(ip)   :: npnew                   ! potential new n of the particle
      logical       :: booms                   ! indicator for the precence of booms
      logical       :: boomseffective          ! indicator for effective precence of booms per particle
      integer(ip)   :: iboomint                ! counter for boom introductions
      integer(ip)   :: iboomtint               ! counter for moving boom introductions
      integer(ip)   :: iboomtot                ! counter for booms plus moving boom introductions
      real(rp)      :: pboomcatch(100)         ! chance for boom to catch particle
      real(sp)      :: xtmp                    !temporary x coordinate of moving boom
      real(sp)      :: xtmp1                   !temporary x coordinate of moving boom
      real(sp)      :: ytmp                    !temporary y coordinate of moving boom
      real(sp)      :: ytmp1                   !temporary y coordinate of moving boom

      real(rp)      :: xaold                   ! old real world x-coordinate
      real(rp)      :: yaold                   ! old real world y-coordinate
      real(rp)      :: xanew                   ! new real world x-coordinate
      real(rp)      :: yanew                   ! new real world y-coordinate
      logical       :: catch
      real(sp)      :: xacatch
      real(sp)      :: yacatch
      logical       :: bounce
      logical       :: bouncefirsttry
      real(sp)      :: xabounce
      real(sp)      :: yabounce
      real(sp)      :: sdir                    ! swimming direction
      real(sp)      :: displacement            ! swimming displacement
      real(sp)      :: dx_swim                 ! dx swimming displacement
      real(sp)      :: dy_swim                 ! dy swimming displacement
      real(sp)      :: xpolb(5,25)            ! polygon showing x new and old moving boom
      real(sp)      :: ypolb(5,25)            ! polygon showing y new and old moving boomF$OMP
      real(sp)      :: xcrossboom(2)          ! x coordinates of points on boom being crossed by a particle
      real(sp)      :: ycrossboom(2)          ! y coordinates of points on boom being crossed by a particle
      real(sp)      :: curangle
      real(sp)      :: boomangle
      real(sp)      :: diffangle
      real(sp)      :: eff_angle
      integer(ip)   :: irow

      integer(ip)   :: nboomtry
      integer(ip)   :: nscreenstry
      logical       :: screensfirsttry
      logical       :: leftside

!**   maintain all values

      save

!**   data assigments

      logical  :: debug  =  .false.
      real(sp) :: cd     =  1.3e-03
      real(sp) :: grav   =  9.81
      real(sp) :: rhoair =  1.25
      real(dp) :: rseed  =  0.5d+00
      real(sp) :: zsurf  =  0.001
      logical  :: first  = .true.
      logical  :: lbott  = .true.
      logical  :: lland  = .true.
      integer(4) ithndl              ! handle to time this subroutine
      data       ithndl / 0 /
      if ( timon ) call timstrt( "part10", ithndl )

!**   initialisations

      if ( first ) then
         first  = .false.
!     Set initial booms flag and booms introduction counter
         booms = .false.
         iboomint = 0
         iboomtint = 0
         sq6    = sqrt( 6.0 )
         twopi  = 8.0 * atan(1.0)
         defang = defang * twopi / 360.0    !  deflection angle oil modelling
         coriol = abs(defang) .ge. 1.0e-6   !  deflection from the equator aparently
         twolay = modtyp .eq. 2
         oilmod = modtyp .eq. 4
         threed = layt .gt. 1
         cdrag  = drand(3) / 100.0          !  wind drag as a fraction
         ptlay  = 1.0 - pblay
         dspmin = 1.0e+10
         dspmax = 0.0
         dsprms = 0.0
         nrms   = 0.0
         wrem=0.0

         if (.not.sedscreen) then
           do irow = 1,nboomint+iboomtint
             boomdepth(irow) = 100000.0
           enddo
         endif
           
         if ( lbott ) then
            write ( lun2, '(6x,a)' ) &
                'Bottom sticking due to diffusion is included '
         else
            write ( lun2, '(6x,a)' )  &
                'Bottom sticking due to diffusion is not included '
         endif
         if ( lland ) then
            write ( lun2, '(6x,a)' ) &
                'Sticking at land due to diffusion is included '
         else
            write ( lun2, '(6x,a)' ) &
                'Sticking at land due to diffusion is not included '
         endif
         if ( .not. ldiffz ) &
            write (lun2,'(6x,a)') 'No vertical diffusivity in part10 '
         if ( .not. ldiffh ) &
            write (lun2, '(6x,a)') 'No horizontal diffusivity in part10'
         if ( chezy .le. 1.0 ) then
            chezy = 50.0
            write(*,*) 'Chezy is missing: assumed is a chezy of: ',chezy
         endif
         if ( rhow .le. 1.0 ) then
            write( *, * ) 'Rho-w is missing ??'
            call stop_exit(1)
         endif
      endif

      !     Check for (update of) boom introductions
      if (iboomint .lt. nboomint) then
         if (iboomset(iboomint + 1) .eq. itime) then
            booms = .true.
            iboomint = iboomint + 1
           if(booms) write( lun2, '(a,i4)' ) ' Booms introduction # ', iboomint
           if(sedscreen) write( lun2, '(a,i4)' ) ' Sediment screen introduction # ', iboomint
            if (tyboom .eq. 1) then
               pboomcatch(1:nfract) = 1 - ((1- min(1.0,efboom(iboomint,1:nfract)))  ** (real(idelt) / 86400.0)) !should this not be just linear ?
!               pboomcatch(1:nfract) = 1 - ((1- min(1.0,efboom(iboomint,1:nfract)))  * (real(idelt) / 86400.0)) !should this not be just linear ?
            else
!              no other function implemented yet!
            endif
         end if
      end if

      !     Check for (update of) moving boom introductions at present only one moving boom is allowed
      ! for the interpolation of the positions we need to read two timers the first one when the time equals itime, but also the second.
      if (iboomtint .lt. nboomtint) then
        
         if (iboomtset(iboomtint + 1) .eq. itime) then
            booms = .true.
            iboomtint = iboomtint +1
            write( lun2, '(a,i4)' ) ' Moving booms introduction # ', iboomtint
            if (tyboomt .eq. 1) then
               pboomcatch(1:nfract) = 1 - ((1- Min(1.0,efboomt(iboomtint,1:nfract)))  ** (real(idelt) / 86400.0)) !linear
            else
!              no other function implemented yet!
            endif
         end if
! interpolate the boom coordinates and generate a polygon needed to "catch" particles. De polygoon moet dan bestaan uit de coordinaten van de
! oorpronkelijke cooridaat en in omgekeerde volgorde die van de nieuwe lijn.
        if(booms.or.sedscreen) then  !also for sediment screens, for sedscreens we need to account for the depth and not just surface
          ! store previous coordinates, for moving booms only
          do irow=2, nrowsboom(iboomtint)
            xpolb(1,irow-1)=xipolboom(irow-1, iboomtint)   ! from previous timestep:this is the polygon per set of two coordinates so in total (n-1) polygons
            ypolb(1,irow-1)=yipolboom(irow-1, iboomtint)
            xpolb(2,irow-1)=xipolboom(irow, iboomtint)   ! this is the polygon per set of two coordinates so in total (n-1) polygons
            ypolb(2,irow-1)=yipolboom(irow, iboomtint)
          enddo
          
          do irow = 1, nrowsboom(iboomtint)
          xtmp = xpolboom(irow, iboomtint)  !original coordinates (not the coordinates of the previous timestep)
          ytmp = ypolboom(irow, iboomtint)
          if (iboomtint.lt.nboomtint)then   !coordinates of the next time stamp in the file
            xtmp1 = xpolboom(irow, iboomtint+1)
            ytmp1 = ypolboom(irow, iboomtint+1)
          else   ! beyond the last time stamp, no further movement
            xtmp1 = xpolboom(irow, iboomtint)
            ytmp1 = ypolboom(irow, iboomtint)
          endif
          if (itime.eq.iboomtset(iboomtint))then  
            xipolboom(irow, iboomtint) = xtmp
            yipolboom(irow, iboomtint) = ytmp
!            write(lun2,*) ' moving boom (X,Y) coordinates:', itime, xipolboom(irow, iboomtint)
          if (irow.gt.1) then
            xpolb(1,irow-1)=xipolboom(irow-1, iboomtint)   ! from previous timestep:this is the polygon per set of two coordinates so in total (n-1) polygons
            ypolb(1,irow-1)=yipolboom(irow-1, iboomtint)
            xpolb(2,irow-1)=xipolboom(irow, iboomtint)   ! this is the polygon per set of two coordinates so in total (n-1) polygons
            ypolb(2,irow-1)=yipolboom(irow, iboomtint)
          endif
          else
            xipolboom(irow, iboomtint) = xtmp + 1.01*(xtmp1 - xtmp)/(iboomtset(iboomtint+1) - &
            iboomtset(iboomtint))*(itime - iboomtset(iboomtint))
            yipolboom(irow, iboomtint) = ytmp + 1.01*(ytmp1 - ytmp)/(iboomtset(iboomtint+1) - &
            iboomtset(iboomtint))*(itime - iboomtset(iboomtint))
!            write(lun2,*) ' interpolated moving boom (X,Y) coordinates:', itime,xipolboom(irow, iboomtint)
          endif
          enddo
          do irow = 1, nrowsboom(iboomtint)-1
            xpolb(4,irow) = xipolboom(irow, iboomtint)
            ypolb(4,irow) = yipolboom(irow, iboomtint)
            xpolb(3,irow) = xipolboom(irow+1, iboomtint)
            ypolb(3,irow) = yipolboom(irow+1, iboomtint)
            xpolb(5,irow) = xpolb(1,irow)
            ypolb(5,irow) = ypolb(1,irow)
          enddo
          
!          do irow = 1, nrowsboom(iboomtint)-1
!            write(lun2,'(10(a16,f20.2, 2x, f20.2))')' polygon boom: ',xpolb(1,irow), ypolb(1,irow)
!            write(lun2,'(10(a16,f20.2, 2x, f20.2))')' polygon boom: ',xpolb(2,irow), ypolb(2,irow)
!            write(lun2,'(10(a16,f20.2, 2x, f20.2))')' polygon boom: ',xpolb(3,irow), ypolb(3,irow)
!            write(lun2,'(10(a16,f20.2, 2x, f20.2))')' polygon boom: ',xpolb(4,irow), ypolb(4,irow)
!            write(lun2,'(10(a16,f20.2, 2x, f20.2))')' polygon boom: ',xpolb(5,irow), ypolb(5,irow)
!          enddo
          
        endif
      elseif (booms.and.(nboomtint.gt.0))then  !no more moving booms because we reached the last time of the boom specs.
        do irow = 1, nrowsboom(iboomtint)
          xipolboom(irow, iboomtint) = xpolboom(irow, iboomtint)
          yipolboom(irow, iboomtint) = ypolboom(irow, iboomtint)
          write(lun2,*) ' no longer moving: moving boom end (X,Y) coordinates:', itime, xipolboom(irow, iboomtint)
        enddo
      end if

!**   to avoid repeated calculations

      c2g     = grav / chezy / chezy
      uscrit  = sqrt( max(0.0,taucs) / rhow )
      uecrit  = sqrt( max(0.0,tauce) / rhow )
      do isub = 1, nosubs
         dfact(isub) = exp( -decays(isub) * idelt / 86400.0 )
         write(lun2,'(a,a,14x,a,es15.7,a,es15.7,a)')  &
            '      Decay ',subst(isub),' : ',dfact(isub), &
            ' [coeff = ',decays(isub),']'
      enddo

!**   store 'old' particles positions, these are also the location with which to check with moving boom

      if ( threed .and. oilmod .and. coriol ) then
         do ipart = 1, nopart
           xpart0(ipart) = xpart(ipart)
           ypart0(ipart) = ypart(ipart)
           npart0(ipart) = npart(ipart)
           mpart0(ipart) = mpart(ipart)
        enddo
      endif

!**   loop over the particles

      ninact = 0
      nstpar = 0

      nopart_sed=0
      nopart_ero=0

      timon = .false.
!$OMP PARALLEL DO PRIVATE ( np, mp, kp, kpp, n0, a, n03d, xp, yp, zp, tp,          &
!$OMP                       itdelt, ddfac, dran1, abuac, deltt, dred, kd, icvis,   &
!$OMP                       icvist, ivisit, lstick, wsum, isub, jsub, ifract,      &
!$OMP                       inside,pstick, wstick, ldispo, trp, t0, dax, day,      &
!$OMP                       n1, n2, dxp, dyp, depth1, idep, vol, vy0, vy1,         &
!$OMP                       vx0, vx1, vvx, vvy, vx, vy, vxr, vyr, ubstar, ubstar_b,&
!$OMP                       vz0, vz1, disp, dvz, depthl, dvzs, dvzt, vzs, icounz,  &
!$OMP                       znew, sangl, zp2, c1, f1, c2, c3, ipc, xnew,           &
!$OMP                       idepm1, vvz, vz, umagp, cf, rtimx, idx, rtimy, idy,    &
!$OMP                       rtimz, idz, rtim, rtim1, ipcgo, xpold, ypold, zpold,   &
!$OMP                       chi0, vxnew, vynew, vznew, chi1, ynew, vxw, vyw,       &
!$OMP                       deppar, yy, n0new, depth2, dstick, xx, n03d2, disp2,   &
!$OMP                       pbounce, xpnew, ypnew, xrest, yrest, dxpnew,           &
!$OMP                       dypnew, mpnew, npnew, ddshift, wdirr, uwstar,          &
!$OMP                       nboomtry, boomseffective, bouncefirsttry, xaold, yaold, xanew, yanew,  &
!$OMP                       xacatch, yacatch, catch, xabounce, yabounce, bounce ,  &
!$OMP                       sdir, displacement, dx_swim, dy_swim ,                 &
!$OMP                       xcrossboom, ycrossboom, ic, h0wav, curangle, boomangle,&
!$OMP                       diffangle, eff_angle, npolbounce,extract,              &
!$OMP                       irow, iboomtot, screensfirsttry, nscreenstry, leftside) , &
!$OMP           REDUCTION ( +   : ninact, nstpar, nopart_ero, nopart_sed,          &
!$OMP                             nrms  , dsprms, wrem ) ,                               &
!$OMP           REDUCTION ( MIN : dspmin ) , REDUCTION ( MAX : dspmax ),           &
!$OMP           SCHEDULE  ( DYNAMIC, max((nopart-npwndw)/100,1)           )

      do 100 ipart = npwndw, nopart

!**     retrieve particle characteristics

         np  = npart(ipart)
         mp  = mpart(ipart)
         kp  = kpart(ipart)
         kpp = kp
         if ( twolay ) kpp = 1             !   two layers are dealt with as one
                                           !   settling particles at the bed come in layer layt + 1
         kpp = min0 ( kp, layt )           !   kpp is a pointer that won't reach the bed..
         n0  = lgrid( np, mp   )
         if ( n0 .lt. 1 ) then             !   not an active grid-cell
            ninact       = ninact + 1
            npart(ipart) = 1
            mpart(ipart) = 1
            kpart(ipart) = 1
            do isub = 1, nosubs
               wpart( isub,ipart ) = 0.0
            enddo
            a = rnd(rseed)                 !   compatibility
            goto 90                        !   next particle
         endif
         n03d = n0 + (kpp-1)*nmax*mmax     !   3d
         xp =       xpart (ipart)
         yp =       ypart (ipart)
         zp =       zpart (ipart)
         tp = float(iptime(ipart))
         itdelt = idelt
         ddfac  = 2.0
         dran1  = drand(1)
         abuac  = 0.0
         if ( twolay ) abuac  = abuoy(ipart)
         if ( tp .lt. 0.0 ) then           !   adaptations because of smooth loading
            tp     = 0.0
            itdelt = idelt + iptime(ipart)
            ddfac  = float(itdelt)/float(idelt)
            dran1  = dran1 * sqrt(ddfac)
            abuac  = abuac * sqrt(ddfac)
         endif

         deltt = float(itdelt)
         dred  = 1.0                       !   depth reduction factor for 2-layer models
         if ( twolay ) then
            if ( kp .eq. 1 ) then
               dred = ptlay
            elseif ( kp .eq. 2 ) then
               dred = pblay
            else
               write (*,*) ' The layer-number is too large for modtyp=2'
               write( lun2,*) ' The layer-number is too large for modtyp=2'
               call stop_exit(1)
            endif
         endif

!**      strange construction to avoid infinite calculation

         icvis     = 1
         icvist    = 1
         ivisit(1) = n03d
         ivisit(2) = 0
         ivisit(3) = 0
         ivisit(4) = 0

!**     look whether the particle sticks

         lstick = .false.
         if ( nstick .gt. 0 ) then         !   nstick is an undocumented parameter
            wsum = 0.0
            do isub = 1 , nosubs
               jsub = mstick(isub)         !   a substance property
               if ( jsub .gt. 0) then      !
                  lstick = .true.          !   this particle might (but need not) stick
                  wsum   = wsum + wpart(jsub,ipart)   ! the weight of the sticked fraction
               endif                  !   the particle sticks if one of its substances sticks
            enddo

            if ( wsum .gt. 0.0 ) then      !   if this is true then
               nstpar = nstpar + 1         !   the particle sticks,
               goto 90                     !   go on with next particle
            endif

            if ( lstick ) then             !   aparently sticked fraction was zero for this particle
               pstick = 0.0                !   determine the net sticking probability
               wstick = 0.0                                       ! ==> to be checked
               do ifract = 1 , nfract      !   these are all oil fractions
                  jsub = 3*(ifract-1)
                  wstick = wstick + wpart(jsub+1,ipart)  &
                                  + wpart(jsub+2,ipart)
                  pstick = pstick + wpart(jsub+1,ipart) * fstick(ifract)  &
                                  + wpart(jsub+2,ipart) * fstick(ifract)
               enddo
               if ( wstick .gt. 1.0e-15 ) then
                  pstick = pstick/wstick
               else
                  pstick = 0.0
               endif
               if ( pstick .gt. rnd(rseed) ) then
                  lstick = .true.
               else
                  lstick = .false.
               endif
            endif
         endif

!**     look whether the (oil)particle floats (version 3.40)

         ldispo = .true.
         if ( oilmod .and. kp==1 .and. zpart(ipart)<=zsurf ) then
            wsum = 0.0                             !   location must be floating
            do isub = 1, nfract                    !   and at least one of the
               jsub = mapsub( (isub-1)*3 + 1 )     !   fractions has floating
               wsum = wsum + wpart( jsub, ipart )  !   mass.
            enddo                                  !   then there is no vertical
            if ( wsum .gt. 0.0 ) ldispo = .false.  !   displacement: ldispo = .false.
         endif

!**      horizontal diffusion ( and if needed buoyancy spreading )

         if ( ldiffh ) then
            trp = dran1 * tp ** drand(2)
            if (twolay ) then
            t0  = t0buoy(ipart)
              if ( t0 .gt. 0.0 .and. kp .eq. 1 ) then
               trp = max( trp , abuac * (tp+t0)**(-0.125) )     ! bouyancy spreading parameter
            endif
            endif
            dax = sq6 * trp * (rnd(rseed) - 0.5) ! This should be changd into a distance and angle rather than x and y direction
            day = sq6 * trp * (rnd(rseed) - 0.5) ! if we want to make dispersion dependent on direction of wind/current
            if (disp_dir) then
              dax = sq6 * disp_xdir * (rnd(rseed) - 0.5) ! only for hor/ver oriented grid (for the moment)
              day = sq6 * disp_ydir * (rnd(rseed) - 0.5) ! 
            endif
         else
            dax = 0.0
            day = 0.0
         endif

!**      bed erosion

         n1     = lgrid2(np - 1, mp    )
         n2     = lgrid2(np    , mp - 1)
         dxp    = dx    (n0)
         dyp    = dy    (n0)
         depth1 = depth (n0)
         if ( depth1 .le. 1.0e-25 ) then
            write(*,*) ' Depth1 too small '
            write(*,*) ' Depth1 = ',depth1,' gridcell n0= ',n0
            write(*,*) ' Error: abort from program  '
            write(lun2,*) ' Depth1 too small '
            write(lun2,*) ' Depth1 = ',depth1,' gridcell n0= ',n0
            write(lun2,*) ' Error: abort from program  '
            call stop_exit(1)
         endif

         ! first calculate or get ubstar_b at the bottom used for sedimentation and erosion

         if ( kpp .ne. layt ) then                  ! tau from file is only defined for bottom layer
            idep   = (layt - 1) * nmax * mmax
            if ( caltau ) then
               vol    = volume(n0 + idep  )
               vy0    = flow  (n1 + idep  ) / vol
               vy1    = flow  (n0 + idep  ) / vol
               vx0    = flow  (n2 + idep + mnmaxk) / vol
               vx1    = flow  (n0 + idep + mnmaxk) / vol
               vvx    = vx1 - vx0
               vvy    = vy1 - vy0
               vx     = vx0 + xp * vvx
               vy     = vy0 + yp * vvy
               vxr    = vx  * dxp
               vyr    = vy  * dyp
               ubstar_b = sqrt(c2g*(vxr*vxr + vyr*vyr))
            else
               ubstar_b = tau(n0 + idep) / rhow
            endif
         endif


         idep   = (kpp - 1) * nmax * mmax
         vol    = volume( n03d )
         vy0    = flow  (n1 + idep  ) / vol
         vy1    = flow  (n0 + idep  ) / vol
         vx0    = flow  (n2 + idep + mnmaxk) / vol
         vx1    = flow  (n0 + idep + mnmaxk) / vol
         vvx    = vx1 - vx0
         vvy    = vy1 - vy0
         vx     = vx0 + xp * vvx
         vy     = vy0 + yp * vvy
         vxr    = vx  * dxp                         ! vxr and vyr must be known
         vyr    = vy  * dyp                         ! from the beginning

         ! get or calculate ubstar for the actual layer

         if ( kpp .eq. layt ) then
            if ( caltau ) then
            ubstar = sqrt(c2g*(vxr*vxr + vyr*vyr))  ! ubstar this is requiered for dispersion
            else
               ubstar = tau(n03d) / rhow
            endif
            ubstar_b   = ubstar                     ! ubstar_bot is required for sedimentation and erosion
         else
            ubstar = sqrt(c2g*(vxr*vxr + vyr*vyr))  ! ubstar this is requiered for dispersion
         endif

         if ( lsettl .and. kp .eq. layt+1 ) then
            if ( ubstar_b .ge. uecrit ) then
               kp = layt
               zp = 0.05
               nopart_ero = nopart_ero + 1
            else
               goto 50
            endif
         endif
!**
!**      vertical diffusion                         ! and for vertical diffusion
!**
         n03d  = n0 + idep
         disp        = 0.0
!         vdiff(n03d) = 0.0
         dvz         = 0.0
         depthl      = volume(n03d)/area(n0)
         if ( ldiffz ) then
            if ( area(n0) .le.  1.0e-25 ) then
               write(*,*) ' Area too small '
               write(*,*) ' n0 = ',n0,' area(n0) = ',area(n0)
               write(*,*) ' Error: abort from program  '
               write(lun2,*) ' Area too small '
               write(lun2,*) ' n0 = ',n0,' area(n0) = ',area(n0)
               write(lun2,*) ' Error: abort from program  '
               call stop_exit(1)
            endif
            depthl = volume(n03d)/area(n0)
            if ( depthl .le.  1.0e-25 ) then
               write(*,*) ' Depthl too small '
               write(*,*) ' Depthl = ',depthl,' n03d = ',n03d
               write(*,*) ' Error: abort from program  '
               write(lun2,*) ' Depthl too small '
               write(lun2,*) ' Depthl = ',depthl,' n03d = ',n03d
               write(lun2,*) ' Error: abort from program  '
               call stop_exit(1)
            endif
            select case ( ioptdv )
               case ( 0 )            !  from version 3.60 also constant vertical diseprsion
                  disp = cdisp*alpha
               case ( 1 )            !  from version 3.43 also required for settling:
                  uwstar = sqrt( cd * rhoair / rhow * wvelo(n0)*wvelo(n0) )
                  disp = 0.83  * 0.41 * depth1 * sqrt(3.0) *     &       !  0.41 is rkappa
                               ( ubstar/6.0 + twopi*uwstar/32.0 )
                  disp = disp*alpha
               case ( 2 )
                  disp = max( cdisp + alpha*vdiff(n03d) , dminim )
               case default
                  write(*,*) ' Ioptdv = ',ioptdv
                  write(*,*) ' this option for vert. dispersion '
                  write(*,*) ' is not implemented               '
                  call stop_exit(1)
            end select
!            vdiff(n03d) = disp
            dvz         = 2.0 * sq6 * sqrt( disp*itdelt ) *   &
                             ( rnd(rseed)-0.5 ) / depthl / dred
         endif
        if (modtyp.ne.7) then
          dvzs = wsettl(ipart)*itdelt/depthl          !  settling      ?? what is the effect of this?? jvb: no bouncing for settling particles
          dvzt = dvzs + dvz                           !  vertical diffusion
        else
          dvzt = dvz                                  !  vertical diffusion
        endif
        
!**      oil: if floating, no vertical dispersion and no settling...

         if ( oilmod .and. .not. ldispo ) then
            disp      = 0.0
!            vdiff(n0) = 0.0
            dvz       = 0.0
            dvzs      = 0.0
            dvzt      = 0.0
         endif

!**      vertically bouncing particles ?

         icounz = 0
         znew   = zp + dvzt
         do while ( znew .gt. 1.0 .or. znew .lt. 0.0 )
!           icounz = icounz + 1
!           if ( icounz .gt. 90 ) then
!              write(*,*) ' Error: particle now',icounz,' times back '
!              write(*,*) ' Vertical displacements too large ?? '
!              write(*,*) ' Particle no.= ',ipart
!              write(lun2,*) ' Error: particle now',icounz,' times back'
!              write(lun2,*) '    Particle number  = ',ipart
!              write(lun2,*) '    Last     layer   = ',kp
!              write(lun2,*) '    Last     z-value = ',dvzt
!              write(lun2,*) '    Last     dvz     = ',dvz
!              write(lun2,*) '    Last     dvzs    = ',dvzs
!              write(lun2,*) '    Last     z-value = ',dvzt
!              write(lun2,*) '    Last     depth   = ',depthl
!              if(icounz > 100) call stop_exit(1)
!           endif

!**      boundary conditions, check here also settling and erosion
!**      of particles with critical velocities at the bed

            if ( znew .gt. 1.0 ) then
               if ( .not. twolay .and. kp .ne. layt ) then
                  n03d2 = n03d + nmax*mmax
                  if ( ioptdv .eq. 2 ) then
                     disp2 = max( cdisp + alpha*vdiff(n03d2) , dminim )
                     if (disp.eq.0.0) then
                       pbounce = 0.0
                     else
                       pbounce = sqrt( disp2/disp )
                     endif
                     
                     if ( disp2 .lt. disp ) then
                        if ( rnd(rseed) .lt. 1.0 - pbounce ) then
                           znew = 2.0 - znew
                           cycle
                        endif
                     endif
                     znew = ( znew - 1.0 ) * depthl * pbounce
                     disp = disp2
                  else
                     znew = ( znew - 1.0 ) * depthl
                  endif
                  kp     = kp + 1
                  n03d   = n03d2
                  depthl = volume(n03d)/area(n0)
                  znew   = znew / depthl
               else                                   ! either twolay or kp .eq. layt
                  if ( lsettl .and. ubstar_b .lt. uscrit   &
                              .and. wsettl(ipart) .gt. 0.0 ) then
                     kpart(ipart) = layt + 1          !  at high vert disp everything settles !!!
                     zpart(ipart) = 0.0
                     nopart_sed = nopart_sed + 1
                     if ( lstick .and. lbott ) then   !  dispersed oil (or other subs) may stick to the bottom
                        lstick = .false.              !  the particle sticks only once in this algorithm
                        do isub = 1, nosubs
                           jsub = mstick(isub)
                           if ( jsub .gt. 0 ) then    !  phase change from floating or dispersed
                              wpart(jsub,ipart) = wpart(isub,ipart)  ! shouldn't this be summed ?
                              wpart(isub,ipart) = 0.0
                           endif
                        enddo                         ! ==> also here
                     endif
                     goto 50
                  else
                     if (vertical_bounce) then
                         znew = 2.0 - znew             !  now it bounces
                     else
                         znew = 0.9990                 !  now it stays near the bottom (no bounce)
                     endif
                  endif
               endif
            elseif ( znew .lt. 0.0 ) then
               if ( .not. twolay .and. kp .ne. 1 ) then
                  n03d2 = n03d - nmax*mmax
                  if ( ioptdv .eq. 2 ) then
                     disp2 = max( cdisp + alpha*vdiff(n03d2) , dminim )
                     if (disp.eq.0.0) then
                       pbounce = 0.0
                     else
                     pbounce = sqrt( disp2/disp )
                     endif
                     if ( disp2 .lt. disp ) then
                        if ( rnd(rseed) .lt. 1.0 - pbounce ) then
                           znew = - znew
                           cycle
                        endif
                     endif
                     znew = znew * depthl * pbounce
                     disp = disp2
                  else
                     znew = znew * depthl
                  endif
                  kp = kp - 1
                  n03d   = n03d2
                  depthl = volume(n03d)/area(n0)
                  znew   = znew / depthl + 1.0
               else
                  if (vertical_bounce) then
                     znew   = -znew                      !  now it bounces
                  else
                     znew   = 0.0001                     !  now it stays near the surface (no bounce)
                  endif
               endif
            endif
         enddo
         if ( znew .eq. 0.0 ) znew = 0.0001
         if ( znew .eq. 1.0 ) znew = 0.9999
         zp = znew                                    !  the final vertical position

!*       for checking results

         dspmin = min(disp,dspmin)
         dspmax = max(disp,dspmax)
         dsprms = dsprms + disp*disp
         nrms   = nrms   + 1.0

!         debugging code
         if ( debug ) vrtdsp(1,ipart) = disp
         if ( debug ) vrtdsp(2,ipart) = dvz  * depthl
         if ( debug ) vrtdsp(3,ipart) = 0.0                    ! jvb moved to advection  = dvzs * depthl
         if ( debug ) vrtdsp(4,ipart) = dvzt * depthl
         if ( debug ) vrtdsp(5,ipart) = depthl
         if ( debug ) vrtdsp(6,ipart) = depth1
         if ( debug ) vrtdsp(7,ipart) = n0
!         vrtdsp(7,ipart) = n0
!**
!**       this is innerloop for particles crossing gridcell borders
!**
   10    n1     = lgrid2(np - 1, mp    )              ! pas op n0 is niet opnieuw gezet hier
         n2     = lgrid2(np    , mp - 1)
         dxp    = dx    (n0)
         dyp    = dy    (n0)
         sangl  = angle (n0)
         depth1 = depth (n0)
         if ( depth1 .le. 1.0e-25 ) then
            write(*,*) ' Depth1 too small '
            write(*,*) ' Depth1 = ',depth1,' gridcell n0= ',n0
            write(*,*) ' Error: abort from program  '
            write(lun2,*) ' Depth1 too small '
            write(lun2,*) ' Depth1 = ',depth1,' gridcell n0= ',n0
            write(lun2,*) ' Error: abort from program  '
            call stop_exit(1)
         endif

! ================a d v e c t i o n=========================================

!** windprofiles and logarithmic profiles
! if winddrag is stochastic sample dragcoefficient from a distribution (for oil only)
! this comes from http://blogs.sas.com/content/iml/2014/06/04/simulate-lognormal-data-with-specified-mean-and-variance.html

         if ( oilmod .and. .not. ldispo .and. dragvar .and. dragsigma.gt.0.0) then   ! floating oil using ! 
!           rnorm = sqrt(-2.*log(rnd(rseed)))*cos(twopi*rnd(rseed)) !normal distribution)/100.
!           cdrag = (dragmean+dragsigma*rnorm)/100.0 !lognormal distribution
!           do while ((cdrag.lt.(dragmean-3*dragsigma)).or.(cdrag.gt.(dragmean+3*dragsigma))) !this is to cut off 3-sigma tails of the distribution
             rnorm = sqrt(-2.*log(rnd(rseed)))*cos(twopi*rnd(rseed)) !normal distribution)/100.
             cdrag = max(0.0,(dragmean+dragsigma*rnorm))/100.0 !normal distribution, but cannot be less than 0.
 !          enddo
        endif
!          if(iptime(ipart).eq.0)then
!            write(lun2,*)rnorm,cdrag
!          else
!            stop
!          endif

         if ( layt .eq. 1 ) then
!              for 2 layers: the pos. in the top layer + thickness bottom layer
            zp2 = 1.0 - zp
            if ( twolay .and. kp .eq. 1 ) zp2 = zp2 - pblay
!              depth curve of the wind component
            c1 = ( 3.0*zp2 - 2.0 ) * zp2 * wvelo(n0) * cdrag * itdelt
!              depth curve of the velocities
            f1 = rough / depth1
            c2 = alog(zp2/f1+1.0) / ( (1.0+f1)*alog(1.0/f1+1.0) - 1.0 )
            c3 = 1.0
         else
            c1 = 0.0  ! no wind profile in 3d
            c2 = 1.0
            c3 = 1.0
         endif

!** determine the velocities in the 3 directions

         kpp  = kp
         if ( twolay ) kpp = 1
         if ( kpp .lt. 1 .or. kpp .gt. layt ) then
            write(*,*) ' program error part10: kpp out of range '
            call stop_exit(1)
         endif
         idep = (kpp - 1) * nmax * mmax
         vol  = volume( n0 + idep          )
         vy0  = flow  ( n1 + idep          ) / vol
         vy1  = flow  ( n0 + idep          ) / vol
         vx0  = flow  ( n2 + idep + mnmaxk ) / vol
         vx1  = flow  ( n0 + idep + mnmaxk ) / vol
         vz0  = 0.0
         vz1  = 0.0

         if ( threed ) then
            if     ( kpp .eq.   1  ) then  ! first layer
               vz1    = flow ( n0 + idep   + 2*mnmaxk) / vol
            elseif ( kpp .eq. layt ) then  ! last  layer
               idepm1 = idep - nmax * mmax
               vz0    = flow ( n0 + idepm1 + 2*mnmaxk) / vol
            else
               idepm1 = idep - nmax * mmax
               vz0    = flow ( n0 + idepm1 + 2*mnmaxk) / vol
               vz1    = flow ( n0 + idep   + 2*mnmaxk) / vol
            endif
         endif

         if (modtyp.eq.7) then
            depthl = volume(n0+idep)/area(n0)
            vzs  = wsettl(ipart)/depthl             !  settling
            if ( kpp .ne.    1 ) vz0  = vz0 + vzs
            if ( kpp .ne. layt ) vz1  = vz1 + vzs
         endif
        
         if ( oilmod .and. .not. ldispo ) then   ! floating oil
            vz0 = 0.0
            vz1 = 0.0
            if ( stickdf .eq. 1 ) then
               if ( vy0 .eq. 0.0 .and. vy1 .eq. 0.0 .and.           &
     &              vx0 .eq. 0.0 .and. vx1 .eq. 0.0        ) then ! the grid cell ran dry
                  do isub = 1, nosubs
                     jsub = mstick(isub)
                     if ( jsub .gt. 0 ) then
                        wpart(jsub,ipart) = wpart(jsub,ipart) + wpart(isub,ipart)
                        wpart(isub,ipart) = 0.0
                     endif
                  enddo
                  goto 50
               endif
            endif
         endif

         vvx = vx1 - vx0
         vvy = vy1 - vy0
         vvz = vz1 - vz0
         vx  = vx0 + xp * vvx
         vy  = vy0 + yp * vvy
         vz  = vz0 + zp * vvz
         vxr = vx  * dxp         !     only used for dispersion
         vyr = vy  * dyp         !     only used for dispersion

!**      velocity correction directly in first step..

         if ( lcorr ) then
            if ( ipc .eq. 5 ) then
               umagp = sqrt(vx*vx + vy*vy + vz*vz)   ! normal magnitude
               if ( abs(umagp) .le. 1.0e-15) then
                  cf = 1.0
               else           ! higher order ( and thus not conserving ! ) velocity
                  cf = umagi ( xp    , yp    , vz    , np    , mp    ,   &
                               kpp   , nmax  , mmax  , layt  , flow  ,   &
                               depth , lgrid , vol   , xcor  , ycor  ,   &
                               lgrid2, mnmaxk, acomp , tcktot)
                  cf = cf/umagp
               endif
            else
               cf = 1.0
            endif
            c2 = c2*cf
         endif

!**      time of flight particle in x direction

         rtimx = float(itdelt)
         idx   = 0
         if ( abs(c2*vvx) .gt. accur ) then
            if ( vx .gt. 0.0 ) then
               if ( vx1 .gt.  1.0e-25 ) then
                  rtimx = alog(vx1/vx) / (c2*vvx)
                  idx   =  1
               endif
            endif
            if ( vx .lt. 0.0 ) then
               if ( vx0 .lt. -1.0e-25 ) then
                  rtimx = alog(vx0/vx) / (c2*vvx)
                  idx   = -1
               endif
            endif
         else
            if     ( vx .gt.  1.0e-25 ) then
               rtimx = (1.0 - xp) / vx / c2
               idx   =  1
            elseif ( vx .lt. -1.0e-25 ) then
               rtimx = -xp / vx / c2
               idx   = -1
            else
               rtimx = 9999*rtimx
            endif
         endif

!**      time of flight particle in y direction

         rtimy = float(itdelt)
         idy   = 0
         if ( abs(c2*vvy) .gt. accur ) then
            if ( vy .gt. 0.0 ) then
               if ( vy1 .gt.  1.0e-25 ) then
                  rtimy = alog(vy1/vy) / (c2*vvy)
                  idy   =  1
               endif
            endif
            if ( vy .lt. 0.0 ) then
               if ( vy0 .lt. -1.0e-25 ) then
                  rtimy = alog(vy0/vy) / (c2*vvy)
                  idy   = -1
               endif
            endif
         else
            if     ( vy .gt.  1.0e-25 ) then
               rtimy = (1.0 - yp) / vy / c2
               idy   =  1
            elseif ( vy .lt. -1.0e-25 ) then
               rtimy = -yp / vy / c2
               idy   = -1
            else
               rtimy = 9999*rtimy
            endif
         endif

!**       time of flight particle in z direction
! 3d..  positive velocity is from z = 0 to z = 1
! 3d..                    is assumed to go downwards
! 3d..                    from small kp to large kp

         rtimz = float(itdelt)
         idz   = 0
         if ( abs(c3*vvz) .gt. accur ) then
            if ( vz .gt. 0.0 ) then
               if ( vz1 .gt.  1.0e-25 ) then
                  rtimz = alog(vz1/vz) / (c3*vvz)
                  idz   =  1
               endif
            endif
            if ( vz .lt. 0.0 ) then
               if ( vz0 .lt. -1.0e-25 ) then
                  rtimz = alog(vz0/vz) / (c3*vvz)
                  idz   = -1
               endif
            endif
         else
            if     ( vz .gt.  1.0e-25 ) then
               rtimz = (1.0 - zp) / vz / c3
               idz   =  1
            elseif ( vz .lt. -1.0e-25 ) then
               rtimz = -zp / vz / c3
               idz   = -1
            else
               rtimz = 9999*rtimz
            endif
         endif

!**      (partial) advection step

         rtim  = min( rtimx, rtimy, rtimz, deltt )
         rtim1 = rtim
         ipcgo = 1

         if ( abs(vvx) .gt. accur ) then
            xp = ( vx/vvx ) * exp( vvx*c2*rtim ) - vx0/vvx
         else
            xp = xp + vx*c2*rtim
         endif

         if ( abs(vvy) .gt. accur ) then
            yp = ( vy/vvy ) * exp( vvy*c2*rtim ) - vy0/vvy
         else
            yp = yp + vy*c2*rtim
         endif

         if ( abs(vvz) .gt. accur ) then
            zp = ( vz/vvz ) * exp( vvz*c3*rtim ) - vz0/vvz
         else
            zp = zp + vz*c3*rtim
         endif

!**
         if ( lcorr ) then
            ipcgo = 1
            xpold = xp
            ypold = yp
            zpold = zp
            umagp = sqrt(vx*vx + vy*vy + vz*vz)   ! normal magnitude
            if ( abs(umagp) .le. 1.0e-15) then
               chi0 = 1.0
            else           ! higher order ( and thus not conserving ! ) velocity
               chi0 = umagi ( xpold , ypold , vz    , np    , mp    ,  &
                              kpp   , nmax  , mmax  , layt  , flow  ,  &
                              depth , lgrid , vol   , xcor  , ycor  ,  &
                              lgrid2, mnmaxk, acomp , tcktot)
               chi0 = chi0/umagp
            endif
            vxnew = vx0 + xp * vvx
            vynew = vy0 + yp * vvy
            vznew = vz0 + zp * vvz
            umagp = sqrt(vxnew*vxnew + vynew*vynew + vznew*vznew)
            if ( abs(umagp) .le. 1.0e-15 ) then
               chi1 = 1.0
            else
               chi1 = umagi ( xp    , yp    , vznew , np    , mp    ,   &
                              kpp   , nmax  , mmax  , layt  , flow  ,   &
                              depth , lgrid , vol   , xcor  , ycor  ,   &
                              lgrid2, mnmaxk, acomp , tcktot)
               chi1 = chi1 / umagp
            endif

            call p10cor ( rtim1 , rtim  , xp    , yp    , zp    ,   &
                          xpold , ypold , zpold , chi0  , chi1  ,   &
                          deltt , vol   , vx    , vy    , vz    ,   &
                          vvx   , vvy   , vvz   , c2    , c3    ,   &
                          ipc   , vx0   , vy0   , vz0   , ipcgo ,   &
                          accur )

         endif

!**      did the time step end ?

         if ( rtim1 .lt. deltt ) then
            if ( rtim  .lt. deltt .and.   &    ! see p10cor for the meaning of these
                 ipcgo .eq.  1         ) then  ! strange things  !lp!

!             logically one of the 3 conditions below should hold,
!             so either rtim=rtimx or rtim=rtimy or rtim=rtimz
!             this follows directly from
!                   rtim = min(rtimx, rtimy, rtimz, deltt)
!             it is however not checked whether the particle is at the border !lp!
!             in order to eliminate round off a relative test is
!             added for this check

               if     ( eql(rtim,rtimx) ) then
                  mp = mp + idx
                  xp = xp - float(idx)
                  ddshift = 1
               elseif ( eql(rtim,rtimy) ) then
                  np = np + idy
                  yp = yp - float(idy)
                  ddshift = 2
               elseif ( eql(rtim,rtimz) ) then
                  kp = kp + idz
                  zp = zp - float(idz)
                  ddshift = 0
                  if (   kp .le. 0 .or.   &
                       ( kp .gt. layt .and. .not. twolay ) ) then
                     write(*,*) ' Particle = ',ipart,  &
                                ' not on an active layer'
                     write(*,*) ' Programming error in rtim in part10'
                     call stop_exit(1)
                  endif
               else
                  write(*,*) ' ipart, rtim, rtimx, rtimy, rtimz '
                  write(*,*) ipart,rtim,rtimx, rtimy, rtimz
                  write(*,*) ' rtim1, rtim, deltt '
                  write(*,*) rtim1,rtim,deltt
                  write(*,*) ' Programming error in rtim in part10'
                  call stop_exit(1)
               endif
            endif

            n0 = lgrid( np, mp)
            if ( n0 .lt. -nbmax ) then
               n0 = - n0 - nbmax
               call p10ddb ( nconn , conn  , n0     , ddshift, np     ,         &
     &                       mp    , xp    , yp     )
               n0 = lgrid( np, mp)
            endif
            if ( n0 .lt. 1 ) then        ! inactive grid cell
               npart(ipart) = 1
               mpart(ipart) = 1
               kpart(ipart) = 1
               do isub = 1, nosubs
                  wpart( isub, ipart ) = 0.0
               enddo
               goto 90                   ! next particle
            endif

            deltt = deltt - rtim1        ! reduce remaining part of time step
            xp = amax1( amin1(xp,1.0), 0.0 ) ! this is not strong
            yp = amax1( amin1(yp,1.0), 0.0 )
            zp = amax1( amin1(zp,1.0), 0.0 )

!**         strange construction to avoid infinite computation

            icvis  = icvis  + 1
            icvist = icvist + 1
            if ( icvis .eq. 5 ) icvis = 1
                                                      !  this is all not so clear
            kpp = kp
            if ( twolay ) kpp = 1
            n03d = n0 + (kpp - 1)*nmax*mmax
            if ( ivisit(icvis) .ne. n03d .and. icvist .lt. 10000 ) then
               ivisit(icvis) = n03d
               goto 10
            endif
         endif

!**      xnew,ynew is inclusive of diffusion (dax,day); c1 = 0.0 for 3d

         if (modtyp.eq.7) then 
            sdir         = d_swim(ipart) * twopi / 360.0
            displacement = v_swim(ipart)*idelt
            dx_swim      = displacement*sin(sdir+sangl)
            dy_swim      = displacement*cos(sdir+sangl)
         else
            dx_swim = 0.0
            dy_swim = 0.0
         endif
         
         wdirr = wdir(n0) * twopi / 360.0
         xnew  = xp + (dax - c1 * sin(wdirr + sangl) + dx_swim) / dxp
         ynew  = yp + (day - c1 * cos(wdirr + sangl) + dy_swim) / dyp

!**      floating oil

         if ( threed .and. oilmod ) then
            if ( ldispo ) then                ! vertical diffusion of oil
               floil( ipart ) = 0
            else
               floil( ipart ) = 1             ! oil is floating
               vxw  = - wvelo(n0) * sin( wdirr + sangl )
               vyw  = - wvelo(n0) * cos( wdirr + sangl )
!             drag on the difference vector: cd * (wind - flow)
               xnew = xnew  + (cdrag*(vxw-vxr)/dxp) * itdelt    !
               ynew = ynew  + (cdrag*(vyw-vyr)/dyp) * itdelt    !
            endif
         endif
         znew = zp

!**      make the depth of the particle from the water surface
!**      this is part of code that aims at bouncing if a particle meets the wall

         deppar = 0.0
         n03d = n0
         do kd = 1, kpp - 1                     ! volume of layers above
            deppar = deppar + volume( n03d )    ! particle
            n03d   = n03d + nmax*mmax
         enddo
         deppar = deppar + volume(n03d)*znew    ! + volume above particle
         deppar = deppar / area(n0)             ! gives depth of particle

         nboomtry = 1
         bouncefirsttry = .true.
         pboomcatch = 1.1
         if((ldispo .and. .not.sedscreen) .or. (.not.booms)) then
            boomseffective = .false.
         else

           if (boom_proc) then 
            ! hier moet efficiency van boom worden aangepast aan de omstandigheden als we afhankelijkheid van wind/golven/stroomsnelheid implementeren. can this be done later?
            ic = lgrid3(npart(ipart), mpart(ipart))
            h0wav = 0.243*wvelo(ic)*wvelo(ic)/grav
!            efboomt(iboomtint,1:nfract) = 1.0-min((max(h0wav-0.5,0))/5.0,1.0) ! just an arbitrary relationship with wave height
            if (maxheigth.gt.minheigth) then
              pboomcatch(1) = (1.0 - min((max(h0wav-minheigth,0.0))/(maxheigth-minheigth),1.0)) ** (idelt / 86400.0)!prob is given per day, scaled for the timestep
            else
              pboomcatch(1) = 1.1
            endif
            
           endif
            a = rnd(rseed)
            boomseffective =  a .lt. pboomcatch(1)
!            if (za .lt.-boomdepth(i)) boomseffective = .false.
         end if

   30    if ( max( xnew-1.0, -xnew ) .lt. max( ynew-1.0, -ynew ) ) then

            idy = 0                             !     y direction first
            if ( ynew .gt. 1.0) then
               idy =  1
               yy  = vy1                        ! velocity at downstream edge
            endif
            if ( ynew .lt. 0.0) then
               idy = -1
               yy  = vy0
            endif
            if ( idy  .ne. 0  ) then
               if (np+idy .gt. nmax) then
                  write(lun2,*) 'particle is now outside grid'
                  write(lun2,*) 'ipart=',ipart
                  write(lun2,*) 'mp   =',mp
                  write(lun2,*) 'np   =',np
                  write(lun2,*) 'kp   =',kp
                  write(lun2,*) 'xp   =',xp
                  write(lun2,*) 'yp   =',yp
                  write(lun2,*) 'zp   =',zp
                  npnew = np
                  mpnew = mp
                  xpnew = xnew
                  goto 49
               endif
               n0new = lgrid( np+idy  , mp    )
               dstick = .false.                 !  logical to determine sticking
               if ( n0new .lt. -nbmax .and. yy .ne. 0.0 ) then
                  xpnew = min(1.0,max(0.0,xnew))
                  xrest = (xnew-xpnew)*dxp
                  n0new = - n0new - nbmax
                  ddshift = 2
                  call p10ddb ( nconn , conn  , n0new  , ddshift, npnew  ,         &
      &                         mpnew , xpnew , ypnew  )
                  n0new = lgrid( npnew, mpnew)
                  if ( n0new .gt. 0 ) then
                     dxpnew = dx( n0new )
                     xpnew  = xpnew + xrest/dxpnew
                  endif
               else
                  npnew = np + idy
                  mpnew = mp
                  xpnew = xnew
               endif
               if ( n0new .gt. 0 .and.     &    !  new cell is active
                    yy    .ne. 0.0     ) then   !  and there is no thin dam
                  depth2 = depth( n0new )
                  if ( deppar .lt. depth2 ) then
                     n0   = n0new
                     np   = npnew
                     mp   = mpnew
                     n1   = lgrid2(np - 1, mp    )
                     n2   = lgrid2(np    , mp - 1)
                     ynew = ( ynew - 0.5*(1+idy) ) * dyp
                     dxp  = dx( n0 )
                     dyp  = dy( n0 )
                     ynew = ynew / dyp + 0.5*(1-idy)
                     vol  = volume( n0 + idep          )
                     vy0  = flow  ( n1 + idep          ) / vol
                     vy1  = flow  ( n0 + idep          ) / vol
                     vx0  = flow  ( n2 + idep + mnmaxk ) / vol
                     vx1  = flow  ( n0 + idep + mnmaxk ) / vol
                  else
                     ynew = (1.0+idy) - ynew    !  it bounces !
                     if ( lstick .and. lbott ) dstick = .true.
                  endif
               else          ! new cell is open boundary or inactive
                  ynew = (1.0+idy) - ynew       ! it bounces !
                  if ( lstick .and. lland .and. n0new .eq. 0 )  &
                                               dstick = .true.
               endif
               if ( dstick ) then
                  ynew = 0.5 + idy*0.4999
                  lstick = .false.
                  do isub = 1, nosubs
                     jsub = mstick(isub)
                     if ( jsub .gt. 0 ) then
                         if (wpart(jsub,ipart) .eq. 0.0 ) then
                            wpart(jsub,ipart) = wpart(isub,ipart)
                            wpart(isub,ipart) = 0.0
                         endif
                     endif
                  enddo
               endif
               goto 30
            endif
         else
            idx = 0                            !     now x direction
            if ( xnew .gt. 1.0) then
               idx =  1
               xx  = vx1
            endif
            if ( xnew .lt. 0.0) then
               idx = -1
               xx  = vx0
            endif
            if ( idx .ne. 0 ) then
               n0new = lgrid( np, mp + idx )
               dstick = .false.                !  logical to determine sticking
               if ( n0new .lt. -nbmax .and. xx .ne. 0.0 ) then
                  ypnew = min(1.0,max(0.0,ynew))
                  yrest = (ynew-ypnew)*dyp
                  n0new = - n0new - nbmax
                  ddshift = 1
                  call p10ddb ( nconn , conn  , n0new  , ddshift, npnew  ,         &
      &                         mpnew , xpnew , ypnew  )
                  n0new = lgrid( npnew, mpnew)
                  if ( n0new .gt. 0 ) then
                     dypnew = dy( n0new )
                     ypnew  = ypnew + yrest/dypnew
                  endif
               else
                  npnew = np
                  mpnew = mp + idx
                  ypnew = ynew
               endif
               if ( n0new .gt. 0 .and.    &    !  new cell is active
                    xx    .ne. 0.0     ) then  !  and there is no thin dam
                  depth2 = depth( n0new )
                  if ( deppar .lt. depth2 ) then
                     n0   = n0new
                     np   = npnew
                     mp   = mpnew
                     n1   = lgrid2(np - 1, mp    )
                     n2   = lgrid2(np    , mp - 1)
                     xnew = ( xnew - 0.5*(1+idx) ) * dxp
                     dxp  = dx( n0 )
                     dyp  = dy( n0 )
                     xnew = xnew / dxp + 0.5*(1-idx)
                     vol  = volume( n0 + idep          )
                     vy0  = flow  ( n1 + idep          ) / vol
                     vy1  = flow  ( n0 + idep          ) / vol
                     vx0  = flow  ( n2 + idep + mnmaxk ) / vol
                     vx1  = flow  ( n0 + idep + mnmaxk ) / vol
                  else
                     xnew = (1.0+idx) - xnew    !  it bounces !
                     if ( lstick .and. lbott ) dstick = .true.
                  endif
               else
                  xnew = (1.0+idx) - xnew       !  it bounces !
                  if ( lstick .and. lland .and. n0new .eq. 0 )  &
                                               dstick = .true.
               endif
               if ( dstick ) then
                  xnew = 0.5 + idx*0.4999
                  lstick = .false.
                  do isub = 1, nosubs
                     jsub = mstick(isub)
                     if ( jsub .gt. 0 ) then   ! ==> sticking
                          if (wpart(jsub,ipart) .eq. 0.0 ) then
                             wpart(jsub,ipart) = wpart(isub,ipart)
                             wpart(isub,ipart) = 0.0
                          endif
                     endif
                  enddo
               endif
               goto 30
            endif
         endif

         if (   kp .le. 0 .or. &
              ( kp .gt. layt .and. .not. twolay ) ) then
            write(*,*) ' Particle = ',ipart, &
                       ' not on an active layer:', kp, layt
            write(*,*) ' Programming error in rtim in part10'
            call stop_exit(1)
         endif

!         if (boomseffective) then
         if (booms) then
            
   !**      compute absolute x's and y's for a single particle end point

            call part11sp ( lgrid , xcor  , ycor  , nmax   , np     , mp    ,    &
                            xnew  , ynew  , xanew , yanew  , lgrid2 , mmax  )
            if (iptime(ipart) .le. 0 .and. nboomtry==1) then  !fmk: can iptime be les than 0? therefore it never executes the following statment, dhy is it here?
!      determine absolute location of starting point as well for new particles
               call part11sp ( lgrid , xcor  , ycor  , nmax   , npart(ipart)     , mpart(ipart)    ,    &  !  new coordinates
                               xpart(ipart)  , ypart(ipart)  , xaold , yaold  , lgrid2 , mmax  )
            else
               xaold = xa(ipart)
               yaold = ya(ipart)

            end if
 !          endif !boomseffective
            npolbounce = 0
            npolinside = 0
!            catch = .false.
!            bounce = .false.
!           in case of sedimentscreens, only check on the boom bouncing if the particle depth is smaller than screen depth 
            if (nboomtint.gt.0.and.boomseffective) then  !in case we deal with a moving boom, we need to check whether the particle is passed by the moving boom.
              do irow=1, nrowsboom(1)-1
               call pinpok(xaold, yaold, 5 , xpolb(1:5,irow), ypolb(1:5,irow), inside)
                if(inside.eq.1) then
                  npolinside = irow
!                  catch = .true.
!                  bounce = .true.
!                  write(lun2,'(a16,i8,2x,i8,2x,2f20.2)')' inside parts: ',irow,ipart, xaold, yaold
                  ! we need to move the particle to outside the polygon, eerste twee zijn van vorige tijdstap.
                  call pmov_boom(xaold, yaold, 5, xpolb(1:5,irow), ypolb(1:5,irow), xanew,yanew)
                  xaold = xanew  ! moved particle locations 
                  yaold = yanew
                  ! make sure that the old particles that are moved are put in the grid system as well
                  call part07nm ( lgrid  , lgrid2 , nmax   , mmax   , xcor  , & ! make relative
                                  ycor   , xaold , yaold, np   , &        ! coordinates
                                  mp, xnew, ynew  , ierror )                    ! again
                     npart(ipart) = np
                     mpart(ipart) = mp
                     xpart(ipart) = xnew
                     ypart(ipart) = ynew
     
                endif
              enddo
                         
              iboomtot = nboomint + 1  ! This implies only one moving boom is allowed
!              if (itime.ge.181500) write(lun2,*)itime, ipart
              call boombounce( xaold, yaold, xanew, yanew, nrowsboom(iboomtot), &
                   xipolboom(1:nrowsboom(iboomtot), iboomtot), &
                   yipolboom(1:nrowsboom(iboomtot), iboomtot), &
                   xacatch, yacatch, catch, xabounce, yabounce, bounce, &
                   xcrossboom, ycrossboom, npolbounce )
              if (npolinside.gt.0) then
                catch =.true.
                bounce=.true.
                xabounce = xanew
                yabounce = yanew   ! to ensure that the particles when moved are placed on the right location (relative) in part07m
                endif
 !             else
 !                  xaold = xanew  ! moved particle locations 
 !                  yaold = yanew
 !              endif
                 

            else
              call boombounce( xaold, yaold, xanew, yanew, nrowsboom(iboomint), &
                   xpolboom(1:nrowsboom(iboomint), iboomint), &
                   ypolboom(1:nrowsboom(iboomint), iboomint), &
                   xacatch, yacatch, catch, xabounce, yabounce, bounce, &
                   xcrossboom, ycrossboom, npolbounce )
            endif
            if (catch .and. (.not.sedscreen .or. (za(ipart) .ge. -boomdepth(iboomint+iboomtint)))) then 
               eff_angle = 1.1
               a = rnd(rseed)
               if (boom_proc) then
                ! here we can adapt the efficiency due to angle of currents with boomsegment being crossed (
                ! here we are dealing with particles being held by the boom, so we need to let a few more through
                ! this is only active when the boom processes have been switched on 
                 if ((ycrossboom(2)-ycrossboom(1))/=0.0) then
                   boomangle = atand(abs((xcrossboom(2)-xcrossboom(1))/(ycrossboom(2)-ycrossboom(1))))
                 else
                   boomangle = 90.0
                 endif
                 if (vyr/=0.0) then
                   curangle = atand(abs(vxr/vyr))
                 else
                   curangle = 90.0
                 endif
                 diffangle =  mod(abs(boomangle - curangle),90.)
!                write (lun2,*), "diffangle", diffangle

                 if (diffangle.gt.minangle)then
                   ! here is the dependency of the efficiency with the current speed and 'angle of attack' 
                   ! eff_angle when 1 does not change the dependency of the waves, but when smaller will let more particles through, 
                   ! all corrected for the timestep (in seconds)
                   eff_angle = max(1.0-max(((sqrt(vxr*vxr+vyr*vyr)*sind(diffangle-minangle)-boom_vlim)/boom_vmax) ** (itdelt/86400.0),0.0),0.0)
                 else
                   eff_angle = 1.1
                 endif
               endif
               if ((a .lt. eff_angle).and.boomseffective) then
               if (bounce .and. bouncefirsttry) then
!                 go back to relative coordinates from the original cell(!), set boomseffective to false (!?) and go back to check for checking of inactive cells etc... (30)
!                 a new bounce may mean that the booms must be effective again... hm, more complicated than I thought... Skip for the first implementation
!                  nm   = lgrid2( n  , m   )
                  np = npart(ipart)
                  mp = mpart(ipart)
                  kp = kpart(ipart)
                  znew = zpart(ipart)
                  call part07nm ( lgrid  , lgrid2 , nmax   , mmax   , xcor  , & ! make relative
                                  ycor   , xabounce , yabounce, np   , &        ! coordinates
                                  mp, xnew, ynew  , ierror )                    ! again
                  if (ierror/=0) then
                        mpart(ipart) = mpart0(ipart)
                        npart(ipart) = npart0(ipart)
                  end if
                  bouncefirsttry = .false.
				  
                  a = rnd(rseed)
                    !check whether boom extract is switched on and the bouncing takes place in one of the defined extraction sections
                    
                    if (a.lt.bstick_prob.or.boom_extract)then
                      ! if the material hits the boom some are removed (eg cleanup), take them out of computation
                      ! if extract keyword is swiched on then look at the lables of the polygon to determine in which 
                      ! section of the polygon the partiles are removed, using the same probability
                      extract = .false.
                      if (a.lt.extract_prob) then
                        if (boom) then
                          do irow =1 ,  nr_extract
                            if (extract_sect(irow,iboomint).eq.npolbounce) extract = .true.  !do the x,y coords need to be reset here?
                          enddo
                        elseif(moving_boom.and.((npolinside.ne.0).or.(npolbounce.ne.0))) then
                          do irow =1 , nr_extract
                            if (extract_sect(irow,iboomtint).eq.npolinside) extract = .true.  !do the x,y coords need to be reset here?
                            if (extract_sect(irow,iboomtint).eq.npolbounce) extract = .true.  !do the x,y coords need to be reset here?
                          enddo
                        endif
                      endif
                      
                     if (a.lt.bstick_prob.and.(.not.extract).and.boomseffective) then
                       extract = .true.
                     endif
                     if (extract) then
                       np=1
                       mp=1
                       kp=1
                       iptime(ipart) = 0
                       xnew = 0.0
                       ynew = 0.0
                       znew = 0.0
                       do ifract =1, nfract
                         wrem(ifract) = wrem(ifract)+wpart((ifract-1)*3+1,ipart)
!                       write(lun2,'(a5,i10,a7,i10,a10,i3,a12, i3,a7,e10.3)')'time: ', itime, ' ipart:', ipart,'nboomtry:',nboomtry,'npolbounce:', npolbounce, 'wpart:',wpart(1,ipart)
                         wpart((ifract-1)*3+1,ipart) = 0.0  !needs to be set to zero because the loop runs past this point several times and there is double counting of wrem otherwise?
                       enddo
                     endif
                      !keep track of removed mass

                  else
                    goto 30   !
                  endif
				  
               else
                  if (nboomtry==1.and.boomseffective) then
                     bouncefirsttry = .true.
!                 go back to the original location and try with dispersion and wind only...
                     np = npart(ipart)
                     mp = mpart(ipart)
                     kp = kpart(ipart)
                     xnew = xpart(ipart)
                     ynew = ypart(ipart)
                     znew = zpart(ipart)
                     !**      xnew,ynew is inclusive of diffusion (dax,day); c1 = 0.0 for 3d
                     !             drag on the difference vector: cd * (wind - flow)
                     xnew  = xnew + (dax - c1 * sin(wdirr + sangl)) / dxp + (cdrag*(vxw-vxr)/dxp) * itdelt
                     ynew  = ynew + (day - c1 * cos(wdirr + sangl)) / dyp + (cdrag*(vyw-vyr)/dyp) * itdelt
                     nboomtry = 2
                     goto 30
                  else if (nboomtry==2.and.boomseffective) then
                     bouncefirsttry = .true.
!                 go back to the original location and try with dispersion only...
                     np = npart(ipart)
                     mp = mpart(ipart)
                     kp = kpart(ipart)
                     xnew = xpart(ipart)
                     ynew = ypart(ipart)
                     znew = zpart(ipart)
                     !**      xnew,ynew is inclusive of diffusion (dax,day); c1 = 0.0 for 3d
                     !             drag on the difference vector: cd * (wind - flow)
                     xnew  = xnew + (dax - c1 * sin(wdirr + sangl)) / dxp
                     ynew  = ynew + (day - c1 * cos(wdirr + sangl)) / dyp
                     nboomtry = 3
                     goto 30
                  else
!                 finally just put it back to the old location!
                     if (boomseffective) then
                     np = npart(ipart)
                     mp = mpart(ipart)
                     kp = kpart(ipart)
                     xnew = xpart(ipart)
                     ynew = ypart(ipart)
                     znew = zpart(ipart)
                     endif
                  end if
               end if
            end if
           end if
         end if

!**      assignment of the coordinates and other properties

  49 continue
!    here we can check on old locations inside the polygon? If so then move it 
!                 but we need to known by how much - calculate distance to the line.
!                  a = rnd(rseed)
                    !check whether boom extract is switched on and the bouncing takes place in one of the defined extraction sections
                    
!                    if (catch.and.boom_extract)then
                      ! if the material hits the boom some are removed (eg cleanup), take them out of computation
                      ! if extract keyword is swiched on then look at the lables of the polygon to determine in which 
                      ! section of the polygon the partiles are removed, using the same probability
!                      extract = .false.
!                      if (a.lt.extract_prob) then
!                        do irow =1 , nrowsboom(iboomint)
!                         if (extract_sect(irow,iboomint).eq.npolbounce) extract = .true.  !do the x,y coords need to be reset here?
!                        enddo
!                      endif
                      
!                     if (a.lt.bstick_prob.and.(.not.extract).and.boomseffective) then
!                       extract = .true.
!                     endif
!                     if (extract) then
!                       np=1
!                       mp=1
!                       kp=1
!                       iptime(ipart) = 0
!                       wrem = wrem+wpart(1,ipart)
!                       write(lun2,*)'1, ipart is :', ipart,'nboomtry: ',nboomtry, 'npolbounce is:', npolbounce
!                     endif
                      !keep track of removed mass
!                    endif
         npart (ipart) = np
         mpart (ipart) = mp
         kpart (ipart) = kp
         xpart (ipart) = xnew
         ypart (ipart) = ynew
         zpart (ipart) = znew

!**      calculate decay wpart for all substances, adapt for smooth loading with ddfac

   50    do isub = 1, nosubs
            if ( ddfac .gt.  1.5 ) then
               wpart(isub, ipart) = wpart(isub, ipart) * dfact(isub)
            else
               wpart(isub, ipart) = wpart(isub, ipart) *  &
                                                     dfact(isub)**ddfac
            endif
         enddo
!         vdiff( n0 ) = disp

!**     end of particle loop

   90    iptime(ipart) = iptime(ipart) + idelt
  100 continue
!$OMP END PARALLEL DO
       timon = .true.
      if (booms.and.(boom_extract.or.boom_stick))then
      write(lun2,'(4x,a,i12,3(2x,es15.7))')' Mass removed at boom(s): ', itime, (wrem (ifract),ifract =1, nfract)
      endif
      
        

      if (lsettl) then
         write(lun2,'(4x,a,i12,a,i4,a)') '  No. of particles settled into bed layer  : ', &
                                          nopart_sed,' (layer ',layt+1,')'
         write(lun2,'(4x,a,i12,a,i4,a)') '  No. of particles eroded from  bed layer  : ', &
                                          nopart_ero,' (layer ',layt+1,')'
      endif

!**   some statistics on vertical diffusion

      if ( nopart .eq. 0 ) then
         dspmin = 0.0
         dspmax = 0.0
         dsprms = 0.0
      else
         dsprms = sqrt(dsprms/max(nrms,1.0))
      endif
      write(lun2,'(6x,a,i12)') &
           'No. of particles in inactive grid cells  : ', ninact
      write(lun2,'(6x,a,i12)') &
           'No. of sticking particles                : ', nstpar
      write(lun2,'(6x,a,es15.7)') &
           'Minimum  (*) for vertical dispersion     : ', dspmin
      write(lun2,'(6x,a,es15.7)') &
           'Maximum  (*) for vertical dispersion     : ', dspmax
      write(lun2,'(6x,a,es15.7)') &
           'Rms value(*) for vertical dispersion     : ', dsprms
      write(lun2,'(15x,a)') '(*) over time and space(=particles)'
!**
!**   apply rotation on (3d) floating oil particles due to effect of waves (coriolis)
!**
      if ( threed .and. oilmod .and. coriol ) then

!**   disable vertically displacing particles

         do ipart = npwndw, nopart
            if ( floil(ipart) .eq. 0 ) then
               npart0(ipart) = -npart0(ipart)
               mpart0(ipart) = -mpart0(ipart)
               npart (ipart) = -npart (ipart)
               mpart (ipart) = -mpart (ipart)
            endif
         enddo

!**      compute absolute x's and y's  ( z's are dummy )

         call part11 ( lgrid , xcor  , ycor  , nmax  , npart0,    &  !  old coordinates
                       mpart0, xpart0, ypart0, xa0   , ya0   ,    &
                       nopart, npwndw, lgrid2, kpart , zpart ,    &
                       za    , locdep, dps   , nolay , mmax  ,    &
                       tcktot)
         call part11 ( lgrid , xcor  , ycor  , nmax  , npart ,    &  !  new coordinates
                       mpart , xpart , ypart , xa    , ya    ,    &
                       nopart, npwndw, lgrid2, kpart , zpart ,    &
                       za    , locdep, dps   , nolay , mmax  ,    &
                       tcktot)

!**      rotate vector by multiplication with rotation matrix (FMK: no  need if the particle has been removed (extract=.true.

         codef = cos(-defang)
         sidef = sin(-defang)
         do ipart = npwndw, nopart
            if ( floil(ipart) .eq. 0 ) then
               npart(ipart)= abs(npart(ipart))
               mpart(ipart)= abs(mpart(ipart))
            else
!               if (npart(ipart).ne.1) then ! FMK: when the particle has been taken (n,m,k = 1) out of the computation then do not calculate these coordinates
                 dxx = xa(ipart)-xa0(ipart)
                 dyy = ya(ipart)-ya0(ipart)
                 xx  = xa0(ipart) + dxx * codef - dyy * sidef
                 yy  = ya0(ipart) + dxx * sidef + dyy * codef
                 call part07 ( lgrid  , lgrid2 , nmax   , mmax   , xcor  , &! make relative
                             ycor   , xx     , yy     , npart(ipart)   , &! coordinates
                             mpart(ipart), xpart(ipart), ypart(ipart)  , &! again
                             ierror )
                 if (ierror/=0) then
                   !
                   ! Location after adding the deflection angle correction is outside the grid,
                   ! so the deflection angle correction will not be performed
                   !
                   mpart(ipart) = mpart0(ipart)
                   npart(ipart) = npart0(ipart)
                 endif
!               endif
            endif
         enddo
      endif

!**   end of subroutine

      if ( timon ) call timstop ( ithndl )
      return
      end subroutine part10

      logical function eql( a, b )

!     function to check if two floating point numbers are equal
!     by checking on the number of equal digits
      implicit none           !   force explicit typing
      real :: a , b
      real, parameter ::  epsrel = 0.0001

      if      ( b  .ne. 0.0 ) then
         eql = abs( a/b - 1.0 ) .lt. epsrel
      else if ( a  .eq. 0.0 ) then
           eql = .true.
      else
           eql = .false.
      endif

      return
      end function
      
end module

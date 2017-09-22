!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017.                                     
!                                                                               
!  This file is part of Delft3D (D-Flow Flexible Mesh component).               
!                                                                               
!  Delft3D is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  Delft3D  is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.             
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D",                  
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting 
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!                                                                               
!-------------------------------------------------------------------------------

! $Id: d3ddimens.f90 52266 2017-09-02 11:24:11Z klecz_ml $
! $HeadURL: https://repos.deltares.nl/repos/ds/branches/dflowfm/20161017_dflowfm_codecleanup/engines_gpl/dflowfm/packages/dflowfm_kernel/src/d3ddimens.f90 $
!-------------------------------------------------------------------------------------------------------
!  Origin:
!     URL: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/flow2d3d/packages/data/include/dimens.igs
!     Revision: 5108

module m_d3ddimens
implicit none
private
    public gd_dimens
    public d3dgrid_dimens

    type gd_dimens
!
       integer :: ncmax
       integer :: nmax       !  Number of gridpoints in the y-dir. (always odd)
       integer :: mmax       !  Number of gridpoints in the x-dir.
       integer :: nlb        !  Minimum index in N direction
       integer :: nub        !  Maximum index in N direction
       integer :: mlb        !  Minimum index in M direction
       integer :: mub        !  Maximum index in M direction
       integer :: nmlb       !  Minimum index when transformed to 1D array
       integer :: nmub       !  Maximum index when transformed to 1D array
       integer :: ddbound
       integer :: nmaxus     !  Description and declaration in esm_alloc_int.f90
       integer :: kmax       !  Description and declaration in esm_alloc_int.f90
       integer :: nmaxd      !  If KMAX > 0 Then := NMAX else := 1
       integer :: mmaxd      !  If KMAX > 0 Then := MMAX else := 1
       integer :: jstart     !  Begin pointer for arrays which have been
                             !  transformed into 1D arrays. Due to the shift in
                             !  the 2nd (M-) index, JSTART = -2*NMAX + 1
       integer :: nmmaxj     !  End   pointer for arrays which have been
                             !  transformed into 1D arrays. Due to the shift in
                             !  the 2nd (M-) index, NMMAXJ = NMMAX + 2 * NMAX
       integer :: nmmax      !  Total number of grid pts. (NMAX*MMAX)
       integer :: lmax       !  Number of constituents (incl. turb. energy
                             !  dissipation and production): LSTSCI+ LTUR
       integer :: lmaxd      !  Maximum of (1,LMAX)
       integer :: lsts       !  Number of Constituents (Salinity, Temperature
                             !  and Suspended Sediment)
       integer :: lstsc      !  Number of Constituents (Salinity, Temperature,
                             !  Suspended Sediment & Conservative Constituents)
       integer :: lstsci     !  Total number of transporable quantities: all
                             !  constituents plus secondary flow intensity.
       integer :: lsal       !  Pointer for Salinity in array R1 for
                             !  constituents (if used, always 1)
       integer :: lsed       !  Number of Suspended Sediment fractions.
       integer :: lsedtot    !  Total number of Sediment fractions (including
                             !  bedload only fractions).
       integer :: ltem       !  Pointer for Temperature in array R1 for
                             !  constituents (if used, LSAL+1)
       integer :: lsecfl     !  Pointer for secondary flow intensity in array
                             !  R1 for constituents
       integer :: lsec       !  Flag for secondary flow
                             !      1 = no equilibrium (advection and diffusion
                             !          effects included)
                             !      2 = equilibrium
       integer :: ltur       !  Description and declaration in esm_alloc_int.f90
       integer :: ltur2d     !  Flag for 2D turbulence model
                             !     (0 = constant      model
                             !      1 = Uittenbogaard model)
       integer :: kmxdt      !  Dimension for IWE equidistant layer array
                             !  Maximum of (1,KMXT)
                             !  Maximum size of interpolation arrays for IWE
       integer :: nbub       !  Total number of 'artificial' discharge points,
                             !  added automatically to model the bubble screens
       integer :: nxbub      !  Number of bubble screens
                             !  Each screen consists of a line of bubble locations
                             !  For each location, for each vertical layer, an 'artificial' discharge point is added
       integer :: npiwe      !  Dimension for IWE frequency arrays
                             !  Maximum of (1,NFREQS)
                             !  Maximum size of frequency arrays for IWE
       integer :: kmxt       !  Nr.of equidistant layers in IWE model
       integer :: nfreqs     !  Number of angular frequency intervals for
                             !  inspecting roots of the TG equation (IWE)
       integer :: nlcest     !  Estimated number of computational rows and
                             !  columns in IROCOL-table = 5 * max(NMAX,MMAX)
       integer :: noroco     !  Number of computational rows & cols
       integer :: norow      !  Number of computational rows in IROCOL-array 
       integer :: nocol      !  Number of computational columns in IROCOL-array 
       integer :: nto        !  Total number of open boundary sections
       integer :: ntof       !  Number of open boundary sections of the Fourier type
                             !  Number of open boundary sections of the Harmonic type
       integer :: ntoq       !  Number of open boundary sections of the QH-rel. type
       integer :: ntot       !  Number of boundary sections of the Time varying type
       integer :: kc         !  Number of Frequencies (incl. mean value) for the Hydrodynamic B.C.
       integer :: kcd        !  Maximum of (KC,KMAX)
       integer :: nopest     !  Estimated number of open boundary
                             !  points = 4 * (NMAX  +  MMAX  )
       integer :: nrob       !  Number of open boundary points
       integer :: nsrc       !  Total number of discharges (source or sink,
                             !  including all automatically added 'artificial' discharge
                             !  points used to model bubble screens
       integer :: nsrcd      !  Number of (ordinary) discharge points plus the number of bubble screens
       integer :: nostat     !  Number of monitoring stations
       integer :: ntruv      !  Total nr. of monit. cross sections
       integer :: ntru       !  Nr. of monitoring cross sections (U)
       integer :: nofou      !  Number of requested fourier analysis
       integer :: ndro       !  Number of released drogues
       integer :: nsluv      !  Number of U- and V-barriers
       integer :: upwsrc     !  Flag for upwind of discharges
                             !  -1 : no upwind for all discharges
                             !   0 : only upwind for momentum discharges
                             !   1 : upwind for all discharges
       integer :: ndx        !  Number of unstructured flownodes (internal + boundary) 
       integer :: lnx        !  Number of unstructured flowlinks (internal + boundary) 
       integer :: ndkx       !  Number of unstructured flownodes (internal + boundary) ! 3D
       integer :: lnkx       !  Number of unstructured flowlinks (internal + boundary) 
    end type gd_dimens

    
    

contains
   subroutine d3dgrid_dimens(gddimens,gddimens_ptr,nmax,mmax,ndx,lnx,ndkx,lnkx)
   !!--description-----------------------------------------------------------------
   !
   ! Initializes a gd_dimens structure (Delft3D) for a simple (1,nmax) and (1,mmax) grid
   ! (Copied from simplegrid_dimens)
   ! Source: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/utils_lgpl/deltares_common/packages/deltares_common/src/grid_dimens_module.f90 
   ! Original licence: LGPL
   !!--declarations----------------------------------------------------------------
       implicit none
   !
   ! Call variables
   !
       type (gd_dimens), target, intent(inout) :: gddimens
       type (gd_dimens), pointer               :: gddimens_ptr 
       integer,                  intent(in)    :: nmax
       integer,                  intent(in)    :: mmax
       integer,                  intent(in)    :: ndx 
       integer,                  intent(in)    :: lnx 
       integer,                  intent(in)    :: ndkx 
       integer,                  intent(in)    :: lnkx 
   !
   ! Local variables
   !
       !integer                   :: istat
   !
   !! executable statements -------------------------------------------------------
   !
       gddimens%nlb    = 1
       gddimens%nub    = nmax
       gddimens%nmax   = nmax
       !
       gddimens%mlb    = 1
       gddimens%mub    = mmax
       gddimens%mmax   = mmax
       !
       gddimens%nmlb   = 1
       gddimens%nmub   = nmax*mmax
       gddimens%nmmax  = nmax*mmax
       !
       ! TODO: Find the counterparts of these fields in the gd_dimens structure and uncomment the lines below  
       !gddimens%nfg    = 1
       !gddimens%nlg    = nmax
       !gddimens%nmaxgl = nmax
       !!
       !gddimens%mfg    = 1
       !gddimens%mlg    = mmax
       !gddimens%mmaxgl = mmax
       !
       !allocate(gddimens%celltype(gddimens%nmmax), stat=istat)
       !if (istat==0) gddimens%celltype = 1
       !gddimens%nmbnd => null()

       gddimens%ndx    = ndx 
       gddimens%lnx    = lnx 
       gddimens%ndkx   = ndkx 
       gddimens%lnkx   = lnkx 

       gddimens_ptr => gddimens
   end subroutine d3dgrid_dimens


end module m_d3ddimens

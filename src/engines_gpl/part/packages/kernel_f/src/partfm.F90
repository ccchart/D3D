!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2021.                                
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

   subroutine partfm(lunpr)

   use MessageHandling
   use partmem, only: hyd, itstrtp, itstopp
   use m_particles, only: Npart, NpartTot
   use alloc_mod
   use fileinfo  , lun=> lunit    ! logical unit numbers for files
   use m_flowtimes
   use timers

   implicit none

   integer(ip), intent(in) :: lunpr

   integer(ip)         :: lunmem, istat, itime, ilay, Ldot

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "partfm", ithndl )

   call SetMessageHandling(lunMessages=lunpr)
   ! initialize the allocation system
   call init_alloc( lunmem , lunpr )

   if (hyd%nolay /= 1) then
      write ( lunpr, * ) ' ERROR: 3D hydrodynamics is not yet supported for unstructured grids!'
      write (   *  , * ) ' ERROR: 3D hydrodynamics is not yet supported for unstructured grids!'
   endif
   dts   = real(hyd%cnv_step_sec, 8)  !idelt
   tzone = 0.0_hp
   refdat = '20200701'
   call setTUDUnitString()

   call ini_part_grid(hyd)

   call rdpart ( lun(1)   , lun(2)   , fname(1) )
   tstart_user = itstrtp
   tstop_user  = itstopp
   Ldot = index(fname(1), '.', .true.)-1
   filebase = ' '
   filebase = fname(1)(1:Ldot)

   call ini_part(partinifile, partrelfile, tstart_user, dts, 0)

   call unc_init_trk()
   call unc_init_map(hyd%crs, hyd%waqgeom)

   time0 = tstart_user
   time1 = time0
   istat = -1 ! skip copying of data durin the first stime stap
   call part_readhydstep(hyd,itime,istat)
   if (istat == 99) then
      write ( lunpr, * ) ' Timing mismatch between input and actual data:', time0, itime
      write (   *  , * ) ' Timing mismatch between input and actual data:', time0, itime
      goto 1234
   else if (istat /= 0) then
      write ( lunpr, * ) ' Error during reading of the hydrodynamic time step data', time0, itime
      write (   *  , * ) ' Error during reading of the hydrodynamic time step data', time0, itime
      goto 1234
   end if
   do while (istat == 0)
      call unc_write_trk()
      call unc_write_map()
      write ( lunpr, '('' Current timestep: '',I15,'' of '',I15)' ) int(time0), int(time1)
      write (   *  , '('' Current timestep: '',I15,'' of '',I15)' ) int(time0), int(time1)
      write ( lunpr, '('' Active particles: '',I15,'' of '',I15)' ) Npart, NpartTot
      write (   *  , '('' Active particles: '',I15,'' of '',I15)' ) Npart, NpartTot
      if (time1 .ge. tstop_user) then
         exit
      endif
      time0 = time1
      time1 = min(tstop_user, time0 + dts)
      call part_readhydstep(hyd,itime,istat)
      call update_part()
   end do
   
1234 continue
   
   call unc_close_trk()
   call unc_close_map()

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine partfm

   subroutine update_particles(q,h0,h1,Dt)
   use m_particles
   use m_flowgeom, only: Ndx, Lnx
   use m_sferic
   use geometry_module, only: Cart3Dtospher
   use MessageHandling
   use timers

   implicit none

   double precision, dimension(Lnx), intent(in) :: q  !< fluxes
   double precision, dimension(Ndx), intent(in) :: h0 !< water deoths at start of time interval
   double precision, dimension(Ndx), intent(in) :: h1 !< water deoths at end of time interval
   double precision,                 intent(in) :: Dt !< time interval

   integer, dimension(1) :: numremaining ! number of remaining particles to be updated

   double precision :: xx, yy

   integer :: i
   integer :: iter
   integer :: ierror

   integer, parameter :: MAXITER = 1000 ! maximum number of substeps

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "update_particles", ithndl )

   ierror = 1

   ! reconstruct velocity field
   call reconst_vel(q, h0, h1, ierror)
   if ( ierror.ne.0 ) goto 1234

   if ( Npart.gt.0 ) then
      ! set remaining time to time step
      dtremaining = Dt
      ! Lpart = 0
      numzero = 0
   end if

   do iter=1,MAXITER
      ! update particles in cells
      call update_particles_in_cells(numremaining(1), ierror)
      if ( ierror.ne.0 ) goto 1234

      write(6,*) 'iter=', iter, 'numremaining=', numremaining(1)
      if ( numremaining(1).eq.0 ) then
         write(6,*) 'iter=', iter
         exit
      end if
   end do


   ! check for remaining particles
   if ( numremaining(1).gt.0 ) then
      ! plot remaining particles
      do i=1,Npart
         if ( dtremaining(i).gt.0d0 .and. kpart(i).gt.0 ) then
            if ( jsferic.eq.0 ) then
               xx = xpart(i)
               yy = ypart(i)
            else
               call Cart3Dtospher(xpart(i),ypart(i),zpart(i),xx,yy,0d0)
            end if
            write(6,"(I0, ':', 2E25.15, ', ', I0)") i, xx, yy, kpart(i)
         end if
      end do
      call mess(LEVEL_WARN, 'update_particles: iter>MAXITER')
      goto 1234
   end if

   ierror = 0
1234 continue

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine update_particles

   !> update positions of particles within triangles
   subroutine update_particles_in_cells(numremaining, ierror)
   use m_particles
   use m_partrecons
   use m_partmesh
   use MessageHandling
   use m_sferic, only: jsferic
   use timers

   implicit none

   integer,        intent(out) :: numremaining !< number of remaining particles to be updated
   integer,        intent(out) :: ierror       !< error (1) or not (0)

   integer                     :: ipart
   integer                     :: i, k, k1, k2, L
   integer                     :: ja
   integer                     :: Lexit

   double precision            :: d, un
   double precision            :: t, tex, dt
   double precision            :: ux0, uy0, uz0, cs, sn
   double precision            :: xn, yn, zn, rl
   double precision            :: dvar, dis, dn

   double precision, dimension(3) :: ddn

   logical                     :: isboundary

   double precision, parameter :: DTOL = 1d-4
   double precision, parameter :: DTOLd  = 1d-4
   double precision, parameter :: DTOLun_rel = 1d-4
   double precision, parameter :: DTOLun = 1e-14

   integer,          parameter :: MAXNUMZERO = 10

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "update_particles_in_cells", ithndl )

   ierror = 0
   numremaining = 0

!$OMP PARALLEL DO PRIVATE (i, k, k1, k2, L, ja, Lexit, d, un, t, tex, dt,          &
!$OMP                      ux0, uy0, uz0, cs, sn, xn, yn, zn, rl, dvar, dis, dn,   &
!$OMP                      ddn, isboundary),                                       &
!$OMP           REDUCTION ( +   : numremaining) ,                                  &
!$OMP           REDUCTION ( MAX : ierror ),                                        &
!$OMP           SCHEDULE  ( DYNAMIC, max(Npart/100,1)           )
   do ipart=1,Npart
      ! check if this particle needs to be updated
      if ( dtremaining(ipart).eq.0d0 .or. kpart(ipart).lt.1 ) cycle
      ! get cell (flownode) particle in in
      k = kpart(ipart)

      ! compute exit time <= dtremaining
      tex = dtremaining(ipart)

      Lexit = 0   ! exit edge (flowlink)

      ! compute velocity at current position
      ux0 = u0x(k) + alpha(k)*(xpart(ipart)-xzwcell(k))
      uy0 = u0y(k) + alpha(k)*(ypart(ipart)-yzwcell(k))
      if ( jsferic.ne.0 ) then
         uz0 = u0z(k) + alpha(k)*(zpart(ipart)-zzwcell(k))
      end if

      ! loop over edges (netlinks) of cells
      do i=jcell2edge(k),jcell2edge(k+1)-1
         L = icell2edge(i)   ! edge

         k1 = edge2node(1,L)
         k2 = edge2node(2,L)

         if ( jsferic.eq.0 ) then
            cs = dnx(1,L)
            sn = dny(1,L)
            if ( edge2cell(2,L).eq.k ) then
               cs = -cs
               sn = -sn
            end if
         else
            if ( edge2cell(1,L).eq.k ) then
               ddn = (/ dnx(1,L), dny(1,L), dnz(1,L) /)
            else
               ddn = (/ dnx(2,L), dny(2,L), dnz(2,L) /)
            end if
         end if

         ! check for boundary edge
         isboundary = ( edge2cell(1,L).eq.0 .or. edge2cell(2,L).eq.0 )

         ! compute normal distance to edge
         if ( jsferic.eq.0 ) then
            if ( isboundary ) then ! boundary: add tolerance
               call dlinedis2(xpart(ipart),ypart(ipart),xnode(k1)+cs*DTOLd,ynode(k1)+sn*DTOLd,xnode(k2)+cs*DTOLd,ynode(k2)+sn*DTOLd,ja,d,xn,yn,rl)
            else
               call dlinedis2(xpart(ipart),ypart(ipart),xnode(k1),ynode(k1),xnode(k2),ynode(k2),ja,d,xn,yn,rl)
            end if
            dis = (xn-xpart(ipart))*cs + (yn-ypart(ipart))*sn
         else
            if ( isboundary ) then ! boundary: add tolerance
               call dlinedis3D(xpart(ipart),ypart(ipart),zpart(ipart),xnode(k1)+DTOLd*ddn(1),  &
                  ynode(k1)+DTOLd*ddn(2),  &
                  znode(k1)+DTOLd*ddn(3),  &
                  xnode(k2)+DTOLd*ddn(1),   &
                  ynode(k2)+DTOLd*ddn(2),   &
                  znode(k2)+DTOLd*ddn(3),   &
                  ja,d,xn,yn,zn,rl)
            else
               call dlinedis3D(xpart(ipart),ypart(ipart),zpart(ipart),xnode(k1),ynode(k1),znode(k1),xnode(k2),ynode(k2),znode(k2),ja,d,xn,yn,zn,rl)
            end if
            dis = (xn-xpart(ipart))*ddn(1) + (yn-ypart(ipart))*ddn(2) + (zn-zpart(ipart))*ddn(3)
         end if

!        BEGIN DEBUG
!         if ( ipart.eq.1 .and. kpart(ipart).eq.5298 ) then
!            write(6,*) i, ':', d, rL, dis
!         end if
!
!         if ( abs(dis-d).gt.1d-1 ) then
!            write(6,*) i, dis, d
!         end if
!        END DEBUG

         ! check inside or outside triangle
         ! if ( dis.lt.-DTOLd .and. rL.ge.0d0 .and. rL.le.1d0 .and. .not.isboundary ) then
         if ( dis.lt.-DTOLd .and. .not.isboundary ) then
            ! outside triangle
            tex = 0d0
            Lexit = L
            exit
         else
            ! inside triangle
            ! compute normal velocity to edge (outward positive)
            if ( jsferic.eq.0 ) then
               un =  ux0*cs + uy0*sn
            else
               un =  ux0*ddn(1) + uy0*ddn(2) + uz0*ddn(3)
            end if

!!           BEGIN DEBUG
!!           check normal velocity at closed boundary
!            if ( edge2cell(1,L).eq.0 .or. edge2cell(2,L).eq.0 ) then
!               dvar = (u0x(k) + alpha(k)* (xn-xzwcell(k)))*ddn(1) + (u0y(k) + alpha(k)*(yn-yzwcell(k)))*ddn(2) + (u0z(k) + alpha(k)*(zn-zzwcell(k)))*ddn(3)
!               if ( abs(dvar) .gt. 1d-4 ) then
!                  continue
!               end if
!            end if
!!           END DEBUG

            if ( un.gt.max(DTOLun_rel*d,DTOLun) ) then   ! normal velocity does not change sign: sufficient to look at u0.n
               ! compute exit time for this edge: ln(1+ d/un alpha) / alpha
               dvar = alpha(k)*dis/un
               if ( dvar.gt.-1d0) then
                  t = dis/un
                  if ( abs(dvar).ge.DTOL ) then
                     t = t * log(1d0+dvar)/dvar
                  end if
               else
                  t = huge(1d0)
               end if

               ! update exit time/edge (flowlink)
               ! if ( t.le.tex .and. t.ge.0d0 ) then
               if ( t.le.tex ) then

                  tex = t
                  Lexit = L
               end if
            else
               continue
            end if

         end if
      end do

      if ( dtremaining(ipart).eq.0d0 ) then
         continue
      end if

      ! compute timestep in cell (flownode)
      dt = min(dtremaining(ipart), tex)

      ! update particle
      if ( abs(alpha(k)).lt.DTOL ) then
         dvar = dt
      else
         dvar = (exp(alpha(k)*dt)-1d0)/alpha(k)
      end if

      xpart(ipart) = xpart(ipart) + dvar * ux0
      ypart(ipart) = ypart(ipart) + dvar * uy0

      if ( jsferic.ne.0 ) then
         zpart(ipart) = zpart(ipart) + dvar * uz0
      end if

!!     BEGIN DEBUG
!      if ( jsferic.eq.1 ) then
!!        project node on triangle
!         dn = (xpart(ipart) - xzwcell(k)) * dnn(1,k) +  &
!              (ypart(ipart) - yzwcell(k)) * dnn(2,k) +  &
!              (zpart(ipart) - zzwcell(k)) * dnn(3,k)
!         xpart(ipart) = xpart(ipart) - dn * dnn(1,k)
!         ypart(ipart) = ypart(ipart) - dn * dnn(2,k)
!         zpart(ipart) = zpart(ipart) - dn * dnn(3,k)
!      end if
!!     END DEBUG

      dtremaining(ipart) = dtremaining(ipart) - dt
      ! Lpart(ipart) = Lexit

      if ( dt.eq.0d0 ) then
         numzero(ipart) = numzero(ipart) + 1
      end if

      if ( numzero(ipart).gt.MAXNUMZERO ) then
         ! disable particle that is not moving
         kpart(ipart) = 0
         dtremaining(ipart) = 0d0

         ! proceed to neighboring cell (if applicable)
      else if ( Lexit.gt.0 ) then
         numremaining = numremaining + 1  ! number of remaining particles for next substep
         if ( edge2cell(1,Lexit).gt.0 .and. edge2cell(2,Lexit).gt.0 ) then   ! internal edge (netlink)
            kpart(ipart) = edge2cell(1,Lexit) + edge2cell(2,Lexit) - k

            if ( kpart(ipart).eq.0 ) then
               continue
            else
               if ( jsferic.eq.1 ) then
                  ! project node on triangle
                  k = kpart(ipart)
                  k1 = edge2node(1,Lexit)
                  k2 = edge2node(2,Lexit)
                  xn = 0.5d0*(xnode(k1)+xnode(k2))
                  yn = 0.5d0*(ynode(k1)+ynode(k2))
                  zn = 0.5d0*(znode(k1)+znode(k2))
                  dn = (xpart(ipart) - xn) * dnn(1,k) +  &
                     (ypart(ipart) - yn) * dnn(2,k) +  &
                     (zpart(ipart) - zn) * dnn(3,k)
                  xpart(ipart) = xpart(ipart) - dn * dnn(1,k)
                  ypart(ipart) = ypart(ipart) - dn * dnn(2,k)
                  zpart(ipart) = zpart(ipart) - dn * dnn(3,k)
               end if
            end if
         else  ! on boundary
            kpart(ipart) = 0
         end if
      else
         ! safety check
         if ( dtremaining(ipart).ne.0d0 ) then
            ierror = 1
         end if
      end if

   end do
!$OMP END PARALLEL DO

   if (ierror == 1) then
      call mess(LEVEL_ERROR, 'update_particles_in_cells: dtremaining <> 0', ' ', ' ')
   endif

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine update_particles_in_cells

   !> reconstruct velocity in cells
   subroutine reconst_vel_coeffs()

   use m_partrecons
   use m_partmesh
   use m_alloc
   use m_sferic
   use geometry_module, only: dbdistance, gaussj
   use timers

   implicit none

   integer,                          parameter   :: N = 4

   double precision, dimension(N,N)              :: Amat ! matrix
   double precision, dimension(N)                :: rhs

   double precision                              :: cs, sn

   integer                                       :: i, icell, j, jj, k, L, NN
   integer                                       :: k1, k2
   integer                                       :: i12, isign
   integer                                       :: ierror

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "reconst_vel_coeffs", ithndl )

   ierror = 1

   ! allocate startpointers
   call realloc(jreconst, numcells+1, keepExisting=.false., fill=0)

   ! set startpointers
   jreconst(1) = 1
   do icell=1,numcells
      jreconst(icell+1) = jcell2edge(icell) + jcell2edge(icell+1)-jcell2edge(icell)
   end do

   ! allocate column indices and matrix entries
   NN = jreconst(numcells+1)-1
   call realloc(ireconst, NN, keepExisting=.false., fill=0)
   if ( jsferic.eq.0 ) then
      call realloc(Areconst, (/3, NN /), keepExisting=.false., fill=0d0)
   else
      call realloc(Areconst, (/4, NN /), keepExisting=.false., fill=0d0)
   end if

   ! dummy rhs
   rhs = 0d0

   ! loop over internal cells
   jj = 0
   do icell=1,numcells
!     check for triangles
!      if ( jcell2edge(k+1)-jcell2edge(k).ne.3 ) then
!         cycle
!      end if

      ! get flownode number (for s, bl)
      k = iabs(cell2nod(icell))

      ! fill system for (ux,uy) = (ux0, uy0) + alpha (x-x0, y-y0)

      ! loop over edges (netlinks)
      i = 0
      do j=jcell2edge(icell),jcell2edge(icell+1)-1
         i = i+1
         L = icell2edge(j) ! netlink

         k1 = edge2node(1,L)
         k2 = edge2node(2,L)

         if ( jsferic.eq.0 ) then
            cs = dnx(1,L)
            sn = dny(1,L)

            ! add to system
            Amat(i,1) = cs
            Amat(i,2) = sn
            Amat(i,3) = (0.5d0*(xnode(k1)+xnode(k2))-xzwcell(icell))*cs + (0.5d0*(ynode(k1)+ynode(k2))-yzwcell(icell))*sn
         else
            i12 = 1
            isign = 1
            if ( edge2cell(2,L).eq.icell ) then
               i12 = 2
               isign = -1d0
            end if

            Amat(i,1) = dnx(i12,L)*isign
            Amat(i,2) = dny(i12,L)*isign
            Amat(i,3) = dnz(i12,L)*isign
            Amat(i,4) = ( (0.5d0*(xnode(k1)+xnode(k2))-xzwcell(icell))*dnx(i12,L) + &
               (0.5d0*(ynode(k1)+ynode(k2))-yzwcell(icell))*dny(i12,L) + &
               (0.5d0*(znode(k1)+znode(k2))-zzwcell(icell))*dnz(i12,L) ) * isign
         end if

         jj = jj+1
         ireconst(jj) = L
      end do

      if ( jsferic.eq.0 ) then
         ! solve system
         call gaussj(Amat,3,N,rhs,1,1)

         do i=1,3
            L = icell2edge(jcell2edge(icell)+i-1)
            Areconst(1:3,jreconst(icell)+i-1) = Amat(1:3,i)
         end do
      else
         ! impose zero cell normal velocity
         Amat(4,1:3) = dnn(:,icell)
         Amat(4,4) = 0d0

         ! solve system
         call gaussj(Amat,4,N,rhs,1,1)

         do i=1,3
            L = icell2edge(jcell2edge(icell)+i-1)
            Areconst(1:4,jreconst(icell)+i-1) = Amat(1:4,i)
         end do

      end if
   end do

   ierror = 0
!  error handling
1234 continue

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine reconst_vel_coeffs


   !> reconstruct velocity in cells
   subroutine reconst_vel(q, h0, h1, ierror)
   use m_flowgeom, only: Ndx, Lnx, bl
   use m_flowparameters, only: epshs
   use m_partrecons
   use m_partmesh
   use m_sferic
   use geometry_module, only: dbdistance
   use timers
   
   implicit none

   double precision, dimension(Lnx), intent(in)  :: q    ! flowlink-based discharge (m3/s)
   double precision, dimension(Ndx), intent(in)  :: h0   ! flownode-based water level (m) at begin of interval
   double precision, dimension(Ndx), intent(in)  :: h1   ! flownode-based water level (m) at end of interval

   integer,                          intent(out) :: ierror

   integer,                          parameter   :: N = 4

   integer                                       :: icell, j, k, L

   double precision                              :: hk0, hk1, h, un

   double precision, parameter                   :: DTOL=1d-10

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "reconst_vel", ithndl )

   ierror = 1

   ! get fluxes at all edges, including internal
   call set_fluxes(Lnx, q, qe)

   ! initialize
   u0x = 0d0
   u0y = 0d0
   if ( jsferic.eq.1 ) then
      u0z = 0d0
   end if
   alpha = 0d0

   do icell=1,numcells
      ! get flownode number (for s, bl)
      k = iabs(cell2nod(icell))

      ! get water depth
      hk0 = h0(k) !s0(k)-bl(k)
      hk1 = h1(k) !s1(k)-bl(k)
      if ( abs(hk1-hk0).gt.DTOL ) then
         if ( hk0.gt.epshs .and. hk1.gt.epshs ) then
            h = (hk1-hk0)/(log(hk1)-log(hk0))
         else if ( hk0.gt.epshs .or. hk1.gt.epshs ) then
            h = 0.5d0*(hk0+hk1)
         else
            h = 0d0
         end if
      else
         h = 0.5d0*(hk0+hk1)
      end if

      if ( h.le.epshs ) cycle

      if ( jsferic.eq.0 ) then
         do j=jreconst(icell),jreconst(icell+1)-1
            L = ireconst(j)
            un = qe(L)/(h*w(L))
            u0x(icell)   = u0x(icell)   + Areconst(1,j)*un
            u0y(icell)   = u0y(icell)   + Areconst(2,j)*un
            alpha(icell) = alpha(icell) + Areconst(3,j)*un
         end do
      else
         do j=jreconst(icell),jreconst(icell+1)-1
            L = ireconst(j)
            un = qe(L)/(h*w(L))
            u0x(icell)   = u0x(icell)   + Areconst(1,j)*un
            u0y(icell)   = u0y(icell)   + Areconst(2,j)*un
            u0z(icell)   = u0z(icell)   + Areconst(3,j)*un
            alpha(icell) = alpha(icell) + Areconst(4,j)*un
         end do
      end if
   end do

   ierror = 0
!  error handling
1234 continue

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine reconst_vel

   !> set pointer
   subroutine part_fill_networkdata(hyd, waqgeom,openbndsect_coll)
   use hydmod
   use io_ugrid
   use network_data, only: kn, xk, yk, zk, xzw, yzw, numk, numL, nump, netcell, lnn, lne
   use m_flowgeom, only: Ndxi, Ndx, Lnx, ba, bl, lne2ln
   use m_flow
   use m_transport
   use m_alloc
   use m_missing
   use m_partmesh
   use m_sferic, only: jsferic, jasfer3D
   use geometry_module, only: dbdistance, sphertocart3D, normaloutchk, comp_masscenter
   use timers

   implicit none

   type(t_hyd)        ,intent(in) :: hyd           !< description of the hydrodynamics
   type(t_ug_meshgeom),intent(in) :: waqgeom
   type(t_openbndsect_coll),intent(in) :: openbndsect_coll       ! collection of openbndsects

   integer :: ierr, i, k, L, ln, N, Nmax
   real(8) :: bottomlevel

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "part_fill_networkdata", ithndl )

   if (hyd%crs%epsg_code == 4326) then ! epsg 4326 is assumed spherical, WGS84 system.
      jsferic = 1
   else
      jsferic = 0
   endif

   numk = waqgeom%numnode
   call realloc(xk, numk, keepExisting=.false., fill = dmiss)
   call realloc(yk, numk, keepExisting=.false., fill = dmiss)
   call realloc(zk, numk, keepExisting=.false., fill = 0.0d0)
   do k = 1, numk
      xk(k) = waqgeom%nodex(k)
      yk(k) = waqgeom%nodey(k)
      if (waqgeom%nodez(k) .ne. dmiss) then
         zk(k) = waqgeom%nodez(k)
      endif
   enddo

   numl = waqgeom%numedge
   call realloc(kn, [3, numl], keepExisting=.false., fill = 0)
   call realloc(lne2ln, numl, keepExisting=.false., fill = 0)
   call realloc(lnn, numl, keepExisting=.false., fill = 0)
   call realloc(lne, [2, numl], keepExisting=.false., fill = 0)
   Lnx = 0
   do L = 1, numl
      kn(1,L) = waqgeom%edge_nodes(1,L)
      kn(2,L) = waqgeom%edge_nodes(2,L)
      if (hyd%edge_type(L) == 1 .or. hyd%edge_type(L) == 2 .or. hyd%edge_type(L) == 21 .or. hyd%edge_type(L) == 22) then
         Lnx = Lnx + 1
         lne2ln(L) = Lnx
         lnn = 2
      else
         lne2ln(L) = 0
         lnn = 1
      endif
      lne(1, L) = max(waqgeom%edge_faces(1,L),0)
      lne(2, L) = max(waqgeom%edge_faces(2,L),0)
   enddo
   call realloc(q0, lnx, keepExisting=.false., fill=0d0)
   call realloc(q1, lnx, keepExisting=.false., fill=0d0)

   nump = waqgeom%numface
   Ndx = waqgeom%numface
   call realloc(xzw, nump, keepExisting=.false., fill = dmiss)
   call realloc(yzw, nump, keepExisting=.false., fill = dmiss)
   call realloc(ba, Ndx, keepExisting=.false., fill = 0.0d0)
   call realloc(bl, Ndx, keepExisting=.false., fill = 0.0d0)
   call realloc(h0, Ndx, keepExisting=.false., fill=0d0)
   call realloc(h1, Ndx, keepExisting=.false., fill=0d0)
   call realloc(vol0, Ndx, keepExisting=.false., fill=0d0)
   call realloc(vol1, Ndx, keepExisting=.false., fill=0d0)

   allocate(netcell(nump), stat = ierr)
   Nmax = waqgeom%maxnumfacenodes
   do k=1,nump
      N = 0
      bottomlevel = 0.0
      do i = 1, Nmax
         if (waqgeom%face_nodes(i, k) > 0) then
            N = i
            bottomlevel = bottomlevel + zk(waqgeom%face_nodes(i, k))
         else
            exit
         endif
      enddo
      netcell(k)%N = N
      allocate(netcell(k)%nod(N))
      netcell(k)%nod = waqgeom%face_nodes(1:N, k)
      allocate(netcell(k)%lin(N))
      netcell(k)%lin = waqgeom%face_edges(1:N, k)
!      allocate(netcell(k)%zlin(N))

      xzw(k) = waqgeom%facex(k)
      yzw(k) = waqgeom%facey(k)
      ba(k) = hyd%surf(k)
      bl(k) = bottomlevel / N
   end do

   if ( timon ) call timstop ( ithndl )

   end subroutine

   !> set pointer
   subroutine part_setmesh()
   use network_data, only: kn, xk, yk, xzw, yzw, numk, numL, nump, netcell, lnn, lne
   use m_flowgeom, only: lne2ln, ba
   use m_alloc
   use m_missing
   use m_partmesh
   use m_sferic, only: jsferic, jasfer3D
   use geometry_module, only: dbdistance, sphertocart3D, normaloutchk, comp_masscenter
   use timers

   implicit none

   integer                        :: i, node1, node2, icell, j, k, L, N
   integer                        :: im1, ip1, L1, L2
   integer                        :: newnode, newedge
   integer                        :: isign, ja, k1, k2, k3, kL, kR
   integer                        :: jaswap

   integer                        :: numnontris
   integer                        :: numaddedges

   double precision, dimension(3) :: xv, yv
   double precision, dimension(3) :: t, t1, t2    ! edge tangential vector

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "part_setmesh", ithndl )

   numnodes = numk
   numedges = numL
   numcells = nump

   ! count number of non-triangles and additional (internal) edges
   numnontris=0
   numaddedges=0
   do k=1,nump
      N = netcell(k)%N
      if ( N.gt.3 ) then
         numnontris = numnontris+1
         numaddedges = numaddedges+N
      end if
   end do

   ! compute sizes
   numnodes = numk + numnontris                 ! add centers of non-triangles
   numedges = numL + numaddedges                ! add internal edges of non-triangles
   numcells = nump - numnontris + numaddedges   ! add internal triangles of non-triangles

   ! (re)allocate
   call realloc_partmesh()

   ! nodes
   if ( jsferic.eq.0 ) then
      do k=1,numk
         xnode(k) = xk(k)
         ynode(k) = yk(k)
      end do
   else
      do k=1,numk
         call sphertocart3D(xk(k),yk(k),xnode(k),ynode(k),znode(k))
      end do
   end if

   ! edges
   do L=1,numL
      edge2node(1,L) = kn(1,L)
      edge2node(2,L) = kn(2,L)
      if ( lne2ln(L).ge.0 ) then
         edge2link(L) = lne2ln(L)
      else
         ! not a flowlink number, see flow_geominit
         continue
      end if
   end do

   ! set jcell2edge startpointers and add new triangles
   icell=0
   jcell2edge(1) = 1
   do k=1,nump
      N = netcell(k)%N
      if ( N.eq.3 ) then
         icell = icell+1
         jcell2edge(icell+1) = jcell2edge(icell) + 3
!         nod2cell(k) = icell
         cell2nod(icell) = k
      else if ( N.gt.3 ) then
         do j=1,N
            icell = icell+1
            jcell2edge(icell+1) = jcell2edge(icell) + 3
!            if ( j.eq.1 ) nod2cell(k) = -icell
            cell2nod(icell) = -k
         end do
      else  ! should not happen
         continue
      end if
   end do

   ! allocate jcell2edge
   N = jcell2edge(numcells+1)-1
   call realloc(icell2edge, N, fill=0, keepExisting=.false.)

   ! copy netcell data into cell2edge and add new triangles (nodes, edges)
   icell = 0
   newnode = numk
   numorigedges = numL
   newedge = numorigedges
   do k=1,nump
      N = netcell(k)%N
      if ( N.eq.3 )then
         icell = icell+1   ! add existing triangle
         i = 0
         do j=jcell2edge(icell),jcell2edge(icell+1)-1
            i = i+1
            L = netcell(k)%lin(i)
            icell2edge(j) = L
         end do
         if ( jsferic.eq.0 ) then
            xzwcell(icell) = xzw(k)
            yzwcell(icell) = yzw(k)
         else
            call sphertocart3D(xzw(k),yzw(k),xzwcell(icell),yzwcell(icell),zzwcell(icell))
         end if
         areacell(icell) = ba(k)
      else
         ! add node
         newnode = newnode+1
         if ( jsferic.eq.0 ) then
            xnode(newnode) = xzw(k)
            ynode(newnode) = yzw(k)
         else
            call sphertocart3D(xzw(k),yzw(k),xnode(newnode),ynode(newnode),znode(newnode))
         end if

         xv(1) = xzw(k)
         yv(1) = yzw(k)

         do i=1,N ! add new edges and triangles
            im1 = i-1; if ( im1.lt.1 ) im1=im1+N
            ip1 = i+1; if ( ip1.gt.N ) ip1=ip1-N
            L = netcell(k)%lin(i)
            call find_common_node(netcell(k)%lin(im1),L,node1)
            call find_common_node(netcell(k)%lin(ip1),L,node2)

            icell = icell+1
            ! add edges
            newedge = newedge+1
            edge2node(1,newedge) = newnode
            edge2node(2,newedge) = node1
            edge2link(newedge) = -1

            ! add edges to cell2edge
            j = jcell2edge(icell)
            icell2edge(j) = newedge
            icell2edge(j+1) = L
            icell2edge(j+2) = newedge+1
            if ( i.eq.N ) then
               icell2edge(j+2) = icell2edge(j+2) - N
            end if

            ! compute mass center
            xv(2) = xk(node1)
            yv(2) = yk(node1)
            xv(3) = xk(node2)
            yv(3) = yk(node2)
            call comp_masscenter(3, xv, yv, xzwcell(icell), yzwcell(icell), areacell(icell), ja, jsferic, jasfer3D, dmiss)

            if ( jsferic.eq.1 ) then
               xv(1) = xzwcell(icell)  ! reuse xv
               yv(1) = yzwcell(icell)   ! reuse yv
               call sphertocart3D(xv(1), yv(1),xzwcell(icell),yzwcell(icell),zzwcell(icell))
            end if
         end do
      end if
   end do

   ! set edge2cell
   edge2cell = 0
   do icell=1,numcells
      do j=jcell2edge(icell),jcell2edge(icell+1)-1
         L = icell2edge(j)
         if ( edge2cell(1,L).eq.0 ) then
            edge2cell(1,L) = icell
         else
            edge2cell(2,L) = icell
         end if
      end do
   end do

   !  fix orientation of edge2cell
   do L=1,numL ! exclude new edges
      k1 = edge2cell(1,L)
      if ( k1.eq.0 ) cycle

      k1 = iabs(cell2nod(k1))

      jaswap = 0

      if ( lnn(L).eq.1 ) then  !     boundaries: inward normal
         jaswap = 1
      else if ( k1.eq.lne(2,L) .and. lnn(L).gt.1 ) then
         jaswap = 1
      end if

      if ( jaswap.eq.1 ) then
         icell = edge2cell(1,L)
         edge2cell(1,L) = edge2cell(2,L)
         edge2cell(2,L) = icell
      end if
   end do

   if ( jsferic.ne.0 ) then
      ! compute cell normal vectors from first two edges
      do k=1,numcells
         L1 = icell2edge(jcell2edge(k))
         L2 = icell2edge(jcell2edge(k)+1)
         k1 = edge2node(1,L1)
         k2 = edge2node(2,L1)
         k3 = edge2node(1,L2)
         if ( k3.eq.k1 .or. k3.eq.k2 ) then
            k3 = edge2node(2,L2)
         end if
         t1 = (/ xnode(k2)-xnode(k1), ynode(k2)-ynode(k1), znode(k2)-znode(k1)/)
         t2 = (/ xnode(k3)-xnode(k2), ynode(k3)-ynode(k2), znode(k3)-znode(k2)/)

         dnn(:,k) = (/ t1(2)*t2(3) - t1(3)*t2(2), t1(3)*t2(1) - t1(1)*t2(3), t1(1)*t2(2) - t1(2)*t2(1) /)
         dnn(:,k) = dnn(:,k) / sqrt( dnn(1,k)**2 + dnn(2,k)**2 + dnn(3,k)**2 )

         ! fix orientation
         if ( dnn(1,k)*xnode(k1) + dnn(2,k)*ynode(k1) + dnn(3,k)*znode(k1) .lt. 0d0 ) then
            dnn(:,k) = -dnn(:,k)
         end if
      end do
   end if

   ! nx, ny, w
   do L=1,numedges
      k1 = edge2node(1,L)
      k2 = edge2node(2,L)

      kL = edge2cell(1,L)
      kR = edge2cell(2,L)

      k = kL   ! outward positive

      isign = 1
      if ( k.le.0 ) then
         k = edge2cell(2,L)   ! inward positive
         isign = -1
      end if

      if ( k.eq.0 ) cycle  ! isolated edge

      ! compute normal vector (outward positive)
      if ( jsferic.eq.0 ) then
         call normaloutchk(xnode(k1),ynode(k1),xnode(k2),ynode(k2),xzwcell(k),yzwcell(k),dnx(1,L),dny(1,L),ja, jsferic, jasfer3D, dmiss, dxymis)

         if ( isign.eq.-1 ) then
            dnx(1,L) = -dnx(1,L)
            dny(1,L) = -dny(1,L)
         end if

         w(L) = dbdistance(xnode(k1),ynode(k1),xnode(k2),ynode(k2),jsferic, jasfer3D, dmiss)
      else
         ! compute outward normal with respect to the left cell
         t = (/ xnode(k2)-xnode(k1), ynode(k2)-ynode(k1), znode(k2)-znode(k1) /)
         w(L) = sqrt( t(1)**2 + t(2)**2 + t(3)**2)

         t = t/w(L)

         if ( kL.gt.0 ) then
            ! left cell normal : nn X t
            dnx(1,L) = dnn(2,kL) * t(3) - dnn(3,kL) * t(2)
            dny(1,L) = dnn(3,kL) * t(1) - dnn(1,kL) * t(3)
            dnz(1,L) = dnn(1,kL) * t(2) - dnn(2,kL) * t(1)

            ! fix orientation
            if ( (xnode(k1)-xzwcell(kL))*dnx(1,L) + (ynode(k1)-yzwcell(kL))*dny(1,L) + (znode(k1)-zzwcell(kL))*dnz(1,L) .lt. 0d0 ) then
               dnx(1,L) = -dnx(1,L)
               dny(1,L) = -dny(1,L)
               dnz(1,L) = -dnz(1,L)
            end if
         end if

         if ( kR.gt.0 ) then
            ! left cell normal : nn X t
            dnx(2,L) = dnn(2,kR) * t(3) - dnn(3,kR) * t(2)
            dny(2,L) = dnn(3,kR) * t(1) - dnn(1,kR) * t(3)
            dnz(2,L) = dnn(1,kR) * t(2) - dnn(2,kR) * t(1)

            ! fix orientation
            if ( (xnode(k1)-xzwcell(kR))*dnx(2,L) + (ynode(k1)-yzwcell(kR))*dny(2,L) + (znode(k1)-zzwcell(kR))*dnz(2,L) .lt. 0d0 ) then
               dnx(2,L) = -dnx(2,L)
               dny(2,L) = -dny(2,L)
               dnz(2,L) = -dnz(2,L)
            end if
         end if
      end if
   end do

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine part_setmesh

   !> set all fluxes, including internal
   subroutine set_fluxes(Lnx,q,qe)
   use m_partmesh
   use m_partfluxes
   use timers

   implicit none

   integer,                               intent(in)  :: Lnx
   double precision, dimension(Lnx),      intent(in)  :: q

   double precision, dimension(numedges), intent(out) :: qe

   integer                                            :: j, L, Lf

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "set_fluxes", ithndl )

   do L=1,numedges
      qe(L) = 0d0
      do j=jflux2link(L),jflux2link(L+1)-1
         Lf = iflux2link(j)
         if ( Lf.gt.0 ) then
            qe(L) = qe(L) + Aflux2link(j)*q(Lf)
         end if
      end do
   end do

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine set_fluxes


   !> compute mapping from prescribed (at flowlinks) to all fluxes (at all edges, including "internal")
   subroutine comp_fluxcoeffs()
   use m_partmesh
   use m_partfluxes
   use m_alloc
   use MessageHandling
   use geometry_module, only: gaussj
   use timers

   implicit none

   integer                                  :: N         ! number of subtriangles

   integer,          dimension(MAXSUBCELLS) :: L1        ! "internal" edges
   integer,          dimension(MAXSUBCELLS) :: L2        ! original "non-internal" edges
   integer,          dimension(MAXSUBCELLS) :: icell     ! subcells
   integer,          dimension(MAXSUBCELLS) :: isign1    ! orientation of "internal: edges L1, all cw or
   integer,          dimension(MAXSUBCELLS) :: isign2    ! orientation of "non-internal: edges L2, outward normal
   double precision, dimension(MAXSUBCELLS) :: circ      ! weight of edge L1 in computation of circulation

   double precision, dimension(MAXSUBCELLS,MAXSUBCELLS) :: A, B

   double precision                         :: areasum, dlen

   integer                                  :: i, im1, ip1, j, j2, k, L, L3, Lf

   integer, external                        :: icommonval

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "comp_fluxcoeffs", ithndl )

   ! allocate
   call realloc_partfluxes()
   ! call realloc(jflux2link, numedges+1, fill=0, keepExisting=.false.)

   ! set startpointers
   jflux2link(1) = 1
   L = 1
   do while ( L.le.numedges )
      Lf = edge2link(L)
      if ( Lf.gt.0 ) then
         ! original flowlink
         jflux2link(L+1) = jflux2link(L) + 1

         ! proceed to next edge
         L = L+1
      else if ( L.gt.numorigedges ) then
         ! internal edge: find number of subtriangels
         N = 0
         do while ( icommonval(edge2cell(:,L+N), edge2cell(:,L+N+1)).ne.0 )
            N = N+1
            if ( L+N.ge.numedges ) exit
         end do

         if ( N.gt.0 ) then
            ! check connection first and last subtriangle
            if ( icommonval(edge2cell(:,L), edge2cell(:,L+N)).ne.0 ) then
               N = N+1

               do i=1,N
                  jflux2link(L+1) = jflux2link(L) + N
                  ! proceed to next edge
                  L = L+1
               end do
            else
               call mess(LEVEL_ERROR, 'comp_fluxcoeffs: numbering error')
            end if

         else  ! should not happen
            call mess(LEVEL_ERROR, 'comp_fluxcoeffs: numbering error')
         end if

         if ( L.ge.numedges ) exit
      else  ! other
         jflux2link(L+1) = jflux2link(L)
         ! proceed to next edge
         L = L+1
      end if
   end do

   ! reallocate
   call realloc(iflux2link, jflux2link(numedges+1)-1, fill=0, keepExisting=.false.)
   call realloc(Aflux2link, jflux2link(numedges+1)-1, fill=0d0, keepExisting=.false.)

   ! fill iflux2link (first in edge-numbers) and compute coefficients
   L = 1
   do while ( L.le.numedges )
      j=jflux2link(L)
      N = jflux2link(L+1)-j

      if ( N.eq.1 ) then
         ! original "non-internal" link
         iflux2link(j) = L
         Aflux2link(j) = 1d0
         ! proceed to next edge
         L = L+1
      else if ( N.gt.1 ) then
         ! "internal" link
         do i=1,N
            ip1 = i+1; if ( ip1.gt.N ) ip1 = ip1-N
            L1(i) = L+i-1
            L3 = L1(i)+ip1

            ! find subcell
            icell(i) = icommonval(edge2cell(:,L+i-1), edge2cell(:,L+ip1-1))

            ! find original edge
            j2 = jcell2edge(icell(i))+1   ! order of cell edges: "left" internal (L1), orginal (L2), "right internal (L3)
            L2(i) = icell2edge(j2)

            ! get orientation of "internal" edges (i-dir positive)
            ! note: will always be (/ -1, 1, 1, ... /), due to generation of edge2cells in part_setmesh
            isign1(i) = 1
            if ( edge2cell(1,L1(i)).eq.icell(i) ) isign1(i) = -1

            ! get orientation of original "non-internal" edges (outward normal)
            isign2(i) = 1
            if ( edge2cell(2,L2(i)).eq.icell(i) ) isign2(i) = -1
         end do

         ! compute summed area and circulation weights
         areasum = 0d0
         do i=1,N
            im1 = i-1; if ( im1.lt.1 ) im1=im1+N
            areasum = areasum + areacell(icell(i))
            dlen = 0.25*(areacell(icell(im1))+areacell(icell(i)))/w(L1(i))
            circ(i) = dlen/w(L1(i))
         end do

         ! build system: A qe = B qlink
         A = 0d0
         B = 0d0
         ! continuity equations
         do i=1,N-1
            A(i,i)   = -isign1(i)   ! outward
            A(i,i+1) =  isign1(i+1) ! outward
            B(i,i)   = -isign2(i)   ! outward, rhs
            do k=1,N
               B(i,k) = B(i,k) + areacell(icell(i))*isign2(k)/areasum
            end do
         end do
         ! circulation
         do k=1,N
            A(N,k) = circ(k) * isign1(k)
         end do

         ! invert matrix: qe = inv(A) B qlink
         call gaussj(A,N,MAXSUBCELLS,B,N,MAXSUBCELLS)

         ! fill data
         do i=1,N
            iflux2link(j:j+N-1) = L2(1:N)
            Aflux2link(j:j+N-1) = B(i,1:N)

            ! proceed to next edge
            j = j+N
            L = L+1
         end do
      else
         ! closed boundary: proceed to next edge
         L = L+1
      end if
   end do

   ! convert to link number
   do j=1,jflux2link(numedges+1)-1
      L  = iflux2link(j)
      Lf = edge2link(L)
      if ( Lf.gt.0 ) then
         iflux2link(j) = Lf
      else  ! closed boundary
         iflux2link(j) = 0
      end if
   end do

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine comp_fluxcoeffs

   !> find common value of two-dimensional arrays i1 and i2
   !> it is assumed there is at most one common value
   integer function icommonval(i1, i2)
   implicit none

   integer, dimension(2), intent(in)  :: i1
   integer, dimension(2), intent(in)  :: i2

   icommonval = 0
   if ( i1(1).eq.i2(1) .or. i1(1).eq.i2(2) ) then
      icommonval = i1(1)
   else if ( i1(2).eq.i2(1) .or. i1(2).eq.i2(2) ) then
      icommonval = i1(2)
   end if

   return
   end function icommonval

   !> add particles
   subroutine add_particles(Nadd, xadd, yadd)
   use m_particles
   use m_partmesh
   use m_alloc
   use m_sferic, only: jsferic
   use geometry_module, only: sphertocart3D
   use timers

   implicit none

   integer,                           intent(in)  :: Nadd       !< number of particles to be added
   double precision, dimension(Nadd), intent(in)  :: xadd       !< x-coordinates of particles to be added
   double precision, dimension(Nadd), intent(in)  :: yadd       !< y-coordinates of particles to be added

   integer,          dimension(:),    allocatable :: kadd       !< cell numbers

   integer                                        :: i, ipoint, Nsize
   integer                                        :: ierror
   integer                                        :: Npartnew
   integer                                        :: Nreplace
   integer                                        :: Nloc

   double precision                               :: xn, yn, zn, dn
   integer                                        :: k, k1

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "add_particles", ithndl )

   ! get new particle cell number
   allocate(kadd(Nadd))
   kadd = 0
   call part_findcell(Nadd,xadd,yadd,kadd,ierror)

   do i=1,Nadd
      if ( kadd(i).eq.0 ) cycle
      Npart = Npart + 1
      if ( jsferic.eq.0 ) then
         xpart(Npart) = xadd(i)
         ypart(Npart) = yadd(i)
      else
         call sphertocart3D(xadd(i),yadd(i),xpart(Npart),ypart(Npart),zpart(Npart))
         if ( jsferic.eq.1 ) then
            ! project particle on triangle
            k = kadd(i)
            if ( k.gt.0 ) then
               k1 = edge2node(1,icell2edge(jcell2edge(k)))
               xn = xnode(k1)
               yn = ynode(k1)
               zn = znode(k1)
               dn = (xpart(Npart) - xn) * dnn(1,k) +  &
                  (ypart(Npart) - yn) * dnn(2,k) +  &
                  (zpart(Npart) - zn) * dnn(3,k)
               xpart(Npart) = xpart(Npart) - dn * dnn(1,k)
               ypart(Npart) = ypart(Npart) - dn * dnn(2,k)
               zpart(Npart) = zpart(Npart) - dn * dnn(3,k)
            end if
         end if

      end if
      kpart(Npart) = kadd(i)
   end do

   !  deallocate
   if ( allocated(kadd) ) deallocate(kadd)

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine add_particles



   !> find in which cells particles are located
   subroutine part_findcell(Npart, xxpart, yypart, kpart, ierror)
   use m_partmesh
   use MessageHandling
   use kdtree2Factory
   use m_sferic, only: jsferic, jasfer3D
   use m_missing, only: jins, dmiss
   use geometry_module, only: pinpok, dbdistance, pinpok3D, cart3Dtospher
   use timers

   implicit none

   integer,                            intent(in)  :: Npart    !< number of particles
   double precision, dimension(Npart), intent(in)  :: xxpart   !< particle x-coordinates, 2D Cartexsion or spherical coordinates (not 3D Cartesian)
   double precision, dimension(Npart), intent(in)  :: yypart   !< particle x-coordinates, 2D Cartexsion or spherical coordinates (not 3D Cartesian)
   integer,          dimension(Npart), intent(out) :: kpart    !< cell numbers

   integer                           , intent(out) :: ierror   !< error (1) or not (0)

   type(kdtree_instance)                           :: kdtree

   double precision, dimension(3)                  :: xv, yv

   double precision                                :: dmaxsize
   double precision                                :: R2search
   double precision                                :: xx, yy

   integer                                         :: i, ip1, j, k, knode, L, Lp1, N, NN
   integer                                         :: inside

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "part_findcell", ithndl )

   ierror = 1

   !  build kdtree
   call build_kdtree(kdtree, Npart, xxpart, yypart, ierror, jsferic, dmiss)
   if ( ierror.ne.0 ) then
      goto 1234
   end if

   if ( Npart.lt.1 ) then
      ierror = 0
      goto 1234
   end if

   kpart = 0

   !  loop over cells
   do k=1,numcells
      !     check cell size
      N = jcell2edge(k+1)-jcell2edge(k)
      if ( N.ne.3 ) then
         call mess(LEVEL_ERROR, 'part_findcell: non-triangle')
         goto 1234
      end if

      !     get cell polygon
      i=0
      do j = jcell2edge(k),jcell2edge(k+1)-1
         i = i+1
         L = icell2edge(j)
         ip1 = i+1; if ( ip1.gt.3 ) ip1=ip1-3
         Lp1 = icell2edge(j-i+ip1)
         !        find common node of L and Lp1
         if ( edge2node(1,L).eq.edge2node(1,Lp1) .or. edge2node(1,L).eq.edge2node(2,Lp1) ) then
            knode = edge2node(1,L)
         else if ( edge2node(2,L).eq.edge2node(1,Lp1) .or. edge2node(2,L).eq.edge2node(2,Lp1) ) then
            knode = edge2node(2,L)
         else  ! should not happen
            continue
         end if
         if ( jsferic.eq.0 ) then
            xv(i) = xnode(knode)
            yv(i) = ynode(knode)
         else
            call Cart3Dtospher(xnode(knode),ynode(knode),znode(knode),xv(i),yv(i),0d0)
         end if
      end do

      !     fill query vector
      if ( jsferic.eq.0 ) then
         call make_queryvector_kdtree(kdtree,xzwcell(k),yzwcell(k), jsferic)
      else
         call cart3Dtospher(xzwcell(k),yzwcell(k),zzwcell(k),xx,yy,0d0)
         call make_queryvector_kdtree(kdtree,xx,yy, jsferic)
      end if

      !     compute maximum flowcell dimension
      dmaxsize = 0d0
      do i=1,N
         ip1=i+1; if ( ip1.gt.N ) ip1=ip1-N
         dmaxsize = max(dmaxsize, dbdistance(xv(i),yv(i),xv(ip1),yv(ip1),jsferic, jasfer3D, dmiss))
      end do

      !     determine square search radius
      R2search = 1.1d0*dmaxsize**2  ! 1.1d0: safety

      !     count number of points in search area
      NN = kdtree2_r_count(kdtree%tree,kdtree%qv,R2search)

      if ( NN.eq.0 ) cycle ! no particles found

      !     reallocate if necessary
      call realloc_results_kdtree(kdtree,NN)

      !     find nearest NN samples
      call kdtree2_n_nearest(kdtree%tree,kdtree%qv,NN,kdtree%results)

      !     check if samples are in cell
      do i=1,NN
         j = kdtree%results(i)%idx
         if ( jsferic.eq.0 ) then
            call pinpok(xxpart(j), yypart(j), 3, xv, yv, inside, jins, dmiss)
         else
            call pinpok3D(xxpart(j), yypart(j), 3, xv, yv, inside, dmiss, jins, jsferic, jasfer3D)
         end if

         if ( inside.eq.1 ) then
            if ( kpart(j).eq.0 ) then
               kpart(j) = k
            end if
         end if
      end do
   end do

   ierror = 0

1234 continue

   !  deallocate
   if ( kdtree%itreestat.ne.ITREE_EMPTY ) call delete_kdtree2(kdtree)

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine part_findcell

   !> (re)allocate
   subroutine realloc_particles(Nsize, LkeepExisting, ierror)
   use m_particles
   use m_alloc
   use m_missing
   use m_sferic, only: jsferic
   implicit none

   integer, intent(in)  :: Nsize          !< array sizes
   logical, intent(in)  :: LkeepExisting  !< keep existing data (1) or not (0)
   integer, intent(out) :: ierror         !< error (1) or not

   ierror = 1

   !  reallocate
   call realloc(xpart, Nsize, keepExisting=LkeepExisting, fill=DMISS)
   call realloc(ypart, Nsize, keepExisting=LkeepExisting, fill=DMISS)
   if ( jsferic.eq.1 ) then
      call realloc(zpart, Nsize, keepExisting=LkeepExisting, fill=DMISS)
   end if
   call realloc(dtremaining, Nsize, keepExisting=LkeepExisting, fill=0d0)
   call realloc(kpart, Nsize, keepExisting=LkeepExisting, fill=0)

   call realloc(numzero, Nsize, keepExisting=LkeepExisting, fill=0)
   numzero = 0

   ierror = 0
1234 continue

   return
   end subroutine realloc_particles

   !> deallocate particle data
   subroutine dealloc_particles()
   use m_particles
   implicit none

   if ( allocated(xpart)       ) deallocate(xpart)
   if ( allocated(ypart)       ) deallocate(ypart)
   if ( allocated(zpart)       ) deallocate(zpart)
   if ( allocated(dtremaining) ) deallocate(dtremaining)
   if ( allocated(kpart)       ) deallocate(kpart)

   if ( allocated(numzero)     ) deallocate(numzero)

   Npart = 0

   if ( allocated(trpart)       ) deallocate(trpart)
   if ( allocated(xrpart)       ) deallocate(xrpart)
   if ( allocated(yrpart)       ) deallocate(yrpart)
   if ( allocated(zrpart)       ) deallocate(zrpart)
   Nrpart = 0

   return
   end subroutine dealloc_particles


   !> (re)allocate partmesh data
   subroutine realloc_partmesh()
   use m_partmesh
   use m_alloc
   use m_missing
   use m_sferic, only: jsferic
   implicit none

   integer :: N

   call realloc(edge2node, (/2, numedges/), fill=0, keepExisting=.false.)
   call realloc(edge2cell, (/2, numedges/), fill=0, keepExisting=.false.)
   call realloc(xnode, numnodes, fill=0d0, keepExisting=.false.)
   call realloc(ynode, numnodes, fill=0d0, keepExisting=.false.)
   if ( jsferic.eq.1 ) then
      call realloc(znode, numnodes, fill=0d0, keepExisting=.false.)
   end if

   call realloc(xzwcell, numcells, fill=DMISS, keepExisting=.false.)
   call realloc(yzwcell, numcells, fill=DMISS, keepExisting=.false.)
   if ( jsferic.eq.1 ) then
      call realloc(zzwcell, numcells, fill=DMISS, keepExisting=.false.)
   end if
   call realloc(areacell, numcells, fill=DMISS, keepExisting=.false.)

   if ( jsferic.eq.0 ) then
      call realloc(dnx, (/1, numedges/), fill=DMISS, keepExisting=.false.)
      call realloc(dny, (/1, numedges/), fill=DMISS, keepExisting=.false.)
   else
      call realloc(dnx, (/2, numedges/), fill=DMISS, keepExisting=.false.)
      call realloc(dny, (/2, numedges/), fill=DMISS, keepExisting=.false.)
      call realloc(dnz, (/2, numedges/), fill=DMISS, keepExisting=.false.)
      call realloc(dnn, (/3, numcells/), fill=DMISS, keepExisting=.false.)
   end if
   call realloc(w, numedges, fill=DMISS, keepExisting=.false.)

   call realloc(edge2link, numedges, fill=0, keepExisting=.false.)
   !   call realloc(nod2cell, numcells, fill=0, keepExisting=.false.)
   call realloc(cell2nod, numcells, fill=0, keepExisting=.false.)

   call realloc(jcell2edge, numcells+1, fill=1, keepExisting=.false.)
   N = jcell2edge(numcells+1)-1
   call realloc(icell2edge, N, fill=0, keepExisting=.false.)

   return
   end subroutine realloc_partmesh

   !> deallocate particle mesh data
   subroutine dealloc_partmesh()
   use m_partmesh
   implicit none

   if ( allocated(edge2node ) ) deallocate(edge2node )
   if ( allocated(edge2cell ) ) deallocate(edge2cell )
   if ( allocated(xnode     ) ) deallocate(xnode     )
   if ( allocated(ynode     ) ) deallocate(ynode     )
   if ( allocated(znode     ) ) deallocate(znode     )

   if ( allocated(xzwcell   ) ) deallocate(xzwcell   )
   if ( allocated(yzwcell   ) ) deallocate(yzwcell   )
   if ( allocated(zzwcell   ) ) deallocate(zzwcell   )
   if ( allocated(areacell  ) ) deallocate(areacell  )

   if ( allocated(dnn       ) ) deallocate(dnn       )
   if ( allocated(dnx       ) ) deallocate(dnx       )
   if ( allocated(dny       ) ) deallocate(dny       )
   if ( allocated(dnz       ) ) deallocate(dnz       )
   if ( allocated(w         ) ) deallocate(w         )

   if ( allocated(edge2link ) ) deallocate(edge2link )
   !  if ( allocated(nod2cell  ) ) deallocate(nod2cell  )
   if ( allocated(cell2nod  ) ) deallocate(cell2nod  )

   if ( allocated(icell2edge) ) deallocate(icell2edge)
   if ( allocated(jcell2edge) ) deallocate(jcell2edge)


   return
   end subroutine dealloc_partmesh


   !> (re)allocate flux coefficients
   subroutine realloc_partfluxes()
   use m_partmesh
   use m_partfluxes
   use m_alloc
   use m_missing
   implicit none

   integer :: N

   call realloc(jflux2link, numedges+1, keepExisting=.false., fill=1)
   N = jflux2link(numedges+1)-1
   call realloc(iflux2link, N,  keepExisting=.false., fill=0)
   call realloc(Aflux2link, N,  keepExisting=.false., fill=0d0)

   return
   end subroutine realloc_partfluxes

   !> deallocate flux_coeffs
   subroutine dealloc_partfluxes()
   use m_partfluxes
   implicit none

   if ( allocated(iflux2link) ) deallocate(iflux2link)
   if ( allocated(jflux2link) ) deallocate(jflux2link)
   if ( allocated(Aflux2link) ) deallocate(Aflux2link)
   end subroutine dealloc_partfluxes

   !> (re)allocate flux coefficients et cetera
   subroutine realloc_partrecons()
   use m_partmesh
   use m_partrecons
   use m_alloc
   use m_missing
   use m_sferic, only: jsferic
   implicit none

   call realloc(qe, numedges, keepExisting=.false., fill=DMISS)

   call realloc(u0x, numcells, keepExisting=.false., fill=DMISS)
   call realloc(u0y, numcells, keepExisting=.false., fill=DMISS)
   if ( jsferic.eq.1 ) then
      call realloc(u0z, numcells, keepExisting=.false., fill=DMISS)
   end if
   call realloc(alpha, numcells, keepExisting=.false., fill=DMISS)

   call realloc(ireconst, numcells+1, keepExisting=.false., fill=0)
   return
   end subroutine realloc_partrecons

   !> deallocate flux_coeffs
   subroutine dealloc_partrecons()
   use m_partrecons
   implicit none

   if ( allocated(qe) ) deallocate(qe)

   if ( allocated(u0x) ) deallocate(u0x)
   if ( allocated(u0y) ) deallocate(u0y)
   if ( allocated(u0z) ) deallocate(u0z)
   if ( allocated(alpha) ) deallocate(alpha)

   if ( allocated(ireconst) ) deallocate(ireconst)
   if ( allocated(jreconst) ) deallocate(jreconst)
   if ( allocated(Areconst) ) deallocate(Areconst)

   return
   end subroutine dealloc_partrecons

!  initialise grid and hydro data   
   subroutine ini_part_grid(hyd)
   use hydmod
   use partmem, only: ihdel, layt, nolayp, tcktot, nmaxp, mmaxp, mnmax2, mnmaxk
   use m_flow                                              
   use m_transport, only: constituents, numconst
   use m_missing
   use m_alloc
   use MessageHandling
   use timers

   implicit none

   type(t_hyd)                    :: hyd           !< description of the hydrodynamics

   integer             :: minp, ilay, lun
   logical             :: lexist

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "ini_part_grid", ithndl )

   !  deallocate
   call dealloc_partmesh()
   call dealloc_partfluxes()
   call dealloc_partrecons()
   call dealloc_particles()
   call dealloc_auxfluxes()

   kmx = 0 ! for now 2D only

   !     fill network_data
   call part_fill_networkdata(hyd, hyd%waqgeom, hyd%openbndsect_coll)

   !     set pointers with mesh connectivity etc.
   call part_setmesh()

   !     set flux coeffecients
   call comp_fluxcoeffs()

   call realloc_partrecons()

   call reconst_vel_coeffs()

   call alloc_auxfluxes()

   ihdel = hyd%cnv_step_sec
   layt = hyd%nolay
   nolayp = hyd%nolay
   allocate ( tcktot(hyd%nolay))
   do ilay = 1, hyd%nolay
      tcktot(ilay) = hyd%waq_layers(ilay)
   enddo
   mnmax2 = hyd%nosegl
   mnmaxk = hyd%nosegl*hyd%nolay

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine ini_part_grid

   !> initialize particles
   subroutine ini_part(partfile, partrelfile, starttime_loc, timestep_loc, threeDtype_loc)
   use hydmod
   use m_particles
   use m_samples
   use m_flow
   use m_flowgeom, only: Ndx
   use m_transport, only: constituents, numconst
   use m_flowtimes, only: tstart_user
   use m_missing
   use m_alloc
   use MessageHandling
   use timers

   implicit none

   character(len=255), intent(in) :: partfile      !< initial particle file
   character(len=255), intent(in) :: partrelfile   !< particle release file
   double precision,   intent(in) :: starttime_loc !< start time (>0) or not (0)
   double precision,   intent(in) :: timestep_loc  !< time step (>0) or every computational time step (0)
   integer,            intent(in) :: threeDtype_loc    !< depth averaged (0) or free surface (1)

   integer             :: minp, ierror
   logical             :: lexist

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "ini_part", ithndl )

   Npart = 0
   Nrpart = 0
   irpart = 0
   NpartTot = 0

!  add particle tracer (when tracers are initialized)
   numconst = 1
   part_iconst = 1
   numconst = 1
   call realloc(constituents, (/ numconst, Ndx /), keepExisting=.false., fill=0d0)

   timenext = 0d0
   timelast = DMISS

   !  start time
   if ( starttime_loc.gt.0d0 ) then
      starttime = starttime_loc
      timenext = starttime
   end if

   !  time step
   if ( timestep_loc.gt.0d0 ) then
      timestep = timestep_loc
   end if

!   !  3D type
!   if ( kmx.gt.0 ) then
!      threeDtype = threeDtype_loc
!   else  ! 2D
!      threeDtype = 0
!   end if

   if ( len_trim(partfile).gt.0 ) then
      !     read initial samples from inputfile
      inquire(FILE = trim(partfile), exist = lexist)
      if ( lexist ) then
         call oldfil(minp, partfile)
         call savesam()
         call reasam(minp, 0)
         NpartTot = NpartTot + Ns
      else
         call mess(LEVEL_ERROR, 'the specified initial particle locations file could not be found: ', trim(partfile))
      end if
   end if
   if ( len_trim(partrelfile).gt.0 ) then
      !     read initial samples from inputfile
      inquire(FILE = trim(partrelfile), exist = lexist)
      if ( lexist ) then
         call read_particles_release_file(partrelfile)
      else
         call mess(LEVEL_ERROR, 'the specified particle release file could not be found: ', trim(partfile))
      end if
   end if

   call realloc_particles(NpartTot, .true., ierror)
   if ( Ns.gt.0 ) then
      call add_particles(Ns, xs, ys)
      call delsam(0)
   end if
   timepart = tstart_user

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine ini_part


   !> read particles release file
   subroutine read_particles_release_file(partrelfile)
   use m_particles
   use m_missing
   use m_alloc
   use MessageHandling
   implicit none

   character(len=255), intent(in) :: partrelfile   !< release particle file
   character(len=1000) :: line
   character(len=1)    :: char
   integer             :: lun, ios, ipart, linenr, ierror
   double precision    :: tr, xr, yr, zr

   call oldfil(lun, partrelfile)

   Nrpart = 0
   linenr = 0
   ios = 0

   do while ( ios==0 )
      read(lun, '(a1000)', iostat=ios) line
      linenr = linenr + 1
      if (ios==0) then
         read(line, '(a)', iostat=ios) char
         if (char.ne.'*'.and.char.ne.'#'.and.char.ne.'!') then
            read(line, *, iostat=ios) tr, xr, yr, zr
            if (ios==0) then
               Nrpart = Nrpart + 1
            endif
         endif
      endif
   end do

   if (Nrpart.gt.0) then
      ipart = 0
      linenr = 0
      ios = 0

      rewind (lun)
      call realloc(trpart, Nrpart)
      call realloc(xrpart, Nrpart)
      call realloc(yrpart, Nrpart)
      call realloc(zrpart, Nrpart)
      call realloc(krpart, Nrpart)

      do while ( ios==0 )
         read(lun, '(a1000)', iostat=ios) line
         linenr = linenr + 1
         if (ios==0) then
            read(line, '(a)', iostat=ios) char
            if (char.ne.'*'.and.char.ne.'#'.and.char.ne.'!') then
               ipart = ipart + 1
               read(line, *, iostat=ios) trpart(ipart), xrpart(ipart), yrpart(ipart), zrpart(ipart)
               if (ios.ne.0 .or. trpart(ipart).eq.dmiss .or. xrpart(ipart).eq.dmiss .or.  yrpart(ipart).eq.dmiss .or. zrpart(ipart).eq.dmiss) then
                  call mess(LEVEL_ERROR, 'error reading particle release file '''//trim(partrelfile)//''' at line', linenr)
               endif
               if (ipart.gt.1) then
                  if (trpart(ipart).lt.trpart(ipart-1)) then
                     call mess(LEVEL_ERROR, 'timing in particle release file '''//trim(partrelfile)//''' is not incremental at line', linenr)
                  endif
               endif
            endif
         endif
      end do
      irpart = 1
      call part_findcell(Nrpart,xrpart,yrpart,krpart,ierror)
      do ipart = 1, Nrpart
         if (krpart(ipart).gt.0) then
            NpartTot = NpartTot + 1
         endif
      enddo
   endif
   close(lun)
   end subroutine read_particles_release_file


   !> add released particles
   subroutine add_particles_from_release_file(time0)
   use m_particles
   use m_partmesh
   use m_sferic, only: jsferic
   use geometry_module, only: sphertocart3D
   use timers
   
   double precision, intent(in) :: time0       !< current   julian (s) of h0

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /

   if (irpart.eq.0 .or. irpart.gt.Nrpart) return

   if ( timon ) call timstrt( "add_particles_from_release_file", ithndl )


   do while (irpart.le.Nrpart)
      if (trpart(irpart)*60.0d0.gt.time0) exit
      if(krpart(irpart).gt.0) then
         Npart = Npart + 1
         if ( jsferic.eq.0 ) then
            xpart(Npart) = xrpart(irpart)
            ypart(Npart) = yrpart(irpart)
         else
            call sphertocart3D(xrpart(irpart),yrpart(irpart),xpart(Npart),ypart(Npart),zpart(Npart))
            ! project particle on triangle
            k = krpart(irpart)
            if ( k.gt.0 ) then
               k1 = edge2node(1,icell2edge(jcell2edge(k)))
               xn = xnode(k1)
               yn = ynode(k1)
               zn = znode(k1)
               dn = (xpart(Npart) - xn) * dnn(1,k) +  &
                  (ypart(Npart) - yn) * dnn(2,k) +  &
                  (zpart(Npart) - zn) * dnn(3,k)
               xpart(Npart) = xpart(Npart) - dn * dnn(1,k)
               ypart(Npart) = ypart(Npart) - dn * dnn(2,k)
               zpart(Npart) = zpart(Npart) - dn * dnn(3,k)
            end if
         end if
         kpart(Npart) = krpart(irpart)
      end if
      irpart = irpart + 1
   end do

   if ( timon ) call timstop ( ithndl )

   end subroutine add_particles_from_release_file


   !> compute concentrations of particles (parts per unit volume) in flownodes
   subroutine comp_concentration(h, nconst, iconst, c)
   use m_particles
   use m_partmesh
   use m_flowgeom, only : Ndx, ba, bl
   use m_flowparameters, only: epshs
   use m_flow, only: Ndkx
   use timers

   implicit none

   double precision, dimension(Ndx),        intent(in)  :: h      !< water depth
   integer,                                 intent(in)  :: nconst !< number of constituents
   integer,                                 intent(in)  :: iconst !< particle tracer constituent number
   double precision, dimension(Nconst,Ndx), intent(out) :: c      !< constituents

   integer :: i, k

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "comp_concentration", ithndl )

   do i=1,Ndx
      c(iconst,i) = 0d0
   end do

   !  count number of particles per cell
   do i=1,Npart
      k = kpart(i)
      if ( k.eq.0 ) cycle

      k = iabs(cell2nod(k))

      c(iconst,k) = c(iconst,k) + 1
   end do

   !  compute concentration (parts per unit volume)
   do k=1,Ndx
      if ( h(k) .gt. epshs ) then
         c(iconst,k) = c(iconst,k) / (ba(k)*(h(k)-bl(k)))
      else
         c(iconst,k) = 0d0
      end if
   end do

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine comp_concentration

   !> allocate auxiliary fluxes
   subroutine alloc_auxfluxes()
   use m_particles
   use m_flowgeom, only: Ndx, Lnx
   use m_alloc
   implicit none

   if ( timestep.gt.0d0 ) then
      call realloc(hbegin, Ndx, fill=0d0, keepExisting=.false.)
      call realloc(qpart, Lnx, fill=0d0, keepExisting=.false.)
   end if

!   if ( threeDtype.eq.1 ) then
!      call realloc(qfreesurf, Lnx, fill=0d0, keepExisting=.false.)
!   end if

   return
   end subroutine alloc_auxfluxes

   !> deallocate auxiliary fluxes
   subroutine dealloc_auxfluxes()
   use m_particles
   implicit none

   if ( allocated(hbegin) ) deallocate(hbegin)
   if ( allocated(qpart)  ) deallocate(qpart)

   if ( allocated(qfreesurf) ) deallocate(qfreesurf)

   return
   end subroutine dealloc_auxfluxes

   subroutine part_readhydstep(hyd,itime,istat)
   use m_flow
   use m_flowgeom, only: ba, ndx, lnx, lne2ln
   use m_flowtimes, only: time1
   use hydmod
   use timers

   implicit none

   type(t_hyd)        ,intent(inout) :: hyd           !< description of the hydrodynamics
   integer            ,intent(inout) :: itime
   integer            ,intent(inout) :: istat

   !  local declarations
   integer iseg, iq, L

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "part_readhydstep", ithndl )

   if (istat==0) then
      do iseg = 1, hyd%noseg
         vol0(iseg) = vol1(iseg)
         h0(iseg) = h1(iseg)
      end do
      do iq = 1, lnx
         q0(iq) = q1(iq)
      end do
   end if

   call read_hyd_step(hyd,itime,istat)
   if (nint(time1) /= itime) then
      istat = 99
   end if
   if (istat == 0) then
      do iseg = 1, ndx !hyd%noseg
         vol1(iseg) = real(hyd%volume(iseg),8)
         h1(iseg) = vol1(iseg)/ba(iseg)
      end do
      do L = 1, lnx !hyd%noq
         iq = lne2ln(L)
         if (L > 0) then
            if (hyd%edge_type(iq) == 2 .or. hyd%edge_type(iq) == 22) then
               q1(iq) = real(hyd%flow(iq),8)
            else
               q1(iq) = real(-hyd%flow(iq),8)
            endif
         endif
      end do
   end if

   if ( timon ) call timstop ( ithndl )

   end subroutine part_readhydstep

   subroutine part_sumfluxes(q1,Dts)
   use m_particles
   use m_partmesh
   use m_flowgeom, only: Lnx
   implicit none

   double precision, dimension(Lnx), intent(in) :: q1  !< fluxes
   double precision,                 intent(in) :: Dts !< time interval

   integer :: L

   do L=1,Lnx
      qpart(L) = qpart(L) + q1(L)*dts
   end do

   return
   end subroutine part_sumfluxes

   !> update particles or add to summed fluxes
   subroutine update_part()
   use m_particles
   use m_flowtimes
   use m_flowgeom, only: Lnx, bl
   use m_flow
   use m_transport, only: numconst, constituents
   use m_missing
   use timers

   implicit none

   integer                                     :: LL, Lb, Lt

   logical                     :: Lsurface

   double precision, parameter :: huni=1d0

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "update_part", ithndl )

   Lsurface = .false. !( threeDtype.eq.1 )

   if ( Lsurface ) then
      do LL=1,Lnx
         !         call getLbotLtop(LL,Lb,Lt)
         qfreesurf(LL) = huni ! u1(Lt)*huni*wu(LL) -> q from top layer!
      end do
   end if

   if ( time0.ge.starttime ) then

      call add_particles_from_release_file(time0)
      if ( timestep.le.0d0 ) then   ! update particles every computational time step
         if ( .not.Lsurface ) then
            call update_particles(q1,h0,h1,dts)
         else
            call update_particles(qfreesurf,bl+huni,bl+huni,dts)
         end if
         timepart = time1
      else

         !        check if timestep has been started
         if ( timelast.eq.DMISS ) then
            !           start particle timestep
            timelast = time0
            timenext = time0+timestep
            hbegin = h0
            qpart = 0d0
         end if

         !        sum fluxes of this computational time step
         !         if ( .not.Lsurface ) then
         !            call part_sumfluxes(q1,Dts)
         !         else
         !            call part_sumfluxes(qfreesurf,Dts)
         !         end if

         if ( time1.ge.timenext ) then
            !           finish particle timestep
            qpart = qpart/(time1-timelast)
            if ( .not.Lsurface ) then
               call update_particles(q0, hbegin, h1, time1-timelast)
               !               call update_particles(qpart, hbegin, h1, time1-timelast)
            else
               call update_particles(qfreesurf, bl+huni, bl+huni, time1-timelast)
            end if
            timepart = time1

            !           start new particle timestep
            timelast = time1
            timenext = time1 + timestep
            hbegin = h1
            qpart = 0d0
         end if
      end if
   end if

   if ( timon ) call timstop ( ithndl )

   return
   end subroutine update_part

   subroutine finalize_part()
   use m_particles
   implicit none

   call dealloc_partmesh()
   call dealloc_partfluxes()
   call dealloc_partrecons()
   call dealloc_particles()
   call dealloc_auxfluxes()

   return
   end subroutine finalize_part

   !> return common node of links L1 and L2
   subroutine find_common_node(L1, L2, node)

   use network_data
   use m_missing

   implicit none

   integer, intent(in)   :: L1, L2           !< links
   integer, intent(out)  :: node             !< common node

   integer, dimension(4) :: a                ! dummy array with nodes of L1 and L2
   ! integer, parameter    :: IMISS = -999999

   a(1:2)    = kn(1:2, L1)
   a(3:4)    = kn(1:2, L2)

   do
      node = IMISS

      if ( a(1).eq.a(3) .or. a(1).eq.a(4) ) node = a(1)
      if ( a(2).eq.a(3) .or. a(2).eq.a(4) ) node = a(2)

      if ( node.ne.IMISS ) exit

      write(6,*) 'find_common_node: no common node found'
      exit
   end do

   end subroutine find_common_node

   subroutine dlinedis2(x3,y3,x1,y1,x2,y2,ja,dis,xn,yn,rl)
   use m_sferic
   use geometry_module, only: getdx, getdy, dbdistance, sphertocart3d, cart3dtospher
   use m_missing, only: dmiss

   implicit none
   integer          :: ja
   double precision :: x1,y1,x2,y2,x3,y3,dis,xn,yn,zn, d2
   double precision :: r2,rl,x21,y21,z21,x31,y31,z31
   double precision :: xx1,xx2,xx3,yy1,yy2,yy3,zz1,zz2,zz3,xxn,yyn,zzn

   ja  = 0

   if (jsferic == 0 .or. jasfer3d == 0) then

      x21 = getdx(x1,y1,x2,y2,jsferic)
      y21 = getdy(x1,y1,x2,y2,jsferic)
      x31 = getdx(x1,y1,x3,y3,jsferic)
      y31 = getdy(x1,y1,x3,y3,jsferic)
      r2  = dbdistance(x2,y2,x1,y1,jsferic, jasfer3d, dmiss)
      r2  = r2*r2
      if (r2 .ne. 0) then
         rl  = (x31*x21 + y31*y21) / r2
         if (0d0 .le. rl .and. rl .le. 1d0) then
            ja = 1
         endif
         xn  = x1 + rl*(x2-x1)
         yn  = y1 + rl*(y2-y1)
         dis = dbdistance(x3,y3,xn,yn,jsferic, jasfer3d, dmiss)
      else
         dis = dbdistance(x3,y3,x1,y1,jsferic, jasfer3d, dmiss)
      endif

   else

      call sphertocart3d(x1,y1,xx1,yy1,zz1)
      call sphertocart3d(x2,y2,xx2,yy2,zz2)
      call sphertocart3d(x3,y3,xx3,yy3,zz3)

      x21 = xx2-xx1
      y21 = yy2-yy1
      z21 = zz2-zz1
      x31 = xx3-xx1
      y31 = yy3-yy1
      z31 = zz3-zz1

      r2  = x21*x21 + y21*y21 + z21*z21
      if (r2 .ne. 0d0) then
         rl = (x31*x21 + y31*y21 + z31*z21) / r2
         if (0d0 .le. rl .and. rl .le. 1d0) then
            ja = 1
         endif
         xxn  = xx1 + rl*x21
         yyn  = yy1 + rl*y21
         zzn  = zz1 + rl*z21
         x31 = xxn-xx3
         y31 = yyn-yy3
         z31 = zzn-zz3
         dis = sqrt(x31*x31 + y31*y31 + z31*z31)

         call cart3dtospher(xxn,yyn,zzn,xn,yn,maxval((/x1,x2,x3/)))
      else
         dis = dbdistance(x3,y3,x1,y1, jsferic, jasfer3d, dmiss)
      endif

   endif

   return
   end subroutine dlinedis2

   subroutine dlinedis3D(xx3,yy3,zz3,xx1,yy1,zz1,xx2,yy2,zz2,JA,DIS,xxn,yyn,zzn,rl)

   implicit none

   integer          :: ja
   double precision :: dis,xn,yn,zn, d2
   double precision :: r2,rl,x21,y21,z21,x31,y31,z31
   double precision :: xx1,xx2,xx3,yy1,yy2,yy3,zz1,zz2,zz3,xxn,yyn,zzn

   !  korste afstand tot lijnelement
   ja  = 0
   x21 = xx2-xx1
   y21 = yy2-yy1
   z21 = zz2-zz1
   x31 = xx3-xx1
   y31 = yy3-yy1
   z31 = zz3-zz1
   r2  = x21*x21 + y21*y21 + z21*z21
   if (r2 .ne. 0d0) then
      rl = (x31*x21 + y31*y21 + z31*z21) / r2
      if (0d0 .le. rl .and. rl .le. 1d0) then
         ja = 1
      endif
      xxn  = xx1 + rl*x21
      yyn  = yy1 + rl*y21
      zzn  = zz1 + rl*z21
      x31 = xxn-xx3
      y31 = yyn-yy3
      z31 = zzn-zz3
      dis = sqrt(x31*x31 + y31*y31 + z31*z31)
   else
      dis = 0d0
   endif

   return

   end subroutine dlinedis3d

   SUBROUTINE REASAM(MSAM, JADOORLADEN)
   USE M_MISSING
   USE M_SAMPLES
   use m_alloc
   use MessageHandling
   implicit none
   integer, intent(inout) :: msam        !< already opened file pointer to sample file
   integer, intent(in)    :: jadoorladen !< whether to append to global set (1) or start empty (0)
   integer :: ierr
   integer :: jflow
   integer :: jqn
   integer :: mcs
   integer :: ncs
   integer :: ndraw
   integer :: nkol
   integer :: nrow
   integer :: ns1
   integer :: nsm
   integer :: num
   integer :: K, K0
   double precision :: x, y, z
   double precision :: XX, YY, ZZ, ZZ2


   COMMON /PHAROSFLOW/  JFLOW
   COMMON /PHAROSLINE/  REC1
   COMMON /SAMPLESADM/  MCS,NCS,NS1
   COMMON /QNRGF/ JQN
   COMMON /DRAWTHIS/ ndraw(50)

   CHARACTER REC*132, TEX*10, REC1*132
   LOGICAL THISISANUMBER

   CALL SAVESAM()
   NSM = 0
   MXSAM = 0
   MYSAM = 0
   IPSTAT = IPSTAT_NOTOK
   nkol = 0
11 READ (MSAM,'()',END = 31)
   NSM = NSM + 1
   GOTO 11
31 NSMAX = 1.2d0*(NSM + JADOORLADEN*NS)
   IF (ALLOCATED (XS) ) DEALLOCATE (XS,YS,ZS)
   ALLOCATE (XS(NSMAX),YS(NSMAX),ZS(NSMAX),STAT=IERR)
   IF (IERR.NE.0) CALL MESS(LEVEL_ERROR, 'XS(NSMAX),YS(NSMAX),ZS(NSMAX)',IERR,NSMAX)
   if ( allocated(ipsam) ) deallocate(ipsam)
   allocate(ipsam(NSMAX),stat=ierr)
   IF (IERR.NE.0) CALL MESS(LEVEL_ERROR, 'ipsam(NSMAX)',ierr,NSMAX)

   REWIND(MSAM)

   WRITE(TEX,'(I10)') NSM
   IF (JADOORLADEN .EQ. 0) THEN
      XS=dmiss
      YS=dmiss
      ZS=dmiss
      K = 0
   ELSE
      CALL RESTORESAM()
      K   = NS
      NS1 = NS
   ENDIF
   K0 = K

   !    check of dit een PHAROS file is
   JFLOW = 1
14 READ (MSAM,'(A)',END = 30) REC1
   if (rec1(1:1) == '*') goto 14

   IF ( .NOT. (THISISANUMBER(REC1)) ) THEN
      READ (MSAM,'(A)',END = 30) REC
      IF ( THISISANUMBER(REC) ) THEN
         READ (REC,*,ERR = 16) NROW,NKOL
         GOTO 15
16       CONTINUE
         READ (MSAM,'(A)',END = 30) REC
         READ (REC,*,ERR = 15) NUM,X,Y,Z
         JFLOW = 3
      ENDIF
   ENDIF
15 CONTINUE

   REWIND (MSAM)


   KMOD = max(1, NSM/100)
10 CONTINUE
   READ (MSAM,'(A)',END = 30) REC
   IF (REC(1:1) .EQ. '*') GOTO 10
   IF ( .NOT. (THISISANUMBER(REC)) ) THEN
      !        we nemen aan dat er net een blokcode is gelezen
      !        en we lezen meteen de nrow ncol regel, maar checken die regel niet
      READ  (MSAM,'(A)',END = 30) REC
   ELSE

      IF (JFLOW .EQ. 3) THEN
         READ (REC,*,ERR = 40) NUM,XX,YY,ZZ
      ELSE IF (NKOL == 4) THEN
         READ (REC,*,ERR = 40) XX,YY, ZZ, ZZ2
         if (zz .ne. -999d0) then
            zz = sqrt(zz*zz + zz2*zz2)
         endif
      ELSE
         READ (REC,*,end = 40) XX,YY,ZZ
         READ (REC,*,ERR = 40) XX,YY,ZZ
      ENDIF

      IF (K  .LE. NSMAX-1 .AND. XX .NE. XYMIS .AND.   &
         ZZ .NE. dmiss .AND. ZZ .NE. 999.999d0 .and. &
         .not.(isnan(XX) .or. isnan(YY) .or. isnan(ZZ)) ) THEN
         K     = K + 1
         NS    = K
         XS(K) = XX
         YS(K) = YY
         ZS(K) = ZZ
      ENDIF
   ENDIF
   GOTO 10

40 CONTINUE
   WRITE(TEX,'(I10)') K
   CALL MESS(LEVEL_ERROR, 'ERROR READING SAMPLES FILE LINE NR ',TEX,REC)

30 CONTINUE
   IF (K .GT. NSMAX) THEN
      WRITE(TEX,'(I8)') NSMAX
      CALL MESS(LEVEL_WARN, 'ONLY',TEX,'SAMPLE POINTS CAN BE LOADED')
      WRITE(TEX,'(I8)') K
      CALL MESS(LEVEL_ERROR, 'YOU TRIED TO LOAD',TEX,'SAMPLE POINTS')
   ENDIF
   WRITE(TEX,'(I10)') NS
   IF (NS .GT. 1) THEN
      CALL TIDYSAMPLES(XS,YS,ZS,IPSAM,NS,MXSAM,MYSAM)
      call get_samples_boundingbox()
      IPSTAT = IPSTAT_OK
   END IF
   call doclose(MSAM)
   RETURN
   END


   !>    delete samples
   !>      jaconfirm=0: do not prompt for confirmation,       keep arrays,        make copy
   !>                1:        prompt for confirmation,       keep arrays,        make copy
   !>               -1: do not prompt for confirmation, deallocate arrays, do not make copy
   SUBROUTINE DELSAM(JACONFIRM)      ! SPvdP: need promptless delsam in orthogonalisenet
   USE M_SAMPLES
   use m_polygon
   USE m_missing
   use geometry_module, only: dbpinpol


   implicit none

   integer, intent(in) :: JACONFIRM  !< prompt for confirmation (1) or not (0)

   integer :: i
   integer :: inhul
   integer :: ja
   integer :: k
   integer :: key
   integer :: nsol
   double precision :: rd
   double precision :: xi
   double precision :: yi

   if (jaconfirm == -1) then
      if (nsmax > 0) then
         nsmax = 0 ; ns = 0
         if ( allocated(xs)    ) deallocate (xs, ys, zs)
         if ( allocated(ipsam) ) deallocate(ipsam)
      endif
      return
   endif

   IF (Npl .LE. 2) THEN
      JA = 1
      IF (JA .EQ. 0) THEN
         KEY = 0
         RETURN
      ENDIF
      CALL SAVESAM()
      DO 5 I = 1,NS
         XS(I) = DMISS
         YS(I) = DMISS
         ZS(I) = DMISS
         ipsam(i) = 0
5     CONTINUE
      NS = 0
      RETURN
   ENDIF
   ! Else: check in polygon
   CALL SAVESAM()
   INHUL = -1
   DO 10 I = 1,NS
      RD = ZS(I)
      XI = XS(I)
      YI = YS(I)
      CALL DBPINPOL(xI, yI, INHUL, dmiss, JINS, NPL, xpl, ypl, zpl)
      IF (INHUL .EQ. 1) ZS(I) = dmiss
10 CONTINUE

   K = 0
   NSOL = NS
   DO 20 I = 1,NS
      IF (ZS(I) .NE. dmiss) THEN
         K     = K + 1
         XS(K) = XS(I)
         YS(K) = YS(I)
         ZS(K) = ZS(I)
         ipsam(k) = ipsam(i)
      ENDIF
20 CONTINUE
   NS = K

   DO 30 I = NS+1,NSOL
      XS(I) = DMISS
      YS(I) = DMISS
      ZS(I) = DMISS
      ipsam(i) = 0
30 CONTINUE

   RETURN
   END SUBROUTINE DELSAM

   function thisisanumber(rec)
   use string_module, only: find_first_char
   !!--description-----------------------------------------------------------------
   ! NONE
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   implicit none
   !
   ! Global variables
   !
   logical                      :: thisisanumber
   character(len=*), intent(in) :: rec
   !
   !
   ! Local variables
   !
   integer                        :: ich
   integer                        :: l
   !
   !
   !! executable statements -------------------------------------------------------
   !
   !
   !     is waar als eerste character van rec een getal is.
   l = find_first_char(rec)
   if (l==0) then
      thisisanumber = .false.
   else
      ich = ichar(rec(l:l))
      if (ich==43 .or. ich==45 .or. ich==46 .or. ich>=48 .and. ich<=57) then
         thisisanumber = .true.
      else
         thisisanumber = .false.
      endif
   endif
   end function thisisanumber

   SUBROUTINE TIDYSAMPLES(XS,YS,ZS,IPSAM,NS,MXSAM,MYSAM)
   use sorting_algorithms, only: indexx
   implicit none
   integer :: ns
   double precision :: XS(NS), YS(NS), ZS(NS)   !< sample coordinates
   integer, dimension(NS), intent(out) :: IPSAM !< permutation array (increasing x-coordinate)
   integer,                intent(in)  :: MXSAM, MYSAM   !< structured sample data dimensions (>0) or unstructured (0)
   !      IF (NS .GT. 1) CALL RSORT3(XS,YS,ZS,NS)

   if ( NS.gt.1 ) then
      call indexx(Ns,xs,IPSAM)
   end if

   !!     remove double/missing samples (non-structured sample data only)
   !      if ( MXSAM*MYSAM.ne.NS ) then
   !         IF (NS .GT. 1) CALL RMDOUBLE(XS,YS,ZS,IPSAM,NS)
   !      end if

   RETURN
   END


   !>    determine sample bounding box
   subroutine get_samples_boundingbox()
   use m_samples
   use m_missing
   implicit none

   integer :: i

   xsammin =  huge(1d0)
   xsammax = -huge(1d0)
   ysammin =  huge(1d0)
   ysammax = -huge(1d0)

   do i=1,NS
      if ( xs(i).ne.DMISS .and. ys(i).ne.DMISS .and. zs(i).ne.DMISS ) then
         xsammin = min(xsammin,xs(i))
         xsammax = max(xsammax,xs(i))
         ysammin = min(ysammin,ys(i))
         ysammax = max(ysammax,ys(i))
      end if
   end do

   return
   end subroutine get_samples_boundingbox

   subroutine unc_init_trk()

   use netcdf
   use m_partfm_trk_netcdf
   use fileinfo, only: filebase
   use m_flowtimes
   use MessageHandling

   implicit none

   integer                      :: ierr

   trkncfilename = trim(filebase)//'_trk.nc'
   ierr = nf90_create(trkncfilename, 0, itrkfile)
   if (ierr /= nf90_noerr) then
      call mess(LEVEL_WARN, 'Could not create tracks file.')
   end if
   call unc_addglobalatts(itrkfile)
   ierr = nf90_def_dim(itrkfile, 'time', nf90_unlimited, id_trk_timedim)
   ierr = nf90_def_var(itrkfile, 'time', nf90_double, id_trk_timedim, id_trk_time)
   ierr = nf90_put_att(itrkfile, id_trk_time,  'units'        , trim(Tudunitstr))
   ierr = nf90_put_att(itrkfile, id_trk_time,  'standard_name', 'time')

   call unc_write_part_header(itrkfile,id_trk_timedim,id_trk_partdim,id_trk_parttime,id_trk_partx,id_trk_party,id_trk_partz)
   ierr = nf90_enddef(itrkfile)
   it_trk = 0
   end subroutine unc_init_trk

   subroutine unc_write_trk()

   use m_flowtimes
   use netcdf
   use m_partfm_trk_netcdf
   use MessageHandling
   use timers

   implicit none

   integer                      :: ierr, time_trk

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "unc_write_trk", ithndl )

   ! Increment output counters in m_flowtimes.
   time_trk = nint(time1)
   it_trk   = it_trk + 1
   ierr = nf90_put_var(itrkfile, id_trk_time, time_trk, (/ it_trk /))

   call unc_write_part(itrkfile,it_trk,id_trk_parttime,id_trk_partx,id_trk_party,id_trk_partz)

   if ( timon ) call timstop ( ithndl )

   return

   end subroutine unc_write_trk

   subroutine unc_close_trk()

   use netcdf
   use m_partfm_trk_netcdf
   use MessageHandling

   implicit none

   integer                      :: ierr

   ierr = nf90_close(itrkfile)
   end subroutine unc_close_trk

   !> Writes the (possibly aggregated) unstructured network and edge type to a netCDF file for DelWAQ.
   !! If file exists, it will be overwritten.
   subroutine unc_init_map(crs, meshgeom)

   use m_partfm_map_netcdf
   use fileinfo, only: filebase
   use netcdf
   use io_ugrid
   use m_flowtimes
   use m_flowgeom
   use m_alloc
   use m_missing

   implicit none

   type(t_crs),         intent(in)  :: crs      !< Optional crs containing metadata of unsupported coordinate reference systems
   type(t_ug_meshgeom), intent(in)  :: meshgeom !< The complete mesh geometry in a single struct.

   character(len=10)                :: cell_method   !< cell_method for this variable (one of 'mean', 'point', see CF for details).
   character(len=50)                :: cell_measures !< cell_measures for this variable (e.g. 'area: mesh2d_ba', see CF for details).
   integer                          :: ierr     !< Result status (UG_NOERR==NF90_NOERR if successful).
   logical                          :: success  !< Helper variable.

   !   character(len=*)                 :: var_name      !< Variable name for in NetCDF variable, will be prefixed with mesh name.
   !   character(len=*)                 :: standard_name !< Standard name (CF-compliant) for 'standard_name' attribute in this variable.
   !   character(len=*)                 :: long_name     !< Long name for 'long_name' attribute in this variable (use empty string if not wanted).
   !   character(len=*)                 :: unit          !< Unit of this variable (CF-compliant) (use empty string for dimensionless quantities).

   mapncfilename = trim(filebase)//'_map.nc'
   ierr = nf90_create(mapncfilename, 0, imapfile)
   if (ierr /= nf90_noerr) then
      call mess(LEVEL_WARN, 'Could not create map file.')
   end if

   ierr = UG_NOERR

   ! Add global attributes to NetCDF file.
   call unc_addglobalatts(imapfile)
   ierr = nf90_def_dim(imapfile, 'time', nf90_unlimited, id_map_timedim)
   ierr = nf90_def_var(imapfile, 'time', nf90_double, id_map_timedim, id_map_time)
   ierr = nf90_put_att(imapfile, id_map_time,  'units'        , trim(Tudunitstr))
   ierr = nf90_put_att(imapfile, id_map_time,  'standard_name', 'time')

   ! Write mesh geometry.
   ierr = ug_write_mesh_struct(imapfile, meshids, networkids, crs, meshgeom)
   if (ierr /= nf90_noerr) then
      call mess(LEVEL_ERROR, 'Could not write geometry to file map')
      return
   end if

   cell_method = 'mean' !< Default cell average.
   cell_measures = ''
   !   cell_measures = 'area: '//trim(mesh2dname)//'_flowelem_ba' ! relies on unc_write_flowgeom_ugrid_filepointer

   !   if (kmx == 1) then
   !               ierr = unc_def_var_map(mapids%ncid, mapids%id_tsp, mapids%id_map_depth_averaged_particle_concentration, nf90_double, UNC_LOC_S, 'depth_averaged_particle_concentration', 'depth_averaged_particle_concentration', 'depth-averaged particle concentration', 'm-3', jabndnd=jabndnd_)

   !              function unc_def_var_map(ncid,id_tsp, id_var, itype, iloc,
   !               var_name = 'depth_averaged_particle_concentration'
   !               standard_name = 'depth_averaged_particle_concentration'
   !               long_name = 'depth_averaged_particle_concentration'
   !               unit = 'm-3'
   !               , is_timedep, dimids, cell_method, which_meshdim, jabndnd) result(ierr)


   ierr = ug_def_var(imapfile, id_map_depth_averaged_particle_concentration, [meshids%dimids(mdim_face), id_map_timedim], nf90_double, UG_LOC_FACE, &
      trim(meshgeom%meshName), 'depth_averaged_particle_concentration', 'depth_averaged_particle_concentration', &
      'depth_averaged_particle_concentration', 'm-3', cell_method, cell_measures, crs, ifill=-999, dfill=dmiss)

   if (ierr /= nf90_noerr) then
      call mess(LEVEL_ERROR, 'Could not create concentration variable in map file')
      return
   end if
   !   else
   !      continue
   !   endif

   ierr = nf90_enddef(imapfile)
   ierr = nf90_sync(imapfile)

   it_map = 0
   call realloc(work, ndx, keepExisting = .false., fill = dmiss)

   end subroutine unc_init_map

   subroutine unc_write_map()

   use m_flowtimes
   use m_flowgeom
   use m_transport
   use m_flow, only: h1
   use m_particles, only: part_iconst
   use netcdf
   use m_partfm_map_netcdf
   use MessageHandling
   use timers

   implicit none

   integer                      :: ierr, time_map, k

   integer(4) ithndl              ! handle to time this subroutine
   data ithndl / 0 /
   if ( timon ) call timstrt( "unc_write_map", ithndl )

   ! Increment output counters in m_flowtimes.
   time_map = nint(time1)
   it_map   = it_map + 1
   ierr = nf90_put_var(imapfile, id_map_time, time_map, (/ it_map /))
   if (ierr /= nf90_noerr) then
      call mess(LEVEL_ERROR, 'Could not write geometry to file map')
      return
   end if

   call comp_concentration(h1,numconst,part_iconst,constituents)
   do k=1,Ndx
      work(k) = constituents(1,k)
   end do
   ierr = nf90_put_var(imapfile, id_map_depth_averaged_particle_concentration, work(1:Ndx), start = (/ 1, it_map /))
   if (ierr /= nf90_noerr) then
      call mess(LEVEL_ERROR, 'Could not write geometry to file map')
      return
   end if
   ierr = nf90_sync(imapfile)

   if ( timon ) call timstop ( ithndl )

   end subroutine unc_write_map

   subroutine unc_close_map()

   use netcdf
   use m_partfm_map_netcdf
   use MessageHandling

   implicit none

   integer                      :: ierr

   ierr = nf90_close(imapfile)
   end subroutine unc_close_map

   !> write particle tracks header to netcdf trk file
   subroutine unc_write_part_header(ifile,id_timedim,id_trk_partdim,id_trk_parttime,id_trk_partx,id_trk_party,id_trk_partz)
   use m_particles
   use netcdf
   use m_flow, only: kmx
   use m_sferic, only: jsferic
   use MessageHandling
   use m_missing
   implicit none

   integer, intent(in)    :: ifile  !< output file identifier
   integer, intent(in)    :: id_timedim
   integer, intent(inout) :: id_trk_partdim, id_trk_parttime, id_trk_partx, id_trk_party, id_trk_partz

   character(len=128)     :: mesg

   integer                :: ierr
   integer                :: jaInDefine

   ! Put dataset in define mode (possibly again) to add dimensions and variables.
   ierr = nf90_redef(ifile)
   if (ierr == nf90_eindefine) jaInDefine = 1 ! Was still in define mode.
   if (ierr /= nf90_noerr .and. ierr /= nf90_eindefine) then
      call mess(LEVEL_ERROR, 'Could not put header in flow geometry file.')
      !       call check_error(ierr)
      return
   end if

   ierr = nf90_def_dim(ifile, 'particles', NpartTot, id_trk_partdim)

   ierr = nf90_def_var(ifile, 'particles_time', nf90_double, id_timedim, id_trk_parttime)
   ierr = nf90_put_att(ifile, id_trk_parttime, 'long_name', 'particles time')

   ierr = nf90_def_var(ifile, 'particles_x_coordinate', nf90_double, (/ id_trk_partdim, id_timedim /), id_trk_partx)
   ierr = nf90_def_var(ifile, 'particles_y_coordinate', nf90_double, (/ id_trk_partdim, id_timedim /), id_trk_party)
   if (jsferic == 0) then
      ierr = nf90_put_att(ifile, id_trk_partx, 'units',         'm')
      ierr = nf90_put_att(ifile, id_trk_party, 'units',         'm')
      ierr = nf90_put_att(ifile, id_trk_partx, 'standard_name', 'projection_x_coordinate')
      ierr = nf90_put_att(ifile, id_trk_party, 'standard_name', 'projection_y_coordinate')
      ierr = nf90_put_att(ifile, id_trk_partx, 'long_name'    , 'x-coordinate')
      ierr = nf90_put_att(ifile, id_trk_party, 'long_name'    , 'y-coordinate')
   else
      ierr = nf90_put_att(ifile, id_trk_partx, 'units',         'degrees_east')
      ierr = nf90_put_att(ifile, id_trk_party, 'units',         'degrees_north')
      ierr = nf90_put_att(ifile, id_trk_partx, 'standard_name', 'longitude')
      ierr = nf90_put_att(ifile, id_trk_party, 'standard_name', 'latitude')
      ierr = nf90_put_att(ifile, id_trk_partx, 'long_name'    , 'longitude')
      ierr = nf90_put_att(ifile, id_trk_party, 'long_name'    , 'latitude')
   end if
   ierr = nf90_put_att(ifile, id_trk_partx, 'long_name', 'x-coordinate of particles')
   ierr = nf90_put_att(ifile, id_trk_party, 'long_name', 'y-coordinate of particles')
   ierr = nf90_put_att(ifile, id_trk_partx, '_FillValue', dmiss)
   ierr = nf90_put_att(ifile, id_trk_party, '_FillValue', dmiss)


   if ( kmx.gt.0 ) then
      ierr = nf90_def_var(ifile, 'particle_z_coordinate', nf90_double, (/ id_trk_partdim, id_timedim /), id_trk_partz)
      ierr = nf90_put_att(ifile, id_trk_partz, 'long_name', 'z-coordinate of particle')
   end if

   ! Leave the dataset in the same mode as we got it.
   if (jaInDefine == 1) then
      ierr = nf90_redef(ifile)
   end if

   return
   end subroutine unc_write_part_header

   !> write particles to netcdf file
   subroutine unc_write_part(ifile,itime,id_trk_parttime,id_trk_partx,id_trk_party,id_trk_partz)
   use m_particles
   use netcdf
   use m_sferic
   use m_flow, only: kmx
   use geometry_module, only: cart3Dtospher
   use m_missing
   use MessageHandling
   use m_alloc
   implicit none

   integer,                        intent(in)  :: ifile  !< output file identifier
   integer,                        intent(in)  :: itime
   integer,                        intent(in)  :: id_trk_parttime, id_trk_partx, id_trk_party, id_trk_partz

   double precision, dimension(:), allocatable :: xx, yy

   double precision                            :: dis2

   integer                                     :: i, i0, ii, iglb
   integer                                     :: ierr, ierror

   double precision,                 parameter :: dtol = 1d-8

   integer, save :: icount=0

   ierror = 1

   icount = icount+1

   if ( icount.ge.24 ) then
      continue
   end if

   !  allocate
   call realloc(xx, NpartTot, keepExisting=.false., fill = dmiss)
   call realloc(yy, NpartTot, keepExisting=.false., fill = dmiss)

   if ( jsferic.eq.1 ) then
      do ii=1,Npart
         call Cart3Dtospher(xpart(ii),ypart(ii),zpart(ii),xx(ii),yy(ii),0d0)
      end do
   else
      do ii=1,NpartTot
         xx(ii) = xpart(ii)
         yy(ii) = ypart(ii)
      end do
   end if

   ierr = nf90_put_var(ifile, id_trk_parttime, timepart, (/ itime /))
   if ( ierr.ne.0 ) goto 1234
   ierr = nf90_put_var(ifile, id_trk_partx, xx, start=(/ 1,itime /), count=(/ NpartTot,1 /) )
   if ( ierr.ne.0 ) goto 1234
   ierr = nf90_put_var(ifile, id_trk_party, yy, start=(/ 1,itime /), count=(/ NpartTot,1 /) )
   if ( ierr.ne.0 ) goto 1234

   if ( kmx.gt.0 ) then
      !     particle vertical coordinate
   end if

   ierror = 0
1234 continue

   !  error handling
   if ( ierror.ne.0 ) then
      call mess(LEVEL_ERROR, 'particles output error')
   end if

   return
   end subroutine unc_write_part

   !> Puts global attributes in NetCDF data set.
   !! This includes: institution, Conventions, etc.
   subroutine unc_addglobalatts(ncid)
   use netcdf
   use MessageHandling
   use part_version_module

   integer, intent(in) :: ncid

   character*8  :: cdate
   character*10 :: ctime
   character*5  :: czone
   integer :: ierr, jaInDefine
   ierr = nf90_noerr
   jaInDefine = 0

   ierr = nf90_redef(ncid)
   if (ierr == nf90_eindefine) jaInDefine = 1 ! Was still in define mode.
   if (ierr /= nf90_noerr .and. ierr /= nf90_eindefine) then
      write (msgbuf, '(a,i0,a,i0,a,a)') 'Could not put global attributes in NetCDF #', ncid, '. Error code ', ierr, ': ', nf90_strerror(ierr)
      call err_flush()
      return
   end if

   ierr = nf90_put_att(ncid, nf90_global,  'institution', trim(part_company))
   ierr = nf90_put_att(ncid, nf90_global,  'references', trim(part_company_url))
   ierr = nf90_put_att(ncid, nf90_global,  'source', part_version_full)

   call date_and_time(cdate, ctime, czone)
   ierr = nf90_put_att(ncid, nf90_global,  'history', &
      'Created on '//cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8)//'T'//ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)//czone(1:5)// ', '//trim(part_program))
   ierr = nf90_put_att(ncid, nf90_global,  'date_created',  cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8)//'T'//ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)//czone(1:5))
   ierr = nf90_put_att(ncid, nf90_global,  'date_modified', cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8)//'T'//ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)//czone(1:5))

   ierr = nf90_put_att(ncid, nf90_global,  'Conventions', 'CF-1.5 Deltares-0.1')

   ! Leave the dataset in the same mode as we got it.
   if (jaInDefine == 0) then
      ierr = nf90_enddef(ncid)
   end if
   end subroutine unc_addglobalatts

   !> Sets the UDUnit timestring based on current model time settings.
   !! Module variable Tudunitstr can the be used in various output routines.
   subroutine setTUDUnitString()
   use m_flowtimes

   implicit none

   integer          :: Tzonehrs
   character(len=1) :: Tzonesgn

   Tzonehrs = int(TZone)
   if (Tzone<0) then
      Tzonesgn = '-'
   else
      Tzonesgn = '+'
   end if
   write(Tudunitstr,'(a,i2.2,a)') 'seconds since '//refdat(1:4)//'-'//refdat(5:6)//'-'//refdat(7:8)//' 00:00:00 '//Tzonesgn, abs(Tzonehrs),':00'

   end subroutine setTUDUnitString

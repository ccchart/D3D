!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2022.                                
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

! $Id$
! $HeadURL$

   subroutine fm_fallve()
   !!--description-----------------------------------------------------------------
   !
   !    Function: Relation between sediment concentration
   !              and vertical fall velocity. Model for
   !              hindered settling.
   !              Fall velocity at layer interfaces.
   !!--declarations----------------------------------------------------------------
   use precision
   use m_physcoef, only: ee, ag, sag, vonkar, frcuni, backgroundsalinity, backgroundwatertemperature
   use m_sediment, only: stmpar, sedtra, stm_included, mtd, vismol
   use m_flowtimes, only: time1
   use m_flowgeom, only: ndx, ln, kfs,bl, wcl, lnx
   use m_flow    , only: ifrctypuni, z0, hs, iturbulencemodel,kbot,ktop,kmx,zws,ucxq,ucyq,sa1,tem1,ucx,ucy,ucz,ndkx,s1,z0urou,ifrcutp,hu,frcu,ucx_mor,ucy_mor
   use m_flowparameters, only: jasal, jatem, jawave, epshs, flowWithoutWaves
   use m_transport, only: constituents, ised1, isalt, itemp
   use m_turbulence, only:turkinepsws, rhowat
   use morphology_data_module
   use message_module, only: write_error
   use unstruc_files, only: mdia
   use m_alloc
   use m_fm_erosed, only: ucxq_mor, ucyq_mor, taub
   use flocculation, only: get_tshear_tdiss
   !
   implicit none
   !
   real(fp), dimension(:)             , pointer :: localpar
   real(fp), dimension(:,:)           , pointer :: ws
   real(fp)                           , pointer :: csoil
   real(fp)           , dimension(:)  , pointer :: rhosol
   real(fp)         , dimension(:,:)  , pointer :: dss
   real(fp)           , dimension(:)  , pointer :: sedd50
   real(fp)           , dimension(:)  , pointer :: sedd50fld

   character(256)     , dimension(:)  , pointer :: dll_usrfil
   character(256)     , dimension(:)  , pointer :: dll_function
   integer(pntrsize)  , dimension(:)  , pointer :: dll_handle
   integer            , dimension(:)  , pointer :: iform_settle
   real(fp)           , dimension(:,:), pointer :: par_settle
   !
   integer                            , pointer :: npar
   integer                            , pointer :: max_integers
   integer                            , pointer :: max_reals
   integer                            , pointer :: max_strings
   integer            , dimension(:)  , pointer :: dll_integers
   real(hp)           , dimension(:)  , pointer :: dll_reals
   character(256)     , dimension(:)  , pointer :: dll_strings
   !
   ! Global variables
   integer                            , pointer :: lsed
   double precision, allocatable  , dimension(:):: ctot
   !
   ! Local variables
   !
   integer                             :: k, kk, k1, k2, L, ll, nm, i, istat, kb, kt
   logical                             :: error

   double precision                    :: chezy
   double precision                    :: cz
   double precision                    :: cz_dum
   double precision                    :: h0
   double precision                    :: rhoint
   double precision                    :: salint
   double precision                    :: temint
   double precision                    :: tka
   double precision                    :: tkb
   double precision                    :: tkt
   double precision                    :: tshear
   double precision                    :: tur_eps
   double precision                    :: tur_k
   double precision                    :: u
   double precision, allocatable, dimension(:)   :: um
   double precision                    :: v
   double precision, allocatable, dimension(:)   :: vm
   double precision                    :: w
   double precision                    :: wsloc
   double precision                    :: z00
   double precision, dimension(:), allocatable    :: z0rou
   double precision                    :: thick

   character(256)              :: errmsg
   !!
   !!! executable statements -------------------------------------------------------
   !!
   call realloc(ctot, ndkx, keepExisting=.false., fill=0d0)
   call realloc(um, ndx, keepExisting=.false., fill=0d0)
   call realloc(vm, ndx, keepExisting=.false., fill=0d0)
   call realloc(z0rou,ndx, keepExisting=.false., fill=0d0)

   csoil               => stmpar%sedpar%csoil
   rhosol              => stmpar%sedpar%rhosol
   dss                 => stmpar%sedpar%dss
   sedd50              => stmpar%sedpar%sedd50
   sedd50fld           => stmpar%sedpar%sedd50fld

   npar                => stmpar%trapar%npar
   dll_usrfil          => stmpar%trapar%dll_usrfil_settle
   dll_function        => stmpar%trapar%dll_function_settle
   dll_handle          => stmpar%trapar%dll_handle_settle
   iform_settle        => stmpar%trapar%iform_settle
   par_settle          => stmpar%trapar%par_settle

   max_integers        => stmpar%trapar%max_integers_settle
   max_reals           => stmpar%trapar%max_reals_settle
   max_strings         => stmpar%trapar%max_strings_settle
   dll_integers        => stmpar%trapar%dll_integers_settle
   dll_reals           => stmpar%trapar%dll_reals_settle
   dll_strings         => stmpar%trapar%dll_strings_settle

   lsed                => stmpar%lsedsus
   ws                  => mtd%ws
   !
   allocate (localpar (npar), stat = istat)
   !
   do k = 1, ndx
      !
      ! bulk mass concentration from previous timestep
      !
!      if (kfs(k) > 0) then
         call getkbotktop(k, kb, kt)
         do kk = kb, kt
            do ll = 1, lsed
               ctot(kk) = ctot(kk) + constituents(ISED1 - 1 +ll, kk)
            end do ! kk
         end do
!      end if
   end do

   ! calculate mean velocity in nodes
!   call setucxucyucxuucyu()

   if (kmx>0) then                       ! 3D
      um = 0d0; vm = 0d0
      do k = 1, ndx
         h0 = max(s1(k)-bl(k), epshs)
         call getkbotktop(k, kb, kt)
         do kk = kb, kt
            thick = zws(kk) - zws(kk-1)
            um(k) = um(k) + thick/h0*ucxq(kk)
            vm(k) = vm(k) + thick/h0*ucyq(kk)
         end do
      end do
   else
      um   = ucxq                       ! discharge based velocities
      vm   = ucyq
   end if

   ! Calculate roughness height
   if (jawave<3 .or. flowWithoutWaves) then  ! current related, potentially including trachy
      do L=1, lnx
         k1 = ln(1,L); k2 = ln(2,L)
         if (frcu(L)>0d0) then
            call getczz0(hu(L), frcu(L), ifrcutp(L), cz_dum, z00)
         else
            call getczz0(hu(L), frcuni, ifrctypuni, cz_dum, z00)
         end if
         z0rou(k1) = z0rou(k1) + wcl(1,L)*z00
         z0rou(k2) = z0rou(k2) + wcl(2,L)*z00
      end do
      !
   else      !  wave enhanced roughness
      do L=1,lnx
         k1=ln(1,L); k2=ln(2,L)
         z0rou(k1) = z0rou(k1) + wcl(1,L)*z0urou(L)
         z0rou(k2) = z0rou(k2) + wcl(2,L)*z0urou(L)
      end do
   end if

   do ll = 1, lsed   
      !
      do i = 1, npar
         localpar(i) = par_settle(i,ll)
      enddo
      !
      do k = 1, ndx
         if (s1(k)-bl(k)<epshs) cycle
         !
         h0 = s1(k)-bl(k)
         chezy = sag * log( 1.0d0 + h0/max(1d-8,ee*z0rou(k)) ) / vonkar
         !
         ! loop over the interfaces in the vertical
         !
         if (kmx > 0) then                ! 3D
            call getkbotktop(k, kb, kt)
            do kk = kb, kt-1
               !
               ! Input parameters are passed via dll_reals/integers/strings-arrays
               !

               ! HK: is this better than first establish fallvelocity in a cell, next interpolate to interfaces?

               tka = zws(kk+1) - zws(kk) ! thickness above
               tkb = zws(kk) - zws(kk-1)   ! thickness below
               tkt = tka+tkb
               if (jasal > 0) then
                  salint = max(0.0_fp, (tka*constituents(isalt,kk+1) + tkb*constituents(isalt,kk)  ) / tkt )
               else
                  salint = backgroundsalinity
               endif
               !                !
               if (jatem > 0) then
                  temint =             (tka*constituents(itemp,kk+1) + tkb*constituents(itemp,kk)  ) / tkt
               else
                  temint = backgroundwatertemperature
               endif
               !
               rhoint = (tka*rhowat(kk+1) + tkb*rhowat(kk)) / tkt
               !
               u = (tka*ucx_mor(kk+1)+tkb*ucx_mor(kk)) / tkt   ! x component
               v = (tka*ucy_mor(kk+1)+tkb*ucy_mor(kk)) / tkt   ! y component
               w = (tka*ucz(kk+1)+tkb*ucz(kk)) / tkt   ! z component

               if (iturbulencemodel == 3) then     ! k-eps
                  tur_k = turkinepsws(1,kk)
               else
                  tur_k = -999.0d0
               endif
               if (iturbulencemodel == 3) then
                  tur_eps = turkinepsws(2,kk)
               else
                  tur_eps = -999.0d0
               endif
               if (iturbulencemodel == 3) then ! k-eps
                  call get_tshear_tdiss( tshear, tur_eps, 3, 2, tke = tur_k )
               else
                  call get_tshear_tdiss( tshear, tur_eps, 3, 0, taub = taub(k), rho_water = rhoint, waterdepth = h0, localdepth = s1(k) - zws(kk))
               endif
               !
               ! Pass to DLL, decoupling because of treatment per layer interface
               !
               if (max_reals < WS_MAX_RP) then
                  write(errmsg,'(a,a,a)') 'Insufficient space to pass real values to settling routine.'
                  call write_error(errmsg, unit=mdia)
                  error = .true.
                  return
               endif
               dll_reals(WS_RP_TIME ) = real(time1  ,hp)
               dll_reals(WS_RP_ULOC ) = real(u      ,hp)
               dll_reals(WS_RP_VLOC ) = real(v      ,hp)
               dll_reals(WS_RP_WLOC ) = real(w      ,hp)
               dll_reals(WS_RP_SALIN) = real(salint ,hp)
               dll_reals(WS_RP_TEMP ) = real(temint ,hp)
               dll_reals(WS_RP_RHOWT) = real(rhoint ,hp)
               dll_reals(WS_RP_CFRCB) = real(constituents(ised1-1+ll, kk),hp)
               dll_reals(WS_RP_CTOT ) = real(ctot(kk),hp)
               dll_reals(WS_RP_KTUR ) = real(tur_k  ,hp)
               dll_reals(WS_RP_EPTUR) = real(tur_eps,hp)
               if (sedd50(ll)<0.0_fp) then
                  dll_reals(WS_RP_D50) = real(sedd50fld(k),hp)
               else
                  dll_reals(WS_RP_D50) = real(sedd50(ll),hp)
               endif
               dll_reals(WS_RP_DSS  ) = real(dss(k,ll) ,hp)
               dll_reals(WS_RP_RHOSL) = real(rhosol(ll) ,hp)
               dll_reals(WS_RP_CSOIL) = real(csoil      ,hp)
               dll_reals(WS_RP_GRAV ) = real(ag         ,hp)
               dll_reals(WS_RP_VICML) = real(vismol     ,hp)
               dll_reals(WS_RP_WDEPT) = real(h0         ,hp)
               dll_reals(WS_RP_UMEAN) = real(um(k)      ,hp)
               dll_reals(WS_RP_VMEAN) = real(vm(k)      ,hp)
               dll_reals(WS_RP_CHEZY) = real(chezy      ,hp)
               dll_reals(WS_RP_SHTUR) = real(tshear     ,hp)
               !
               if (max_integers < WS_MAX_IP) then
                  write(errmsg,'(a,a,a)') 'Insufficient space to pass integer values to settling routine.'
                  call write_error(errmsg, unit=mdia)
                  error = .true.
                  return
               endif
               dll_integers(WS_IP_NM  ) = kk
               dll_integers(WS_IP_ISED) = ll
               !
               if (max_strings < WS_MAX_SP) then
                  write(errmsg,'(a,a,a)') 'Insufficient space to pass strings to settling routine.'
                  call write_error(errmsg, unit=mdia)
                  error = .true.
                  return
               endif
               dll_strings(WS_SP_USRFL) = dll_usrfil(ll)
               !             !
               call eqsettle(dll_function(ll), dll_handle(ll), max_integers, max_reals, max_strings, &
                           & dll_integers, dll_reals, dll_strings, mdia, iform_settle(ll),  &
                           & localpar, npar, wsloc, error)
               if (error) then
                  write(errmsg,'(a,a,a)') 'fm_fallve:: Error from eqsettle.'
                  call write_error(errmsg, unit=mdia)
                  error = .true.
                  return
               endif
               !
               ws(kk, ll) = wsloc
            end do
            if (kmx>1) then
               ws(kb-1,ll) = ws(kb,ll)     ! to check
            endif
         else                           ! 2D
            if (jasal > 0) then
               salint = max(0d0, constituents(isalt,k))
            else
               salint = backgroundsalinity
            endif
            !             !
            if (jatem > 0) then
               temint = constituents(itemp,k)
            else
               temint = backgroundwatertemperature
            endif
            !
            rhoint = rhowat(k)   
            !
            u       = ucx_mor(k)   ! x component
            v       = ucy_mor(k)   ! y component
            w       = -999d0   ! z component
            tur_k   = -999d0
            tur_eps = -999d0
            call get_tshear_tdiss( tshear, tur_eps, 2, taub = taub(k), rho_water = rhoint, waterdepth = h0 )
            !
            ! Pass to DLL
            !
            if (max_reals < WS_MAX_RP) then
               write(errmsg,'(a,a,a)') 'Insufficient space to pass real values to settling routine.'
               call write_error(errmsg, unit=mdia)
               error = .true.
               return
            endif
            dll_reals(WS_RP_TIME ) = real(time1  ,hp)
            dll_reals(WS_RP_ULOC ) = real(u      ,hp)
            dll_reals(WS_RP_VLOC ) = real(v      ,hp)
            dll_reals(WS_RP_WLOC ) = real(w      ,hp)
            dll_reals(WS_RP_SALIN) = real(salint ,hp)
            dll_reals(WS_RP_TEMP ) = real(temint ,hp)
            dll_reals(WS_RP_RHOWT) = real(rhoint ,hp)
            dll_reals(WS_RP_CFRCB) = real(constituents(ised1-1+ll, k),hp)
            dll_reals(WS_RP_CTOT ) = real(ctot(k),hp)
            dll_reals(WS_RP_KTUR ) = real(tur_k  ,hp)
            dll_reals(WS_RP_EPTUR) = real(tur_eps,hp)
            if (sedd50(ll)<0.0_fp) then
               dll_reals(WS_RP_D50) = real(sedd50fld(k),hp)
            else
               dll_reals(WS_RP_D50) = real(sedd50(ll),hp)
            endif
            dll_reals(WS_RP_DSS  ) = real(dss(k,ll) ,hp)
            dll_reals(WS_RP_RHOSL) = real(rhosol(ll) ,hp)
            dll_reals(WS_RP_CSOIL) = real(csoil      ,hp)
            dll_reals(WS_RP_GRAV ) = real(ag         ,hp)
            dll_reals(WS_RP_VICML) = real(vismol     ,hp)
            dll_reals(WS_RP_WDEPT) = real(h0         ,hp)
            dll_reals(WS_RP_UMEAN) = real(um(k)      ,hp)
            dll_reals(WS_RP_VMEAN) = real(vm(k)      ,hp)
            dll_reals(WS_RP_CHEZY) = real(chezy      ,hp)
            dll_reals(WS_RP_SHTUR) = real(tshear     ,hp)
            !
            if (max_integers < WS_MAX_IP) then
               write(errmsg,'(a,a,a)') 'Insufficient space to pass integer values to settling routine.'
               call write_error(errmsg, unit=mdia)
               error = .true.
               return
            endif
            dll_integers(WS_IP_NM  ) = k
            dll_integers(WS_IP_ISED) = ll
            !
            if (max_strings < WS_MAX_SP) then
               write(errmsg,'(a,a,a)') 'Insufficient space to pass strings to settling routine.'
               call write_error(errmsg, unit=mdia)
               error = .true.
               return
            endif
            dll_strings(WS_SP_USRFL) = dll_usrfil(ll)
            !             !
            call eqsettle(dll_function(ll), dll_handle(ll), max_integers, max_reals, max_strings, &
                        & dll_integers, dll_reals, dll_strings, mdia, iform_settle(ll),  &
                        & localpar, npar, wsloc, error)
            if (error) then
               write(errmsg,'(a,a,a)') 'fm_fallve:: Error from eqsettle.'
               call write_error(errmsg, unit=mdia)
               error = .true.
               return
            endif
            !
            ws(k, ll) = wsloc
         end if    ! kmx>0
      enddo        ! k
   enddo           ! ll

   deallocate (localpar, stat = istat)
   end subroutine fm_fallve

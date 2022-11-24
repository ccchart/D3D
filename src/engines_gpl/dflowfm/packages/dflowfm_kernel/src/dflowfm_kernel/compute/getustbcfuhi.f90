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

   subroutine getustbcfuhi( LL,Lb,ustbLL,cfuhiLL,hdzb, z00,cfuhi3D)                ! see Uittenbogaard's subroutine USTAR
   use m_flow
   use m_flowgeom  , only : ln, dxi, csu, snu, acL, lnxi
   use m_flowtimes , only : dti
   use m_waves     , only : ustokes, vstokes, wblt, hwav
   use m_sediment  , only : stm_included
   use m_turbulence, only : tkepro
   use m_flowtimes, only: dts

   implicit none
   integer,          intent (in)  :: LL, Lb
   double precision, intent (out) :: ustbLL, cfuhiLL, hdzb, z00
   double precision, intent (out) :: cfuhi3D                                       ! 3D bedfriction coeffient, advi(Lb) = advi(Lb) + cfuhi3D

   integer          :: ifrctyp, L
   double precision :: frcn, sqcf, cz, umod, u1Lb, gsx, ustw2, ustc2, fw, cdrag, abscos, dfuc, costu
   double precision :: taubpuLL                                ! taubpu = umod*ag/C2 or ypar*(taucur+tauwav)/rho/umod or ustar*ustar/u
   double precision :: taubxuLL                                ! taubxu = ymxpar*(taucur+tauwav)

   double precision :: csw, snw                                ! wave direction cosines
   double precision :: Dfu, Dfu0, Dfu1, htop, dzu              ! wave dissipation by bed friction, / (rhomean*c*deltau)
   double precision :: deltau                                  ! wave dissipation layer thickness
   double precision :: zbot, ztop, u2dh, frac
   double precision :: z0urouL, cf, ust, rz, umod1, rhoL, dzuu, uorbu
   double precision :: cwall
   double precision :: umodeps

   integer          :: nit, nitm = 100
   double precision :: r, rv = 123.8d0, e = 8.84d0 , eps = 1d-2
   double precision :: s, sd, er, ers, dzb, uu, vv, dzw, alin
   double precision :: cphi, sphi
   double precision :: fsqrtt = sqrt(2d0)
   double precision :: hul1,hul0,hul

   cfuhi3D = 0d0
   ustbLL = 0d0;  cfuhiLL = 0d0;  hdzb = 0d0; z00 = 0d0; cz = 0d0; nit = 0

   umodeps = 1d-4

   frcn = frcu(LL)
   if (frcn == 0d0 ) return
   ifrctyp = ifrcutp(LL)

   if ( hu(LL) < trsh_u1Lb) then
      gsx = ag*( s1(ln(2,LL)) - s1(ln(1,LL)) ) * dxi(LL)
   endif

   if (ifrctyp < 10) then
      if (frcn > 0d0 ) then
         call getczz0(hu(LL), frcn, ifrctyp, cz, z00)

         hdzb  = 0.5d0*hu(Lb)     + c9of1*z00                ! half bottom layer plus 9z0

         if (z00 > 0d0) then

            if (jaustarint == 0) then
               ! sqcf = vonkar/log(c9of1 + hdzb/z00)            ! till 012015
               sqcf = vonkar/log(hdzb/z00)
            else if (jaustarint == 1) then                      ! Yoeri 2014 long time default for jaustarint == 1
               dzb  = hu(Lb) + c9of1*z00
               sqcf = vonkar / ( log(dzb/z00)-1d0 )
            else if (jaustarint == 2) then                      ! remobilised through jaustarint == 2, good convergence
               dzb  = hu(Lb)/ee + c9of1*z00
               sqcf = vonkar / ( log(dzb/z00) )
            else if (jaustarint == 3) then                      ! Delft3D
               hdzb  = 0.5d0*hu(Lb)     + z00
               sqcf = vonkar / ( log(1d0+0.5d0*hu(Lb)/z00) )
            else if(jaustarint == 4) then
               !hdzb  = 0.5d0*hu(Lb)     + c9of1*z00/0.65d0
               dzb  = hu(Lb)/ee + c9of1*z00 *0.66d0
               sqcf = vonkar / ( log(dzb/z00) )
            else if (jaustarint == 5) then
               dzb  = hu(Lb)
               sqcf = vonkar / ( ( 1d0 + c9of1 * z00 / dzb ) * log(dzb/z00+c9of1) - c9of1 * z00/dzb * log(c9of1) - 1d0 )
            endif
            z0ucur(LL) = z00
         else
            sqcf = 0d0
         endif
      else
         hdzb = 0.5d0*hu(Lb)
         sqcf = 0d0
      endif

      u1Lb = u1(Lb)

10    continue

      umod = sqrt( u1Lb*u1Lb + v(Lb)*v(Lb) )

      if (umod == 0d0) then            ! from dry to wet
         umod = max(umodeps, dts*ag*dxi(LL)*min( abs( s1(ln(1,LL)) - s1(ln(2,LL)) ), 0.333333d0*hu(LL) ) )
      else
         umod = max(umod, umodeps)     ! 1d-6 for klopman     ! until 3D handled like 2D iterative loop , solves Roses problem: ust=1.1e-104 to the power 3 is underflow
      endif

      ustbLL = sqcf*umod                                   ! ustar based upon bottom layer/layer integral velocity

      if (jawave>0 .and. .not. flowWithoutWaves) then
         !
         ! get stokes drift and some wave parameters
         call getustwav(LL, z00, umod, uorbu)
         ac1 = acl(LL); ac2 = 1d0-ac1
         k1  = ln(1,LL); k2 = ln(2,LL)
         ! Overwrite ustbLL with Nguyen version
         hrmsLL = ac1*hwav(k1) + ac2*hwav(k2)
         twavLL = ac1*twav(k1) + ac2*twav(k2)
         omeg   = twopi/twavLL
         aorb   = uorbu/omeg
         kn     = 30d0*z00
         !
         if (hrmsLL<1d-3) then                              !current only
            deltau = z00*ee
         elseif (abs(umod).lt.(0.1*uorbu)) then             !wave only
            deltau = 0.072*aorb*(aorb/kn)**(-0.25)          !Johnsen and Carlsen (1976)
         else    !wave + current
            !
            deltau = 0.2*aorb*(aorb/kn)**(-0.25) * (1.+abs(umod/uorbu))                 ! Modified van Rijn 2011 by Nguyen (2021)
         endif
         ka       = 30d0*deltau/ee
         fc       = 0.242/(log10(12.*huLL/ka))**2
         fc0      = 0.242/(log10(12.*huLL/kn))**2
         tauw     = 0.5d0*rhomean*fw*uorbu**2
         !     According to Feddersen (2000) (for random waves)
         tauwc    = 0.5*rhomean*fc*sqrt(umod**2+0.5*(1.16*uorbu)**2)*u1Lb               ! Feddersen (2000) (for random waves)
         tauc0    = 0.5*rhomean*fc0*umod*u1Lb
         !
         ustar_c  = sqrt(abs(tauc0)/rhomean)                                            ! should correspond to ustbLL above 
         !
         taub     = tauwc*(1d0+1.2d0*(tauw/(tauw+abs(tauwc))**3.2))                     ! Soulsby(1997)
         ustar_cw = sqrt(abs(taub)/rhomean)-ustar_c
      endif

      cfuhiLL   = sqcf*sqcf/hu(Lb)                              ! cfuhiLL   = g / (H.C.C) = (g.K.K) / (A.A)
      cfuhi3D   = cfuhiLL*umod                                  ! cfuhi3D = frc. contr. to diagonal

    if (jawave==0 .or. flowWithoutWaves) then
         z0urou(LL) = z0ucur(LL)                                ! morfo, bedforms, trachytopes
    endif

    if (jawave>0 .and. jawaveStokes >= 1 .and. .not. flowWithoutWaves) then                               ! Ustokes correction at bed
         adve(Lb)  = adve(Lb) - cfuhi3D*ustokes(Lb)
      endif

    else if (ifrctyp == 10) then                                 ! Hydraulically smooth, glass etc
      nit = 0
      u1Lb = u1(Lb)
      umod  = sqrt( u1Lb*u1Lb + v(Lb)*v(Lb) )
      if (jawave>0) then
         call getustwav(LL, z00, umod, fw, ustw2, csw, snw, Dfu, Dfuc, deltau, costu, uorbu) ! get ustar wave squared, fw and wavedirection cosines based upon Swart, ustokes
         !
         if (jawaveStokes >= 1) then
            umod  = sqrt( (u1Lb-ustokes(Lb))*(u1Lb-ustokes(Lb)) + (v(Lb)-vstokes(Lb))*(v(Lb)-vstokes(Lb)) )   ! was ustokes(LL)
         endif
      endif

      r   = umod*hu(Lb)/viskin                                  ! Local re-number:
      r   = max(r,0.001d0)
      er  = e*r
      if (r.lt.rv) then                                         ! Viscous sublayer:
         s   = sqrt(r)
      else

         s   = 12d0                                             ! In log-layer; initial trial for s:
100      continue
         nit = nit+1
         sd  = s
         ers = max(er/sd, 1.0001d0)
         s   = log(ers)/vonkar

         if (nit.ge.nitm) then
            call error ('***ERROR in USTAR: no convergence.', ' ', ' ' )
         endif
         if (s.gt.r) then
            call error ('***ERROR in USTAR: S too large.', ' ', ' ' )
         endif


         if (abs(sd-s).gt.(eps*s)) then
            go to 100                                          ! Convergence criterium:
         endif
      endif

      if (s > 0d0) then
         sqcf = 1d0/s
      else
         sqcf = 0d0
      endif
      ustbLL = sqcf*umod                                        ! ustar based upon bottom layer velocity
      cfuhiLL  = sqcf*sqcf/hu(Lb)
      hdzb   = 0.5d0*hu(Lb)

      if (cfuhiLL > 100d0) then
         nit = nit + 1
      endif

      !     advi(Lb) = advi(Lb) +  cfuhiLL*umod                        ! g / (H.C.C) = (g.K.K) / (A.A) travels in cfuhi
      cfuhi3D = cfuhiLL*umod

   else if (ifrctyp == 11) then                                    ! Noslip

      !    advi(Lb) = advi(Lb) +  2d0*(vicwwu(Lb)+vicouv)/hu(Lb)**2
      cfuhi3D = 2d0*(vicwwu(Lb)+vicoww)/hu(Lb)**2

   endif

   if ( hu(LL) < trsh_u1Lb .and. abs(gsx) > 1d-3 .and. nit <= 3) then
      ! u1Lb = ( u1(Lb)*dti - adve(Lb) - gsx ) / (cfuhi3D + dti)
      u1Lb = ( u1(Lb)*dti            - gsx ) / (cfuhi3D + dti)
      nit  = nit + 1
      goto 10
   endif

   if (jafrculin > 0) then
      cfuhi3D = cfuhi3D + frculin(LL)/hu(Lb)
   endif

   end subroutine getustbcfuhi

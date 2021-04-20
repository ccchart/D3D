!>    transform local spherical coordinates (xloc,yloc) around reference point (xref,yref) to global spherical coordinates (xglob,yglob)
      subroutine loc2spher(xref,yref,N,xloc,yloc,xglob,yglob)
         use m_sferic
         use m_missing, only: dmiss
         use geometry_module, only: sphertocart3D, cart3Dtospher
         implicit none

         double precision,               intent(in)  :: xref,  yref    !< global coordinates of reference point (longitude, latitude)
         integer,                        intent(in)  :: N              !< number of global coordinates
         double precision, dimension(N), intent(in)  :: xloc,  yloc    !< local coordinates
         double precision, dimension(N), intent(out) :: xglob, yglob   !< global coordinates, (longitude, latitude)

         double precision, dimension(3)              :: exxp, eyyp, ezzp   ! base vectors of rotated 3D Cartesian reference frame

         double precision                            :: xx, yy, zz     !  3D Cartesian coordinates
         double precision                            :: xxp, yyp, zzp  !  3D Cartesian coordinates in rotated frame

         double precision                            :: phi0, lambda0

         integer                                     :: i

         if ( jsferic.eq.0 .or. jasfer3D.eq.0 ) then
            do i=1,N
               xglob(i) = xloc(i)+xref
               yglob(i) = yloc(i)+yref
            end do

         else
            phi0 = yref*dg2rd
            lambda0 = xref*dg2rd

!           compute base vectors
            exxp = (/  cos(phi0) * cos(lambda0),  cos(phi0) * sin(lambda0), sin(phi0) /)
            eyyp = (/             -sin(lambda0),              cos(lambda0), 0d0       /)
            ezzp = (/ -sin(phi0) * cos(lambda0), -sin(phi0) * sin(lambda0), cos(phi0) /)

            do i=1,N
!              get 3D-coordinates in rotated frame
               call sphertocart3D(xloc(i),yloc(i),xxp,yyp,zzp)

!              project to fixed frame
!               xxp = exxp(1) * xx + exxp(2) * yy + exxp(3) * zz
!               yyp = eyyp(1) * xx + eyyp(2) * yy + eyyp(3) * zz
!               zzp = ezzp(1) * xx + ezzp(2) * yy + ezzp(3) * zz

               xx = exxp(1) * xxp + eyyp(1) * yyp + ezzp(1) * zzp
               yy = exxp(2) * xxp + eyyp(2) * yyp + ezzp(2) * zzp
               zz = exxp(3) * xxp + eyyp(3) * yyp + ezzp(3) * zzp

!              tranform to spherical coordinates
               call Cart3Dtospher(xx,yy,zz,xglob(i),yglob(i),xref)
            end do

         end if

         return
      end subroutine loc2spher

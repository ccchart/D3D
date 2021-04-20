      SUBROUTINE dLINEDIS2(X3,Y3,X1,Y1,X2,Y2,JA,DIS,XN,YN,rl)
      use m_sferic
      use geometry_module, only: getdx, getdy, dbdistance, sphertoCart3D, cart3Dtospher
      use m_missing, only: dmiss

      implicit none
      integer          :: ja
      DOUBLE PRECISION :: X1,Y1,X2,Y2,X3,Y3,DIS,XN,YN,ZN, d2
      DOUBLE PRECISION :: R2,RL,X21,Y21,Z21,X31,Y31,Z31
      DOUBLE PRECISION :: xx1,xx2,xx3,yy1,yy2,yy3,zz1,zz2,zz3,xxn,yyn,zzn
      !     korste afstand tot lijnelement

      JA  = 0

      if (jsferic == 0 .or. jasfer3D == 0) then

         X21 = getdx(x1,y1,x2,y2,jsferic)
         Y21 = getdy(x1,y1,x2,y2,jsferic)
         X31 = getdx(x1,y1,x3,y3,jsferic)
         Y31 = getdy(x1,y1,x3,y3,jsferic)
         R2  = dbdistance(x2,y2,x1,y1,jsferic, jasfer3D, dmiss)
         R2  = R2*R2
         IF (R2 .NE. 0) THEN
            RL  = (X31*X21 + Y31*Y21) / R2
            IF (0d0 .LE. RL .AND. RL .LE. 1d0) then
               JA = 1
            endif
            XN  = X1 + RL*(x2-x1)
            YN  = Y1 + RL*(y2-y1)
            DIS = dbdistance(x3,y3,xn,yn,jsferic, jasfer3D, dmiss)
         else
            DIS = dbdistance(x3,y3,x1,y1,jsferic, jasfer3D, dmiss)
         ENDIF

      else

         call sphertocart3D(x1,y1,xx1,yy1,zz1)
         call sphertocart3D(x2,y2,xx2,yy2,zz2)
         call sphertocart3D(x3,y3,xx3,yy3,zz3)

         x21 = xx2-xx1
         y21 = yy2-yy1
         z21 = zz2-zz1
         x31 = xx3-xx1
         y31 = yy3-yy1
         z31 = zz3-zz1

         r2  = x21*x21 + y21*y21 + z21*z21
         if (r2 .ne. 0d0) then
            RL = (X31*X21 + Y31*Y21 + Z31*Z21) / R2
            IF (0d0 .LE. RL .AND. RL .LE. 1d0) then
               JA = 1
            endif
            XXN  = xx1 + RL*x21
            YYN  = yy1 + RL*y21
            ZZN  = zz1 + RL*z21
            x31 = xxn-xx3
            y31 = yyn-yy3
            z31 = zzn-zz3
            DIS = sqrt(x31*x31 + y31*y31 + z31*z31)

            call Cart3Dtospher(xxn,yyn,zzn,xn,yn,maxval((/x1,x2,x3/)))
         else
            DIS = dbdistance(x3,y3,x1,y1, jsferic, jasfer3D, dmiss)
         endif

      endif

      RETURN
      END subroutine DLINEDIS2

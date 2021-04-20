    SUBROUTINE ISPOIN(      X,      Y,     mmax, nmax, MC,     NC,   RD1, &
                            XL,     YL,     MV,     NV)
      use m_missing
      use m_wearelt
      implicit none
      integer, intent(in)  :: mmax, nmax, mc, nc
      integer, intent(out) :: mv, nv
      double precision ::   X(MMAX,NMAX), Y(MMAX,NMAX), RD1(MMAX,NMAX)
      double precision :: xl, yl


      integer :: m1, n1, m2, n2, ishot, mvol, nvol, i, j

      DATA MVOL /0/, NVOL /0/
      MV    = 0
      NV    = 0
      ISHOT = 0

  666 CONTINUE
      IF (ISHOT .EQ. 0 .AND. MVOL .NE. 0) THEN
         M1    = MAX(1,MVOL - 3)
         N1    = MAX(1,NVOL - 3)
         M2    = MIN(MC,MVOL + 3)
         N2    = MIN(NC,NVOL + 3)
         ISHOT = 1
      ELSE
         M1    = 1
         N1    = 1
         M2    = MC
         N2    = NC
         ISHOT = 0
      ENDIF

      DO 10 I = M1,M2
         DO 10 J = N1,N2
            IF (X(I,J) .NE. XYMIS) THEN
               IF (ABS(XL - X(I,J)) .LT. RCIR) THEN
                  IF (ABS(YL - Y(I,J)) .LT. RCIR) THEN
                     MV   = I
                     NV   = J
                     XL   = X(I,J)
                     YL   = Y(I,J)
                     MVOL = MV
                     NVOL = NV
                     CALL DISVAL(MV,NV,RD1(MV,NV))
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
   10 CONTINUE
      IF (ISHOT .EQ. 1) GOTO 666
      MVOL = 0
      CALL DISVAL(0,0,0d0)
      RETURN
      END subroutine ispoin

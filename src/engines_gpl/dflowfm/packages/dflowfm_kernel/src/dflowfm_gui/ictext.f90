      SUBROUTINE ICTEXT(TEX,NX,NY,NCOL)
      use unstruc_colors
      implicit none
      integer :: l
      integer :: ncol
      integer :: nx
      integer :: ny
      double precision :: x
      double precision :: y
!     grafische tekst op normale text posities
      CHARACTER TEX*(*)
      X = X1 + (X2-X1)*dble(NX)/dble(IWS)
      Y = Y2 + (Y1-Y2)*dble(NY)/dble(IHS)
      CALL SETCOL(NCOL)
      L = len_trim(TEX)
      CALL DRAWTEXT(real(X),real(Y),TEX(1:L))
      RETURN
      END

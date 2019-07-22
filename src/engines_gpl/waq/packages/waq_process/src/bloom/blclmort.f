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

!     This file contains the subroutines concerning the chlorinity
!     dependend mortality rates:
!     BLCLST - adapt the rates
!     BLCLRS - reset the rates
!     And a routine to set ppmax for a specific alg
!     BLSPPM - set ppmax

      subroutine blclst (mrtm1,mrtm2,mrtb1,mrtb2,ntyp_a,cl)

      use bloom_data_dim
      use bloom_data_size 

      implicit none

!     Arguments
!
!     Name    Type  Length   I/O  Description
!
!     MRTM1   R     NTYP_A   O    Original mortality rates
!     MRTM2   R     NTYP_A   I    M2 mort rate coeff
!     MRTB1   R     NTYP_A   I    B1 mort rate sal stress coeff
!     MRTB2   R     NTYP_A   I    B2 mort rate sal stress coeff
!
      INTEGER NTYP_A
      REAL    MRTM1(NTYP_A),MRTM2(NTYP_A),MRTB1(NTYP_A),MRTB2(NTYP_A),CL
      INTEGER IALG

!     Loop over algae types
      DO 10 IALG = 1, NTYP_A
!       Store the original value
        MRTM1(IALG) = RMORT1(IALG)
!       Salinity dep. mortality ??
        IF (MRTM2(IALG).GT.0.) THEN
          CL = MIN(CL,35000.)
          RMORT1(IALG) =  (MRTM2(IALG)-MRTM1(IALG))/
     1      (1.+EXP(MRTB1(IALG)*(CL-MRTB2(IALG))))+MRTM1(IALG)
        ENDIF
   10 CONTINUE

      RETURN

      END
      SUBROUTINE BLCLRS (MRTM1,NTYP_A)

      use bloom_data_dim
      use bloom_data_size 

      implicit none

!     Arguments
!
!     Name    Type  Length   I/O  Description
!
!     MRTM1   R     NTYP_A   O    Original mortality rates
!
      INTEGER NTYP_A
      REAL    MRTM1(NTYP_A)
!
!     Local variables
!
!     Name    Type  Length   I/O  Description
!
!     IALG    I     1             counter over algae types

      INTEGER      IALG

!     Loop over algae types
      DO 20 IALG = 1, NTYP_A
!       Store the original value
        RMORT1(IALG)= MRTM1(IALG)
   20 CONTINUE

      RETURN

      END
      SUBROUTINE BLSPPM (IALG  , PPMAX )

      use bloom_data_dim
      use bloom_data_size 
      
      implicit none

!
!     Set PPMAX in BLOOM array for specific alg
!
!     Arguments
!
!     Name    Type  Length   I/O  Description
!
!     IALG    1     1        I    index alg involved
!     PPMAX   R     1        I    PPMAX value to be set
!
      INTEGER IALG
      REAL    PPMAX
!
!      INCLUDE 'blmdim.inc'
!      INCLUDE 'size.inc'
!
      PMAX1(IALG) = PPMAX
!
      RETURN
      END
      SUBROUTINE BLSSDM (IALG  , SDMIXN )

      use bloom_data_dim
      use bloom_data_size 
      
      implicit none

!
!     Set SDMIX in BLOOM array for specific alg
!
!     Arguments
!
!     Name    Type  Length   I/O  Description
!
!     IALG    1     1        I    index alg involved
!     SDMIXN  R     1        I    SDMIX value to be set
!
      INTEGER IALG
      REAL    SDMIXN
!
!      INCLUDE 'blmdim.inc'
!      INCLUDE 'size.inc'
!
      SDMIX(IALG) = SDMIXN
!
      RETURN
      END
      SUBROUTINE BLSAEF (IGROUP  , EFFIN )

      use bloom_data_dim
      use bloom_data_size 
      use bloom_data_phyt
      
      implicit none
!
!     Set AVEFFI in BLOOM array for specific alg
!
!     Arguments
!
!     Name    Type  Length   I/O  Description
!
!     IALG    1     1        I    index alg involved
!     EFFIN   R     1        I    EFFI value to be set
!
      INTEGER IGROUP, ITYPE
      REAL    EFFIN
!
!      include 'blmdim.inc'
!      include 'size.inc'
!      include 'phyt2.inc'
!
      DO ITYPE = IT2(IGROUP,1),IT2(IGROUP,2)
         AVEFFI(ITYPE) = EFFIN
      END DO
!
      RETURN
      END

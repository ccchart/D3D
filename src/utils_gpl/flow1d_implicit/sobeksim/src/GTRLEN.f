      FUNCTION GTRLEN (PNTNAM)

c=======================================================================
c            Rijkswaterstaat/RIZA and DELFT HYDRAULICS
c                One Dimensional Modelling System
c                           S O B E K
c-----------------------------------------------------------------------
c Subsystem:          Memory Management Module
c
c Programmer:         A. Hoekstra / S.L. van der Woude
c
c Module:             GTRLEN (GeT Real variable LENgth)
c
c Module description: Report length of variable in RPool
c
c                     This function reports the length of a real vari-
c                     able allocated in the RPool
c
c-----------------------------------------------------------------------
c Parameters:
c NR NAME              IO DESCRIPTION
c  0 gtrlen            O  Length of real variable
c  1 pntnam            P  -
c-----------------------------------------------------------------------
c Subprogram calls:
c NAME    DESCRIPTION
c zzgtln  ZZ GeT LeNgth
c=======================================================================
c
c***********************************************************************
c CVS log information:
c
c $Id$
c
c History:
c $Log: gtrlen.pf,v $
c Revision 1.2  1995/05/30  07:03:32  hoeks_a
c file changed from dos to ux
c
c Revision 1.1  1995/04/13  07:08:45  hoeks_a
c Initial check-in
c
c Revision 1.1.1.1  1993/07/21  14:44:02  kuipe_j
c Initial version
c
c
c***********************************************************************
c
      INTEGER  GTRLEN

      CHARACTER PNTNAM*(*)

      include '..\include\pointrs.i'

      INTEGER   ZZGTLN
      EXTERNAL  ZZGTLN

      GTRLEN = ZZGTLN (PNTNAM, RPNTRS, RLENGT, NRPNTR, ROLDPT)

      END

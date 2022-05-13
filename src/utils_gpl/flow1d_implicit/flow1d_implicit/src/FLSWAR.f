      subroutine FLSWAR (istru  ,strpar ,sign   ,zs     ,wstr   ,
     +                   cw     ,slim   ,itab   )

c=======================================================================
c            Rijkswaterstaat/RIZA and DELFT HYDRAULICS
c                One Dimensional Modelling System
c                           S O B E K
c-----------------------------------------------------------------------
c Subsystem:          Flow Module
c
c Programmer:         J.Brouwer
c
c Module:             FLSWAR (FLow get Simple Weir ARguments)
c
c Module description: Parameters for a simple weir are extracted from
c                     the packed array strpar.
c
c-----------------------------------------------------------------------
c Parameters:
c NR NAME              IO DESCRIPTION
c  6 cw                O  Correction coefficient for weir flow.
c  1 istru             I  Number of structure.
c  8 itab              O  Table number.
c  3 sign              I  Flow direction (+/-).
c  7 slim              O  Submergence limit.
c  2 strpar(21,nstru)  I  Each structure is characterized by a number of
c                         specific parameters. strpar (i,j) = parameter
c                         i of structure j:
c                         - Simple weir:
c                         (1,j) = Crest height Zs.
c                         (2,j) = Crest width Ws.
c                              Positive flow:
c                         (3,j) = Correction coefficient cw.
c                         (4,j) = Submergence limit Slim.
c                         (5,j) = Table pointer for drowned reduction
c                                 curve f(h2/h1).
c                              Negative flow:
c                         (6,j) = Correction coefficient cw.
c                         (7,j) = Submergence limit Slim.
c                         (8,j) = Table pointer for drowned reduction
c                                 curve f(h2/h1).
c                         - Advanced weir:
c                         (1,j) = Crest height Zs.
c                         (2,j) = Total net width Wn.
c                         (3,j) = Number of piers N.
c                              Positive flow:
c                         (4,j) = Heigth of upstream face P.
c                         (5,j) = Design head H0 of the weir.
c                         (6,j) = Pier contraction coefficient Kp.
c                         (7,j) = Abutment contraction coefficient Ka.
c                              Negative flow:
c                         (8,j) = Heigth of upstream face P.
c                         (9,j) = Design head H0 of the weir.
c                         (10,j)= Pier contraction coefficient Kp.
c                         (11,j)= Abutment contraction coefficient Ka.
c                         - Pump:
c                         (1,j) = Control direction:
c                                 cpmpup (-1) : upward control
c                                 cpmpdw (+1) : downward control
c                         (2,j) = Table pointer for pump capacitity re-
c                                 duction factor.
c                         (3,j) = Capacity.
c                         (4,j) = Water level which starts pump.
c                         (5,j) = Water level which stops pump.
c                         - General structure:
c                         (1,j) = Width left side of structure W1.
c                         (2,j) = Bed level left side of structure Zb1.
c                         (3,j) = Width structure left side Wsdl.
c                         (4,j) = Bed left side of structure Zbsl.
c                         (5,j) = Width structure centre Ws.
c                         (6,j) = Bed level centre Zs.
c                         (7,j) = Width structure right side Wsdr.
c                         (8,j) = Bed right side of structure Zbsr.
c                         (9,j) = Width right side of structure W2.
c                         (10,j)= Bed level right side of structure Zb2.
c                         (11,j)= Gate opening heigth dg.
c                              Positive flow:
c                         (12,j)= Correction coefficient for free gate
c                                 flow cgf.
c                         (13,j)= Correction coefficient for drowned
c                                 gate flow cgd.
c                         (14,j)= Correction coefficient for free weir
c                                 flow cwf.
c                         (15,j)= Correction coefficient for drowned
c                                 weir flow cwd.
c                         (16,j)= Contraction coefficient for free gate
c                                 flow MU-gf.
c                              Negative flow:
c                         (17,j)= Correction coefficient for free gate
c                                 flow cgf.
c                         (18,j)= Correction coefficient for drowned
c                                 gate flow cgd.
c                         (19,j)= Correction coefficient for free weir
c                                 flow cwf.
c                         (20,j)= Correction coefficient for drowned
c                                 weir flow cwd.
c                         (21,j)= Contraction coefficient for free gate
c                                 flow MU-gf.
c  5 wstr              O  Width at centre of structure.
c  4 zs                O  Bed level at centre of structure.
c=======================================================================
c
c
c***********************************************************************
c CVS log information:
c
c $Id$
c
c History:
c $Log: flswar.pf,v $
c Revision 1.5  1999/03/15  15:50:52  kuipe_j
c tabs removed
c
c Revision 1.4  1995/09/22  10:02:19  kuipe_j
c variable dimensions, new headers
c
c Revision 1.3  1995/05/30  09:55:31  hoeks_a
c Minor changes
c
c Revision 1.2  1995/05/30  06:59:32  hoeks_a
c file converted from dos to ux
c
c Revision 1.1  1995/04/13  07:08:11  hoeks_a
c Initial check-in
c
c Revision 1.2  1993/11/26  15:31:39  kuipe_j
c Update after finishing Sobeksel.
c
c Revision 1.1.1.1  1993/07/21  14:43:55  kuipe_j
c Initial version
c
c
c***********************************************************************
c
c     Include constants for array dimensions
c
      include '..\include\sobdim.i'
c
c     Declaration of parameters:
c
      integer istru, itab
      real    strpar(dmstrpar,*), sign, zs, wstr, cw, slim
c
c     Determine crest height and crest width
c
      zs   = strpar(1,istru)
      wstr = strpar(2,istru)
c
c     Determine Cw, Slim and table pointer
c     (flow direction dependent)
c
      if ( sign .gt. 0.0 ) then
         cw   = strpar(3,istru)
         slim = strpar(4,istru)
         itab = int(strpar(5,istru))
      else
         cw   = strpar(6,istru)
         slim = strpar(7,istru)
         itab = int(strpar(8,istru))
      endif
c
      end

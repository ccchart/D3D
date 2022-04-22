      subroutine wqexad (ngrid , nbran , nnode , nbrnod,
     j                   nsegmt, npntr , nsegtb, nexdef,
     j                   branch, brnode, exdef , pntr  ,
     j                   segmnt, segtab, segcon, grdcon,
     j                   dircon)
c=======================================================================
c            Rijkswaterstaat/RIZA and DELFT HYDRAULICS
c                One Dimensional Modelling System
c                           S O B E K
c-----------------------------------------------------------------------
c Subsystem:          Water Quality Interface Module
c
c Programmer:         J.A.G. van Gils
c
c Module:             WQEXAD (Water Quality EXchanges ADministration)
c
c Module description: This module creates the administration of exchanges
c                     for use in de water quality model DELWAQ. This
c                     module sets NPNTR, PNTR, NEXDEF and EXDEF. This
c                     task was originally executed by the UI. The input
c                     information is the mapping of water quality
c                     segments on Working Units: in SEGMNT and SEGTAB.
c                     These arrays are constructed by the network editor
c                     in the user interface. The relation between
c                     boundary numbers and gridpoints is stored in the
c                     array SEGMNT as well.
c
c-----------------------------------------------------------------------
c Parameters:
c NR NAME              IO DESCRIPTION
c  1 ngrid             I  nr of gridpoints
c  2 nbran             I  nr of branches
c  3 nnode             I  nr of nodes
c  4 nbrnod            I  max nr of branches in a node
c  5 nsegmt            I  nr of entries in segmnt (segments + boundaries)
c  6 npntr             I  nr of exchanges
c  7 nsegtb            I  nr of entries in segtab
c  8 nexdef            I  nr of entries in exdef
c  9 branch(4,nbran)   I  Branch information:
c                         (1,i) = Node number n1 at begin of branch i.
c                         (2,i) = Node number n2 at end of branch i.
c                         (3,i) = Grid point i1 at begin of branch i.
c                         (4,i) = Grid point i2 at end of branch i.
c 10 brnode(nbrnod+1,  I  Node-Branch relation table. The first index
c        ,nnode)          contains the number of connected branches
c                         (index 1) for each node. The second index
c                         contains the first connected branch number
c                         etc.
c 11 exdef(6,nexdef)   O  This table contains the elementary exchange
c                         definitions. Each exchange between segments
c                         can be described by one or more elementary
c                         definitions in this table.
c                         (1,i) = Exchange type:
c                                 cexigp (1) : Exchange in gridpoint.
c                                 cexbgp (2) : Exchange between gpoints.
c                                 cexind (3) : Exchange in a node.
c                                 cexqlt (4) : Exchange from Qlat to seg
c                         - Exchange in a gridpoint (Type 1):
c                         (2,i) = Gridpoint
c                         (3,i) = Section from
c                         (4,i) = Section to
c                         (5,i) = Direction of "from" --> "to":
c                                 cpopnt (+1) : pos branch direction.
c                                 cnepnt (-1) : neg branch direction.
c                         - Exchange between gridpoints (Type 2):
c                         (2,i) = Gridpoint from
c                         (3,i) = Gridpoint to
c                         (4,i) = Section
c                         (5,i) = Length factor (0< length factor <1)
c                         - Exchange in a node (Type 3):
c                         (2,i) = Node number
c                         (3,i) = Gridpoint from
c                         (4,i) = Gridpoint to
c                         (5,i) = Section from
c                         (6,i) = Section to
c                         - Exchange from Qlat stat to segment (Type 4)
c                         (2,i) = Gridpoint
c                         (3,i) = Length factor (0 < factor <= 1).
c                         (4,i) = Lateral station number.
c 12 pntr(4,npntr)     O  Definition of the pointer table. From here it
c                         is possible to find the exchanges between the
c                         segments and the starting location in the
c                         exdef array.
c                         (1,j) = From segment number.
c                         (2,j) = To segment number.
c                         (3,j) = Pointer to exchange table exdef.
c                         (4,j) = Number of exchange definitions in ex-
c                                 def.
c 13 segmnt(3,nsegmt)  I  Segment definition:
c                         (1,i) = segment number (if > 0)
c                         (2,i) = first entry in segtab for this segment
c                         (3,i) = nr.of entries in segtab for this
c                                 segment
c                         (1,i) = negative boundary number (if < 0)
c                         (2,i) = grdpoint where boundary is situated
c                         (3,i) = should have value 0!
c 14 segtab            I
c 15 segtab(5,nsegtb)  I  This table contains for each segment the en-
c                         closed grid-cells together with the length
c                         factors and section indication.
c                         (1,j) = Gridpoint 1
c                         (2,j) = Gridpoint 2
c                         (3,j) = Length factor Lb
c                         (4,j) = Length factor Le
c                         (5,j) = Section
c                                 cnopar (0) : No parallel sections
c                                 cmainc (1) : Main channel
c                                 csub1  (2) : Sub section 1
c                                 csub2  (3) : Sub section 2
c 16 segcon            O  Auxilliary array with all segments at a node.
c 17 grdcon            O  Auxilliary array with all gridpoints at a node
c 18 dircon            O  Auxilliary array with all directions at a node
c-----------------------------------------------------------------------
c Subprogram calls:
c NAME    DESCRIPTION
c wqtnod  get Node number for a gridpoint
c wqtcns  get Connected Segments for a grid point
c wqtbou  get Boundary gridpoint
c=======================================================================
c
c
c
c***********************************************************************
c CVS log information:
c
c $Id$
c
c History:
c $Log: wqexad.pf,v $
c Revision 1.2  1999/03/12  12:34:08  kuipe_j
c parallel segments added
c
c Revision 1.1  1996/10/31  09:51:43  kuipe_j
c Calculation of exchanges added
c
c
c***********************************************************************
c
c       Parameters
c

      integer      ngrid ,
     j             nbran ,
     j             nnode ,
     j             nbrnod,
     j             nsegmt,
     j             npntr ,
     j             nsegtb,
     j             nexdef
      integer      branch(4,nbran),
     j             segmnt(3,nsegmt),
     j             pntr(4,npntr),
     j             brnode(nbrnod+1,nnode)
      real         exdef(6,nexdef),
     j             segtab(5,nsegtb)

c     Locals

      integer      igrid , inode , isegm , isgtb , jsegm , jsgtb ,
     j             nsgtb , msgtb , isgtb0, jsgtb0, nsegto, nsegfr,
     j             igrdfr, igrdto, iconn , iboun , icon2 , nconn ,
     j             iexdef, ipntr , ibran , ksecfr, ksecto, mainsc,
     j             i     , j     , ksec  , ksec1 , ksec2 , pntr1 ,
     j             pntr2 , nsegno, isegto, isegfr
      integer      segcon(nbrnod),
     j             grdcon(nbrnod),
     j             dircon(nbrnod)
      real         dista

      include '..\include\sobcon.i'


c**********************************************************************
c     Exchanges of type 1: exchanges in gridpoints (NOT KNODES)
c**********************************************************************

      iexdef = 0
      ipntr  = 0

      do 100 igrid = 1,ngrid

c        Only if this is not a knode!

         call wqtnod ( igrid , inode , branch , nbran )
         if ( inode .le. 0 ) then

c           find segments connected to igrid:
c             nsegfr is segment which includes WU ending in igrid
c             nsegto is segment which includes WU starting in igrid

            call wqtcns ( igrid , nsegmt, isegto, isegfr,
     j                    nsegtb, segmnt, segtab)

            nsegfr = segmnt(1,isegfr)
            nsegto = segmnt(1,isegto)

c           two different connected segments found??

            if ( nsegfr .ne. 0 .and. nsegto .ne. 0 .and.
     j           nsegfr .ne. nsegto ) then

c              yes! create exchange of type 1
c
c              Parallel segments have adjacent numbers
c              (Main=j, Sub1=j+1, Sub2=j+2)
c              This will be used by making exchange definitions
c
               ksecfr = nint(segtab(5,segmnt(2,isegfr)))
               ksecto = nint(segtab(5,segmnt(2,isegto)))
               mainsc = max(ksecfr,ksecto)
               ksec   = max(mainsc*3,1)
               do i=1,ksec
                  if (i .gt. 1) then
                    if (ksecfr.gt.0) then
                      isegfr = isegfr + 1
                      ksecfr = nint(segtab(5,segmnt(2,isegfr)))
                    endif
                    if (ksecto.gt.0) then
                      isegto = isegto + 1
                      ksecto = nint(segtab(5,segmnt(2,isegto)))
                    endif
                  endif

                  ipntr = ipntr + 1
                  pntr(1,ipntr) = segmnt(1,isegfr)
                  pntr(2,ipntr) = segmnt(1,isegto)
                  pntr(3,ipntr) = ipntr
                  pntr(4,ipntr) = 1

                  iexdef = iexdef + 1
                  exdef(1,iexdef) = real(cexigp)
                  exdef(2,iexdef) = real(igrid)
c                  exdef(3,iexdef) = real(max(ksecfr,mainsc))
c                  exdef(4,iexdef) = real(max(ksecto,mainsc))
                  exdef(3,iexdef) = real(ksecfr)
                  exdef(4,iexdef) = real(ksecto)

                  exdef(5,iexdef) = 1.
               enddo
            endif

         endif
  100 continue

c**********************************************************************
c     Exchanges of type 2: exchanges between gridpoints
c**********************************************************************

c     Loop over segments

      do 200 isegm = 1,nsegmt
        isgtb0 = segmnt(2,isegm)
        nsegto = segmnt(1,isegm)
c
c       Search only for the main segment, as sub segments
c       do have adjacent numbers
c

        if (nsegto .gt. 0 .and. nint(segtab(5,isgtb0)) .le. 1) then
          nsgtb  = segmnt(3,isegm)

c         Loop over grid cells (working units) within the segment

          do 150 isgtb = isgtb0,isgtb0+nsgtb-1
            dista = segtab(3,isgtb)
            igrdfr = nint(segtab(1,isgtb))
            igrdto = nint(segtab(2,isgtb))

c           Is this a grid cell cut in two pieces?

            if ( dista .gt. 0.0 ) then

c              Yes! We will look for the other piece!
c              Second loop over segments

               do 130 jsegm = 1,nsegmt
                 jsgtb0 = segmnt(2,jsegm)
                 nsegfr = segmnt(1,jsegm)
c
c                Search only for the main segment, as sub segments
c                do have adjacent numbers
c
                 if (nsegfr .gt. 0 .and.
     +               nint(segtab(5,jsgtb0)) .le. 1) then
                   msgtb  = segmnt(3,jsegm)

c                  Second loop over grid cells within the segment

                   do 120 jsgtb = jsgtb0,jsgtb0+msgtb-1
                     if ( abs(segtab(4,jsgtb)-dista) .lt. 1e-10 .and.
     j               nint(segtab(1,jsgtb)).eq. igrdfr .and.
     j               nint(segtab(2,jsgtb)).eq. igrdto ) then

c                       This is the other piece! create exchange of type 2

c                       Parallel segments have adjacent numbers
c                       (Main=j, Sub1=j+1, Sub2=j+2)
c                       This will be used by making exchange definitions
c
                        isegfr = jsegm
                        isegto = isegm
                        ksecfr = nint(segtab(5,segmnt(2,isegfr)))
                        ksecto = nint(segtab(5,segmnt(2,isegto)))
                        mainsc = max(ksecfr,ksecto)
                        ksec   = max(mainsc*3,1)
                        do i=1,ksec
                          if (i .gt. 1) then
                             if (ksecfr.gt.0) then
                               isegfr = isegfr + 1
                               ksecfr = nint(segtab(5,segmnt(2,isegfr)))
                             endif
                             if (ksecto.gt.0) then
                               isegto = isegto + 1
                               ksecto = nint(segtab(5,segmnt(2,isegto)))
                             endif
                           endif
                           ipntr = ipntr + 1
                           pntr(1,ipntr) = segmnt(1,isegfr)
                           pntr(2,ipntr) = segmnt(1,isegto)
                           pntr(3,ipntr) = ipntr
                           pntr(4,ipntr) = 1

                           iexdef = iexdef + 1
                           exdef(1,iexdef)= real(cexbgp)
                           exdef(2,iexdef)= real(igrdfr)
                           exdef(3,iexdef)= real(igrdto)
                           exdef(4,iexdef)= real(max(ksecfr,ksecto))
                           exdef(5,iexdef)= dista
                        enddo
                        goto 150
                     endif
  120              continue
                 endif
  130          continue
            endif
  150     continue
        endif
  200 continue

c**********************************************************************
c     Exchanges of type 3: exchanges in knodes
c     Includes exchanges in boundaries, which will be defined as type 1
c**********************************************************************

c     we will determine for every node which segments and gridcells
c     are connected to this knode, and whether the corresponding branch
c     starts or ends in the node
c     for this purpose we use the local arrays segcon, dircon and grdcon
c     with dimension nbrnod

c     loop over nodes

      do 300 inode = 1,nnode

c        zero local arrays

         do 240 iconn = 1,nbrnod
            segcon(iconn) = 0
            dircon(iconn) = 0
            grdcon(iconn) = 0
  240    continue

c        loop over branches connected to this node

         nconn = brnode(1,inode)
         do 250 iconn = 1,nconn
            ibran = brnode(1+iconn,inode)

c           find grid point on this branch in the node
c           and set dircon: = 1: branch starts, = -1: branch stops

            if ( branch(1,ibran) .eq. inode ) then
               grdcon(iconn) = branch(3,ibran)
               dircon(iconn) = 1
            else
               grdcon(iconn) = branch(4,ibran)
               dircon(iconn) = -1
            endif
            igrid = grdcon(iconn)

c           find connected segment

            call wqtcns ( igrid , nsegmt, isegto, isegfr,
     j                    nsegtb, segmnt, segtab)

            if ( isegfr .ne. 0 ) then
                segcon(iconn) = isegfr
            elseif ( isegto .ne. 0 ) then
                segcon(iconn) = isegto
            endif
  250    continue

         if ( nconn .eq. 1 ) then

c           create boundary: exchange of type 1

c           find boundary number at the gridpoint grdcon(1)

            call wqtbou ( iboun , segmnt , nsegmt , grdcon(1) )

c           Parallel segments have adjacent numbers
c           (Main=j, Sub1=j+1, Sub2=j+2)
c           This will be used by making exchange definitions

            mainsc = nint(segtab(5,segmnt(2,segcon(1))))
            ksec  = max(mainsc*3,1)
            if ( dircon(1) .eq. 1 ) then
               pntr1  = iboun
               pntr2  = segmnt(1,segcon(1))
               ksecfr = 0
               ksecto = mainsc
            else
               pntr1  = segmnt(1,segcon(1))
               pntr2  = iboun
               ksecfr = mainsc
               ksecto = 0
            endif

            do i=1,ksec
               ipntr = ipntr + 1
               pntr(1,ipntr) = pntr1
               pntr(2,ipntr) = pntr2
               pntr(3,ipntr) = ipntr
               pntr(4,ipntr) = 1

               iexdef = iexdef + 1
               exdef(1,iexdef)= real(cexigp)
               exdef(2,iexdef)= real(grdcon(1))
               exdef(3,iexdef)= real(ksecfr)
               exdef(4,iexdef)= real(ksecto)
               exdef(5,iexdef)= 1.0
               if ( dircon(1) .eq. 1 ) then
                  ksecto = i + 1
                  pntr2  = pntr2 + 1
               else
                  ksecfr = i + 1
                  pntr1  = pntr1 + 1
               endif
            enddo
         else

c           loop over every combination of branches

            do 270 iconn = 1,nconn
               do 260 icon2 = iconn+1,nconn

                  isegfr = segcon(iconn)
                  isegto = segcon(icon2)
                  igrdfr = grdcon(iconn)
                  igrdto = grdcon(icon2)
                  nsegfr = segmnt(1,isegfr)
                  nsegto = segmnt(1,isegto)

c                 are they in the same segment?

                  if ( nsegfr .ne. nsegto ) then

c                    No! Create exchange of type 3

c                    Parallel segments have adjacent numbers
c                   (Main=j, Sub1=j+1, Sub2=j+2)
c                    This will be used by making exchange definitions
c
                     ksecfr = nint(segtab(5,segmnt(2,isegfr)))
                     ksecto = nint(segtab(5,segmnt(2,isegto)))
                     mainsc = max(ksecfr,ksecto)
                     ksec1  = max(ksecfr*3,1)
                     ksec2  = max(ksecto*3,1)
                     do i=1,ksec1
                        if (ksecfr.gt.0 .and. i.gt.1) then
                           isegfr = isegfr + 1
                           ksecfr = nint(segtab(5,segmnt(2,isegfr)))
                        endif
                        isegto = segcon(icon2)
                        ksecto = nint(segtab(5,segmnt(2,isegto)))
                        do j=1,ksec2
                           if (ksecto.gt.0 .and. j.gt.1) then
                              isegto = isegto + 1
                              ksecto = nint(segtab(5,segmnt(2,isegto)))
                           endif
                           ipntr = ipntr + 1
                           pntr(1,ipntr) = segmnt(1,isegfr)
                           pntr(2,ipntr) = segmnt(1,isegto)
                           pntr(3,ipntr) = ipntr
                           pntr(4,ipntr) = 1

                           iexdef = iexdef + 1
                           exdef(1,iexdef)= real(cexind)
                           exdef(2,iexdef)= real(inode)
                           exdef(3,iexdef)= real(igrdfr)
                           exdef(4,iexdef)= real(igrdto)
                           exdef(5,iexdef)= real(ksecfr)
                           exdef(6,iexdef)= real(ksecto)
                        enddo
                     enddo
                  endif
  260          continue
  270       continue

         endif

  300 continue

c**********************************************************************
c     Exchanges of type 5: exchanges between parallel segments
c     Parallel segments have adjacent numbers
c     (Main=j, Sub1=j+1, Sub2=j+2)
c     This will be used by making exchange definitions
c**********************************************************************

      do 400 isegm = 1,nsegmt
        isgtb0 = segmnt(2,isegm)
        nsegno = segmnt(1,isegm)
c
c       Search only for the main segment, as sub segments
c       do have adjacent numbers
c
        if (nsegno .gt. 0) then
           ksec = nint(segtab(5,isgtb0))
           if (ksec .gt. 1) then
              do i=2,ksec
                 ipntr = ipntr + 1
c                segment numbers
                 pntr(1,ipntr) = nsegno - ksec - 1 + i
                 pntr(2,ipntr) = nsegno
                 pntr(3,ipntr) = ipntr
                 pntr(4,ipntr) = 1

                 iexdef = iexdef + 1
                 exdef(1,iexdef)= real(cexbsc)
c                segment indices
                 exdef(2,iexdef)= real(isegm - ksec - 1 + i)
                 exdef(3,iexdef)= real(isegm)
                 exdef(4,iexdef)= real(i - 1)
                 exdef(5,iexdef)= real(ksec)
               enddo
           endif
        endif
 400  continue
c
c      WRITE (11,*)  'EXDEF'
c      DO j=1,nexdef
c         WRITE (11,'(i4,6f8.0)') j,(exdef(i,j),i=1,6)
c      ENDDO
c      WRITE (11,*)  'PNTR'
c      DO j=1,npntr
c         WRITE (11,'(5i6)') j,(pntr(i,j),i=1,4)
c      ENDDO


      end

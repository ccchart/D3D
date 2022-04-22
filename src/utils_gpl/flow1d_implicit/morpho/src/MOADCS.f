      subroutine moadcs ( igp,     isec,   ibrtyp,
     +                    deltaa,  time,   moropt,  
     +                    nboun,   ngrid,  nnode,
     +                    branch,  ibr,    node,   mbdpar,
     +                    maxtab,  ntabm,  ntab,   table,
     +                    maxlev,  nlev,   hlev,
     +                    wft,    ws,
     +                    flwdir
     +                  )

c=======================================================================
c            Rijkswaterstaat/RIZA and DELFT HYDRAULICS
c                One Dimensional Modelling System
c                           S O B E K
c-----------------------------------------------------------------------
c Subsystem:          Morphology module
c
c Programmer:         S.L. van der Woude
c
c Module:             MOADCS (MORPHology Adapt Cross Sections)
c
c Module description: Adapt cross section dimensions for a cross secti-
c                     on. The following cross section types are possi-
c                     ble:
c
c                     o     sedredge cross sections;
c                     o     tabulated cross sections.
c
c                     Cross sections of type circle can not be proces-
c                     sed.
c
c                     For tabulated cross sections this routine will
c                     call a cross section adaption routine depending on
c                     the users choice. Two choices are possible for
c                     tabulated cross sections:
c
c                     o   erosion/sedimentation spread equally over
c                         transport width;
c                     o   option = erosion/sedimentation proportional to
c                         local depth.
c
c                     If a grid point is a boundary and the boundary
c                     condition is a bed level and sediment is moving in
c                     no delta z is calculated. In that case the bounda-
c                     ry condition is used to adapt the cross section.
c                     In all other cases the calculated delta A is used
c                     to calculate a delta z. For sedredge cross secti-
c                     ons this routine will calculate a delta z and
c                     update the bottom of the section involved. In case
c                     the grid point is a boundary and the boundary
c                     conditions is a bed level and sediment is moving
c                     on no delta z is calculated. The new bottom is
c                     interpolated from the boundary conditions table
c                     and will be assigned to the section involved.
c
c-----------------------------------------------------------------------
c Parameters:
c NR NAME              IO DESCRIPTION
c 11 branch(4,nbran)   I  Branch information:
c                         (1,i) = Node number n1 at begin of branch i.
c                         (2,i) = Node number n2 at end of branch i.
c                         (3,i) = Grid point i1 at begin of branch i.
c                         (4,i) = Grid point i2 at end of branch i.
c  4 deltaa            I  Calculated change in area
c 20 hlev(ngrid,       IO (i,j) = H at level j in cross section i.
c       maxlev)           - For a circle cross section:
c                         (i,1) = Reference level.
c                         - For a sedredge cross section:
c                         (i,1) = Bed level of main section (i.e. left
c                                 channel).
c                         (i,2) = Bed level of sub section 1 (i.e. right
c                                 channel).
c  3 ibrtyp            I  Type of branch
c                           ccrtab (1) : tabulated branch
c                           ccrcir (2) : circle branch
c                           ccrsed (3) : sedredge branch
c  1 igp               I  Gridpoint number
c  2 isec              I  Section number (1 or 2)
c 18 maxlev            I  Maximum+1 over nlev(1:ngrid); needed for
c                         declarations
c 14 maxtab            I  Maximum number of defined tables.
c 13 mbdpar(5,nboun)   I  Morphodynamic boundary conditions:
c                         (1,i) = Type of boundary condition:
c                                 cmbsft (1) : Condition S=f(t).
c                                 cmbsfq (2) : Condition S=f(Q).
c                                 cmbzft (3) : Condition z=f(t).
c                         (2,i) = Location (node number).
c                         (3,i) = Branch number that is connected.
c                         (4,i) = Table pointer for boundary table. In
c                                 case of a connected sedredge branch
c                                 the pointer will be assigned to the
c                                 left channel.
c                         (5,i) = Table pointer for boundary table. In
c                                 case of a connected sedredge branch
c                                 the pointer will be assigned to the
c                                 right channel. In other cases undefi-
c                                 ned.
c  6 moropt            I  Method of adapting cross sections
c                         ceqows (1) : Equally over the actual sediment
c                                      transp. width of cross section
c                         cprodp (2) : Proportional to the local water
c                                      depth across the cross section
c  8 nboun             I  Number of boundary nodes.
c  9 ngrid             I  Number of grid points in network.
c 19 nlev              P  -
c 10 nnode             I  Number of nodes.
c 12 node(4,nnode)     I  Definition of nodes:
c                         (1,i) = Type of node i:
c                                 cintnd (1) : Internal node
c                                 chbou  (2) : H-boundary
c                                 cqbou  (3) : Q-boundary
c                                 cqhbou (4) : QH-boundary
c                                 chqbou (5) : HQ-boundary
c                         (2,i) = Gridpoint in case of boundary, else
c                                 undefined.
c                         (3,i) = Station number for boundary, undefined
c                                 for internal nodes:
c                                 HQ, H-boundary: station nr H-station.
c                                 QH, Q-boundary: station nr Q-station.
c                         (4,i) = Boundary number in case of boundary.
c 16 ntab(4,maxtab)    I  Table descriptor. Tables are numbererd from 1
c                         to maxtab. For a specific table number k and
c                         function Y = f (X) the following definitions
c                         exist:
c                         (1,k) = Length of table k.
c                         (2,k) = Start address X in table.
c                         (3,k) = Start address Y in table.
c                         (4,k) = Access method and period control: xy
c                                 x = ctbnpf (0) : No period defined
c                                 x = ctbpfu (1) : Period defined
c                                 y = ctbico (0) : Continue interpltn
c                                 y = ctbidi (1) : Discrete interpltn
c 15 ntabm             I  Maximum size of table (Used for dimensioning
c                         table).
c 24 sedtr(ngrid,*)    I  Sediment transport results for each gridpoint.
c               1|2       (At first transports per unit width, finaly
c                         total transports)
c                         - Normal branches:
c                         (i,1) = Sediment transport in gridpoint i of
c                                 main section.
c                         - Sedredge branches:
c                         (i,1) = Sediment transport in gridpoint i,
c                                 left channel.
c                         (i,2) = Sediment transport in gridpoint i,
c                                 right channel.
c                         (i,3) = Sediment exchange in gridpoint i.
c 17 table             P  -
c  5 time              P  -
c 22 wft(ngrid,maxlev) I  (i,j) = flow width at h = hlev(i,j) for grid
c                                 point i.
c                         - For a circle cross section:
c                         (i,1) = Radius of the circle.
c                         - For a sedredge cross section:
c                         (i,1) = Width of main section (i.e. left chan-
c                                 nel).
c                         (i,2) = Width of sub section 1 (i.e. right
c                                 channel).
c 23 ws                P  -
c-----------------------------------------------------------------------
c Subprogram calls:
c NAME    DESCRIPTION
c inttab  INTerpolate in TABle
c moeqdz  MORPHology EQually over Transport Width Delta Z calculation
c moeqtw  MORPHology erosion/sedimentation spread EQually over Transport Width
c momlev  MORPHology Morphodynamic LEVel
c moprdz  MORPHology PRoportional Delta Z calculation
c moprld  MORPHology erosion/sedimentation PRoportional to Local Depth
c moseci  MOrphology SECtion I
c=======================================================================
c
c
c***********************************************************************
c CVS log information:
c
c $Id$
c
c History:
c $Log: moadcs.pf,v $
c Revision 1.7  1999/03/15  15:52:44  kuipe_j
c tabs removed
c
c Revision 1.6  1998/06/11  11:47:05  kuipe_j
c Estuary special integrated
c
c Revision 1.5  1997/06/17  11:18:27  kuipe_j
c Remove undefined vars
c
c Revision 1.4  1996/09/03  14:48:46  kuipe_j
c frequency time hist,Improved sed distribution at nodes
c
c Revision 1.3  1995/05/30  09:55:45  hoeks_a
c Minor changes
c
c Revision 1.2  1995/05/30  07:04:31  hoeks_a
c file changed from dos to ux
c
c Revision 1.1  1995/04/13  07:09:01  hoeks_a
c Initial check-in
c
c Revision 1.4  1994/12/02  13:27:18  kuipe_j
c At a boundary always 'equal' distribution.
c
c Revision 1.3  1994/11/28  08:52:18  kuipe_j
c Time and timestep in double precision.
c
c Revision 1.2  1993/11/26  15:32:24  kuipe_j
c Update after finishing Sobeksel.
c
c Revision 1.1.1.1  1993/07/21  14:44:05  kuipe_j
c Initial version
c
c
c***********************************************************************
c
c     Parameters
c
      integer   igp, ibr ,
     +          ibrtyp,
     +          isec,
     +          maxlev,
     +          maxtab,
     +          moropt,
     +          ntabm,
     +          nboun,
     +          ngrid,
     +          nnode
c Aanpassing Kees Sloff 8-8-1997: branch dimensies 
      integer   branch (4,ibr),
     +          mbdpar (5,nboun),
     +          nlev   (ngrid),
     +          ntab   (4,maxtab),
     +          node   (4,nnode),
     +          flwdir (ngrid)

      real      ws(ngrid)

      real      table  (ntabm),
     +          wft    (ngrid,maxlev)

      double precision    time, hlev   (ngrid,maxlev),
     +                    deltaa
c
c     Local variables
c
      integer   i1,  i2,  iboun,  itab,  k,   n1,  n2
c
      real      z, wsact
c
      double precision hws, deltaz 
c
      logical   entbou
c
c     Include sobek constants
c
      include '..\include\sobcon.i'
c
c     Read parameters from branch
c
      n1 = branch(1,ibr)
      n2 = branch(2,ibr)
      i1 = branch(3,ibr)
      i2 = branch(4,ibr)
c
c     Check if grid point is entering boundary:
c
c     This is the case if sediment is moving in at a node and the
c     boundary condition is given as a delta Z.
c
      if     (igp .eq. i1) then
c
c        Begin of branch (S > 0: entering)
c
         if (flwdir(igp) .gt. 0 ) then
c
c           Fetch boundary number if not internal node
c
            if (node (1,n1) .ne. cintnd) then
               iboun = node (4,n1)
c
c              Read morphodynamic boundary condition table number itab
c
               if (mbdpar(1,iboun) .eq. cmbzft) then
c
c                 z = f(t)
c
                  entbou = .true.
c
c                 Boundary table number (index 4 or 5)
c
                  itab = mbdpar(isec+3,iboun)
               else
                  entbou = .false.
               endif
            else
               entbou = .false.
            endif
         else
            entbou = .false.
         endif

      elseif (igp .eq. i2) then
c
c        End of branch (S < 0: entering boundary)
c
         if (flwdir(igp) .lt. 0 ) then
c
c           Fetch boundary number if not internal node
c
            if (node(1,n2) .ne. cintnd) then
               iboun = node (4,n2)
c
c              Read morphodynamic boundary condition table number itab
c
               if (mbdpar(1,iboun) .eq. cmbzft) then
c
c                 z = f(t)
c
                  entbou = .true.
c
c                 Boundary table number (index 4 or 5)
c
                  itab = mbdpar(isec+3,iboun)
               else
                  entbou = .false.
               endif
            else
               entbou = .false.
            endif
         else
            entbou = .false.
         endif

      else
c
c        Not a boundary gridpoint
c
         entbou = .false.
      endif


      if (ibrtyp .eq. ccrsed) then
c
c        Sedredge branch is processed differently
c
         if (.not. entbou) then
c
c           Calculate delta z for current section (left or right)
c
            deltaz = deltaa / wft(igp,isec)
c
c           Update bottom of section
c
            hlev(igp,isec) = hlev(igp,isec) - deltaz

         else
c
c           Determine Z from boundary condition z = f(t)
c
            call inttab ( int (ntab(1,itab)),
     +                    int (ntab(4,itab)),
     +                    table (int (ntab(2,itab))),
     +                    table (int (ntab(3,itab))),
     +                    time,
     +                    z
     +                  )
c
c           Store new bottom level
c
            hlev(igp,isec) = z

         endif
      else
c
c        Code for tabulated cross sections
c

c
c        Determine maximum level k in cross sectional table
c
         call momlev ( igp,   ngrid,  maxlev,
     +                 nlev,  wft,    ws,
     +                 k,     wsact
     +               )
c
c        Check for morphology option (equally or proportional)
c
         if (moropt .eq. ceqows) then
c
c           Erosion/sedimentation spread equally over transport width
c
            if (.not. entbou) then
c
c              Determine delta Z by calculation
c
               call moeqdz ( igp,    k,    ngrid,   maxlev,
     +                       nlev,   wft,  deltaa,  deltaz
     +                     )

            else
c
c              Determine delta Z from boundary condition z = f(t)
c
               call inttab ( int (ntab(1,itab)),
     +                       int (ntab(4,itab)),
     +                       table (int (ntab(2,itab))),
     +                       table (int (ntab(3,itab))),
     +                       time,
     +                       z
     +                     )
               deltaz = hlev(igp,1) - z
            endif
c
c           Adapt cross sectional table
c
            call moeqtw ( igp,   k,      ngrid,    maxlev,
     +                    hlev,  deltaz
     +                  )

         elseif (moropt .eq. cprodp) then
c
c           Erosion/sedimentation proportional to local depth
c
c           Calculate level of sediment transporting width
c
            call moseci ( ngrid  ,igp    ,maxlev ,nlev   ,
     +                    wft    ,hlev   ,wsact  ,hws    )
c
            if (.not. entbou) then
c
c              Determine delta Z by calculation
c
               call moprdz ( igp,     hws,    k,     ngrid,
     +                       maxlev,  nlev,   hlev,  wft,
     +                       deltaa,  deltaz
     +                     )
        
c
c              Adapt cross sectional table
c
               call moprld ( igp,     hws,    k,      ngrid,
     +                       maxlev,  hlev,   deltaz
     +                     )
            else
c
c              Determine delta Z from boundary condition z = f(t)
c
               call inttab ( int (ntab(1,itab)),
     +                       int (ntab(4,itab)),
     +                       table (int (ntab(2,itab))),
     +                       table (int (ntab(3,itab))),
     +                       time,
     +                       z
     +                     )
               deltaz = hlev(igp,1) - z
c
c              Adapt cross sectional table
c
               call moeqtw ( igp,   k,      ngrid,    maxlev,
     +                       hlev,  deltaz
     +                     )
            endif
         endif
      endif
      return
      end

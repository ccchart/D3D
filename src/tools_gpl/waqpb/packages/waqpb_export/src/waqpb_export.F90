!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2023.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  
!  

!
!     Program to compose a PROCES.ASC file from tables
!
!     This program consists of the following parts:
!     - Reading of tables holding PROCLIB data structure
!     - Loop over the processes:
!       - Empty PDF structure
!       - Construct and write PDF structure
!       - Dump structure to PROCES.ASC

!     Include data structures for tables and PDF-file
      include 'data_ff.inc'
      include 'pdf_ff.inc'
      integer      jndex , iproc , iinpu , iitem , ioutp , idisp , &
                   ioutf , isubs , naanta, ioffse, ioffs2, ivelo , &
                   istoc , iconf , naant2, serial, &
                   ierror, icnsb , imodv , i
      logical      itmswi(nitemm)
      logical      newfrm
      character*10 c10   
      character*20 c20
      character*50 adduni
      character*255 ArgumentString
      real         actdef, versio
      integer      lu_inp, lu_mes, status, lunfil
      
!     Defaults for command line arguments

      versio = 5.00
      serial = 20130101
      newfrm = .true.

      do i=1,9999
            call getarg (i,ArgumentString)
            if (ArgumentString.eq.'') exit
            if (index(ArgumentString,'-version').gt.0) then
                c20 = ArgumentString(9:28)
                read (c20,'(f20.0)',iostat=status) versio
            endif
            if (index(ArgumentString,'-serial').gt.0) then
                c20 = ArgumentString(8:27)
                read (c20,'(i20)',iostat=status) serial
            endif
            if (index(ArgumentString,'-newfrm').gt.0) newfrm = .true.
            if (index(ArgumentString,'-oldfrm').gt.0) newfrm = .false.
      enddo


      itmswi = .false.
      open ( newunit=lu_mes , file = 'waqpb_export.log' )
      if (newfrm) then
        write (lu_mes,'(''Using NEW format'')')
	else
        write (lu_mes,'(''Using OLD format'')')
	endif

      write (*,'('' Reading data......'')')

!----------------------------------------------------------------------c
!     READ DATABASE
!----------------------------------------------------------------------c

      call readdb ( lu_inp, lu_mes )

!     Check validity of table R9

      do imodv = 1,nmodv
          call zoek (modvci(imodv),nconf,confid,10,iconf)
          if ( iconf .le. 0 ) then
              write ( lu_mes , '(''Unknown config in TABLE5: '',a10,1x, a10)') modvci(imodv),modvit(imodv)
          end if
          call zoek (modvit(imodv),nitem,itemid,10,iitem)
          if ( iitem .le. 0 ) then
              write ( lu_mes , '(''Unknown item in TABLE5: '',a10,1x, a10)') modvci(imodv),modvit(imodv)
          end if
      end do

!     Create auxiliary table of substances

      nsubs = 0
      do icnsb = 1,ncnsb
          c10 = r2_sid(icnsb)

!         Lookup substance in item array
          call zoek (c10,nitem,itemid,10,iitem)
          if ( iitem .le. 0 ) then
              write (*,*) ' ITEM: ',c10
              STOP 'Unknown substance in R2 table'
          endif

!         Add to substances array
          call zoek (c10,nsubs,subsid,10,isubs)
          if ( isubs .le. 0 ) then
              if ( nsubs+1 .gt. nsubsm ) STOP 'Dimension NSUBSM'
              nsubs = nsubs+1
              subsid(nsubs) = c10
          endif
      enddo

!     Dump TRM tables

!      write (*,'('' Writing TRM tables......'')')
!      call writrm
      write (*,'('' Writing TRM tables for LaTeX......'')')
      call writex
      
!----------------------------------------------------------------------c
!     SET VERSION, SERIAL AND WRITE NEFIS FILE
!----------------------------------------------------------------------c

      write (lu_mes,'(''Writing NEFIS process definition file'')')
      call makind()
      call pdfnef(lu_mes    , serial, versio, ierror)
      if ( ierror .ne. 0 ) then
         write (lu_mes,'(''ERROR writing NEFIS file'')')
         write (*,'(''ERROR writing NEFIS file, see report file'')')
      endif

!----------------------------------------------------------------------c
!     LOOP OVER PROCESSES
!----------------------------------------------------------------------c

      write (*,'('' Making PROCES.ASC......'')')
      write (*,*)
      open ( newunit=lunfil , file = 'procesm.asc' )
      write ( lunfil , '(i10,50x,f8.2,2x,i10)' ) nproc,versio,serial

      do 800 iproc=1,nproc

          write (*,'(''+Process: '',a10)') procid(iproc)

!----------------------------------------------------------------------c
!         CONSTRUCT PROCESS
!----------------------------------------------------------------------c

!         Clear PDF structure

          ins = 0
          ine = 0
          ous = 0
          oue = 0
          flu = 0
          sto = 0
          dis = 0
          vel = 0

!         Fill PDF structure

!         INPUT ITEMS ON SEGMENT LEVEL/EXCHANGE LEVEL

!         scan input items table for FIRST occurence of proces
          call zoek ( procid(iproc), ninpu, inpupr, 10, ioffse )
          naanta = 0
          if ( ioffse .gt. 0 ) then

!             loop over all INPU rows related to this process

  410         continue
              naanta = naanta + 1

!             Process current row

!             Lookup item in items table
              iinpu = ioffse + naanta-1
              call zoek ( inpuit(iinpu), nitem, itemid, 10, iitem)
              if ( iitem .le. 0 ) stop 'unknown ITEM'

!             Documented items are marked for COEFEDIT.DAT
              if ( inpudo(iinpu) .eq. 'x' ) itmswi(iitem) = .true.

!             Find item properties and store in PDF structure
              if ( inpude(iinpu) .eq. 'Y' ) then
                  actdef = itemde(iitem)
              elseif ( inpude(iinpu) .eq. 'G' ) then
                  actdef = -888.
              elseif ( inpude(iinpu) .eq. 'B' ) then
                  actdef = -101.
              elseif ( inpude(iinpu) .eq. 'M' ) then
                  actdef = -11.
              elseif ( inpude(iinpu) .eq. 'O' ) then
                  actdef = -1.
              else
                  actdef = -999.
              endif
              if ( inpusx(iinpu) .eq. 1 ) then
                  ins = ins + 1
                  if ( ins .gt. insmax ) stop 'DIMENSION insmax'
                  ins_id(ins) = itemid(iitem)
                  if (newfrm) then
                  ins_nm(ins) = itemnm(iitem)
                  ins_un(ins) = itemun(iitem)
                  else
                  ins_nm(ins) = adduni(itemnm(iitem),itemun(iitem))
                  endif
                  ins_va(ins) = actdef
                  ins_do(ins) = inpudo(iinpu)
              else
                  ine = ine + 1
                  if ( ine .gt. inemax ) stop 'DIMENSION inemax'
                  ine_id(ine) = itemid(iitem)
                  if (newfrm) then
                  ine_nm(ine) = itemnm(iitem)
                  ine_un(ine) = itemun(iitem)
                  else
                  ine_nm(ine) = adduni(itemnm(iitem),itemun(iitem))
                  endif
                  ine_va(ine) = actdef
                  ine_do(ine) = inpudo(iinpu)
              endif

!             Back for next row in table INPU,
!             if it still matches current proces

              if ( (iinpu+1) .le. ninpu ) then
                  call zoek ( procid(iproc), 1, inpupr(iinpu+1), 10, jndex )
                  if ( jndex .gt. 0 ) goto 410
              endif
          endif

!         OUTPUT ITEMS ON SEGMENT LEVEL/EXCHANGE LEVEL

!         scan output items table for FIRST occurence of proces
          call zoek ( procid(iproc), noutp, outppr, 10, ioffse )
          naanta = 0
          if ( ioffse .gt. 0 ) then

!             loop over all OUTP rows related to this process

  440         continue
              naanta = naanta + 1

!             Process current row

!             Lookup item in items table
              ioutp = ioffse + naanta-1
              call zoek ( outpit(ioutp), nitem, itemid, 10, iitem)
              if ( iitem .le. 0 ) stop 'unknown ITEM'

!             Find item properties and store in PDF structure
              if ( outpsx(ioutp) .eq. 1 ) then
                  ous = ous + 1
                  if ( ous .gt. ousmax ) stop 'DIMENSION ousmax'
                  ous_id(ous) = itemid(iitem)
                  if (newfrm) then
                  ous_nm(ous) = itemnm(iitem)
                  ous_un(ous) = itemun(iitem)
                  else
                  ous_nm(ous) = adduni(itemnm(iitem),itemun(iitem))
                  endif
                  ous_do(ous) = outpdo(ioutp)
              else
                  oue = oue + 1
                  if ( oue .gt. ouemax ) stop 'DIMENSION ouemax'
                  oue_id(oue) = itemid(iitem)
                  if (newfrm) then
                  oue_nm(oue) = itemnm(iitem)
                  oue_un(oue) = itemun(iitem)
                  else
                  oue_nm(oue) = adduni(itemnm(iitem),itemun(iitem))
                  endif
                  oue_do(oue) = outpdo(ioutp)

!                 SCAN VELO and DISP TABLES FOR LINES ASSOCIATED WITH
!                 CURRENT OUTPUT ITEM ON EXCHANGE LEVEL

!                 scan dispersion lines table for FIRST occurence of item
                  call zoek ( itemid(iitem), ndisp, dispit, 10, ioffs2)
                  naant2 = 0
                  if ( ioffs2 .gt. 0 ) then

!                     loop over all DISP rows related to this item

  450                 continue
                      naant2 = naant2+1
                      dis = dis + 1
                      if ( dis .gt. dismax ) stop 'dimension DISMAX'

!                     Process current row

                      idisp = ioffs2 + naant2-1
                      dis_su(dis) = dispsu(idisp)
                      dis_it(dis) = dispit(idisp)
                      dis_sc(dis) = dispsc(idisp)

!                     Back for next row in table DISP,
!                     if it still matches current item

                      if ( (idisp+1) .le. ndisp ) then
                          call zoek ( itemid(iitem), 1, dispit(idisp+1), 10, jndex )
                          if ( jndex .gt. 0 ) goto 450
                      endif
                  endif

!                 scan velocity lines table for FIRST occurence of item
                  call zoek ( itemid(iitem), nvelo, veloit, 10, ioffs2)
                  naant2 = 0
                  if ( ioffs2 .gt. 0 ) then

!                     loop over all VELO rows related to this item

  460                 continue
                      naant2 = naant2+1
                      vel = vel + 1
                      if ( vel .gt. velmax ) stop 'dimension VELMAX'

!                     Process current row

                      ivelo = ioffs2 + naant2-1
                      vel_su(vel) = velosu(ivelo)
                      vel_it(vel) = veloit(ivelo)
                      vel_sc(vel) = velosc(ivelo)

!                     Back for next row in table VELO,
!                     if it still matches current item

                      if ( (ivelo+1) .le. nvelo ) then
                          call zoek ( itemid(iitem), 1, veloit(ivelo+1), 10, jndex )
                          if ( jndex .gt. 0 ) goto 460
                      endif
                  endif

!                 END of processing output item on exchange level!

              endif

!             Back for next row in table OUTP,
!             if it still matches current proces

              if ( (ioutp+1) .le. noutp ) then
                  call zoek ( procid(iproc), 1, outppr(ioutp+1), 10, jndex )
                  if ( jndex .gt. 0 ) goto 440
              endif
          endif

!         FLUXES

!         scan output fluxes table for FIRST occurence of proces
          call zoek ( procid(iproc), noutf, outfpr, 10, ioffse )
          if ( ioffse .gt. 0 ) then

!             loop over all FLUX rows related to this process

  470         continue
              flu = flu + 1
              if ( flu .gt. flumax ) stop 'dimension FLUMAX'

!             Process current row

!             Lookup flux in items table
              ioutf = ioffse + flu-1
!             write (lu_mes,*) ' flu ',flu,' ioutf ', ioutf
              call zoek ( outffl(ioutf), nitem, itemid, 10, iitem)
              if ( iitem .le. 0 ) stop 'unknown FLUX'

!             Find and store flux properties
              flu_id(flu) = itemid(iitem)
              if (newfrm) then
              flu_nm(flu) = itemnm(iitem)
              flu_un(flu) = itemun(iitem)
              else
              flu_nm(flu) = adduni(itemnm(iitem),itemun(iitem))
              endif
              flu_do(flu) = outfdo(ioutf)

!             SCAN STOCHI TABLE FOR LINES ASSOCIATED WITH PRESENT FLUX

!             scan stochi lines table for FIRST occurence of flux
              call zoek ( itemid(iitem), nstoc, stocfl, 10, ioffs2)
              naant2 = 0
              if ( ioffs2 .gt. 0 ) then

!                 loop over all STOC rows related to this flux

  480             continue
                  naant2 = naant2+1
                  sto = sto + 1
                  if ( sto .gt. stomax ) stop 'dimension STOMAX'

!                 Process current row

                  istoc = ioffs2 + naant2-1
                  sto_su(sto) = stocsu(istoc)
                  sto_fl(sto) = stocfl(istoc)
                  sto_sc(sto) = stocsc(istoc)

!                 Back for next row in table STOC,
!                 if it still matches current flux

                  if ( (istoc+1) .le. nstoc ) then
                      call zoek ( itemid(iitem), 1, stocfl(istoc+1), &
                            10, jndex )
                      if ( jndex .gt. 0 ) goto 480
                  endif
              endif

!             Back for next row in table OUTF,
!             if it still matches current proces

              if ( (ioutf+1) .le. noutf ) then
                  call zoek ( procid(iproc), 1, outfpr(ioutf+1), 10, jndex )
                  if ( jndex .gt. 0 ) goto 470
              endif
          endif

!----------------------------------------------------------------------c
!         WRITE PROCESS
!----------------------------------------------------------------------c

!         Write PDF file (formats as in HARMONIZE to allow comparison)

          if (newfrm) then
          call wripdn ( procid(iproc), procnm(iproc), procco(iproc), procfo(iproc), lunfil )
          else
          call wripdf ( procid(iproc), procnm(iproc), procco(iproc), procfo(iproc), lunfil )
          endif
  800 continue
      close (lunfil)

!     Write all active coefficients to COEFEDIT.DAT in the Sobek-format
      call coefed(serial,itmswi)

  900 continue
      close (lu_mes)

      stop 'Normal end'
      end
      function adduni(name,unit)
      character*50 adduni, name
      character*20 unit

      integer      lennam, lenuni, i

!     find length of name and unit

      lennam = -1
      do 10 i = 50,1,-1
          if ( name(i:i) .ne. ' ' ) then
              lennam = i
              goto 11
          endif
   10 continue
   11 continue
      if ( lennam .le. 1 ) then
          write (*,*) ' ',name
          stop 'ADDUNI Fatal Error 1'
      endif

      lenuni = 0
      if ( unit(2:3) .eq. 'no' .and. unit(5:8) .eq. 'unit' ) then
          lenuni = 0
      else
          do 20 i = 20,1,-1
              if ( unit(i:i) .ne. ' ' ) then
                  lenuni = i
                  goto 21
              endif
   20     continue
   21     continue
      endif
      if ( lenuni .lt. 0 ) then
          write (*,*) ' ',unit
          stop 'ADDUNI Fatal Error 1'
      endif

      if ( lennam + lenuni .gt. 50 ) then
          lennam = 50-lenuni
      endif

      write ( adduni(1          :lennam) , '(a)' ) name(1:lennam)
      do 30 i = 1,50-lennam-lenuni
   30 adduni(lennam+i:lennam+i) = ' '
      if (lenuni.gt.0) then
         write ( adduni(50-lenuni+1:50    ) , '(a)' ) unit(1:lenuni)
      end if

      return
      end

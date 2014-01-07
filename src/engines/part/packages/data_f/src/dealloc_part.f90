!!  Copyright(C) Stichting Deltares, 2012-2014.
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

!!  Note: The "part" engine is not yet Open Source, but still under
!!  development. This package serves as a temporary dummy interface for
!!  the references in the "waq" engine to the "part" engine.

module dealloc_part_mod

use global_pointers

use precision

use dealloc_mod

implicit none

contains
      subroutine dealloc_part()

      call dealloc(abuoy      )
      call dealloc(acf        )
      call dealloc(aconc      )
      call dealloc(aconud     )
      call dealloc(adepth     )
      call dealloc(amap       )
      call dealloc(amapsett   )
      call dealloc(amassc     )
      call dealloc(amassd     )
      call dealloc(angle      )
      call dealloc(apeak      )
      call dealloc(area       )
      call dealloc(atotal     )
      call dealloc(atrack     )
      call dealloc(cbuff      )
      call dealloc(cdelv      )
      call dealloc(chismp     )
      call dealloc(chispl     )
      call dealloc(conc       )
      call dealloc(const      )
      call dealloc(constev    )
      call dealloc(c2         )
      call dealloc(decay      )
      call dealloc(decays     )
      call dealloc(depth      )
      call dealloc(dfact      )
      call dealloc(dps        )
      call dealloc(drand      )
      call dealloc(dx         )
      call dealloc(dy         )
      call dealloc(elt_names  )
      call dealloc(elt_types  )
      call dealloc(elt_dims   )
      call dealloc(elt_bytes  )
      call dealloc(finud      )
      call dealloc(floil      )
      call dealloc(flow       )
      call dealloc(flres      )
      call dealloc(fractd     )
      call dealloc(fracte     )
      call dealloc(fstick     )
      call dealloc(ftime      )
      call dealloc(fwatoil    )
      call dealloc(ibuff      )
      call dealloc(ictime     )
      call dealloc(ictmax     )
      call dealloc(idtime     )
      call dealloc(ifopt      )
      call dealloc(iftime     )
      call dealloc(ihplot     )
      call dealloc(imap       )
      call dealloc(imask      )
      call dealloc(ioptrad    )
      call dealloc(ipnt       )
      call dealloc(ipset      )
      call dealloc(iptime     )
      call dealloc(isfile     )
      call dealloc(isfud      )
      call dealloc(isub       )
      call dealloc(iutime     )
      call dealloc(ivtime     )
      call dealloc(iwndtm     )
      call dealloc(iwtime     )
      call dealloc(kpart      )
      call dealloc(kwaste     )
      call dealloc(lgrid      )
      call dealloc(lgrid2     )
      call dealloc(linear     )
      call dealloc(locdep     )
      call dealloc(mapsub     )
      call dealloc(mcell      )
      call dealloc(mpart0     )
      call dealloc(mpart      )
      call dealloc(mplsta     )
      call dealloc(mstat      )
      call dealloc(mstick     )
      call dealloc(mwaste     )
      call dealloc(nbin       )
      call dealloc(ncell      )
      call dealloc(ncheck     )
      call dealloc(ndprt      )
      call dealloc(nmconr     )
      call dealloc(nmdyer     )
      call dealloc(nmstat     )
      call dealloc(nosud      )
      call dealloc(nosys      )
      call dealloc(npart0     )
      call dealloc(npart      )
      call dealloc(nplay      )
      call dealloc(nplot      )
      call dealloc(nplsta     )
      call dealloc(nstat      )
      call dealloc(nwaste     )
      call dealloc(radius     )
      call dealloc(rbuff      )
      call dealloc(rbuffr     )
      call dealloc(recovr     )
      call dealloc(rem        )
      call dealloc(rhooilv    )
      call dealloc(stoch      )
      call dealloc(stoil      )
      call dealloc(subst      )
      call dealloc(substi     )
      call dealloc(subst2     )
      call dealloc(subsud     )
      call dealloc(t0buoy     )
      call dealloc(t0cf       )
      call dealloc(tcktot     )
      call dealloc(title      )
      call dealloc(tmass      )
      call dealloc(tmassc     )
      call dealloc(tmassu     )
      call dealloc(tmasud     )
      call dealloc(totfe      )
      call dealloc(track      )
      call dealloc(qentr      )
      call dealloc(uscal      )
      call dealloc(vdiff      )
      call dealloc(velo       )
      call dealloc(viso       )
      call dealloc(visowat    )
      call dealloc(volfracw   )
      call dealloc(vol1       )
      call dealloc(vol2       )
      call dealloc(volume     )
      call dealloc(vrtdsp     )
      call dealloc(vsfact     )
      call dealloc(vsfour     )
      call dealloc(wdira      )
      call dealloc(wevap      )
      call dealloc(window     )
      call dealloc(wparm      )
      call dealloc(wpart      )
      call dealloc(wsettl     )
      call dealloc(wveloa     )
      call dealloc(xa0        )
      call dealloc(xa         )
      call dealloc(xb         )
      call dealloc(xpart0     )
      call dealloc(xpart      )
      call dealloc(xpol       )
      call dealloc(xstat      )
      call dealloc(xwaste     )
      call dealloc(xyztrk     )
      call dealloc(ya0        )
      call dealloc(ya         )
      call dealloc(yb         )
      call dealloc(ypart0     )
      call dealloc(ypart      )
      call dealloc(ypol       )
      call dealloc(ystat      )
      call dealloc(ywaste     )
      call dealloc(za         )
      call dealloc(zlevel     )
      call dealloc(zpart      )
      call dealloc(zwaste     )
      return
      end subroutine dealloc_part
end module dealloc_part_mod

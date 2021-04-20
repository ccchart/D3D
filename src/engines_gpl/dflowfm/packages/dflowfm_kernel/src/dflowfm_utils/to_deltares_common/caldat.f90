 subroutine caldat(julian,mm,id,iyyy)
 implicit none
 integer :: julian,mm,id,iyyy
 integer :: igreg
 parameter (igreg=2299161)
 integer :: jalpha, ja, jb, jc, jd, je
 if(julian.ge.igreg)then
   jalpha=int(((julian-1867216)-0.25)/36524.25)
   ja=julian+1+jalpha-int(0.25*jalpha)
 else
   ja=julian
 endif
 jb=ja+1524
 jc=int(6680.+((jb-2439870)-122.1)/365.25)
 jd=365*jc+int(0.25*jc)
 je=int((jb-jd)/30.6001)
 id=jb-jd-int(30.6001*je)
 mm=je-1
 if(mm.gt.12)mm=mm-12
 iyyy=jc-4715
 if(mm.gt.2)iyyy=iyyy-1
 if(iyyy.le.0)iyyy=iyyy-1
 return
 end subroutine caldat

module larvm2_mod
!
!  data definition module(s)
!
use precision_part               ! single/double precision
use timers
!
!  module procedure(s)
!
!
implicit none

contains
      subroutine larvm2 ( lunrep   , itime    , nosegl   , nolay    , nosubs   ,    &
                          conc     )

      ! function            : calculate larvae per m2 over the water column, output to map file

      ! arguments :

      integer(ip), intent(in)    :: lunrep              ! report file
      integer(ip), intent(in)    :: itime               ! time in seconds
      integer(ip), intent(in)    :: nosegl              ! number segments per layer
      integer(ip), intent(in)    :: nolay               ! number of layers in calculation
      integer(ip), intent(in)    :: nosubs              ! number of substances per particle
      real   (sp), pointer       :: conc  ( : , : )     ! concentration array in transport grid

      ! function            : calculate larvae per m2 over the water column

      integer, parameter                :: lunmm2 = 845
      character(len=256)                :: filmm2
      real   , allocatable, save        :: c_m2(:,:)
      character(len=20), allocatable    :: syname(:)
      character(len=40)                 :: moname(4)
      integer, save                     :: notot
      integer                           :: iseg
      integer                           :: isegl
      integer                           :: ilay
      integer                           :: isys
      integer                           :: istage
      integer, save                     :: ifirst = 1

      if ( ifirst .eq. 1 ) then
         ifirst = 0
         filmm2    = 'larvae_m2.map'
         moname(1) = 'larvae per m2'
         moname(2) = ' '
         moname(3) = ' '
         moname(4) = ' '
         notot     = nosubs + 1
         allocate(syname(nosubs-3))
         allocate(c_m2(nosubs-3,nosegl))
         syname(1) = 'sole'
         do isys = 2, nosubs-3
            istage = isys - 1
            write(syname(isys),'(''sole'',i2.2,''stage'')') istage
         enddo
         open (lunmm2, file=filmm2, form='binary')
         write(lunmm2) moname
         write(lunmm2) nosubs-3,nosegl*nolay
         write(lunmm2) syname
      endif

      do isegl = 1, nosegl

         ! if seperate stages accumulate stages (thus skipping stage 0), also seperate per stage, notot = localdepth

         if ( nosubs .gt. 4 ) then

            c_m2(1,isegl) = 0.0
            do istage = 1, nosubs-4
               c_m2(1,isegl)        = c_m2(1,isegl) + conc(4+istage,isegl)*conc(notot,isegl)
               c_m2(1+istage,isegl) = conc(4+istage,isegl)*conc(notot,isegl)
               do ilay = 2, nolay
                  iseg = isegl + (ilay-1)*nosegl
                  c_m2(1,isegl)        = c_m2(1,isegl) + conc(4+istage,iseg)*(conc(notot,iseg)-conc(notot,iseg-nosegl))
                  c_m2(1+istage,isegl) = c_m2(1+istage,isegl) + conc(4+istage,iseg)*(conc(notot,iseg)-conc(notot,iseg-nosegl))
               enddo
               iseg = isegl + ilay*nosegl
               c_m2(1,isegl)        = c_m2(1,isegl) + conc(4+istage,iseg)
               c_m2(1+istage,isegl) = c_m2(1+istage,isegl) + conc(4+istage,iseg)
            enddo

         else

            ! just cummulate substance 1 per layer

            c_m2(1,isegl) = conc(1,isegl)*conc(notot,isegl)
            do ilay = 2, nolay
               iseg = isegl + (ilay-1)*nosegl
               c_m2(1,isegl) = c_m2(1,isegl) + conc(1,iseg)*(conc(notot,iseg)-conc(notot,iseg-nosegl))
            enddo
            iseg = isegl + ilay*nosegl
            c_m2(1,isegl) = c_m2(1,isegl) + conc(1,iseg)

         endif
      enddo
      write(lunmm2) itime, (c_m2, ilay = 1, nolay) ! for consistency with the used lga/cco

      return
      end subroutine
end module

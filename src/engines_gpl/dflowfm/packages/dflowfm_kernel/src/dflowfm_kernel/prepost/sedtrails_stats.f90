module m_sedtrails_stats
   use m_sedtrails_data
   integer                                  :: is_numndvals         !< Number of variables on flow nodes for which statistics are recorded.
   integer, parameter                       :: IDX_BL         = 1   !< Index for avg bottom level
   integer, parameter                       :: IDX_HS         = 2   !< Index for avg water depth
   integer, parameter                       :: IDX_UCX        = 3   !< Index for avg x velocity
   integer, parameter                       :: IDX_UCY        = 4   !< Index for avg y velocity
   integer, parameter                       :: IDX_TAUS       = 5   !< Index for avg x bed shear stress
   integer, parameter                       :: IDX_TAUSMAX    = 6   !< Index for avg y shear stress
   integer, parameter                       :: IDX_SBX        = 7    
   integer, parameter                       :: IDX_SBY        = 8   
   integer, parameter                       :: IDX_SSX        = 9   
   integer, parameter                       :: IDX_SSY        = 10  
   integer, parameter                       :: IDX_SSC        = 11  !< Index for avg sediment concentration  
                                            
   double precision, allocatable, target    :: is_sumvalsnd(:,:,:) 
   
   character(len=1024), allocatable, target :: is_valnamesnd(:) 
   double precision, target                 :: is_dtint             !< [s] total time interval since last statistics reset.
   
   contains
   
   !> Sets ALL (scalar) variables in this module to their default values.
   !! For a reinit prior to flow computation, only call reset_integralstats() instead.
   subroutine default_sedtrails_stats()
   
      is_numndvals = 0
      is_valnamesnd(:)  = ''
      is_valnamesnd(1)  = 'bl'
      is_valnamesnd(2)  = 'hs'
      is_valnamesnd(3)  = 'ucx'
      is_valnamesnd(4)  = 'ucy'
      is_valnamesnd(5)  = 'taus'      ! change to vector comps after sedmor merge
      is_valnamesnd(6)  = 'tausmax'
      is_valnamesnd(7)  = 'sbx'
      is_valnamesnd(8)  = 'sby'
      is_valnamesnd(9)  = 'ssx'
      is_valnamesnd(10) = 'ssy'
      is_valnamesnd(11) = 'ssc'
   
      ! Remaining of variables is handled in reset_integralstats()
      call reset_sedtrails_stats()
   end subroutine default_sedtrails_stats
   
   !> Resets only stats variables intended for a restart of flow simulation.
   !! Upon loading of new model/MDU, call default_sedtrails_stats() instead.
   subroutine reset_sedtrails_stats()
       is_sumvalsnd(1:is_numndvals,:,:) = 0d0   
       is_dtint = 0d0
   end subroutine reset_sedtrails_stats
   
   subroutine alloc_sedtrails_stats()
      use m_alloc
      use m_fm_erosed, only: lsedtot
      use m_flowgeom, only: ndx
      
      implicit none
      
      if (is_numndvals > 0) then
         call realloc(is_sumvalsnd, (/ is_numndvals, ndx, lsedtot /), keepExisting = .false., fill = 0d0)
         call realloc(is_valnamesnd, is_numndvals, keepExisting = .false., fill = '')
      end if
   
   end subroutine alloc_sedtrails_stats
   
   !> Update the (time-)integral statistics for all flow nodes, typically after each time step.
   subroutine update_sedtrails_stats()
      use m_flowtimes, only: dts
      use m_flow, only: hs, ucx, ucy, taus, tausmax
      use m_flowgeom, only: ndx, bl
      use m_fm_erosed
      use m_transport, only: constituents, ISED1
      use m_sediment, only: sedtot2sedsus, stm_included
   
      integer :: k
   
      if (is_numndvals <= 0) then
         return
      end if
      
      if (jawave<3) then      ! do not overwrite current+wave induced bed shear stresses from tauwave
         call gettaus(1, 1)   ! tausmax known from taubxu_nowave
      endif
   
      do k=1,ndx
         is_sumvalsnd(IDX_BL  , k, 1) = is_sumvalsnd(IDX_BL      ,k, 1) + dts * bl(k)
         is_sumvalsnd(IDX_HS  , k, 1) = is_sumvalsnd(IDX_HS      ,k, 1) + dts * hs(k)
         is_sumvalsnd(IDX_UCX , k, 1) = is_sumvalsnd(IDX_UCX     ,k, 1) + dts * ucx(k)     ! assumes depth-averaged value in base node index
         is_sumvalsnd(IDX_UCY , k, 1) = is_sumvalsnd(IDX_UCY     ,k, 1) + dts * ucy(k)
         is_sumvalsnd(IDX_TAUS, k, 1) = is_sumvalsnd(IDX_TAUS    ,k, 1) + dts * taus(k)
      enddo
      
      if (stm_included) then
         do k=1,ndx
             is_sumvalsnd(IDX_TAUSMAX, k, 1) = is_sumvalsnd(IDX_TAUSMAX ,k, 1) + dts * tausmax(k)  
         enddo
         
         do l=1,lsedtot 
            do k=1, ndx
               is_sumvalsnd(IDX_SBX , k, l) = is_sumvalsnd(IDX_SBX  ,k, l) + dts * (sbcx(k,l) + sbwx(k,l))
               is_sumvalsnd(IDX_SBY , k, l) = is_sumvalsnd(IDX_SBY  ,k, l) + dts * (sbcy(k,l) + sbwy(k,l))
               is_sumvalsnd(IDX_SSX , k, l) = is_sumvalsnd(IDX_SSX  ,k, l) + dts * (sscx(k,l) + sswx(k,l))
               is_sumvalsnd(IDX_SSY , k, l) = is_sumvalsnd(IDX_SSY  ,k, l) + dts * (sscy(k,l) + sswy(k,l))
            enddo
         end do 
         
         do l=1, lsedsus
            do k=1,ndx
               is_sumvalsnd(IDX_SSC , k, sedtot2sedsus(l)) = is_sumvalsnd(IDX_SSC  ,k, sedtot2sedsus(l)) + dts * constituents(ISED1+l-1,k)   ! this works just for 2D now
            enddo   
         enddo   
      endif
      
      is_dtint = is_dtint + dts         
                                           
   end subroutine update_sedtrails_stats
                                        
end module m_sedtrails_stats         
                                        
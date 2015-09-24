module m_hash_search
   
   use m_alloc
   use m_netw, only: lenc
   use messageHandling
   
   implicit none
   
   public hashfill
   public hashsearch

   
   type, public :: t_hashlist
      integer :: hashcon = 1009
      integer :: id_count = 0
      character(len=idLen), allocatable, dimension(:) :: id_list
      integer, allocatable, dimension(:) :: hashfirst
      integer, allocatable, dimension(:) :: hashnext
   end type
   
contains   
   
   integer function hashfun(string, hashcon)

      ! Hashing function
      ! Original by: Geert Prinsen
      ! Module description: Modified version of the hashing system used
      !                     in SOBEK_LITE/PLUVIUS. 

      character(*), intent(in)               :: string
      integer, intent(in)                    :: hashcon
      
      integer                                :: ires
      integer                                :: length
      integer                                :: i

      ires = 0
      length = len_trim(string)

      do i = 1, length
         ires = ires + iachar(string(i:i))
      enddo

      ires = mod (ires, hashcon)
      
      if (ires == 0) ires = hashcon
      
      hashfun = ires

   end function hashfun

   subroutine hashfill(hashlist)
 
      ! Module description: Fill hashing arrays
 
      ! Global Variables
      type(t_hashlist), pointer, intent(inout) :: hashlist

      ! Local Variables
      integer                                  :: icount
      integer                                  :: hashcode
      integer                                  :: inr
      integer                                  :: next
      integer                                  :: hashfun
      integer                                  :: ierr
      character(len=lenc)                      :: locid
      
      call realloc(hashlist%hashfirst, hashlist%hashcon - 1, lindex = 0, stat = ierr)
      call aerr('hashfirst(0:hashcon - 1)', ierr, hashlist%hashcon)
    
      call realloc(hashlist%hashnext, hashlist%id_count, stat = ierr)
      call aerr('hashnext(id_count)', ierr, hashlist%id_count)
      
      hashlist%hashfirst = 0
      hashlist%hashnext  = 0

      do icount = 1, hashlist%id_count
      
         locid = hashlist%id_list(icount)
         call str_upper(locid)
         
         hashcode = hashfun(locid)
  
         !      write(*,*) ' Hashfill ', id,' ', hashcode
  
         if (hashlist%hashfirst(hashcode) .eq. 0) then
            
            hashlist%hashfirst(hashcode) = icount
            hashlist%hashnext(icount)    = 0
           
         else
            
            inr  = hashlist%hashfirst(hashcode)
            next = hashlist%hashnext(inr)
           
            do while (next .ne. 0)
               inr  = next
               next = hashlist%hashnext (inr)
            enddo
           
            hashlist%hashnext(inr) = icount
            
         endif
 
      end do
      
   end subroutine hashfill
 
 !  inode = hashsearch(node_local_id, nodefirst, nodenext, node_id)
   
   integer function hashsearch(hashlist, id)
 
      ! Module description: Search in hashing arrays
 
      ! Global Variables
      character(len=*), intent(in)                    :: id
      type(t_hashlist), pointer, intent(inout)        :: hashlist
      
      ! Local Variables
      character(len=lenc)                             :: locid
      character(len=lenc)                             :: idtest
      integer                                         :: hashcode
      integer                                         :: inr
      integer                                         :: next
      integer                                         :: ifound
 
      ifound = -1
      locid  = id
      call str_upper(locid)
      
      hashcode = hashfun(locid, hashlist%hashcon)
 
      !   write(*,*) ' Hashsearch', id, hashcode
 
      if (hashlist%hashfirst(hashcode) > 0) then
 
        inr  = hashlist%hashfirst(hashcode)
        next = inr
        
        do while (next .ne. 0)
           
          idtest = hashlist%id_list(next)
          call upperc (idtest)
          
          if (locid .ne. idtest) then
            inr  = next
            next = hashlist%hashnext (inr)
          else
            ifound = next
            exit
          endif
          
        enddo
        
      else
         ! 'Hash search failed'
      endif

      hashsearch = ifound
 
   end function hashsearch

end module m_hash_search

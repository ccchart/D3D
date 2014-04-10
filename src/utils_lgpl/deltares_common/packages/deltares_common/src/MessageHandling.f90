!----- LGPL --------------------------------------------------------------------
!
!  Copyright (C)  Stichting Deltares, 2011-2014.
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation version 2.1.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, see <http://www.gnu.org/licenses/>.
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
!  $Id$
!  $HeadURL$

!> Specifies the interface for MessageHandling's callback functionality.
! (A bit awkward, but including it in MessageHandling's module header
! did not make it visible in subprograms when using ifort 11.0.)
module MHCallBack
   abstract interface
      subroutine mh_callbackiface(level, message)
         integer, intent(in) :: level !< The severity level
         character(len=*), intent(in) :: message !< log message
      end subroutine mh_callbackiface
   end interface

   abstract interface
      subroutine c_callbackiface(level, msg)
        use iso_c_binding
        use iso_c_utils
        integer(c_int), intent(in) :: level !< severity
        character(c_char), intent(in) :: msg(MAXSTRINGLEN) !< c message null terminated
      end subroutine c_callbackiface
   end interface
   
   abstract interface
      subroutine progress_iface(msg, progress)
        character(len=*), intent(in) :: msg !< c message null terminated
        double precision,  intent(in) :: progress
      end subroutine progress_iface
   end interface


   abstract interface
      subroutine progress_c_iface(msg, progress)
        use iso_c_binding
        use iso_c_utils
        character(c_char), intent(in) :: msg(MAXSTRINGLEN) !< c message null terminated
        real(c_double), intent(in) :: progress !< progress in fraction
      end subroutine progress_c_iface
   end interface


end module MHCallBack

!> Diagnostics output module.
!! Prints and/or logs messages from an application.
!! Three variants:
!! -# write2screen: writes directly to stdout
!! -# useLog: log all messages in a message buffer (a queue), can be read out by any other application.
!! -# lunMessages: writes directly to an already opened specified file pointer.
!!
!! See MessageHandling::SetMessageHandling for more details.
!!
!! Messages have a severity level: LEVEL_(DEBUG|INFO|WARN|ERROR|FATAL).
module MessageHandling
   use MHCallBack
   use iso_c_utils
   implicit none

   procedure(progress_iface), pointer :: progress_callback => null()
   procedure(progress_c_iface), pointer :: progress_c_callback => null()

   integer, parameter, public    :: BUFLEN = 1024
   !> The message buffer allows you to write any number of variables in any
   !! order to a character string. Call msg_flush or err_flush to output
   !! the actual message or error.
   character(len=BUFLEN), public :: msg_mh
   character(len=MAXSTRINGLEN), public :: errmsg

   public SetMessage
   public GetMessageCount
   public SetMessageHandling
   public mess
   public err
   public GetMessage_MH
   public getOldestMessage
   public resetMessageCount_MH
   public getMaxErrorLevel
   public resetMaxerrorLevel
   public set_mh_c_callback
   public set_mh_callback
   public aerr
   public msg_flush
   public dbg_flush
   public warn_flush
   public err_flush
   public set_progress_callback
   public progress
   
   integer,parameter, public     :: LEVEL_DEBUG = 1
   integer,parameter, public     :: LEVEL_INFO  = 2
   integer,parameter, public     :: LEVEL_WARN  = 3
   integer,parameter, public     :: LEVEL_ERROR = 4
   integer,parameter, public     :: LEVEL_FATAL = 5
   integer,parameter, public     :: LEVEL_NONE  = 6
   integer,parameter, public     :: Charln = 256
   integer,parameter, public     :: Idlen = 40
   integer,parameter, public     :: max_level = 5
   character(len=12), dimension(max_level), private    :: level_prefix = (/'** DEBUG  : ',  &
                                                                           '** INFO   : ',  &
                                                                           '** WARNING: ',  &
                                                                           '** ERROR  : ',  &
                                                                           '** FATAL  : '/)
   character(len=12), public              :: space12 = ' '
   integer, dimension(max_level), public  :: mess_level_count


   interface mess
   module procedure message1string
   module procedure message2string
   module procedure message3string
   module procedure message4string
   module procedure message1char1real
   module procedure message1char2real
   module procedure message2char1real
   module procedure message2char2real
   module procedure message1char1int
   module procedure message1char2int
   module procedure message1char3int
   module procedure message1char1double
   module procedure message2int1char
   module procedure message1char1int1double
   module procedure message1double1int1char
   end interface

   interface err
   module procedure error1char
   module procedure error2char
   module procedure error3char
   module procedure error4char
   module procedure error1char1real
   module procedure error1char2real
   module procedure error2char1real
   module procedure error2char2real
   module procedure error1char1double
   module procedure error1char1int
   module procedure error1char2int
   module procedure error1char1int1double
   end interface

private

   integer               , parameter,              private :: maxMessages = 3000
   integer               ,                         private :: messagecount = 0 !< Number of messages currently in message buffer (queue).
   character(len=charln) , dimension(maxMessages), private :: Messages
   integer               , dimension(maxMessages), private :: Levels
   integer               ,                         private :: ibuffertail  = 0 !< Index of newest message in message buffer.

   integer,                                    private :: maxErrorLevel = 0
   integer,                                    public  :: thresholdLvl = 0

   integer, save                  :: lunMess          = 0
   logical, save                  :: writeMessage2Screen = .false.
   logical, save                  :: useLogging = .true.
   logical, save                  :: alreadyInCallback=.false.                   !< flag for preventing recursive calls to callback subroutine
   !> Callback routine invoked upon any mess/err (i.e. SetMessage)
   procedure(mh_callbackiface), pointer :: mh_callback => null()
   procedure(c_callbackiface), pointer :: f_callback => null()

contains

!> Sets up the output of messages. All three formats are optional
!! and can be used in any combination.
subroutine SetMessageHandling(write2screen, useLog, lunMessages, callback, thresholdLevel, reset_counters)
   logical, optional, intent(in)       :: write2screen !< Print messages to stdout.
   logical, optional, intent(in)       :: useLog       !< Store messages in buffer.
   integer, optional, intent(in)       :: lunMessages  !< File pointer whereto messages can be written.
   integer, optional, intent(in)       :: thresholdLevel  !< Messages with level lower than the thresholdlevel
                                                          !< will be discarded.
   logical, optional, intent(in)       :: reset_counters  !< If present and True then reset message counters.
                                                          !< SetMessageHandling is called more than once.

   procedure(mh_callbackiface), optional :: callback

   if (present(write2screen) ) writeMessage2Screen = write2screen
   if (present(lunMessages) )  lunMess             = lunMessages
   if (present(useLog) )       useLogging          = useLog
   if (present(callback) ) then
      call set_mh_callback(callback)
   endif
   if (present(thresholdLevel) )  thresholdLvl     = thresholdLevel

   if (present(reset_counters)) then
     if (reset_counters) then
        mess_level_count = 0
        messagecount = 0
        ibuffertail  = 0
     end if
   endif

   alreadyInCallback = .false.

end subroutine SetMessageHandling

subroutine set_mh_callback(callback)
  procedure(mh_callbackiface) :: callback
  mh_callback => callback
end subroutine set_mh_callback


subroutine set_mh_c_callback(c_callback) bind(C, name="set_mh_c_callback")
  !DEC$ ATTRIBUTES DLLEXPORT::set_mh_c_callback

  use iso_c_binding
  implicit none
  type(c_funptr) :: c_callback

  ! Set a callback that will be cauled with new messages

  call c_f_procpointer(c_callback, f_callback)
end subroutine set_mh_c_callback


!> The main message routine. Puts the message string to all output
!! channels previously set up by SetMessageHandling
recursive subroutine SetMessage(level, string)
  use iso_c_binding
  use iso_c_utils


  integer, intent(in)           :: level  !< One of: LEVEL_(DEBUG|INFO|WARN|ERROR|FATAL).
  character(len=*), intent(in)  :: string !< Complete message string.
  character(c_char)             :: c_string(MAXSTRINGLEN)

  integer :: levelact


   levelact = max(1,min(max_level, level))

   if (level >= thresholdLvl) then

      if (writeMessage2Screen) then
         write (*, '(a)') level_prefix(levelact)//trim(string)
      endif

      if (lunMess > 0) then

         write (lunMess, '(a)') level_prefix(levelact)//trim(string)

         ! Only count for Log-File, otherwise confusing.....
         mess_level_count(levelact) = mess_level_count(levelact) + 1

      end if

      if (level > maxErrorLevel) then
         maxErrorLevel = level
      endif

      if (useLogging) then
         call pushMessage(levelact, string)
      endif

   elseif (level < 0) then

      ! If negative level just put string to all output channels without prefix and counting
      if (writeMessage2Screen) then
         write (*, '(a)') trim(string)
      endif

      if (lunMess > 0) then
         write (lunMess, '(a)') trim(string)
      end if

   endif

   ! Optional callback routine for any user-specified actions (e.g., upon error)
   if (associated(mh_callback).and. .not. alreadyInCallback) then
      alreadyInCallback = .true.
      call mh_callback(level, trim(string)) !In future, possibly also error #ID
      alreadyInCallback = .false.
   end if

   if (associated(f_callback).and. .not. alreadyInCallback) then
      alreadyInCallback = .true.
      c_string = string_to_char_array(trim(string))
      call f_callback(level, c_string)
      alreadyInCallback = .false.
   end if

end subroutine SetMessage


!> Pushes a new message at the tail of the message queue.
subroutine pushMessage(level, string)
   integer,          intent(in)  :: level  !< One of: LEVEL_(DEBUG|INFO|WARN|ERROR|FATAL).
   character(len=*), intent(in)  :: string !< Complete message string.


   ibuffertail = mod(ibuffertail, maxmessages) + 1
   messagecount = min(maxmessages, messagecount+1)

   messages(ibuffertail)  = string
   levels(ibuffertail)    = level
end subroutine pushMessage


!> Pops the oldest message from the head of the message queue.
!! When the message queue is empty, level=LEVEL_NONE is returned.
subroutine getOldestMessage(level, msg)
   integer,          intent(out) :: level !< Set to the level of the newest message.
   character(len=*), intent(out) :: msg   !< Set to the newest message text.

   integer :: ibufferhead, itrimlen

   if (messagecount == 0) then
      level = LEVEL_NONE
      return
   else
      ! ibufferhead: front element in the message queue, is tail minus count, but notice the cyclic list/queue (mod) and +1 for 1-based index.
      ibufferhead = mod(ibuffertail - messagecount + maxmessages, maxmessages) + 1
      msg = ' '
      itrimlen = min(len(msg), len_trim(messages(ibuffertail))) ! Shortest of actual message and the available space in output variable.

      msg   = messages(ibuffertail)(1:itrimlen)
      level = levels(ibufferhead)
      messagecount = messagecount-1
   end if
end subroutine getOldestMessage


!> Returns the number of messages that are still in the message buffer queue.
!! Note: it is advised to use getOldestMessage to pop messages from the queue.
integer function getMessageCount()
   getMessageCount = messagecount
end function


!> Gets a message from the message buffer queue, without any checks on whether it's there.
!! Returns the integer level as the function's return value, and stores the accompanying message in the message dummy argument.
integer function GetMessage_MH(imessage, message)
   integer,            intent(in)  :: imessage  !< Position in the message queue.
   character(len=200), intent(out) :: message   !< The message text.

   message=messages(imessage)(1:200)
   GetMessage_MH = levels(imessage)
end function


subroutine resetMessageCount_MH()

   messageCount = 0
end subroutine

integer function getMaxErrorLevel()
   getMaxErrorLevel = maxErrorLevel
end function

subroutine resetMaxerrorLevel()
   maxErrorLevel = 0
end subroutine

subroutine message1string(level, w1)
  use iso_c_utils
    character(*)    :: w1
    integer         :: level

    integer                        :: l1
    character(MAXSTRINGLEN)                 :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    write (rec(1:), '(a)') w1(:l1)

    call setMessage(level, rec)
end subroutine message1string

subroutine message2string(level, w1, w2)
    character(*) :: w1, w2
    integer         :: level

    integer :: l1, l2
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    l2 = max(1, len_trim(w2))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(a)') w2(:l2)

    call SetMessage(level, rec)
end subroutine message2string

subroutine message3string(level, w1, w2, w3)
    character(*) :: w1, w2, w3
    integer         :: level

    integer :: l1, l2, l3
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    l2 = max(1, len_trim(w2))
    l3 = max(1, len_trim(w3))
    write (rec(1:), '(a)') w1(1:l1)
    write (rec(2 + l1:), '(a)') w2(1:l2)
    write (rec(3 + l1 + l2:), '(a)') w3(1:l3)

    call SetMessage(level, rec)
 end subroutine message3string

subroutine message4string(level, w1, w2, w3, w4)
    character(*) :: w1, w2, w3, w4
    integer         :: level

    integer :: l1, l2, l3, l4
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    l2 = max(1, len_trim(w2))
    l3 = max(1, len_trim(w3))
    l4 = max(1, len_trim(w4))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(a)') w2(:l2)
    write (rec(3 + l1 + l2:), '(a)') w3(:l3)
    write (rec(4 + l1 + l2 + l3:), '(a)') w4(:l4)

    call SetMessage(level, rec)
end subroutine message4string

subroutine message2char1real(level, w1, w2, r3)

    real :: r3
    character(*) :: w1, w2
    intent (in) r3
    integer         :: level

    integer :: l1, l2
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    l2 = max(1, len_trim(w2))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(a)') w2(:l2)
    write (rec(3 + l1 + l2:), '(f14.6)') r3

    call SetMessage(level, rec)
end subroutine message2char1real

subroutine message2char2real(level, w1, w2, r3, r4)
    real :: r3, r4
    character(*) :: w1, w2
    intent (in) r3, r4
    integer         :: level

    integer :: l1, l2
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    l2 = max(1, len_trim(w2))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(a)') w2(:l2)
    write (rec(3 + l1 + l2:), '(2f14.6)') r3, r4

    call SetMessage(level, rec)
end subroutine message2char2real

subroutine message1char1real(level, w1, r2)
    real :: r2
    character(*) :: w1
    intent (in) r2
    integer         :: level

    integer :: l1
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(F14.6)') r2

    call SetMessage(level, rec)
end subroutine message1char1real

subroutine message1char1double(level, w1, d2)
    double precision, intent(in) :: d2
    character(*) :: w1
    integer         :: level

    integer :: l1
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(F14.6)') d2

    call SetMessage(level, rec)
end subroutine message1char1double

subroutine message1char1int(level, w1, i2)
    integer :: i2
    character(*) :: w1
    intent (in) i2
    integer         :: level

    integer :: l1
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(I14)') i2

    call SetMessage(level, rec)
end subroutine message1char1int

subroutine message1char2int(level, w1, i2, i3)
    integer :: i2, i3
    character(*) :: w1
    intent (in) i2, i3
    integer         :: level

    integer :: l1
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(2I14)') i2, i3

    call SetMessage(level, rec)

end subroutine message1char2int

subroutine message2int1char(level, i1, i2, w3)
    integer :: i1, i2
    character(*) :: w3
    intent (in) i1, i2
    integer         :: level

    integer :: l3
    character(600) :: rec

    rec = ' '
    l3 = max(1, len_trim(w3))
    write (rec( 1:28), '(2I14)') i1, i2
    write (rec(30:)  , '(a)'   ) w3(:l3)


    call SetMessage(level, rec)

end subroutine message2int1char

subroutine message1char3int(level, w1, i2, i3, i4)
    integer :: i2, i3, i4
    character(*) :: w1
    intent (in) i2, i3, i4
    integer        :: level
    integer        :: l1
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(3I14)') i2, i3, i4

    call SetMessage(level, rec)
end subroutine message1char3int

subroutine message1char2real(level, w1, r2, r3)
    real         :: r2, r3
    character(*) :: w1
    integer      :: level
    intent(in)   :: level, w1, r2, r3

    integer        :: l1
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(2F14.6)') r2, r3

    call SetMessage(level, rec)
end subroutine message1char2real

subroutine message1char1int1double(level, w1, i2, d3)
    integer          :: level
    character(*)     :: w1
    integer          :: i2
    double precision :: d3

    integer :: l1
    character(600) :: rec

    rec = ' '
    l1 = max(1, len_trim(w1))
    write (rec(1:), '(a)') w1(:l1)
    write (rec(2 + l1:), '(i14)') i2
    write (rec(16 + l1:), '(F14.6)') d3

    call SetMessage(level, rec)
end subroutine message1char1int1double

subroutine message1double1int1char(level, d1, i2, w3)
    integer          :: level
    character(*)     :: w3
    integer          :: i2
    double precision :: d1

    integer :: l3
    character(600) :: rec

    rec = ' '
    l3 = max(1, len_trim(w3))
    write (rec(1 :16), '(F16.6)') d1
    write (rec(18:31), '(i14)'  ) i2
    write (rec(33:  ), '(a)'    ) w3(:l3)

    call SetMessage(level, rec)
end subroutine message1double1int1char
!-- Error interfaces ----------------------------
subroutine error4char(w1, w2, w3, w4)
    character(*) :: w1, w2, w3, w4

    call mess(LEVEL_ERROR, w1, w2, w3, w4)
end subroutine error4char

subroutine error3char(w1, w2, w3)
    character(*) :: w1, w2, w3

    call mess(LEVEL_ERROR, w1, w2, w3)
end subroutine error3char

subroutine error2char(w1, w2)
    character(*) :: w1, w2

    call mess(LEVEL_ERROR, w1, w2)
end subroutine error2char

subroutine error1char(w1)
    character(*) :: w1

    call mess(LEVEL_ERROR, w1)
end subroutine error1char

subroutine error2char1real(w1, w2, r3)
    real :: r3
    character(*) :: w1, w2

    call mess(LEVEL_ERROR, w1, w2, r3)
end subroutine error2char1real

subroutine error2char2real(w1, w2, r3, r4)
    real :: r3, r4
    character(*) :: w1, w2

    call mess(LEVEL_ERROR, w1, w2, r3, r4)
end subroutine error2char2real

subroutine error1char1real(w1, r2)
    real :: r2
    character(*) :: w1

    call mess(LEVEL_ERROR, w1, r2)
end subroutine error1char1real

subroutine error1char1double(w1, d2)
    double precision, intent(in) :: d2
    character(*) :: w1

    call mess(LEVEL_ERROR, w1, d2)
end subroutine error1char1double

subroutine error1char1int(w1, i2)
    integer :: i2
    character(*) :: w1

    call mess(LEVEL_ERROR, w1, i2)
end subroutine error1char1int

subroutine error1char2real(w1, r2, r3)
    real :: r2, r3
    character(*) :: w1

    call mess(LEVEL_ERROR, w1, r2, r3)
end subroutine error1char2real

subroutine error1char2int(w1, i2, i3)
    integer :: i2, i3
    character(*) :: w1

    call mess(LEVEL_ERROR, w1, i2, i3)
end subroutine error1char2int

subroutine error1char1int1double(w1, i2, d3)
    character(*)     :: w1
    integer          :: i2
    double precision :: d3

    call mess(LEVEL_ERROR, w1, i2, d3)
end subroutine error1char1int1double

subroutine aerr(w1, iostat, isize, errmsg) ! Allocation error          ,continue if iostat = 0   )
  character(*) :: w1
  integer      :: iostat
  integer      :: isize
  integer      :: i3
  character(len=*), optional :: errmsg
  double precision :: rmemtot
  data rmemtot/0/ ! AvD: causes access problems with multiple OpenMP threads

  if (iostat==0) then
!$OMP CRITICAL
     rmemtot = rmemtot + isize
     i3 = rmemtot*1e-6
     if (abs(isize) > 1000) then
        write (msg_mh,*) i3, isize*1e-6, ' ', w1
        call dbg_flush()
     endif
!$OMP END CRITICAL
  else

     if (present(errmsg)) then
        write (msg_mh,*) ' Allocation Error: ', w1, ', Allocate status = ', iostat, ', Integer parameter = ', isize, '=>', trim(errmsg(1:500))
     else
        write (msg_mh,*) ' Allocation Error: ', w1, ', Allocate status = ', iostat, ', Integer parameter = ', isize
     end if
     call err_flush()
     !call err ('Allocation Error  :',w1)
     !call err ('Allocate status   =',i1)
     !call err ('Integer parameter =',i2)
  endif

end subroutine aerr

!> Output the current message buffer as a 'debug' message.
subroutine dbg_flush()
! We could check on empty buffer, but we omit this to stay lightweight. [AvD]
    call mess(LEVEL_DEBUG,  msg_mh)
end subroutine dbg_flush

!!> Output the current message buffer as an 'info' message.
subroutine msg_flush()
! We could check on empty buffer, but we omit this to stay lightweight. [AvD]
    call mess(LEVEL_INFO,  msg_mh)
end subroutine msg_flush

!> Output the current message buffer as a 'warning' message.
subroutine warn_flush()
! We could check on empty buffer, but we omit this to stay lightweight. [AvD]
    call mess(LEVEL_WARN,  msg_mh)
end subroutine warn_flush

!> Output the current message buffer as an 'error' message.
subroutine err_flush()
    call mess(LEVEL_ERROR, msg_mh)

end subroutine err_flush

!subroutine fewsdiag_init()
!end subroutine fewsdiag_init
!
!
!!--------------------------------------------------------------------------------------------------
!
!
subroutine progress(message, fraction)
  use iso_c_binding
  use iso_c_utils
  implicit none

  character(len=*), intent(in) :: message
  double precision, intent(in) :: fraction

  ! call the registered progress bar
  if (associated(progress_callback)) then
     call progress_callback(message, fraction)
  end if

  ! call the c callback if registered
  if (associated(progress_c_callback)) then
     call progress_c_callback(string_to_char_array(message), fraction)
  end if

end subroutine progress
!
subroutine set_progress_callback(callback)
  ! set the progress handler
  procedure(progress_iface) :: callback

  ! TODO check if we need cptr2fptr
  progress_callback => callback

end subroutine set_progress_callback

subroutine set_progress_c_callback(c_callback) bind(C, name="set_progress_c_callback")
  !DEC$ ATTRIBUTES DLLEXPORT:: set_progress_c_callback

  use iso_c_binding
  use iso_c_utils
  implicit none
  type(c_funptr) :: c_callback

  ! Set a callback that will be cauled with new messages

  call c_f_procpointer(c_callback, progress_c_callback)
end subroutine set_progress_c_callback




end module MessageHandling

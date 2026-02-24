module ModEM_logger

    use utilities
    use Declaration_MPI

implicit none

    character (len=512), private ::  log_fname
    integer, private :: log_fid = 0

contains

subroutine ModEM_log_init(mainOnly)

      implicit none

      logical, optional, intent(in) :: mainOnly 
      character(len=*), parameter :: log_str_fmt = '(A,I4.4,A)'

      logical :: mainOnly_lcl

      if (present(mainOnly)) then
        mainOnly_lcl = mainOnly
      else
        mainOnly_lcl = .true.
      end if

      if ((taskid == 0 .and. mainOnly_lcl) .or. (.not. mainOnly_lcl)) then
          write(log_fname, log_str_fmt) 'log.', taskid, '.modem.out'
          open(newunit=log_fid, file=log_fname, status='replace')
          call ModEM_log("Log Initalized $i - "//trim(log_fname), intArgs=(/taskid/), mainOnly=mainOnly)
      end if

end Subroutine ModEM_log_init

! The following License applies to ModEM_Log:
!
! Copyright (c) 2013-2019,  Los Alamos National Security, LLC (LANS) (Ocean: LA-CC-13-047;
! Land Ice: LA-CC-13-117) and the University Corporation for Atmospheric Research (UCAR).
!
! All rights reserved.
!
! LANS is the operator of the Los Alamos National Laboratory under Contract No.
! DE-AC52-06NA25396 with the U.S. Department of Energy.  UCAR manages the National
! Center for Atmospheric Research under Cooperative Agreement ATM-0753581 with the
! National Science Foundation.  The U.S. Government has rights to use, reproduce,
! and distribute this software.  NO WARRANTY, EXPRESS OR IMPLIED IS OFFERED BY
! LANS, UCAR OR THE GOVERNMENT AND NONE OF THEM ASSUME ANY LIABILITY FOR THE USE
! OF THIS SOFTWARE.  If software is modified to produce derivative works, such
! modified software should be clearly marked, so as not to confuse it with the
! version available from LANS and UCAR.
!
! Additionally, redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1) Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2) Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation and/or
! other materials provided with the distribution.
!
! 3) None of the names of LANS, UCAR or the names of its contributors, if any, may
! be used to endorse or promote products derived from this software without
! specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!***********************************************************************
!
!  routine ModEM_log
!
!> \brief   Writes a message to the log file
!> \author  Matt Hoffman - Modified for ModEM by Miles Curry
!> \date    14 February 2017, 5 Feburary 2026
!> \details
!>  This routine writes a message to the log file.  The message is the only
!>  required argument.  A number of optional arguments control additional
!>  behavior:
!>    intArgs, int8Args, realArgs, logicArgs:  arrays of variable values to be inserted into the
!>                                   message to replace the following characters: $i, $j, $r, $l
!>                                   See routine expand_string below for details.
!>    fid: (Optional) Instead of using the log created/opened in ModEM_log_init, write to the unit
!>         in found in fid
!>    mainOnly: flag indicating only the master task should write this message,
!>                regardless of if all tasks have open log files.
!>    flush_log: flag indicating that the log should be flushed after writing
!>                
!
!-----------------------------------------------------------------------
subroutine ModEM_log(msg, intArgs, realArgs, logicArgs, fid, mainOnly, flush_log)

   implicit none

   character(len=*), intent(in) :: msg
   logical, intent(in), optional :: mainOnly
   integer, dimension(:), intent(in), optional :: intArgs  !< Input: integer variable values to insert into message
   real(kind=prec), dimension(:), intent(in), optional :: realArgs  !< Input: real variable values to insert into message
        !< Input: exponential notation variable values to insert into message
   logical, dimension(:), intent(in), optional :: logicArgs  !< Input: logical variable values to insert into message
   integer, intent(in), optional :: fid
   logical, intent(in), optional :: flush_log

   logical :: mainOnly_lcl
   logical :: flush_lcl
   integer :: fid_lcl

   character(len=512) :: messageExpanded !< message after expansion of $ variable insertions

   if (present(mainOnly)) then
       mainOnly_lcl = mainOnly
   else
       mainOnly_lcl = .true.
   end if

   if (present(fid)) then
        fid_lcl = fid
   else
        fid_lcl = log_fid
   end if


   if (present(flush_log)) then
       flush_lcl = flush_log
   else
        flush_lcl = .false.
   end if

   call expand_string(msg, messageExpanded, intArgs, logicArgs, realArgs)


   if ((mainOnly_lcl .and. taskid == 0) .or. (.true.)) then
       write(0,*) fid_lcl
       write(fid_lcl,*) trim(messageExpanded)

       if (.true.) then
           call flush(fid_lcl)
       end if
    end if

end subroutine ModEM_log

end module ModEM_logger


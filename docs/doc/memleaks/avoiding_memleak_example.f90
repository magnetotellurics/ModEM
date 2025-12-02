! test_chain --
! By Arjen Markus <arjen.markus@wldelft.nl>
!    Test program to see if memory leaks originating from derived-types
!    can be circumvented
!    The idea:
!    - The user is responsible for cleaning up his/her own variables
!    - The module is responsible for cleaning up its intermediate
!      results (flagged as "temporary")
!
!    Note:
!    - seqno and alloc_seq are only used for debugging purposes
!
module chains
   type chain
      integer, dimension(:), pointer :: values => null()
      logical                        :: tmp    =  .false.
      integer                        :: seqno  =  0
   end type chain

   integer :: seqno = 0
   logical, dimension(1:100) :: alloc_seq = .false.

   interface assignment(=)
      module procedure assign_chain
      module procedure assign_array
   end interface

   interface operator(.concat.)
      module procedure concat_chain
   end interface

contains

subroutine assign_array( ic, jc )
   type(chain),intent(out)           :: ic
   integer, dimension(:), intent(in) :: jc

   call cleanup( ic, .false. )

   allocate( ic%values(1:size(jc)) )
   seqno = seqno + 1
   alloc_seq(seqno) = .true.
   ic%values = jc
   ic%seqno  = seqno
   ic%tmp    = .false.

end subroutine assign_array

subroutine assign_chain( ic, jc )
   type(chain), intent(inout) :: ic
   type(chain), intent(in)    :: jc

   call cleanup( ic, .false. )

   allocate( ic%values(1:size(jc%values)) )
   seqno = seqno + 1
   alloc_seq(seqno) = .true.
   ic%values = jc%values
   ic%seqno  = seqno
   ic%tmp    = .false.

   call cleanup( jc, .true. )

end subroutine assign_chain

function concat_chain( ic, jc )
   type(chain), intent(in) :: ic
   type(chain), intent(in) :: jc
   type(chain)             :: concat_chain

   integer                 :: nic
   integer                 :: njc

   nic = size(ic%values)
   njc = size(jc%values)

   allocate( concat_chain%values(1:nic+njc) )
   seqno = seqno + 1
   alloc_seq(seqno) = .true.

   concat_chain%values(1:nic)         = ic%values(1:nic)
   concat_chain%values(nic+1:nic+njc) = jc%values(1:njc)
   concat_chain%seqno  = seqno
   concat_chain%tmp    = .true.

   call cleanup( ic, .true. )
   call cleanup( jc, .true. )

end function concat_chain

subroutine cleanup( ic, only_tmp )
   type(chain) :: ic
   logical, optional :: only_tmp

   logical :: clean_tmp

   clean_tmp = .false.
   if ( present(only_tmp) ) clean_tmp = only_tmp

   if ( .not. clean_tmp .or. ic%tmp ) then
      if ( associated( ic%values) ) deallocate( ic%values )
      if ( ic%seqno .gt. 0 ) alloc_seq(ic%seqno) = .false.
   endif

end subroutine cleanup

subroutine report_chain

   integer :: i
   integer :: count

   count = 0
   do i = 1,size(alloc_seq)
      if ( alloc_seq(i) ) then
         write(*,*) 'Allocated item', i
         count = count + 1
      endif
   enddo
   write(*,*) 'Number of allocated items:', count

end subroutine report_chain

end module chains

program test_chain
   use chains

   type(chain) :: ic
   type(chain) :: jc
   type(chain) :: kc

   ic = (/1,2,3/)
   jc = (/4,5/)
   kc = ic .concat. jc
   call report_chain
   kc = jc .concat. ic
   call report_chain
   call cleanup(ic)
   call cleanup(jc)
   call cleanup(kc)
   call report_chain

end program test_chain

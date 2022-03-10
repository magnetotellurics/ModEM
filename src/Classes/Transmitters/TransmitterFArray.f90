!*************
!
! Class to provide a dynamic and polymorphic transmitters of Transmitter_t objects
!
! Last modified at 03/2021 by Paulo Werdt
!
!*************
!
module TransmitterFArray
    !
	use Constants
	!
	use TransmitterMT
    use TransmitterCSEM
	!
    implicit none
	!
	! Allocatable Transmitter element of the array, for Old Fortran polymorphism !!!
    type, public :: Tx_t
        !
        class( Transmitter_t ), allocatable :: Tx
        !
    end type Tx_t
	!
	!
	!type( Tx_t ), pointer, save, public, dimension(:)	:: transmitters
	type( Tx_t ), pointer, dimension(:)	:: transmitters
    !
	integer                             	:: index_of_tx_array
    !
	public :: getTransmitter, addTransmitter, showTransmitterArray
    public :: initTransmitterArray, updateTransmitterArray
    public :: IsTxExist
    !
contains
    !
    ! Initializes and allocates the transmitters
	subroutine initTransmitterArray()
        implicit none
		!
		allocate( transmitters( 1 ) )
		!
        index_of_tx_array = 0;
        !
    end subroutine initTransmitterArray
    !
    ! Add a new Transmitter_t and initialize it if necessary
    subroutine updateTransmitterArray( new_tx )
        implicit none
		!
        class( Transmitter_t ), intent( in )					:: new_tx
		!
        integer													:: iTx, nTx
        !
		if( index_of_tx_array == 0 ) then
            call initTransmitterArray()
            call addTransmitter( new_tx )
            return
        end if
		!
        nTx = size( transmitters )
        do iTx = 1, nTx
            if ( IsTxExist( new_tx, transmitters( iTx )%Tx ) ) then
                return
            end if
        end do
		!
        call addTransmitter( new_tx )
        !
    end subroutine updateTransmitterArray
    !
    ! Compare two transmitters: returns true if they are the same 
    function IsTxExist( Txa, Txb ) result ( YESNO )
        class( Transmitter_t ), intent( in )	:: Txa, Txb
        logical                         		:: YESNO
		!
        real( kind=prec )                 		:: a_xyz(3), b_xyz(3)
		!
        YESNO = .FALSE.
		!
        select type( Txa )
            !
            class is ( TransmitterCSEM_t )
			    !
				select type( Txb )
				!
                class is ( TransmitterCSEM_t )
			       !
				   a_xyz = Txa%location
				   b_xyz = Txb%location
				   if( ABS( Txa%period - Txb%period ) < TOL6 .AND.	&
					   ABS( a_xyz(1) - b_xyz(1) ) < TOL6 .AND.   				&
					   ABS( a_xyz(2) - b_xyz(2) ) < TOL6 .AND.   				&
					   ABS( a_xyz(3) - b_xyz(3) ) < TOL6 .AND.   				&
					   ABS( Txa%azimuth - Txb%azimuth ) < TOL6 ) then
						   YESNO = .true.
				   end if
				class is ( TransmitterMT_t )
				   YESNO = .false.
				class default
				   YESNO = .false.
                !
                end select
			class is ( TransmitterMT_t )
			    !
				select type( Txb )
				!
                class is ( TransmitterMT_t )
			       !
				   if( ABS( Txa%period - Txb%period ) < TOL6 ) then
					   YESNO = .true.
				   end if
			       !
				class is ( TransmitterCSEM_t )
				   YESNO = .false.
				class default
				   YESNO = .false.
                !
                end select
        end select
		!
    end function IsTxExist
	!
    function getTransmitter( iTx ) result( tx )
        !
        integer										:: iTx
		!
		class( Transmitter_t ), pointer			:: tx
		!
		tx => transmitters( iTx )%Tx
		!
    end function getTransmitter
	!
    ! Add one Transmitter_t, increasing the size of the transmitters
    subroutine addTransmitter( new_tx )
        implicit none
		!
        class( Transmitter_t ), intent( in )        	:: new_tx
		!
		type( Tx_t ), allocatable, dimension(:)			:: temp
        integer                                         :: istat
        type( Tx_t ), allocatable                 		:: temp_tx
        !
        if( index_of_tx_array == 0 ) then
            allocate( transmitters(1), STAT = istat )
            transmitters(1)%Tx = new_tx
			transmitters(1)%Tx%id = 1
        else
            allocate( temp( index_of_tx_array + 1 ), STAT=istat )
            temp( 1 : size( transmitters ) ) = transmitters
            allocate( Tx_t::temp_tx )
            temp_tx%Tx = new_tx
            temp( index_of_tx_array + 1 ) = temp_tx
			temp( index_of_tx_array + 1 )%Tx%id = index_of_tx_array + 1
            deallocate( transmitters,STAT=istat )
            allocate( transmitters( index_of_tx_array + 1 ), STAT=istat )
            transmitters = temp
            deallocate( temp, STAT=istat)
        end if
        !
        index_of_tx_array = index_of_tx_array + 1
        !
    end subroutine addTransmitter
    !
    ! Prints the content of the transmitters on screen
    subroutine showTransmitterArray()
        integer                                     :: index_of_tx_array
		class( Tx_t ), allocatable					:: alloc_tx
        !
        print *, size( transmitters ), " TransmitterFArray_t"
        !
        do index_of_tx_array=1, size( transmitters )
			alloc_tx = transmitters( index_of_tx_array )
            call alloc_tx%Tx%write()
        end do
		!
    end subroutine showTransmitterArray
	!
end module TransmitterFArray

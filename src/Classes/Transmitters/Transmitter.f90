! *************
! 
! Base class to define a Transmitter
!
! *************
! 
module Transmitter
    ! 
    use Constants
    use Source
    use ForwardSolver
    use ModelParameter
    !
    ! Global name for e_solution file
    character(:), allocatable :: e_solution_file_name
    !
    type, abstract :: Transmitter_t
        !
        integer :: id, n_pol, fwd_key(8)
        !
        real( kind=prec ) :: period
        !
        class( ForwardSolver_t ), pointer :: forward_solver
        !
        class( Source_t ), allocatable :: source
        !
        class( Vector_t ), allocatable, dimension(:) :: e_all
        !
        integer, allocatable, dimension(:) :: receiver_indexes
        ! !
        ! procedure( interface_p_mult_tx ), pointer, nopass :: pMult_ptr
        ! !
        ! procedure( interface_p_mult_t_tx ), pointer, nopass :: pMult_t_ptr
        ! !
    contains
        !
        procedure, public :: init     => initializeTx
        procedure, public :: dealloc  => deallocateTx
        !
        procedure, public :: updateFwdKey => updateFwdKeyTx
        !
        procedure, public :: updateReceiverIndexesArray
        !
        procedure( interface_solve_fwd_tx ), deferred, public :: solveFWD
        !
        procedure( interface_is_equal_tx ), deferred, public :: isEqual
        !
        procedure( interface_print_tx ), deferred, public :: print
        !
        procedure, public :: pMult => pMultTx
        !
        procedure, public :: pMult_t => pMult_t_Tx
        !
    end type Transmitter_t
    !
    abstract interface
        !
        subroutine interface_solve_fwd_tx( self )
            import :: Transmitter_t
            class( Transmitter_t ), intent( inout ) :: self
        end subroutine interface_solve_fwd_tx
        !        !
        function interface_is_equal_tx( self, other ) result( equal )
            import :: Transmitter_t
            class( Transmitter_t ), intent( in ) :: self, other
            logical                              :: equal
        end function interface_is_equal_tx
        !        
        subroutine interface_print_tx( self )
            import :: Transmitter_t
            class( Transmitter_t ), intent(in) :: self
        end subroutine interface_print_tx
        ! !
        ! pure subroutine interface_p_mult_tx( m0, dm, bSrc )
            ! import :: ModelParameter_t, Source_t
            ! class( ModelParameter_t ), intent( in ) :: m0, dm
            ! class( Source_t ), intent( inout )      :: bSrc
        ! end subroutine interface_p_mult_tx
        ! !
        ! pure subroutine interface_p_mult_t_tx( m0, eSens, d_m )
            ! import :: ModelParameter_t, Vector_t
            ! class( ModelParameter_t ), intent( in )    :: m0
            ! class( Vector_t ), intent( in )           :: eSens
            ! class( ModelParameter_t ), intent( inout ) :: d_m
        ! end subroutine interface_p_mult_t_tx
        ! !
    end interface
    !
    contains
        !
        subroutine initializeTx( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            self%id = 0
            self%n_pol = 0
            call self%updateFwdKey()
            !
            self%period = 0.0
            !
            self%forward_solver => null()
            !
        end subroutine initializeTx
        !
        subroutine deallocateTx( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            if( allocated( self%source ) ) deallocate( self%source )
            !
            if( allocated( self%e_all ) ) deallocate( self%e_all )
            !
            if( allocated( self%receiver_indexes ) ) deallocate( self%receiver_indexes )
            !
        end subroutine deallocateTx
        !
        subroutine updateFwdKeyTx( self )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            !
            call date_and_time( values=self%fwd_key )
            !
        end subroutine updateFwdKeyTx
        !
        subroutine updateReceiverIndexesArray( self, new_int )
            implicit none
            !
            class( Transmitter_t ), intent( inout ) :: self
            integer, intent( in )                   :: new_int
            !
            integer, allocatable, dimension(:)      :: temp_array
            integer                                 :: nRi, idx
            !
            if( .NOT. allocated( self%receiver_indexes ) ) then
                allocate( self%receiver_indexes(1) )
                self%receiver_indexes(1) = new_int
            else
                !
                nRi = size( self%receiver_indexes )
                !
                do idx = 1, nRi
                    if ( new_int == self%receiver_indexes( idx ) ) then
                        return
                    end if
                end do
                !
                allocate( temp_array( nRi + 1 ) )
                temp_array( 1 : nRi ) = self%receiver_indexes(:)
                temp_array( nRi + 1 ) = new_int
                !
                if( allocated( self%receiver_indexes ) ) deallocate( self%receiver_indexes )
                self%receiver_indexes = temp_array
                !
                deallocate( temp_array )
                !
            endif
            !
        end subroutine updateReceiverIndexesArray
        !
        subroutine pMultTx( self, m0, dm, bSrc )
            implicit none
            !
            class( Transmitter_t ), intent( in )    :: self
            class( ModelParameter_t ), intent( in ) :: m0
            class( ModelParameter_t ), intent( in ) :: dm
            !
            class( Source_t ), intent( inout ) :: bSrc
            !
            complex( kind=prec ) :: miwm
            class( ModelParameter_t ), allocatable :: temp
            logical :: adjt
            integer :: k
            class( Vector_t ), allocatable:: eVec
            !
            ! WHERE THE HELL AM I GOING TO GET PERIOD ????
            miwm = -ONE_I * MU_0 * isign * cmplx( 0.0, 1./ ( 2.0 * PI / self%period ), kind=prec )
            !
            allocate( temp, source = m0 )
            !
            ! WHAT TO DO WITH eVec????
            call temp%dPDEmapping( dm, eVec )
            !
            adjt = .FALSE.
            !
            ! AND NPOL ????
            do k = 1, self%n_pol
                !
                ! ADJOINT SOURCE ????
                !
            enddo
            !
            ! MATLAB IMPLEMENTATION
            !
            !miwm = -1i*Tx.fwd.modOp.mu0*Tx.fwd.isign*Tx.omega;
            !temp = m0.dPDEmapping(dm);
            !
            !adjt = false;
            !
            !bSrc(2) = TSourceInteriorForce(Tx.fwd.modOp,adjt);
            !bSrc(2).SetSourceParams(miwm*Tx.e(nPol).*temp);
            !for k = 1:Tx.nPol-1
                !bSrc(k) = TSourceInteriorForce(Ts.fwd.modOp,adjt);
                !bSrc(k).SetSourceParams(miwm.*Tx.e(k).*temp);
            !end
            !
        end subroutine pMultTx
        !
        subroutine pMult_t_Tx( self, m0, eSens, d_m )
            implicit none
            !
            class( Transmitter_t ), intent( in )                    :: self
            class( ModelParameter_t ), intent( in )                 :: m0
            class( Vector_t ), intent( inout )                     :: eSens(:)
            class( ModelParameter_t ), allocatable, intent( inout ) :: d_m
            !
            complex( kind=prec ) :: miwm
            class( ModelParameter_t ), allocatable :: temp
            logical :: adjt
            integer :: k
            !
            !
            miwm = -ONE_I * MU_0 * isign * cmplx( 0.0, 1./ ( 2.0 * PI / self%period ), kind=prec )
            !
            allocate( d_m, source = m0%dPDEmappingT( eSens(1) ) )
            !
            ! MATLAB IMPLEMENTATION
            !
            !Tx = dTx.Tx;  %  transmitter for this DataVectorTx object;  As for Pmult_E
            !
            !miwm = -1i*Tx.fwd.modOp.mu0*Tx.fwd.isign*Tx.omega;
            !eSens(1) = miwm.*Tx.e(1).*eSens(1);
            !
            !for k = 2:Tx.nPol
                !eSens(1) = eSens(1) + miwm.*Tx.e(k).*eSens(k);
                !end
            !d_m = m0.dPDEmappingT(eSens(1));
            !
        end subroutine pMult_t_Tx
        !
        ! ! PMult
        ! elemental subroutine pMultTx( self, m0, dm, bSrc )
            ! implicit none
            ! !
            ! class( Transmitter_t ), intent( in )    :: self
            ! class( ModelParameter_t ), intent( in ) :: m0, dm
            ! class( Source_t ), intent( inout ) :: bSrc
            ! !
            ! call self%pMult_ptr( m0, dm, bSrc )
            ! !
        ! end subroutine pMultTx
        ! !
        ! pure subroutine pMult_E( m0, dm, bSrc )
            ! implicit none
            ! !
            ! class( ModelParameter_t ), intent( in ) :: m0
            ! class( ModelParameter_t ), intent( in ) :: dm
            ! !
            ! class( Source_t ), intent( inout ) :: bSrc
            ! !
            ! complex( kind=prec ) :: miwm
            ! class( ModelParameter_t ), allocatable :: temp
            ! logical :: adjt
            ! integer :: k
            ! class( Vector_t ), allocatable:: eVec
            ! !
            ! ! WHERE THE HELL AM I GOING TO GET PERIOD ????
            ! !miwm = -ON_I * MU_0 * isign * cmplx( 0.0, 1./ ( 2.0 * PI / self%period ), kind=prec )
            ! !
            ! allocate( temp, source = m0 )
            ! !
            ! ! WHAT TO DO WITH eVec????
            ! call temp%dPDEmapping( dm, eVec )
            ! !
            ! adjt = .FALSE.
            ! !
            ! ! AND NPOL ????
            ! !do k = 1, self%npol
                ! !
                ! ! ADJOINT SOURCE ????
                ! !
            ! !enddo
            ! !
            ! ! MATLAB IMPLEMENTATION
            ! !
            ! !miwm = -1i*Tx.fwd.modOp.mu0*Tx.fwd.isign*Tx.omega;
            ! !temp = m0.dPDEmapping(dm);
            ! !
            ! !adjt = false;
            ! !
            ! !bSrc(2) = TSourceInteriorForce(Tx.fwd.modOp,adjt);
            ! !bSrc(2).SetSourceParams(miwm*Tx.e(nPol).*temp);
            ! !for k = 1:Tx.nPol-1
                ! !bSrc(k) = TSourceInteriorForce(Ts.fwd.modOp,adjt);
                ! !bSrc(k).SetSourceParams(miwm.*Tx.e(k).*temp);
            ! !end
            ! !
        ! end subroutine pMult_E
        ! !
        ! ! PMult_t
        ! elemental subroutine pMult_t_Tx( self, m0, eSens, d_m )
            ! implicit none
            ! !
            ! class( Transmitter_t ), intent( in )    :: self
            ! class( ModelParameter_t ), intent( in ) :: m0
            ! class( Vector_t ), intent( in )        :: eSens
            ! !
            ! class( ModelParameter_t ), intent( inout ) :: d_m
            ! !
            ! call self%pMult_t_ptr( m0, eSens, d_m )
            ! !
        ! end subroutine pMult_t_Tx
        ! !
        ! pure subroutine pMult_t_E( m0, eSens, d_m )
            ! implicit none
            ! !
            ! class( ModelParameter_t ), intent( in )    :: m0
            ! class( Vector_t ), intent( in )           :: eSens
            ! !
            ! class( ModelParameter_t ), intent( inout ) :: d_m
            ! !
            ! complex( kind=prec ) :: miwm
            ! class( ModelParameter_t ), allocatable :: temp
            ! logical :: adjt
            ! integer :: k
            ! !
            ! ! WHERE THE HELL AM I GOING TO GET PERIOD ????
            ! !miwm = -ON_I * MU_0 * isign * cmplx( 0.0, 1./ ( 2.0 * PI / self%period ), kind=prec )
            ! !
            ! d_m = m0%dPDEmappingT( eSens(1) )
            ! !
            ! ! MATLAB IMPLEMENTATION
            ! !
            ! !Tx = dTx.Tx;  %  transmitter for this DataVectorTx object;  As for Pmult_E
            ! !
            ! !miwm = -1i*Tx.fwd.modOp.mu0*Tx.fwd.isign*Tx.omega;
            ! !eSens(1) = miwm.*Tx.e(1).*eSens(1);
            ! !
            ! !for k = 2:Tx.nPol
                ! !eSens(1) = eSens(1) + miwm.*Tx.e(k).*eSens(k);
                ! !end
            ! !d_m = m0.dPDEmappingT(eSens(1));
            ! !
        ! end subroutine pMult_t_E
        ! !
end module Transmitter

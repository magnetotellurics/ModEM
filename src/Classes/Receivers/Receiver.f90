!*************
!
! Base class to define a Receiver
!
! Last modified at 20/10/2021 by Paulo Werdt
!
!*************
!
module Receiver
    !
    use FileUnits
    use Transmitter
    use cVector3D_SG
    use ModelOperator
    use DataGroupArray
    use PredictedDataHandle
    use Grid3D_SG
    !
    type, abstract :: Receiver_t
        !
        integer                               :: id, n_comp
        !
        character(:), allocatable             :: code
        !
        real( kind=prec )                     :: location(3)
        !
        logical                               :: is_complex = .FALSE., interpolation_set = .FALSE.
        !
        class( Grid_t ), pointer              :: grid
        !
        character(2), allocatable             :: EHxy(:)
        !
        character(3), allocatable             :: comp_names(:)
        !
        complex( kind=prec ), allocatable  :: I_BB(:,:), EE(:,:), Z(:)
        !
        class( cVector_t ), allocatable     :: Lex, Ley, Lez, Lbx, Lby, Lbz
        !
        class( DataGroupArray_t ), pointer :: data_groups
        !
        type( PredictedDataHandle_t ), allocatable :: predicted_data_entries(:)
        !
    contains
        !
        ! DEFERRED INTERFACES
        procedure( interface_predicted_data ), deferred, public :: predictedData
        !
        procedure( interface_write_rx ), deferred, public                     :: write
        !
        ! CLASS PROCEDURES
        procedure, public :: evaluationFunction => evaluationFunctionRx
        !
        procedure, public :: init => initializeRx
        procedure, public :: dealloc => deallocateRx
        !
        procedure, public :: isEqual => isEqualRx
        !
        procedure, public :: has => hasDataGroupRx
        procedure, public :: add => addDataGroupRx
        procedure, public :: get => getDataGroupRx
        procedure, public :: getNdg => getNumberOfDataGroupRx
        !
        procedure, public :: updatePredictedDataArray
        procedure, public :: savePredictedData
        procedure, public :: writePredictedData
        !
    end type Receiver_t
    !
    abstract interface
        !
        subroutine interface_predicted_data( self, model_operator, transmitter )
            !
            import :: Receiver_t, ModelOperator_t, Transmitter_t
            !
            class( Receiver_t ), intent( inout )  :: self
            class( ModelOperator_t ),intent( in ) :: model_operator
            class( Transmitter_t ), intent( in )  :: transmitter
            !
        end subroutine interface_predicted_data
        !
        subroutine interface_save_predicted_data_rx( self, tx )
            !
            import :: Receiver_t, Transmitter_t
            !
            class( Receiver_t ), intent( in )     :: self
            class( Transmitter_t ), intent( in ) :: tx
            !
        end subroutine interface_save_predicted_data_rx
        !
        subroutine interface_write_predicted_data_rx( self )
            !
            import :: Receiver_t
            !
            class( Receiver_t ), intent( in ) :: self
            !
        end subroutine interface_write_predicted_data_rx
        !
        subroutine interface_write_rx( self )
            !
            import :: Receiver_t
            !
            class( Receiver_t ), intent(in) :: self
            !
        end subroutine interface_write_rx
        !
    end interface
    !
contains
    !
    subroutine initializeRx( self )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        !
        self%id = 0
        !
        self%n_comp = 0
        !
        self%location = 0.0
        !
        self%grid => null()
        !
        allocate( self%data_groups, source = DataGroupArray_t() )
        !
    end subroutine initializeRx
    !
    subroutine deallocateRx( self )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        !
        if( associated( self%grid ) ) deallocate( self%grid )
        !
        if( allocated( self%EHxy ) ) deallocate( self%EHxy )
        !
        if( allocated( self%I_BB ) ) deallocate( self%I_BB )
        !
        if( allocated( self%EE ) ) deallocate( self%EE )
        !
        if( allocated( self%Z ) ) deallocate( self%Z )
        !
        if( allocated( self%Lex ) ) deallocate( self%Lex )
        !
        if( allocated( self%Ley ) ) deallocate( self%Ley )
        !
        if( allocated( self%Lez ) ) deallocate( self%Lez )
        !
        if( allocated( self%Lbx ) ) deallocate( self%Lbx )
        !
        if( allocated( self%Lby ) ) deallocate( self%Lby )
        !
        if( allocated( self%Lbz ) ) deallocate( self%Lbz )
        !
        !deallocate( self%data_groups )
        !
    end subroutine deallocateRx
    !
    subroutine evaluationFunctionRx( self, model_operator, omega )
        implicit none
        !
        class( Receiver_t ), intent( inout )    :: self
        class( ModelOperator_t ), intent( in ) :: model_operator
        real( kind=prec ), intent( in )          :: omega
        !
        integer                  :: k
        complex( kind=prec ) :: comega
        !
        ! THESE SHOULD BE rVECTORS ????
        class( cVector_t ), allocatable :: e, h
        class( cVector_t ), allocatable :: lh
        !
        !comega = cmplx( 1./omega, 0.0, kind=prec )
        comega = cmplx( 0.0, 1./omega, kind=prec )
        !
        do k = 1, self%n_comp
            !
            selectcase( self%EHxy(k) )
                !
                case( "Ex", "Ey" )
                    !
                     select type( grid => model_operator%metric%grid )
                          class is( Grid3D_SG_t )
                                if( .not. allocated( e ) ) allocate( e, source = cVector3D_SG_t( grid, EDGE ) )
                          class default
                                stop "Receiver: Unclassified model_operator%metric%grid for e"
                     end select
                    !
                case( "Bx", "By", "Bz" )
                    !
                     select type( grid => model_operator%metric%grid )
                          class is( Grid3D_SG_t )
                                if( .not. allocated(h) ) allocate( h, source = cVector3D_SG_t( grid, FACE ) )
                          class default
                                stop "Receiver: Unclassified model_operator%metric%grid for h"
                     end select
                    !
            end select
            !
        end do
        !
        do k = 1, self%n_comp
            !
            if( allocated( lh ) ) deallocate( lh )
            !
            selectcase( self%EHxy(k) )
                !
                case( "Ex" )
                    call e%interpFunc( self%location, "x", self%Lex )
                !
                case( "Ey" )
                    call e%interpFunc( self%location, "y", self%Ley )
                !
                case( "Ez" )
                    call e%interpFunc( self%location, "z", self%Lez )
                !
                case( "Bx" )
                    !
                    call h%interpFunc( self%location, "x", lh )
                    !
                    select type( lh )
                        class is(cVector3D_SG_t)
                            if( allocated( self%Lbx ) ) deallocate( self%Lbx )
                            allocate( self%Lbx, source = cVector3D_SG_t( lh%grid, EDGE ) )
                            !
                        class default
                            write(*, *) 'ERROR:Receiver::evaluationFunction:'
                            stop          '            Unkonow lh type'
                    end select
                    !
                    call model_operator%multCurlT( lh, self%Lbx )
                    !call self%Lbx%mults( isign * ONE_I / comega )
                call self%Lbx%mults( isign * comega )
                    !
                case( "By" )
                    ! 
                    call h%interpFunc( self%location, "y", lh )
                    !
                    select type( lh )
                        class is(cVector3D_SG_t)
                            if( allocated( self%Lby ) ) deallocate( self%Lby )
                            allocate( self%Lby, source = cVector3D_SG_t( lh%grid, EDGE ) )
                            !
                        class default
                            write(*, *) 'ERROR:Receiver::evaluationFunction:'
                            stop          '            Unkonow lh type'
                    end select
                    !
                    call model_operator%multCurlT( lh, self%Lby )
                    !call self%Lby%mults( isign * ONE_I / comega )
                call self%Lby%mults( isign * comega )
                !
                case( "Bz" )
                    !
                    call h%interpFunc( self%location, "z", lh )
                    !
                    select type( lh )
                        class is(cVector3D_SG_t)
                            if( allocated( self%Lbz ) ) deallocate( self%Lbz )
                            allocate( self%Lbz, source = cVector3D_SG_t( lh%grid, EDGE ) )
                            !
                        class default
                            write(*, *) 'ERROR:Receiver::evaluationFunction:'
                            stop          '            Unkonow lh type'
                    end select
                    !
                    call model_operator%multCurlT( lh, self%Lbz )
                    !call self%Lbz%mults( isign * ONE_I / comega )
                call self%Lbz%mults( isign * comega )
                !
            end select
            !
        end do
        !
        if( allocated( e ) ) deallocate( e )
        if( allocated( h ) ) deallocate( h )
        !
    end subroutine evaluationFunctionRx
    !
    function isEqualRx( self, other ) result( equal )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: self
        class( Receiver_t ), intent( in ) :: other
        !
        logical :: equal
        !
        equal = .FALSE.
        !
        if( self%location(1) == other%location(1) .AND.    &
             self%location(2) == other%location(2) .AND.    &
             self%location(3) == other%location(3) ) then
            equal = .TRUE.
        endif
        !
    end function isEqualRx
    !
    function hasDataGroupRx( self, data_group ) result( found )
        implicit none
        !
        class( Receiver_t ), intent( in )  :: self
        class( DataGroup_t ), intent( in ) :: data_group
        !
        logical :: found
        integer :: iDg, nDg
        !
        found = .FALSE.
        !
        nDg = self%data_groups%size()
        !
        do iDg = 1, nDg
            !
            if( data_group%isEqual( self%data_groups%get( iDg ) ) ) then
                found = .TRUE.
            end if
        end do
        !
    end function hasDataGroupRx
    !
    subroutine addDataGroupRx( self, data_group )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        class( DataGroup_t ), intent( in )    :: data_group
        !
        call self%data_groups%add( data_group )
        !
    end subroutine addDataGroupRx
    !
    function getDataGroupRx( self, index ) result( data_group )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: self
        integer, intent( in )                 :: index
        class( DataGroup_t ), allocatable :: data_group
        !
        data_group = self%data_groups%get( index )
        !
    end function getDataGroupRx
    !
    function getNumberOfDataGroupRx( self ) result( counter )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: self
        integer                           :: counter
        !
        counter = self%data_groups%size()
        !
    end function getNumberOfDataGroupRx
    !
    !
    subroutine savePredictedData( self, tx )
        implicit none
        !
        class( Receiver_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: tx
        !
        character(:), allocatable :: code, component
        real( kind=prec )         :: period, real_part, imaginary
        real( kind=prec )         :: xyz(3)
        integer                   :: i
        !#Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
        !
        if( allocated( self%predicted_data_entries ) ) deallocate( self%predicted_data_entries )
        !
        do i = 1, self%n_comp
             !
             period = real( tx%period, kind=prec )
             code = trim( self%code )
             xyz = (/real( self%location( 1 ), kind=prec ), real( self%location( 2 ), kind=prec ), real( self%location( 3 ), kind=prec )/)
             component = trim( self%comp_names( i ) )
             real_part = real( self%Z( i ), kind=prec )
             imaginary = real( aimag( self%Z( i ) ), kind=prec )
             !
             call self%updatePredictedDataArray( buildPredictedDataHandle( code, component, period, xyz, real_part, imaginary ) )
             !
        enddo
        !
    end subroutine savePredictedData
    !
    subroutine writePredictedData( self )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: self
        !
        integer :: i
        !
        open( ioPredData, file = "predicted_data.dat", action="write", position="append" )
        !
        !
        do i = 1, size( self%predicted_data_entries )
            !
            write( ioPredData, "(es12.6, A20, f15.3, f15.3, f15.3, f15.3, f15.3, A20, es16.6, es16.6, es16.6)" ) self%predicted_data_entries(i)%period, self%predicted_data_entries(i)%code, R_ZERO, R_ZERO, self%predicted_data_entries(i)%xyz(1), self%predicted_data_entries(i)%xyz(2), self%predicted_data_entries(i)%xyz(3), self%predicted_data_entries(i)%component, self%predicted_data_entries(i)%real, self%predicted_data_entries(i)%imaginary, 1.0
            !
        enddo
        !
        close( ioPredData )
        !
    end subroutine writePredictedData
    !
    subroutine updatePredictedDataArray( self, new_data )
          implicit none
          !
          class( Receiver_t ), intent( inout )        :: self
          type( PredictedDataHandle_t ), intent( in ) :: new_data
          !
          !
          type( PredictedDataHandle_t ), allocatable, dimension(:) :: temp_array
          integer :: istat
          !
          if( .NOT. allocated( self%predicted_data_entries )  ) then
                allocate( self%predicted_data_entries(1) )
                self%predicted_data_entries(1) = new_data
          else
                !
                allocate( temp_array( size( self%predicted_data_entries ) + 1 ), STAT = istat )
                temp_array( 1 : size( self%predicted_data_entries ) ) = self%predicted_data_entries
                temp_array( size( self%predicted_data_entries ) + 1 ) = new_data
                self%predicted_data_entries = temp_array
                !
          endif
          !
     end subroutine updatePredictedDataArray
     !
end module Receiver

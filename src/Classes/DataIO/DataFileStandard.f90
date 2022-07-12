!*************
!
! Class to read a data file and create an array with all data entries(lines)
!
!*************
!
module DataFileStandard
    !
    use DataFile
    use DataEntryMT
    use DataEntryMT_REF
    use DataEntryCSEM
    !
    type, extends( DataFile_t ) :: DataFileStandard_t
        !
        ! No derived properties
        !
    contains
        !
        final :: DataFileStandard_dtor
        !
    end type DataFileStandard_t
    !
    interface DataFileStandard_t
        module procedure DataFileStandard_ctor
    end interface DataFileStandard_t
    !
contains
    !
    ! Read line by line of the data file, create Data Entry objects (MT, MT_REF or CSEM)
    function DataFileStandard_ctor( funit, fname ) result( self )
        implicit none
        !
        integer, intent( in )                   :: funit
        character(:), allocatable, intent( in ) :: fname
        !
        type( DataFileStandard_t ) :: self
        !
        character(1000)                   :: full_line_text
        character(len=200), dimension(20) :: args
        !
        character(:), allocatable         :: line_text, actual_type, code, code_ref, component, dipole
        integer                           :: iDe, io_stat, p_nargs, nRx
        integer                           :: header_counter, header_line_counter, mt_counter, csem_counter
        real( kind=prec )                 :: period, rvalue, imaginary, error
        real( kind=prec )                 :: xyz_ref(3), latitude_ref, longitude_ref
        real( kind=prec )                 :: latitude, longitude, xyz(3), tx_xyz(3), moment, azimuth, dip
        !
        !write( *, * ) "Constructor DataFileStandard_t"
        !
        call self%init()
        !
        self%fine_name = fname
        !
        call Compact( fname )
        !
        open( unit = funit, file = fname, iostat = io_stat, status = "old" )
        !
        if( io_stat == 0 ) then
            !
            header_counter = 0
            header_line_counter = 0
            mt_counter = 0
            csem_counter = 0
            !
            do
                read( funit, "(a)", END = 10 ) full_line_text
                line_text = adjustl( full_line_text )
                line_text = trim( line_text )
                !
                call Parse( line_text, " ", args, p_nargs )
                !
                if( index( line_text, "#" ) == 0 .and. index( line_text, ">" ) == 0 ) then
                     !
                     iDe = size( self%data_entries ) + 1
                     !
                     selectcase( actual_type )
                          !
                          ! MT file line
                          case( "Full_Impedance", "Off_Diagonal_Impedance", "Full_Vertical_Components", &
                          "Off_Diagonal_Rho_Phase", "Phase_Tensor" )
                                !
                                !# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
                                !
                                read( args(1), "(f16.6)" )  period
                                code = trim( args(2) )
                                read( args(3), "(f16.6)" )  latitude
                                read( args(4), "(f16.6)" )  longitude
                                read( args(5), "(f16.6)" )  xyz(1)
                                read( args(6), "(f16.6)" )  xyz(2)
                                read( args(7), "(f16.6)" )  xyz(3)
                                component = trim( args(8) )
                                read( args(9), "(f16.6)" )  rvalue
                                read( args(10), "(f16.6)" ) imaginary
                                read( args(11), "(f16.6)" ) error
                                !
                                call self%loadReceiversAndTransmitters( DataEntryMT_t( iDe, actual_type, period, code, &
                                latitude, longitude, xyz, component, rvalue, imaginary, error ) )
                                !
                                mt_counter = mt_counter + 1
                                !
                          ! MT REF file line
                          case( "Full_Interstation_TF" )
                                !
                                !# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Code_REF GG_Lat_REF GG_Lon_REF X(m)_REF Y(m)_REF Z(m)_REF Component Real Imag Error
                                !
                                read( args(1), "(f16.6)" )  period
                                code = trim( args(2) )
                                read( args(3), "(f16.6)" )  latitude
                                read( args(4), "(f16.6)" )  longitude
                                read( args(5), "(f16.6)" )  xyz(1)
                                read( args(6), "(f16.6)" )  xyz(2)
                                read( args(7), "(f16.6)" )  xyz(3)
                                code_ref = trim( args(8) )
                                read( args(9), "(f16.6)" )  latitude_ref
                                read( args(10), "(f16.6)" ) longitude_ref
                                read( args(11), "(f16.6)" ) xyz_ref(1)
                                read( args(12), "(f16.6)" ) xyz_ref(2)
                                read( args(13), "(f16.6)" ) xyz_ref(3)
                                component = trim( args(14) )
                                read( args(15), "(f16.6)" ) rvalue
                                read( args(16), "(f16.6)" ) imaginary
                                read( args(17), "(f16.6)" ) error
                                !
                                call self%loadReceiversAndTransmitters( DataEntryMT_REF_t( iDe, actual_type,    &
                                period, code, latitude, longitude, xyz, code_ref,    &
                                latitude_ref, longitude_ref, xyz_ref, component, rvalue, imaginary, error ) )
                                !
                                mt_counter = mt_counter + 1
                                !
                          ! CSEM file line
                          case( "Ex_Field", "Ey_Field", "Bx_Field", "By_Field", "Bz_Field" )
                                !
                                !# Dipole Period(s) Moment(Am) Azi Dip Tx_X(m) Tx_Y(x) Tx_Z(m) Code X(m) Y(x) Z(m) Component Real Imag, Error
                                !
                                dipole = args(1)
                                read( args(2), "(f16.6)" ) period
                                read( args(3), "(f16.6)" ) moment
                                read( args(4), "(f16.6)" ) azimuth
                                read( args(5), "(f16.6)" ) dip
                                read( args(6), "(f16.6)" ) tx_xyz(1)
                                read( args(7), "(f16.6)" ) tx_xyz(2)
                                read( args(8), "(f16.6)" ) tx_xyz(3)
                                code = trim( args(9) )
                                read( args(10), "(f16.6)" ) xyz(1)
                                read( args(11), "(f16.6)" ) xyz(2)
                                read( args(12), "(f16.6)" ) xyz(3)
                                component = trim( args(13) )
                                read( args(14), "(f16.6)" ) rvalue
                                read( args(15), "(f16.6)" ) imaginary
                                read( args(16), "(f16.6)" ) error
                                !
                                call self%loadReceiversAndTransmitters( DataEntryCSEM_t( iDe, actual_type,    &
                                dipole, period, moment, azimuth, dip, tx_xyz,    &
                                code, xyz, component, rvalue, imaginary, error ) )
                                !
                                csem_counter = csem_counter + 1
                                !
                          case default
                                !
                                write( *, * ) "Unknown type :[", actual_type, "]"
                                stop "DataFileStandard.f08: DataFileStandard_ctor()"
                                !
                     end select
                     !
                     header_line_counter = 0
                     !
                else
                     !# Synthetic 3D MT data written in Matlab
                     !# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
                     !> Full_Impedance
                     !> exp(-i\omega t)
                     !> [V/m]/[T]
                     !> 0.00
                     !> 0.000 0.000
                     !> 4 2
                     header_line_counter = header_line_counter + 1
                     selectcase( header_line_counter )
                         !
                         ! Main Header, Data Fields
                         case( 1, 2 )
                         !
                         ! Data Type
                         case( 3 )
                             actual_type = args(2)
                         !
                         ! exp(-i\omega t) ????, [V/m]/[T] ????, 0.00 ????, 0.000 0.000 ????
                         case( 4, 5, 6, 7 )
                         !
                         ! nTx, nRx
                         case( 8 )
                             read( args(2), "(I8)" ) self%nTx
                             read( args(3), "(I8)" ) nRx
                             !
                             header_counter = header_counter + 1
                             write( *, "(A17, I8, A5, A30, A3, I8, A9, I8, A5)" ) "Header", header_counter, " -> [", trim(actual_type), "]: ", self%nTx, " Txs and ", nRx, " Rxs."
                             !
                             self%nRx = self%nRx + nRx
                             !
                         case default
                             !
                             write( *, * ) "Unknown header format in line :[", header_line_counter, "]"
                             stop "DataFileStandard.f08: DataFileStandard_ctor()"
                             !
                     end select
                     !
                end if
                !
            end do
            !
10          close( unit = funit )
            !
            if( mt_counter > 0 )   write( *, * ) "          Read ", mt_counter,   " MT Entries"
            if( csem_counter > 0 ) write( *, * ) "          Read ", csem_counter, " CSEM Entries"
            !
        else
            write( *, * ) "Error opening [", fname, "] in DataFileStandard_ctor"
            stop
        end if
        !
    end function DataFileStandard_ctor
    !
    subroutine DataFileStandard_dtor( self )
        implicit none
        !
        type( DataFileStandard_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor DataFileStandard_t"
        !
        call self%dealloc()
        !
    end subroutine DataFileStandard_dtor
    !
end module DataFileStandard

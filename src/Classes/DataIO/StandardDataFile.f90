!*************
!
! Class to read a data file and create an array with all data entries(lines)
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
!
module StandardDataFile
   !
   use DataFile
   use DataEntryMT
   use DataEntryMT_REF
   use DataEntryCSEM
   !
   type, extends( DataFile_t ) :: StandardDataFile_t
      !
      integer :: line_counter = 0
      !
   contains
      !
      final :: StandardDataFile_dtor
      !
   end type StandardDataFile_t
   !
   interface StandardDataFile_t
      module procedure StandardDataFile_ctor
   end interface StandardDataFile_t
   !
   public :: getLineNumber
   !
contains
   !
   ! Read line by line of the data file, create Data Entry objects (MT, MT_REF or CSEM)
   function StandardDataFile_ctor( funit, fname ) result( self )
      implicit none
      !
      integer, intent( in )                   :: funit
      character(:), allocatable, intent( in ) :: fname
	  !
	  class( StandardDataFile_t ), pointer :: self
      !
      character(1000)                   :: full_line_text
      character(len=200), dimension(20) :: args
      !1111
      class( DataEntry_t ), pointer     :: data_entry
      character(:), allocatable         :: line_text, actual_type, code, code_ref, component, dipole
      integer                           :: iDe, io_stat, p_nargs, header_counter, mt_counter, csem_counter
      real( kind=prec )                 :: period, real, imaginary, error
      real( kind=prec )                 :: xyz_ref(3), latitude_ref, longitude_ref
      real( kind=prec )                 :: latitude, longitude, xyz(3), tx_xyz(3), moment, azimuth, dip
      !
      ! write(*,*) "Constructor StandardDataFile_t"
      !
      allocate( StandardDataFile_t :: self )
      !
      call self%init()
      !
      call Compact( fname )
      !
      open( unit = funit, file = fname, iostat = io_stat, status = 'old' )
      !
      if( io_stat /= 0 ) then
         write(*,*) 'Unable to open [', fname, '], Stat: ', io_stat
      else
         self%line_counter = getLineNumber( funit )
         !
         rewind( funit )
         !
         header_counter = 0
         mt_counter = 0
         csem_counter = 0
         !
         do
            read( funit, '(a)', END = 10 ) full_line_text
            line_text = adjustl( full_line_text )
            line_text = trim( line_text )
            !
            call Parse( line_text, " ", args, p_nargs )
            !
            if( index( line_text, "#" ) == 0 .and. index( line_text, ">" ) == 0 ) then
                !
                iDe = self%data_entries%size() + 1
                !
                selectcase( actual_type )
                    !
                    ! MT file line
                    case( "Full_Impedance", "Off_Diagonal_Impedance", "Full_Vertical_Components", &
                    "Off_Diagonal_Rho_Phase", "Phase_Tensor" )
                        !
                        !# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
                        !
                        read( args(1), '(f15.5)' )   period
                        code = trim( args(2) )
                        read( args(3), '(f15.5)' )   latitude
                        read( args(4), '(f15.5)' )   longitude
                        read( args(5), '(f15.5)' )   xyz(1)
                        read( args(6), '(f15.5)' )   xyz(2)
                        read( args(7), '(f15.5)' )   xyz(3)
                        component = trim( args(8) )
                        read( args(9), '(f15.5)' )   real
                        read( args(10), '(f15.5)' )   imaginary
                        read( args(11), '(f15.5)' )   error
                        !
                        data_entry => DataEntryMT_t( iDe, actual_type, period, code, &
                        latitude, longitude, xyz, component, real, imaginary, error )
                        !
                        !if( .NOT. self%data_entries%has( data_entry ) ) then 
                            call self%data_entries%add( data_entry )
                        !end if
                        !
                        mt_counter = mt_counter + 1
                        !
                    ! MT REF file line
                    case( "Full_Interstation_TF" )
                        !
                        !# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Code_REF GG_Lat_REF GG_Lon_REF X(m)_REF Y(m)_REF Z(m)_REF Component Real Imag Error
                        !
                        read( args(1), '(f15.5)' )   period
                        code = trim( args(2) )
                        read( args(3), '(f15.5)' )   latitude
                        read( args(4), '(f15.5)' )   longitude
                        read( args(5), '(f15.5)' )   xyz(1)
                        read( args(6), '(f15.5)' )   xyz(2)
                        read( args(7), '(f15.5)' )   xyz(3)
                        code_ref = trim( args(8) )
                        read( args(9), '(f15.5)' )   latitude_ref
                        read( args(10), '(f15.5)' )   longitude_ref
                        read( args(11), '(f15.5)' )   xyz_ref(1)
                        read( args(12), '(f15.5)' )   xyz_ref(2)
                        read( args(13), '(f15.5)' )   xyz_ref(3)
                        component = trim( args(14) )
                        read( args(15), '(f15.5)' )   real
                        read( args(16), '(f15.5)' )   imaginary
                        read( args(17), '(f15.5)' )   error
                        !
                        data_entry => DataEntryMT_REF_t( iDe, actual_type,   &
                        period, code, latitude, longitude, xyz, code_ref,   &
                        latitude_ref, longitude_ref, xyz_ref, component, real, imaginary, error )
                        !
                        !if( .NOT. self%data_entries%has( data_entry ) ) then 
                            call self%data_entries%add( data_entry )
                        !end if
                        !
                        mt_counter = mt_counter + 1
                        !
                    ! CSEM file line
                    case( "Ex_Field", "Ey_Field", "Bx_Field", "By_Field", "Bz_Field" )
                        !
                        !# Dipole Period(s) Moment(Am) Azi Dip Tx_X(m) Tx_Y(x) Tx_Z(m) Code X(m) Y(x) Z(m) Component Real Imag, Error
                        !
                        dipole = args(1)
                        read( args(2), '(f15.5)' )  period
                        read( args(3), '(f15.5)' )  moment
                        read( args(4), '(f15.5)' )  azimuth
                        read( args(5), '(f15.5)' )  dip
                        read( args(6), '(f15.5)' )  tx_xyz(1)
                        read( args(7), '(f15.5)' )  tx_xyz(2)
                        read( args(8), '(f15.5)' )  tx_xyz(3)
                        code = trim( args(9) )
                        read( args(10), '(f15.5)' )   xyz(1)
                        read( args(11), '(f15.5)' )   xyz(2)
                        read( args(12), '(f15.5)' )   xyz(3)
                        component = trim( args(13) )
                        read( args(14), '(f15.5)' )   real
                        read( args(15), '(f15.5)' )   imaginary
                        read( args(16), '(f15.5)' )   error
                        !
                        data_entry => DataEntryCSEM_t( iDe, actual_type,   &
                        dipole, period, moment, azimuth, dip, tx_xyz,   &
                        code, xyz, component, real, imaginary, error )
                        !
                        !if( .NOT. self%data_entries%has( data_entry ) ) then 
                            call self%data_entries%add( data_entry )
                        !end if
                        !
                        csem_counter = csem_counter + 1
                        !
                    case default
                        !
                        write(*,*) "unknow type :[", actual_type, "]"
                        STOP "StandardDataFile.f08: StandardDataFile_ctor()"
                        !
                end select
                !
                header_counter = 0
                !
            else
                header_counter = header_counter + 1
                if( header_counter == 3 ) then
                    actual_type = args(2)
                end if
            end if
         end do
         !
10        close( unit = funit )
         !
         write(*,*) 'Finish read file [', fname, ']:'
         if( mt_counter > 0 )   write(*,*) mt_counter, ' MT Entries'
         if( csem_counter > 0 )   write(*,*) csem_counter, ' CSEM Entries'
         !
      end if
      !
   end function StandardDataFile_ctor
   !
   subroutine StandardDataFile_dtor( self )
      implicit none
      !
      type( StandardDataFile_t ), intent( in out ) :: self
      !
      ! write(*,*) "Destructor StandardDataFile_t"
      !
      deallocate( self%data_entries )
      !
   end subroutine StandardDataFile_dtor
   !
   ! Return the number of lines of a given file
   function getLineNumber( funit ) result( line_counter )
      implicit none
      !
      integer, intent( in )      :: funit
      !
      integer               :: line_counter
      character(80)         :: line_text
      !
      line_counter = 0
      !
      rewind( funit )
      !
      do
         read( funit, "(a)", END = 10 ) line_text
         line_text = adjustl( line_text )
         if( index( line_text, "#" ).eq.0 ) then
            !
            line_counter = line_counter + 1
         end if
      end do
      !
10      return
      !
   end function getLineNumber
   !
end module StandardDataFile

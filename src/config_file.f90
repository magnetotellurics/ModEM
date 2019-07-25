module config_file
    !
    use UserCtrl
	use file_units
    !
    private
    !
    public :: configs
    !
    !######################################
    ! BASE CLASS configs
    !######################################
    !
    type configs
        !
        private
        integer                     :: nlines
        character(50), pointer      :: names(:), values(:)
        ! Character-based information specified by the user
        type (userdef_control)      :: cUserDef
    contains
        procedure, public           :: init => initConfigs
        procedure, public           :: print => printConfigs
        procedure, public           :: validate => validateConfigs
		procedure, public           :: config => configUserDef
        procedure, public           :: usage => printUsage
        procedure, public           :: getCtrl => getcUserDef
        procedure, public           :: getNameId, getValidJob
        !

    end type configs
    !
contains
    !
    ! INIT READ CONFIGURATION FILE
    subroutine initConfigs( this, cUserDef )
        implicit none
        !
        class( configs ), intent( inout )   		:: this
        type ( userdef_control ), intent( in )    	:: cUserDef
        character(50), pointer              		:: aux_names(:), aux_values(:)
        character(50)                       		:: new_name, new_value
        integer                             		:: io_stat
        !
        this%cUserDef = cUserDef
        !
        open( unit=ioConfig, file=cUserDef%rFile_Config, iostat=io_stat, status='old' )
        !
        if ( io_stat /= 0 ) then
            write(*,*) 'STAT1: ', io_stat
        end if
        !
        ! READ FIRST LINE
        read( ioConfig, *, iostat=io_stat ) new_name, new_value
        !
        allocate( this%names(1) )
        this%names(1) = trim(new_name)
        !
        allocate( this%values(1) )
        this%values(1) = trim(new_value)
        !
        ! READ THE REST OF THE LINES
        this%nlines = 1
        do while( io_stat == 0 )
            !
            read( ioConfig, *, iostat=io_stat ) new_name, new_value
            !
            if( io_stat == 0 ) then
                allocate( aux_names( this%nlines + 1 ) )
                aux_names( 1 : this%nlines ) = this%names
                aux_names( this%nlines + 1 ) = trim(new_name)
                deallocate( this%names )
                allocate( this%names( this%nlines + 1 ) )
                this%names = aux_names
                deallocate( aux_names )
                !
                allocate( aux_values( this%nlines + 1 ) )
                aux_values( 1 : this%nlines ) = this%values
                aux_values( this%nlines + 1 ) = trim(new_value)
                deallocate( this%values )
                allocate( this%values( this%nlines + 1 ) )
                this%values = aux_values
                deallocate( aux_values )
                !
                this%nlines = this%nlines + 1
            end if
        end do
		!
		close( unit=ioConfig )
		!
		call this%validate()
		!
    end subroutine initConfigs
    !
    function getcUserDef( this ) result( cUserDef )
        class( configs ), intent( inout )   :: this
        type (userdef_control)              :: cUserDef
        !
        cUserDef = this%cUserDef
    end function getcUserDef
    !
    function getValidJob( this ) result( job )
        class( configs ), intent( inout )   :: this
        character(50)                       :: str_job
        character(1)                        :: job
        !
        ifield_name = this%getNameId( 'job' )
        str_job = this%values( ifield_name )
        !
        if( index( trim(str_job), 'CONF_FILE' ) /= 0 ) then
            job = 'W'
        end if
        !
        if( index( trim(str_job), 'READ_WRITE' ) /= 0 ) then
            job = 'R'
        end if
        !
        if( index( trim(str_job), 'FORWARD' ) /= 0 ) then
            job = 'F'
        end if
        !
        if( index( trim(str_job), 'COMPUTE_J' ) /= 0 ) then
            job = 'J'
        end if
        !
        if( index( trim(str_job), 'MULT_BY_J' ) /= 0 ) then
            job = 'M'
        end if
        !
        if( index( trim(str_job), 'MULT_BY_J_T' ) /= 0 ) then
            job = 'T'
        end if
        !
        if( index( trim(str_job), 'MULT_BY_J_T_multi_Tx' ) /= 0 ) then
            job = 'x'
        end if
        !
        if( index( trim(str_job), 'INVERSE' ) /= 0 ) then
            job = 'I'
        end if
        !
        if( index( trim(str_job), 'APPLY_COV' ) /= 0 ) then
            job = 'C'
        end if
        !
        if( index( trim(str_job), 'EXTRACT_BC' ) /= 0 ) then
            job = 'b'
        end if
        !
        if( index( trim(str_job), 'TEST_GRAD' ) /= 0 ) then
            job = 'g'
        end if
        !
        if( index( trim(str_job), 'TEST_ADJ' ) /= 0 ) then
            job = 'A'
        end if
        !
        if( index( trim(str_job), 'TEST_SENS' ) /= 0 ) then
            job = 'S'
        end if
        !
    end function getValidJob
    !
    function getNameId( this, field_name ) result( name_id )
        class( configs ), intent( inout )   :: this
        character(*), intent( in )          :: field_name
        integer                             :: name_id
        !
        ! LOOKS FOR field_name TAG
        name_id = 0
        do i = 1, this%nlines
            if( index( this%names(i), trim(field_name) ) /= 0 ) then
                name_id = i
            end if
        end do
        !
    end function getNameId
    !
    subroutine validateConfigs( this )
        class( configs ), intent( inout ) :: this
        integer                           :: ifield_name
        !
        !TEST JOB
        ifield_name = this%getNameId( 'job' )
        this%cUserDef%job = this%values( ifield_name )
        !
        if( ifield_name == 0 ) then
            write(*,*) 'There is no tag job!'
            stop
        else
            !
            select case ( this%cUserDef%job )

                case ( 'READ_WRITE' ) !R
					write(6,*) 'CONFIG: READ_WRITE'
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                    end if
                    !TEST rFile_Model OPTIONAL
                    ifield_name = this%getNameId( 'wFile_Model' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%wFile_Model = this%values( ifield_name )
                    end if
                    !TEST rFile_Data OPTIONAL
                    ifield_name = this%getNameId( 'wFile_Data' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%wFile_Data = this%values( ifield_name )
                    end if

                case ( 'FORWARD' ) !F
					write(6,*) 'CONFIG: FORWARD'
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                    end if
                    !TEST wFile_Data
                    ifield_name = this%getNameId( 'wFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag wFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%wFile_Data = this%values( ifield_name )
                    end if
                    !TEST wFile_EMsoln OPTIONAL
                    ifield_name = this%getNameId( 'wFile_EMsoln' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%wFile_EMsoln = this%values( ifield_name )
                    end if
                    !TEST rFile_fwdCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_fwdCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_fwdCtrl = this%values( ifield_name )
                    end if
                    !TEST rFile_EMrhs OPTIONAL
                    ifield_name = this%getNameId( 'rFile_EMrhs' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_EMrhs = this%values( ifield_name )
                    end if

                case ( 'COMPUTE_J' ) ! J
					write(6,*) 'CONFIG: COMPUTE_J'
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                    end if
                    !TEST wFile_Sens
                    ifield_name = this%getNameId( 'wFile_Sens' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag wFile_Sens!'
                        call this%usage()
                    else
                        this%cUserDef%wFile_Sens = this%values( ifield_name )
                    end if
                    !TEST rFile_fwdCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_fwdCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_fwdCtrl = this%values( ifield_name )
                    end if

                case ( 'MULT_BY_J' ) ! M
					write(6,*) 'CONFIG: MULT_BY_J'
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                    end if
                    !TEST rFile_dModel
                    ifield_name = this%getNameId( 'rFile_dModel' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_dModel!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_dModel = this%values( ifield_name )
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                    end if
                    !TEST wFile_Data
                    ifield_name = this%getNameId( 'wFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag wFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%wFile_Data = this%values( ifield_name )
                    end if
                    !TEST rFile_fwdCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_fwdCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_fwdCtrl = this%values( ifield_name )
                    end if

                case ( 'MULT_BY_J_T', 'MULT_BY_J_T_multi_Tx' ) ! J
					write(6,*) 'CONFIG: MULT_BY_J_T, MULT_BY_J_T_multi_Tx'
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                    end if
                    !TEST wFile_Sens
                    ifield_name = this%getNameId( 'wFile_Sens' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag wFile_Sens!'
                        call this%usage()
                    else
                        this%cUserDef%wFile_Sens = this%values( ifield_name )
                    end if
                    !TEST rFile_fwdCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_fwdCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_fwdCtrl = this%values( ifield_name )
                    end if

                case ( 'INVERSE' ) ! I
                    write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data [lambda eps]'

                case ( 'APPLY_COV' ) ! C
                    write(0,*) 'Usage: -C  [FWD|INV] rFile_Model wFile_Model [rFile_Cov rFile_Prior]'

                case ( 'EXTRACT_BC' ) ! b
                    write(0,*) 'Usage: -b  rFile_Model rFile_Data wFile_EMrhs [rFile_fwdCtrl]'


                case ( 'TEST_GRAD' ) !g
                    write(0,*) 'Usage: -g  rFile_Model rFile_Data rFile_dModel [rFile_fwdCtrl rFile_EMrhs]'
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                    end if
                    !TEST rFile_dModel
                    ifield_name = this%getNameId( 'rFile_dModel' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_dModel!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_dModel = this%values( ifield_name )
                    end if
                    !TEST rFile_fwdCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_fwdCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_fwdCtrl = this%values( ifield_name )
                    end if
                    !TEST rFile_EMrhs OPTIONAL
                    ifield_name = this%getNameId( 'rFile_EMrhs' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_EMrhs = this%values( ifield_name )
                    end if

                case ( 'TEST_ADJ' ) ! A
                    write(0,*) 'Usage: Test the adjoint implementation for each of the critical'

                case ( 'TEST_SENS' ) ! S
                    write(0,*) 'Usage: -S  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl wFile_Sens]'

                case default
                    write(0,*) 'Unknown job. Please check your command line options'
                    stop

            end select
        !
        this%cUserDef%job = this%getValidjob()
        !
        end if
    end subroutine validateConfigs
    !
	subroutine configUserDef( this, cUserDef )
        class( configs ), intent( inout ) :: this
		type ( userdef_control ), intent( in )    	:: cUserDef
		!
		 ! special case W - read from configuration file
		 if( cUserDef%job .eq. CONF_FILE ) then !W
			!
			! IF CONFIG FILE EXIST...
			inquire( file=cUserDef%rFile_Config, exist=config_file_exists )
			!
			if ( config_file_exists ) then
				!
				write(6,*) '#### INIT CONFIGS'
				call this%init( cUserDef )
				!
				write(6,*) '#### PRINT'
				call this%print()
				!
			else
				write(6,*) 'File [', cUserDef%rFile_Config, '] does not exists!'
				stop
			end if
		end if
	end subroutine configUserDef
    !
    subroutine printUsage( this )
        class( configs ), intent( in ) :: this
        !
        select case ( this%cUserDef%job )

            case ('READ_WRITE') !R
                write(0,*) 'Usage: -R  rFile_Model rFile_Data [wFile_Model wFile_Data]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             READ_WRITE'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_Model     <name> (OPTIONAL)'
                write(0,*) 'wFile_Data      <name> (OPTIONAL)'

            case ('FORWARD') !F
                write(0,*) 'Usage: -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl rFile_EMrhs]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             FORWARD'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_Data      <name>'
                write(0,*) 'wFile_EMsoln    <name> (OPTIONAL)'
                write(0,*) 'rFile_fwdCtrl   <path> (OPTIONAL)'
                write(0,*) 'rFile_EMrhs     <path> (OPTIONAL)'

            case ('COMPUTE_J') ! J
                write(0,*) 'Usage: -J  rFile_Model rFile_Data wFile_Sens [rFile_fwd!this%cUserDef]'

            case (MULT_BY_J_T) ! T
                write(0,*) 'Usage: -T  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'

            case (MULT_BY_J_T_multi_Tx) ! x
                write(0,*) 'Usage: -x  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'

            case (INVERSE) ! I
                write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data [lambda eps]'

            case (APPLY_COV) ! C
                write(0,*) 'Usage: -C  [FWD|INV] rFile_Model wFile_Model [rFile_Cov rFile_Prior]'

            case (EXTRACT_BC) ! b
                write(0,*) 'Usage: -b  rFile_Model rFile_Data wFile_EMrhs [rFile_fwdCtrl]'

            case (TEST_GRAD) !g
                write(0,*) 'Usage: -g  rFile_Model rFile_Data rFile_dModel [rFile_fwdCtrl rFile_EMrhs]'

            case (TEST_ADJ) ! A
                write(0,*) 'Usage: Test the adjoint implementation for each of the critical'

            case (TEST_SENS) ! S
                write(0,*) 'Usage: -S  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl wFile_Sens]'

            case default
                write(0,*) 'Unknown job. Please check your command line options'
        end select
        !
        stop
            !
    end subroutine printUsage
    !
    !
    subroutine printConfigs( this )
        class( configs ), intent( in ) :: this
        !
        write(*,*) 'The configuration under [', this%cUserDef%rFile_Config, '] has ', this%nlines, ' lines :'
        do i = 1, this%nlines
            write(*,*) this%names(i), ' = ', this%values(i)
        end do
        !
        write(*,*) 'USER CONTROL JOB: ', this%cUserDef%job
        !
    end subroutine printConfigs
    !
end module config_file

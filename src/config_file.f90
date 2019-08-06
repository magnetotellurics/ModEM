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
        procedure, public           :: init     => initConfigs
        procedure, public           :: print    => printConfigs
        procedure, public           :: validate => validateConfigs
        procedure, public           :: config   => configUserDef
        procedure, public           :: usage    => printUsage
        procedure, public           :: getCtrl  => getcUserDef
        procedure, public           :: getNameId
        procedure, public           :: getValidJob
        !
    end type configs
    !
contains
    !
    ! INIT READ CONFIGURATION FILE
    subroutine initConfigs( this, cUserDef )
        implicit none
        !
        class( configs ), intent( inout )           :: this
        type ( userdef_control ), intent( in )      :: cUserDef
        character(50), pointer                      :: aux_names(:), aux_values(:)
        character(50)                               :: new_name, new_value
        integer                                     :: io_stat
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
        !
        class( configs ), intent( inout )   :: this
        integer                             :: ifield_name
        logical                             :: exists
        !
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
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
                    end if
                    !TEST wFile_Model OPTIONAL
                    ifield_name = this%getNameId( 'wFile_Model' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%wFile_Model = this%values( ifield_name )
                    end if
                    !TEST wFile_Data OPTIONAL
                    ifield_name = this%getNameId( 'wFile_Data' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%wFile_Data = this%values( ifield_name )
                    end if

                case ( 'FORWARD' ) !F
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
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
                        inquire(FILE=this%cUserDef%rFile_fwdCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid rFile_fwdCtrl'
                            stop
                        end if
                    end if
                    !TEST rFile_EMrhs OPTIONAL
                    ifield_name = this%getNameId( 'rFile_EMrhs' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_EMrhs = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_EMrhs,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid rFile_EMrhs'
                            stop
                        end if
                    end if

                case ( 'COMPUTE_J' ) ! J
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
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
                        inquire(FILE=this%cUserDef%rFile_fwdCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid rFile_fwdCtrl'
                            stop
                        end if
                    end if

                case ( 'MULT_BY_J' ) ! M
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST rFile_dModel
                    ifield_name = this%getNameId( 'rFile_dModel' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_dModel!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_dModel = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_dModel,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_dModel'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
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
                        inquire(FILE=this%cUserDef%rFile_fwdCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid rFile_fwdCtrl'
                            stop
                        end if
                    end if

                case ( 'MULT_BY_J_T', 'MULT_BY_J_T_multi_Tx' ) ! J
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
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
                        inquire(FILE=this%cUserDef%rFile_fwdCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid rFile_fwdCtrl'
                            stop
                        end if
                    end if

                case ( 'INVERSE' ) ! I
                    !TEST search
                    ifield_name = this%getNameId( 'search' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag search!'
                        call this%usage()
                    else
                        this%cUserDef%search = this%values( ifield_name )
                        select case ( this%cUserDef%search )
                            case ( 'NLCG','DCG','Hybrid','LBFGS' )
                                write(*,*) 'Inverse search ',trim(this%cUserDef%search),' selected.'
                            case default
                                write(*,*) 'Unknown inverse search. Usage: -I [NLCG | DCG | Hybrid | LBFGS]'
                                stop
                        end select
                    end if
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
                    end if
                    !TEST lambda OPTIONAL
                    ifield_name = this%getNameId( 'lambda' )
                    if( ifield_name /= 0 ) then
                        read( this%values( ifield_name ), '(f10.8)' )  this%cUserDef%lambda
                    end if
                    !TEST eps OPTIONAL
                    ifield_name = this%getNameId( 'eps' )
                    if( ifield_name /= 0 ) then
                        read( this%values( ifield_name ), '(f10.8)' )  this%cUserDef%eps
                    end if
                    !TEST rFile_invCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_invCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_invCtrl = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_invCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid inverse control file or damping parameter'
                            stop
                        end if
                    end if
                    !TEST rFile_fwdCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_fwdCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_fwdCtrl = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_fwdCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid fwd control file or damping parameter'
                            stop
                        end if
                    end if


                case ( 'APPLY_COV' ) ! C
                    !TEST option
                    ifield_name = this%getNameId( 'option' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag option!'
                        call this%usage()
                    else
                        this%cUserDef%option = this%values( ifield_name )
                        select case ( this%cUserDef%option )
                            case ( 'FWD','INV' )
                                write(*,*) 'Covariance option ',trim(this%cUserDef%option),' selected.'
                            case default
                                write(*,*) 'Unknown Covariance option. Usage: -C  [FWD|INV]'
                                stop
                        end select
                    end if
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST wFile_Model
                    ifield_name = this%getNameId( 'wFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag wFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%wFile_Model = this%values( ifield_name )
                    end if
                    !TEST rFile_Cov OPTIONAL
                    ifield_name = this%getNameId( 'rFile_Cov' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_Cov = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Cov,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Cov'
                            stop
                        end if
                    end if
                    !TEST rFile_Prior OPTIONAL
                    ifield_name = this%getNameId( 'rFile_Prior' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_Prior = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Prior,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Prior'
                            stop
                        end if
                    end if

                case ( 'EXTRACT_BC' ) ! b
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
                    end if
                    !TEST wFile_EMrhs
                    ifield_name = this%getNameId( 'wFile_EMrhs' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag wFile_EMrhs!'
                        call this%usage()
                    else
                        this%cUserDef%wFile_EMrhs = this%values( ifield_name )
                    end if
                    !TEST rFile_fwdCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_fwdCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_fwdCtrl = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_fwdCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid fwd control file or damping parameter'
                            stop
                        end if
                    end if

                case ( 'TEST_GRAD' ) !g
                    !TEST rFile_Model
                    ifield_name = this%getNameId( 'rFile_Model' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Model!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Model = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Model'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
                    end if
                    !TEST rFile_dModel
                    ifield_name = this%getNameId( 'rFile_dModel' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_dModel!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_dModel = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_dModel,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_dModel'
                            stop
                        end if
                    end if
                    !TEST rFile_fwdCtrl OPTIONAL
                    ifield_name = this%getNameId( 'rFile_fwdCtrl' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_fwdCtrl = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_fwdCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid rFile_fwdCtrl'
                            stop
                        end if
                    end if
                    !TEST rFile_EMrhs OPTIONAL
                    ifield_name = this%getNameId( 'rFile_EMrhs' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%rFile_EMrhs = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_EMrhs,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid rFile_EMrhs'
                            stop
                        end if
                    end if

                case ( 'TEST_ADJ' ) ! A
                    !TEST option
                    ifield_name = this%getNameId( 'option' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag option!'
                        call this%usage()
                    else
                        this%cUserDef%option = this%values( ifield_name )
                        !
                        select case ( this%cUserDef%option )
                            case ( 'J' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST rFile_dModel
                                ifield_name = this%getNameId( 'rFile_dModel' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_dModel!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_dModel = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_dModel,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_dModel'
                                        stop
                                    end if
                                end if
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if
                                !TEST wFile_Model OPTIONAL
                                ifield_name = this%getNameId( 'wFile_Model' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_Model = this%values( ifield_name )
                                end if
                                !TEST wFile_Data OPTIONAL
                                ifield_name = this%getNameId( 'wFile_Data' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_Data = this%values( ifield_name )
                                end if

                            case ( 'L' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST rFile_EMsoln
                                ifield_name = this%getNameId( 'rFile_EMsoln' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_EMsoln!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_EMsoln = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_EMsoln,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_EMsoln'
                                        stop
                                    end if
                                end if
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if
                                !TEST wFile_EMrhs OPTIONAL
                                ifield_name = this%getNameId( 'wFile_EMrhs' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_EMrhs = this%values( ifield_name )
                                end if
                                !TEST wFile_Data OPTIONAL
                                ifield_name = this%getNameId( 'wFile_Data' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_Data = this%values( ifield_name )
                                end if

                            case ( 'S' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST rFile_EMrhs
                                ifield_name = this%getNameId( 'rFile_EMrhs' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_EMrhs!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_EMrhs = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_EMrhs,EXIST=exists)
                                    if ( .not. exists ) then
                                        ! problem - invalid argument
                                        write(0,*) 'Please specify a valid rFile_EMrhs'
                                        stop
                                    end if
                                end if
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if
                                !TEST wFile_EMsoln OPTIONAL
                                ifield_name = this%getNameId( 'wFile_EMsoln' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_EMsoln = this%values( ifield_name )
                                end if

                            case ( 'P' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST rFile_dModel
                                ifield_name = this%getNameId( 'rFile_dModel' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_dModel!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_dModel = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_dModel,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_dModel'
                                        stop
                                    end if
                                end if
                                !TEST rFile_EMsoln
                                ifield_name = this%getNameId( 'rFile_EMsoln' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_EMsoln!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_EMsoln = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_EMsoln,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_EMsoln'
                                        stop
                                    end if
                                end if
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if
                                !TEST wFile_Model OPTIONAL
                                ifield_name = this%getNameId( 'wFile_Model' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_Model = this%values( ifield_name )
                                end if
                                !TEST wFile_EMrhs OPTIONAL
                                ifield_name = this%getNameId( 'wFile_EMrhs' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_EMrhs = this%values( ifield_name )
                                end if

                            case ( 'Q' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST rFile_dModel
                                ifield_name = this%getNameId( 'rFile_dModel' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_dModel!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_dModel = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_dModel,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_dModel'
                                        stop
                                    end if
                                end if
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if
                                !TEST wFile_Model OPTIONAL
                                ifield_name = this%getNameId( 'wFile_Model' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_Model = this%values( ifield_name )
                                end if
                                !TEST wFile_Data OPTIONAL
                                ifield_name = this%getNameId( 'wFile_Data' )
                                if( ifield_name /= 0 ) then
                                    this%cUserDef%wFile_Data = this%values( ifield_name )
                                end if

                            case ( 'O' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if

                            case ( 'm' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST wFile_Model
                                ifield_name = this%getNameId( 'wFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag wFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%wFile_Model = this%values( ifield_name )
                                end if
                                !TEST delta OPTIONAL
                                ifield_name = this%getNameId( 'delta' )
                                if( ifield_name /= 0 ) then
                                    read( this%values( ifield_name ), '(f10.8)' )  this%cUserDef%delta
                                end if

                            case ( 'd' )
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if
                                !TEST wFile_Data
                                ifield_name = this%getNameId( 'wFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag wFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%wFile_Data = this%values( ifield_name )
                                end if
                                !TEST delta OPTIONAL
                                ifield_name = this%getNameId( 'delta' )
                                if( ifield_name /= 0 ) then
                                    read( this%values( ifield_name ), '(f10.8)' )  this%cUserDef%delta
                                end if

                            case ( 'e' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if
                                !TEST rFile_EMsoln
                                ifield_name = this%getNameId( 'rFile_EMsoln' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_EMsoln!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_EMsoln = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_EMsoln,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_EMsoln'
                                        stop
                                    end if
                                end if
                                !TEST wFile_EMsoln
                                ifield_name = this%getNameId( 'wFile_EMsoln' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag wFile_EMsoln!'
                                    call this%usage()
                                else
                                    this%cUserDef%wFile_EMsoln = this%values( ifield_name )
                                end if
                                !TEST delta OPTIONAL
                                ifield_name = this%getNameId( 'delta' )
                                if( ifield_name /= 0 ) then
                                    read( this%values( ifield_name ), '(f10.8)' )  this%cUserDef%delta
                                end if

                            case ( 'b' )
                                !TEST rFile_Model
                                ifield_name = this%getNameId( 'rFile_Model' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Model!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Model = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Model,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Model'
                                        stop
                                    end if
                                end if
                                !TEST rFile_Data
                                ifield_name = this%getNameId( 'rFile_Data' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_Data!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_Data = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                                    if ( .not. exists ) then
                                        write(0,*) 'Please specify a valid rFile_Data'
                                        stop
                                    end if
                                end if
                                !TEST rFile_EMrhs
                                ifield_name = this%getNameId( 'rFile_EMrhs' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag rFile_EMrhs!'
                                    call this%usage()
                                else
                                    this%cUserDef%rFile_EMrhs = this%values( ifield_name )
                                    inquire(FILE=this%cUserDef%rFile_EMrhs,EXIST=exists)
                                    if ( .not. exists ) then
                                        ! problem - invalid argument
                                        write(0,*) 'Please specify a valid rFile_EMrhs'
                                        stop
                                    end if
                                end if
                                !TEST wFile_EMrhs
                                ifield_name = this%getNameId( 'wFile_EMrhs' )
                                if( ifield_name == 0 ) then
                                    write(*,*) 'There is no tag wFile_EMrhs!'
                                    call this%usage()
                                else
                                    this%cUserDef%wFile_EMrhs = this%values( ifield_name )
                                end if
                                !TEST delta OPTIONAL
                                ifield_name = this%getNameId( 'delta' )
                                if( ifield_name /= 0 ) then
                                    read( this%values( ifield_name ), '(f10.8)' )  this%cUserDef%delta
                                end if

                            case default
                                write(0,*) 'Unknown option. Please check your configuration file'
                                write(0,*) 'Valid options: [J, L, S, P, Q, O, m, d, e, b]'
                                stop

                        end select
                    end if

                case ( 'TEST_SENS' ) ! S
                    !TEST rFile_dModel
                    ifield_name = this%getNameId( 'rFile_dModel' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_dModel!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_dModel = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_dModel,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_dModel'
                            stop
                        end if
                    end if
                    !TEST rFile_Data
                    ifield_name = this%getNameId( 'rFile_Data' )
                    if( ifield_name == 0 ) then
                        write(*,*) 'There is no tag rFile_Data!'
                        call this%usage()
                    else
                        this%cUserDef%rFile_Data = this%values( ifield_name )
                        inquire(FILE=this%cUserDef%rFile_Data,EXIST=exists)
                        if ( .not. exists ) then
                            write(0,*) 'Please specify a valid rFile_Data'
                            stop
                        end if
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
                        inquire(FILE=this%cUserDef%rFile_fwdCtrl,EXIST=exists)
                        if ( .not. exists ) then
                            ! problem - invalid argument
                            write(0,*) 'Please specify a valid rFile_fwdCtrl'
                            stop
                        end if
                    end if
                    !TEST wFile_Sens OPTIONAL
                    ifield_name = this%getNameId( 'wFile_Sens' )
                    if( ifield_name /= 0 ) then
                        this%cUserDef%wFile_Sens = this%values( ifield_name )
                    end if

                case default
                    write(0,*) 'Unknown job. Please check your configuration file'
                    write(0,*) ' '
                    write(0,*) 'Valis job names:'
                    write(0,*) ' '
                    write(0,*) 'READ_WRITE, FORWARD, COMPUTE_J, MULT_BY_J_T'
                    write(0,*) 'MULT_BY_J_T_multi_Tx, INVERSE, APPLY_COV'
                    write(0,*) 'EXTRACT_BC, TEST_GRAD, TEST_ADJ and TEST_SENS'
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
        type ( userdef_control ), intent( in )      :: cUserDef
        !
        ! special case W - read from configuration file
        if( cUserDef%job .eq. CONF_FILE ) then !W
            !
            ! IF CONFIG FILE EXIST...
            inquire( file=cUserDef%rFile_Config, exist=config_file_exists )
            !
            if ( config_file_exists ) then
                !
                call this%init( cUserDef )
                !
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
                write(0,*) 'job         READ_WRITE'
                write(0,*) 'rFile_Model <path>'
                write(0,*) 'rFile_Data  <path>'
                write(0,*) 'wFile_Model <name> [OPTIONAL]'
                write(0,*) 'wFile_Data  <name> [OPTIONAL]'

            case ('FORWARD') !F
                write(0,*) 'Usage: -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl rFile_EMrhs]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             FORWARD'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_Data      <name>'
                write(0,*) 'wFile_EMsoln    <name> [OPTIONAL]'
                write(0,*) 'rFile_fwdCtrl   <path> [OPTIONAL]'
                write(0,*) 'rFile_EMrhs     <path> [OPTIONAL]'

            case ('COMPUTE_J') ! J
                write(0,*) 'Usage: -J  rFile_Model rFile_Data wFile_Sens [rFile_fwdCtrl]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             COMPUTE_J'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_Sens      <name>'
                write(0,*) 'rFile_fwdCtrl   <path> [OPTIONAL]'

            case ( 'MULT_BY_J' ) ! M
                write(0,*) 'Usage: -M  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             MULT_BY_J'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_dModel    <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_Data      <name>'
                write(0,*) 'rFile_fwdCtrl   <path> [OPTIONAL]'

            case ('MULT_BY_J_T') ! T
                write(0,*) 'Usage: -T  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             MULT_BY_J_T'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_dModel    <name>'
                write(0,*) 'rFile_fwdCtrl   <path> [OPTIONAL]'

            case ('MULT_BY_J_T_multi_Tx') ! x
                write(0,*) 'Usage: -x  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             MULT_BY_J_T_multi_Tx'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_dModel    <name>'
                write(0,*) 'rFile_fwdCtrl   <path> [OPTIONAL]'

            case ('INVERSE') ! I
                write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data [lambda eps]'
                write(0,*) 'or'
                write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data [rFile_invCtrl rFile_fwdCtrl]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             INVERSE'
                write(0,*) 'search          <NLCG, DCG, Hybrid, LBFGS>'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'lambda          <value> [OPTIONAL]'
                write(0,*) 'eps             <value> [OPTIONAL]'
                write(0,*) 'rFile_invCtrl   <path>  [OPTIONAL]'
                write(0,*) 'rFile_fwdCtrl   <path>  [OPTIONAL]'

            case ('APPLY_COV') ! C
                write(0,*) 'Usage: -C [FWD|INV] rFile_Model wFile_Model [rFile_Cov rFile_Prior]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             APPLY_COV'
                write(0,*) 'option          FWD or INV'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'wFile_Model     <name>'
                write(0,*) 'rFile_Cov       <path> [OPTIONAL]'
                write(0,*) 'rFile_Prior     <path> [OPTIONAL]'

            case ('EXTRACT_BC') ! b
                write(0,*) 'Usage: -b  rFile_Model rFile_Data wFile_EMrhs [rFile_fwdCtrl]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             EXTRACT_BC'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_EMrhs     <name>'
                write(0,*) 'rFile_fwdCtrl   <path> [OPTIONAL]'

            case ('TEST_GRAD') !g
                write(0,*) 'Usage: -g  rFile_Model rFile_Data rFile_dModel [rFile_fwdCtrl rFile_EMrhs]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             TEST_GRAD'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'rFile_dModel    <path>'
                write(0,*) 'rFile_fwdCtrl   <path> [OPTIONAL]'
                write(0,*) 'rFile_EMrhs     <path> [OPTIONAL]'

            case ('TEST_ADJ') ! A

                select case ( this%cUserDef%option )
                    case ( 'J' )
                        write(0,*) 'Usage: -A  J rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          J'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'rFile_dModel    <path>'
                        write(0,*) 'rFile_Data      <path>'
                        write(0,*) 'wFile_Model     <name> [OPTIONAL]'
                        write(0,*) 'wFile_Data      <name> [OPTIONAL]'

                    case ( 'L' )
                        write(0,*) 'Usage: -A  L rFile_Model rFile_EMsoln rFile_Data [wFile_EMrhs wFile_Data]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          L'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'rFile_EMsoln    <path>'
                        write(0,*) 'rFile_Data      <path>'
                        write(0,*) 'wFile_EMrhs     <name> [OPTIONAL]'
                        write(0,*) 'wFile_Data      <name> [OPTIONAL]'

                    case ( 'S' )
                        write(0,*) 'Usage: -A  S rFile_Model rFile_EMrhs rFile_Data [wFile_EMsoln]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          S'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'rFile_EMrhs     <path>'
                        write(0,*) 'rFile_Data      <path>'
                        write(0,*) 'wFile_EMsoln    <name> [OPTIONAL]'

                    case ( 'P' )
                        write(0,*) 'Usage: -A  P rFile_Model rFile_dModel rFile_EMsoln rFile_Data [wFile_Model wFile_EMrhs]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          P'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'rFile_dModel    <path>'
                        write(0,*) 'rFile_EMsoln    <path>'
                        write(0,*) 'rFile_Data      <path>'
                        write(0,*) 'wFile_Model     <name> [OPTIONAL]'
                        write(0,*) 'wFile_EMrhs     <name> [OPTIONAL]'

                    case ( 'Q' )
                        write(0,*) 'Usage: -A  Q rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          Q'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'rFile_dModel    <path>'
                        write(0,*) 'rFile_Data      <path>'
                        write(0,*) 'wFile_Model     <name> [OPTIONAL]'
                        write(0,*) 'wFile_Data      <name> [OPTIONAL]'

                    case ( 'O' )
                        write(0,*) 'Usage: -A  O rFile_Model rFile_Data'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          O'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'rFile_Data      <path>'

                    case ( 'm' )
                        write(0,*) 'Usage: -A  m rFile_Model wFile_Model [delta]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          m'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'wFile_Model     <name>'
                        write(0,*) 'delta           <value> [OPTIONAL]'

                    case ( 'd' )
                        write(0,*) 'Usage: -A  d rFile_Data wFile_Data [delta]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          d'
                        write(0,*) 'rFile_Data      <path>'
                        write(0,*) 'wFile_Data      <name>'
                        write(0,*) 'delta           <value> [OPTIONAL]'

                    case ( 'e' )
                        write(0,*) 'Usage: -A  e rFile_Model rFile_Data rFile_EMsoln wFile_EMsoln [delta]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          e'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'rFile_Data      <path>'
                        write(0,*) 'rFile_EMsoln    <path>'
                        write(0,*) 'wFile_EMsoln    <name>'
                        write(0,*) 'delta           <value> [OPTIONAL]'

                    case ( 'b' )
                        write(0,*) 'Usage: -A  b rFile_Model rFile_Data rFile_EMrhs wFile_EMrhs [delta]'
                        write(0,*) ' '
                        write(0,*) 'Your configuration file must have these fields:'
                        write(0,*) ' '
                        write(0,*) 'job             TEST_ADJ'
                        write(0,*) 'option          b'
                        write(0,*) 'rFile_Model     <path>'
                        write(0,*) 'rFile_Data      <path>'
                        write(0,*) 'rFile_EMrhs     <path>'
                        write(0,*) 'wFile_EMrhs     <name>'
                        write(0,*) 'delta           <value> [OPTIONAL]'

                    case default
                        write(0,*) 'Unknown option. Please check your configuration file'
                        write(0,*) 'Valid options: [J, L, S, P, Q, O, m, d, e, b]'
                        stop

                end select

            case ('TEST_SENS') ! S
                write(0,*) 'Usage: -S  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl wFile_Sens]'
                write(0,*) ' '
                write(0,*) 'Your configuration file must have these fields:'
                write(0,*) ' '
                write(0,*) 'job             TEST_SENS'
                write(0,*) 'rFile_Model     <path>'
                write(0,*) 'rFile_dModel    <path>'
                write(0,*) 'rFile_Data      <path>'
                write(0,*) 'wFile_Data      <name>'
                write(0,*) 'rFile_fwdCtrl   <path> [OPTIONAL]'
                write(0,*) 'wFile_Sens      <path> [OPTIONAL]'

            case default
                write(0,*) 'Unknown job. Please check your configuration file'
                write(0,*) ' '
                write(0,*) 'Valis job names:'
                write(0,*) ' '
                write(0,*) 'READ_WRITE, FORWARD, COMPUTE_J, MULT_BY_J_T'
                write(0,*) 'MULT_BY_J_T_multi_Tx, INVERSE, APPLY_COV'
                write(0,*) 'EXTRACT_BC, TEST_GRAD, TEST_ADJ and TEST_SENS'
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

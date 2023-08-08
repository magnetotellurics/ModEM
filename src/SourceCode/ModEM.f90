!
!> Prototype of a full ModEM serial|parallel program
!
program ModEM
    !
#ifdef MPI
    !
    use WorkerMPI
    !
    call constructorMPI
    !
    ! MPI Master process
    if( mpi_rank == 0 ) then
        !
        call runProgram
        !
        call MPI_Finalize( ierr )
        !
    ! MPI Worker process
    else
        !
        call workerMainLoop
        !
    endif
    !
#else
    !
    use InversionJob
    !
    call runProgram
    !
#endif
    !
contains
    !
    !> Main program protocol
    !
    subroutine runProgram()
        implicit none
        !
        character( len=20 ) :: str_msg
        real( kind=prec ) :: t_start, t_finish
        integer :: int_time
        !
        call date_and_time( str_date, str_time )
        !
        !> Start runtime countdown
        call cpu_time( t_start )
        !
        modem_job = "unknown"
        !
        !> Set default value of program variables
        call setupDefaultParameters
        !
        !> Handle arguments passed by the user on the command line
        call handleArguments
        !
        write( *, * )
        write( *, "(a18, a8, a1, a6, a1)" ) "Start ModEM-OO at ", str_date, "_", str_time, "."
        write( *, * )
        !
        !> If it was passed by argument,
        !> Check parameters at the forward control file
        if( has_fwd_control_file ) call handleForwardControlFile
        !
        !> If it was passed by argument,
        !> Check parameters at the inversion control file
        if( has_inv_control_file ) call handleInversionControlFile
        !
        !> Execute the job specified in the arguments
        call handleJob
        !
        !> Deallocate remaining main program memory
        call garbageCollector
        !
        !> End runtime countdown
        call cpu_time( t_finish )
        !
        int_time = int( t_finish - t_start )
        !
        write( *, * )
        write( *, * ) "Finish ModEM-OO: ", getLiteralTime( int_time ), aux_counter
        !
        if( warning_counter .GT. 0 ) then
            !
            write( str_msg, "(I2)") warning_counter
            !
            call warning( str_msg )
            !
        endif
        !
        write( *, * )
        !
    end subroutine runProgram
    !
    !> No subroutine briefing
    !
    subroutine handleJob()
        implicit none
        !
        !> Except for the case of Inversion,
        !> Free the memory used by the global control file, which is no longer useful
        if( index( modem_job, "Inversion" ) .LE. 0 .AND. allocated( inv_control_file ) ) then
            deallocate( inv_control_file )
        endif
        !
        select case( modem_job )
            !
            case( "JobForwardModeling" )
                !
                call jobForwardModeling
                !
            case( "JobJMult" )
                !
                call jobJMult
                !
            case( "JobJMult_T" )
                !
                call jobJMult_T
                !
            case( "JobInversion" )
                !
                call jobInversion
                !
            case( "JobSplitModel" )
                !
                call jobSplitModel
                !
            case default
                !
                call warning( "Unknown job: ["//modem_job//"]!" )
                call printHelp
                stop
            !
        end select
        !
    end subroutine handleJob
    !
end program ModEM
!
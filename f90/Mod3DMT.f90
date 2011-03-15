! *****************************************************************************
program Mod3DMT
! Program for running 3D MT forward, sensitivity and inverse modelling
! Copyright (c) 2004-2010 Oregon State University
!              AUTHORS  Gary Egbert, Anna Kelbert & Naser Meqbel
!              College of Atmospheric and Oceanic Sciences

     use SensComp
     use SymmetryTest
     use Main
     use NLCG
     use DCG
     !use mtinvsetup

#ifdef MPI
     Use MPI_main
#endif

     implicit none

     ! Character-based information specified by the user
     type (userdef_control)	:: cUserDef

     ! Variable required for storing the date and time
     type (timer_t)         :: timer



#ifdef MPI
              call  MPI_constructor
			  if (taskid==0) then
			      call parseArgs('Mod3DMT',cUserDef) ! OR readStartup(rFile_Startup,cUserDef)
			      write(6,*)'I am a PARALLEL version'
			      call Master_job_Distribute_userdef_control(cUserDef)
	              open(ioMPI,file=cUserDef%wFile_MPI)
	              write(ioMPI,*) 'Total Number of nodes= ', number_of_workers
			  else
			       call RECV_cUserDef(cUserDef)
			 end if

#else
             call parseArgs('Mod3DMT',cUserDef) ! OR readStartup(rFile_Startup,cUserDef)
			 write(6,*)'I am a SERIAL version'
#endif




      call initGlobalData(cUserDef)
      ! set the grid for the numerical computations

#ifdef MPI
      call setGrid_MPI(grid)
#else
      call setGrid(grid)
#endif


#ifdef MPI
    if (taskid.gt.0) then
			    call Worker_job(sigma0,allData)
	            if (trim(worker_job_task%what_to_do) .eq. 'Job Completed')  then
	               	 call deallGlobalData()
		             call cleanUp_MPI()
	                 call MPI_destructor
	              stop
	            end if
    end if
#endif



	 ! Start the (portable) clock
	 call reset_time(timer)

     select case (cUserDef%job)

     case (READ_WRITE)
        if (write_model .and. write_data) then
        	write(*,*) 'Writing model and data files and exiting...'
        	call write_modelParam(sigma0,cUserDef%wFile_Model)
        	call write_dataVectorMTX(allData,cUserDef%wFile_Data)
		else if (write_model) then
        	write(*,*) 'Writing model and exiting...'
        	call write_modelParam(sigma0,cUserDef%wFile_Model)
		end if

     case (FORWARD)

     write(*,*) 'Calculating predicted data...'

#ifdef MPI
        call Master_job_fwdPred(sigma0,allData,eAll)
        call Master_job_STOP_MESSAGE
#else
        call fwdPred(sigma0,allData,eAll)
#endif

        ! write out all impedances
        call write_dataVectorMTX(allData,cUserDef%wFile_Data)

        if (write_EMsoln) then
        	! write out EM solutions
        	write(*,*) 'Saving the EM solution...'
        	call write_solnVectorMTX(fidWrite,cUserDef%wFile_EMsoln,eAll)
        end if

     case (COMPUTE_J)
        write(*,*) 'Calculating the full sensitivity matrix...'
#ifdef MPI
        call Master_job_fwdPred(sigma0,allData,eAll)
        call Master_job_calcJ(allData,sigma0,sens,eAll)
        call Master_job_STOP_MESSAGE
#else
        call calcJ(allData,sigma0,sens)
#endif
        call write_sensMatrixMTX(sens,cUserDef%wFile_Sens)

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'

#ifdef MPI
            !call Master_job_Jmult(dsigma,sigma0,allData)
            call Master_job_STOP_MESSAGE
#else
            call Jmult(dsigma,sigma0,allData)
#endif


        call write_dataVectorMTX(allData,cUserDef%wFile_Data)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
#ifdef MPI
         call Master_job_fwdPred(sigma0,allData,eAll)
         call Master_job_JmultT(sigma0,allData,dsigma,eAll)
         call Master_job_STOP_MESSAGE
#else
         call JmultT(sigma0,allData,dsigma)
#endif

         call write_modelParam(dsigma,cUserDef%wFile_dModel)

     case (INVERSE)
     	if (trim(cUserDef%search) == 'NLCG') then
            ! sigma1 contains mHat on input (zero = starting from the prior)
        	write(*,*) 'Starting the NLCG search...'
            sigma1 = dsigma
           	call NLCGsolver(allData,cUserDef%lambda,sigma0,sigma1,cUserDef%rFile_invCtrl)

#ifdef MPI
        	call Master_job_STOP_MESSAGE
#endif
    	elseif (trim(cUserDef%search) == 'DCG') then
        	write(*,*) 'Starting the DCG search...'
        	sigma1 = dsigma
        	call DCGsolver(allData,sigma0,sigma1,cUserDef%lambda)
            !call Marquardt_M_space(allData,sigma0,sigma1,cUserDef%lambda)
#ifdef MPI
        	call Master_job_STOP_MESSAGE
#endif

       else
        	write(*,*) 'Inverse search ',trim(cUserDef%search),' not yet implemented. Exiting...'
        	stop
        end if
        if (write_model) then
            call write_modelParam(sigma1,cUserDef%wFile_Model)
        end if
        if (write_data) then
        	call write_dataVectorMTX(allData,cUserDef%wFile_Data)
        end if

     case (APPLY_COV)
        select case (cUserDef%option)
            case('FWD')
                write(*,*) 'Multiplying input model parameter by square root of the covariance ...'
                sigma1 = multBy_CmSqrt(dsigma)
                call linComb(ONE,sigma1,ONE,sigma0,sigma1)
            case('INV')
                write(*,*) 'Multiplying input model parameter by inverse square root of the covariance ...'
                call linComb(ONE,dsigma,MinusONE,sigma0,dsigma)
                sigma1 = multBy_CmSqrtInv(dsigma)
            case default
               write(0,*) 'Unknown covariance option ',trim(cUserDef%option),'; please use FWD or INV.'
               stop
        end select
        call write_modelParam(sigma1,cUserDef%wFile_Model)

     case (TEST_ADJ)
       select case (cUserDef%option)
           case('J')
               call Jtest(sigma0,dsigma,allData)
           case('Q')
               call Qtest(sigma0,dsigma,allData)

           case default
               write(0,*) 'Symmetry test for operator ',trim(cUserDef%option),' not yet implemented.'
       end select
         if (write_model .and. write_data) then
            write(*,*) 'Writing model and data files and exiting...'
            call write_modelParam(dsigma,cUserDef%wFile_Model)
            call write_dataVectorMTX(allData,cUserDef%wFile_Data)
         else if (write_model) then
            write(*,*) 'Writing model and exiting...'
            call write_modelParam(dsigma,cUserDef%wFile_Model)
         end if

     case default

        write(0,*) 'No job ',trim(cUserDef%job),' defined.'

     end select
	 ! cleaning up
	 call deallGlobalData()

#ifdef MPI
            close(ioMPI)
	    call cleanUp_MPI()
#else
            call cleanUp()
#endif





#ifdef MPI
	 write(0,*) ' elapsed time (mins) ',elapsed_time(timer)/60.0
	 call MPI_destructor
#else
	 write(0,*) ' elapsed time (mins) ',elapsed_time(timer)/60.0
#endif



end program


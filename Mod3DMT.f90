! *****************************************************************************
program Mod3DMT
! Program for running 3D MT forward, sensitivity and inverse modelling
! Copyright (c) 2004-2014 Oregon State University
!              AUTHORS  Gary Egbert, Anna Kelbert & Naser Meqbel
!              College of Earth, Ocean and Atmospheric Sciences

     use SensComp
     use SymmetryTest
     use Main
     use NLCG
     use DCG
     !use mtinvsetup

#ifdef MPI
     Use Main_MPI
#endif

     implicit none

     ! Character-based information specified by the user
     type (userdef_control) :: cUserDef

     ! Variable required for storing the date and time
     type (timer_t)         :: timer

     ! Output variable
     character(80)          :: header
     integer                :: ios

#ifdef MPI
      call  constructor_MPI
      if (taskid==0) then
          call parseArgs('Mod3DMT',cUserDef)  
          ! OR readStartup(rFile_Startup,cUserDef)
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
    ! Check if a large grid file with E field is defined:
    ! NOTE: right now both grids share the same transmitters.
    ! This why, reading and setting the large grid and its E solution comes after setting the trasnmitters Dictionary.
    if (NESTED_BC) then
      if (taskid==0) then
            call read_Efield_from_file(solverParams)
            call Interpolate_BC(grid)
            call Master_job_Distribute_nTx_nPol(nTx_nPol)
      else
            call RECV_nTx_nPol
            call ini_BC_from_file(grid)
            call RECV_BC_form_Master
      end if    
    end if    
#else
      call setGrid(grid)
#endif

#ifdef MPI
    if (BC_FROM_RHS_FILE) then
        if (taskid==0) then
            call read_rhsVectorMTX(grid,bAll,cUserDef%rFile_EMrhs)
        else
            ! need to logic to fetch the BCs from the master node
        end if
    end if
#else
    if (BC_FROM_RHS_FILE) then
        call read_rhsVectorMTX(grid,bAll,cUserDef%rFile_EMrhs)
    end if
#endif

#ifdef MPI
      ! for debug
      ! write(6,*)'Reporting from Node #', taskid
      if (taskid.gt.0) then
          call Worker_job(sigma0,allData)
          if (trim(worker_job_task%what_to_do) .eq. 'Job Completed')  then
              call deallGlobalData()
              call cleanUp_MPI()
              call destructor_MPI
              stop
          end if
      end if
#endif
      ! Start the (portable) clock
      call reset_time(timer)

      select case (cUserDef%job)
      case (READ_WRITE)
        if (output_level > 3) then
            call print_txDict()
            call print_rxDict()
        end if
        if (write_model .and. write_data) then
            write(*,*) 'Writing model and data files and exiting...'
            call write_modelParam(sigma0,cUserDef%wFile_Model)
            call write_dataVectorMTX(allData,cUserDef%wFile_Data)
        else if (write_model) then
            write(*,*) 'Writing model and exiting...'
            call write_modelParam(sigma0,cUserDef%wFile_Model)
        end if
      case (FORWARD)
        write(6,*) 'Calculating predicted data...'
#ifdef MPI
        call Master_job_fwdPred(sigma0,allData,eAll)
#else
        call fwdPred(sigma0,allData,eAll)
#endif


        ! write out all impedances
        call write_dataVectorMTX(allData,cUserDef%wFile_Data)

     case (COMPUTE_J)
        write(*,*) 'Calculating the full sensitivity matrix...'
#ifdef MPI
        call Master_job_fwdPred(sigma0,allData,eAll)
        call Master_job_calcJ(allData,sigma0,sens,eAll)
#else
        call calcJ(allData,sigma0,sens)
#endif
        call write_sensMatrixMTX(sens,cUserDef%wFile_Sens)

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'

#ifdef MPI
        call Master_job_Jmult(dsigma,sigma0,allData)
#else
        call Jmult(dsigma,sigma0,allData)
#endif


        call write_dataVectorMTX(allData,cUserDef%wFile_Data)

     case (MULT_BY_J_T)
         write(*,*) 'Multiplying by J^T...'
#ifdef MPI
         !call Master_job_fwdPred(sigma0,allData,eAll)
         call Master_job_JmultT(sigma0,allData,dsigma)
#else
         call JmultT(sigma0,allData,dsigma)
#endif

         call write_modelParam(dsigma,cUserDef%wFile_dModel)
         
     case (MULT_BY_J_T_multi_Tx)
         write(*,*) 'Multiplying by J^T...output multi-Tx model vectors'
#ifdef MPI
         call Master_job_fwdPred(sigma0,allData,eAll)
         call Master_job_JmultT(sigma0,allData,dsigma,eAll,JT_multi_Tx_vec)
#else
         call fwdPred(sigma0,allData,eAll)
         call JmultT(sigma0,allData,dsigma,eAll,JT_multi_Tx_vec)
#endif
         open(unit=ioSens, file=cUserDef%wFile_dModel, form='unformatted', iostat=ios)
         write(0,*) 'Output JT_multi_Tx_vec...'
         write(header,*) 'JT multi_Tx vectors'
         write(ioSens) header
         call writeVec_modelParam_binary(size(JT_multi_Tx_vec),JT_multi_Tx_vec,header,cUserDef%wFile_dModel)
         close(ioSens)

     case (INVERSE)
         if (trim(cUserDef%search) == 'NLCG') then
            ! sigma1 contains mHat on input (zero = starting from the prior)
             write(*,*) 'Starting the NLCG search...'
             sigma1 = dsigma
             call NLCGsolver(allData,cUserDef%lambda,sigma0,sigma1,       &
     &            cUserDef%rFile_invCtrl)

         elseif (trim(cUserDef%search) == 'DCG') then
             write(*,*) 'Starting the DCG search...'
             sigma1 = dsigma
             call DCGsolver(allData,sigma0,sigma1,cUserDef%lambda,        &
     &            cUserDef%rFile_invCtrl)
            !call Marquardt_M_space(allData,sigma0,sigma1,cUserDef%lambda)

         else
           write(*,*) 'Inverse search ',trim(cUserDef%search),            &
     &          ' not yet implemented. Exiting...'
           stop
         end if
         if (write_model) then
             call write_modelParam(sigma1,cUserDef%wFile_Model)
         end if
         if (write_data) then
             call write_dataVectorMTX(allData,cUserDef%wFile_Data)
         end if
 
      case (APPLY_COV)
#ifdef MPI
       write(0,*) 'Covariance options can be run in serial only.'
#else
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
#endif

     case (EXTRACT_BC)
        ! no need to run the forward solution to extract the boundary
        ! conditions from the initial electric field
        call dryRun(sigma0,allData,bAll,eAll)

     case (TEST_GRAD)
        ! note that on input, dsigma is the non-smoothed model parameter
        sigma1 = sigma0
        call zero_modelParam(sigma1)
        alpha = 0.0d0

        write(0,*) '  Compute f(m0)'
        call func(cUserDef%lambda,allData,sigma0,sigma1,f1,mNorm,predData,eAll,rms)
        !call printf('f(m0)',lambda,alpha,f1,mNorm,rms)
        write(0,'(a10)',advance='no') 'f(m0)'//':'
        write(0,'(a3,es13.6,a4,es13.6,a5,f11.6)') ' f=',f1,' m2=',mNorm,' rms=',rms

        write(0,*) '  Compute df/dm|_{m0}'
        call gradient(cUserDef%lambda,allData,sigma0,sigma1,sigmaGrad,predData,eAll)
        write(0,'(a20,g15.7)') '| df/dm|_{m0} |: ',dotProd(sigmaGrad,sigmaGrad)

        write(0,*) '  Compute f(m0+dm)'
        call func(cUserDef%lambda,allData,sigma0,dsigma,f2,mNorm,predData,eAll,rms)
        !call printf('f(m0+dm)',lambda,alpha,f2,mNorm,rms)
        write(0,'(a10)',advance='no') 'f(m0+dm)'//':'
        write(0,'(a3,es13.6,a4,es13.6,a5,f11.6)') ' f=',f2,' m2=',mNorm,' rms=',rms

        ! sometimes, it makes more sense to use the gradient at the second point...
        ! which would be equivalent to taking a step from -dm to zero.
        ! Note that in the limit rms->0 such as when the artifically large errors
        ! are used, Taylor series breaks for the model norm - can't ignore the
        ! second derivative. Then analytically, this approximation yields
        ! f(dm) - f(0) = df/dm|_0 x dm = 0, or
        ! f(-dm) - f(0) = df/dm|_{-dm} x {-dm} = 2 \lambda dm^2 = 2 * f(-dm)
        ! In both of these cases, f(0) = 0 and f(dm) = f(-dm) = \lambda * dm^2.
        !write(0,*) '  Compute df/dm|_{m0+dm}'
        !call gradient(cUserDef%lambda,allData,sigma0,dsigma,sigmaGrad,predData,eAll)
        !write(0,'(a20,g15.7)') '| df/dm|_{m0+dm} |: ',dotProd(sigmaGrad,sigmaGrad)

        !sigma1 = multBy_CmSqrt(dsigma)
        write(0,'(a82)') 'Total derivative test based on Taylor series [f(m0+dm) ~ f(m0) + df/dm|_{m0} x dm]'
        write(0,'(a20,g15.7)') 'f(m0+dm) - f(m0): ',f2-f1
        write(0,'(a20,g15.7)') 'df/dm|_{m0} x dm: ',dotProd(dsigma,sigmaGrad)

     case (TEST_SENS)
        ! compute d = J m row by row using the full sensitivity matrix
        write(*,*) 'Calculating the full sensitivity matrix...'
#ifdef MPI
        call Master_job_fwdPred(sigma0,allData,eAll)
        call Master_job_calcJ(allData,sigma0,sens,eAll)
#else
        call calcJ(allData,sigma0,sens)
#endif

        call multBy_sensMatrixMTX(sens,dsigma,predData)
        allData = predData

        call write_dataVectorMTX(allData,cUserDef%wFile_Data)

        ! now, compute d = J m using Jmult
        write(*,*) 'Multiplying by J...'

#ifdef MPI
        call Master_job_Jmult(dsigma,sigma0,predData,eAll)
#else
        call Jmult(dsigma,sigma0,predData)
#endif

        write(0,'(a82)') 'Comparison between d = J m using full Jacobian, row by row (calcJ) and using Jmult'
        write(0,'(a20,g15.7)') '|d| using calcJ: ',dotProd(allData,allData)
        write(0,'(a20,g15.7)') '|d| using Jmult: ',dotProd(predData,predData)

     case (TEST_ADJ)
#ifdef MPI
       write(0,*) 'Adjoint symmetry tests are implemented to be run in serial only.'
#else
       select case (cUserDef%option)
           case('J')
               call Jtest(sigma0,dsigma,allData)
           case('P')
               call Ptest(sigma0,allData,dsigma,eAll)
           case('L')
               call Ltest(sigma0,eAll,allData)
           case('Q')
               call Qtest(sigma0,dsigma,allData)
           case('S')
               call Stest(sigma0,bAll,eAll)
           case('m')
               dsigma = sigma0
               call random_modelParam(dsigma,cUserDef%delta)
           case('d')
               call random_dataVectorMTX(allData,cUserDef%delta)
           case('e')
               call random_solnVectorMTX(eAll,cUserDef%delta)
           case('b')
               call random_rhsVectorMTX(bAll,cUserDef%delta)
           case default
               write(0,*) 'Adjoint symmetry test for operator ',trim(cUserDef%option),' not yet implemented.'
       end select

       ! file writing...
       if (write_model) then
           write(*,*) 'Writing the model file...'
           call write_modelParam(dsigma,cUserDef%wFile_Model)
       end if

       if (write_data) then
           write(*,*) 'Writing the data file...'
           call write_dataVectorMTX(allData,cUserDef%wFile_Data)
       end if
#endif

     case default

        write(0,*) 'No job ',trim(cUserDef%job),' defined.'

     end select

     if (write_EMsoln) then
         ! write out EM solutions
         write(*,*) 'Saving the EM solution...'
         call write_solnVectorMTX(eAll,cUserDef%wFile_EMsoln)
     end if

     if (write_EMrhs) then
         ! write out the RHS
         write(*,*) 'Saving the RHS...'
         call write_rhsVectorMTX(bAll,cUserDef%wFile_EMrhs,'sparse')
     end if

     write(*,*) 'Exiting...'


#ifdef MPI
      ! cleaning up already completed by worker jobs... no "cleanup" needed on master node
      if (taskid==0) then
         call deallGlobalData()
      end if
      call Master_job_STOP_MESSAGE
      close(ioMPI)
#else
      ! cleaning up
      call deallGlobalData()
      call cleanUp()
#endif


#ifdef MPI
      write(0,*) ' elapsed time (mins) ',elapsed_time(timer)/60.0
      call destructor_MPI
#else
      write(0,*) ' elapsed time (mins) ',elapsed_time(timer)/60.0
#endif

end program

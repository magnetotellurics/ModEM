! *****************************************************************************
module ModelSpace
    ! This module contains the type definitions for the structures that contain
    ! the full information about the model parametrization

  use model_operators
  use modeldef
  use modelmap
  implicit none

Contains

  ! ***************************************************************************
  ! * read_modelParam reads the parametrization info to store in the derived data
  ! * type variables (modelCoeff_t, modelLayer_t) in the special case when the
  ! * parametrization is in terms of spherical harmonics of various degree/order
  ! * per layer, each of the parameters being either variable (keyword 'range') or
  ! * constant (keyword 'const'). We keep this information internally in such
  ! * a format, that generalizations of this case would be easy to implement.
  ! * In the parametrization type variables, each layer has the same (maximum)
  ! * degree and order of spherical harmonics, and hence identical numbers of
  ! * coefficients to store. Out of these coefficients, those that have been
  ! * defined in the script bear logical keyword 'exists'. Out of those, variable
  ! * coefficients have logical 'frozen==.FALSE.'
  ! * Obviously, the information on the range is not required for the forward
  ! * solver to operate. This is provided for the inversion, which will share
  ! * the same input format for now.
  ! * The optional output Pimag allows to use the same file format for the complex
  ! * spherical harmonic source vector storage...

  subroutine read_modelParam(P,cfile,Pimag)

    character(*), intent(in)            :: cfile
    type (modelParam_t), intent(inout)  :: P
    type (modelParam_t), intent(inout), optional  :: Pimag ! used for sources only
    integer                             :: ilayer,i,j,k,n,l,m,w
    integer                             :: nF,nL
    integer                             :: sum,sum0,degree
    integer                             :: ios,istat
    real(8)                             :: upperb,lowerb,width,depth,alpha,beta,weight
    character(6)                        :: if_log_char,if_var_char
    logical                             :: if_log, if_tan, if_fixed, exists, isComplex
    character(200)                      :: prmname, string
    real(8)                             :: v,vimag,vmin,vmax
    real(8)                             :: period ! read in place of depth for sources

    lowerb = EARTH_R
    depth = 0.0d0
    width = 0.0d0
    sum = 0
    sum0 = 0
    period = 0.0d0

    if(present(Pimag)) then
      isComplex = .true.
    else
      isComplex = .false.
    end if

    inquire(FILE=trim(cfile),EXIST=exists)
    if(exists) then
      open(ioPrm,file=cfile,status='old',form='formatted',iostat=ios)
      write(6,*) node_info,'Reading from the parametrization file ',trim(cfile)
    else
      write(0,*) node_info,'Error: (read_modelParam) input file ',trim(cfile),' does not exist'
      stop
    end if

    read(ioPrm,'(a8,a80)') string,prmname

    if (index(prmname,'harmonic')==0) then
       write(0, *) node_info,'Error: (read_modelParam) not a spherical harmonic parametrization'
       stop
    else
      i = max(index(prmname,'layers'),index(prmname,'periods'))
      read(prmname(i+7:len(prmname)),*) nL
      i = index(prmname,'degree')
      read(prmname(i+7:len(prmname)),*) degree
      write(6,*) node_info,trim(prmname)
    end if


    call create_modelParam(P,nL,degree)
    if (isComplex) then
       call create_modelParam(Pimag,nL,degree)
    end if

    write(6,'(a60,i3)') 'Number of layers/periods in script: ',P%nL


    ! Unwind spherical harmonic parametrization

    do n=1,nL

       read(ioPrm, *)

       upperb = lowerb

       ! If you wish to read the line in sections, use advance='no' specifier
     !  read(ioPrm,'(a3,1x,a6,1x,i2,a6,g15.7)',iostat=ios) &
            !if_log_char,string,degree,string,depth

            !read(ioPrm,'(a3)',iostat=ios,advance='no') if_log_char

            read(ioPrm,'(a200)',iostat=ios) string
            i = index(string,'degree')
            j = index(string,'layer')
            k = index(string,'reg')
            w = index(string,'weight')

            if (w==0) then
                    ! no weighting specified for this layer
                    weight = 1.0d0
                    w = len(string)
            else
                    read(string(w+7:len(string)),*) weight
            end if

            if (k==0) then
                    ! no regularisation specified for this layer
                    alpha = 0.0d0
                    beta  = 1.0d0
                    k = len(string)
            else
                    read(string(k+4:w),*) alpha,beta
            end if

            read(string(1:i-1),*) if_log_char
            if (j>0) then
                read(string(i+7:j),*) degree
                read(string(j+6:k),*) depth
                !write(*,*) 'DEBUG 1: ', string(j+6:k)
            else ! try source structure
                j = index(string,'period')
                read(string(i+7:j),*) degree
                read(string(j+6:k),*) period
                !write(*,*) 'DEBUG 2: ', string(j+6:k)
            end if



       read(ioPrm,'(a80)',iostat=ios) string

       if (if_log_char == 'log') then  ! log means log_{10}
          if_log = .TRUE.
       else
          if_log = .FALSE.
       end if

       if (if_log_char == 'tan') then
          if_tan = .TRUE.
       else
          if_tan = .FALSE.
       end if

       lowerb = EARTH_R - depth


       if (period > 0.0d0) then
          call setLayer_modelParam(P,n,upperb,lowerb,alpha,beta,weight,if_log,if_tan,period)
          call setLayer_modelParam(Pimag,n,upperb,lowerb,alpha,beta,weight,if_log,if_tan,period)
       else
          call setLayer_modelParam(P,n,upperb,lowerb,alpha,beta,weight,if_log,if_tan)
       end if

       sum0=sum0+sum
       sum=0
       do l=0,degree
          sum = sum + (2*l+1)
       end do
       !write(6,'(a46,i2,a2,i3)') 'Number of coefficients in layer ',n,': ',sum

       do i=1,sum
          if (isComplex) then
            read(ioPrm,*,iostat=ios) l,m,v,vimag,vmin,vmax,if_var_char
            !write(*,'(i4,2e14.6)') i,v,vimag
          else
            read(ioPrm,*,iostat=ios) l,m,v,vmin,vmax,if_var_char
          end if
          if (if_var_char == 'range') then
            if_fixed = .FALSE.
          else if (if_var_char == 'const') then
            if_fixed = .TRUE.
          else
            write(0, *) 'Error: (read_modelParam) wrong character constant for ',n,l,m,': ',if_var_char
            write(0, *) 'Error: (read_modelParam) if file looks OK, try running dos2unix on it. Exiting...'
            stop
          end if
          call setCoeffValue_modelParam(P,n,l,m,v,vmin,vmax,if_fixed)
          if (isComplex) then
            call setCoeffValue_modelParam(Pimag,n,l,m,vimag,vmin,vmax,if_fixed)
          end if
          !print *,'Values: ',l,m,v,vmin,vmax,if_fixed
       end do

    end do

    ! this check is not necessary, but helps to debug parametrization scripts
    !write(6,'(a50,i3)') 'Number of variable parameters in script: ',count(.not.P%c%frozen)
    write(6,*)

    close(ioPrm)

    P%zeroValued = .FALSE.
    P%smoothed = .FALSE.

    if (isComplex) then
        Pimag%zeroValued = .FALSE.
        Pimag%smoothed = .FALSE.
    end if

    return

  end subroutine read_modelParam

  ! **********************************************************************
  ! * This is an essential subroutine that is used to output the model
  ! * solution in e.g. NLCG. This solution is then used as the prior/initial
  ! * model, or to compute the fields.
  ! * We output both the final "smooth" model m = Cm^{1/2} mhat + m0,
  ! * and the "rough" model mhat, which is the result of the inverse search.
  ! * The smooth model m is best used to save and to plot the model; while
  ! * the rough model mhat is best used in conjunction with the prior to
  ! * restart the inversion, or compute the responses / derivatives.
  ! * Therefore, to avoid confusion we DO NOT write the regularization
  ! * parameters for a model if the smoothing has already been applied.
  ! * Regularization parameters are only needed to smooth a rough model,
  ! * so we leave them out for m (and leave them in for mhat).
  ! * Does not at present support spherical harmonic source output.
  ! * Here, we choose to output model in spherical harmonics for the model
  ! * parameter and on the grid for the final (smoothed) model.
  subroutine write_modelParam(P,cfile)

    implicit none
    type (modelParam_t), intent(in)         :: P
    character(*), intent(in)                :: cfile

    type (modelParam_t)                     :: P_grid
    integer                                 :: lmax,i,j,k,istat
    character(6)                            :: if_log_char,if_var_char

    if(.not.P%allocated) then
       call warning('(write_modelParam) parametrization not allocated yet; writing to file skipped')
       return
    end if

    if(index(cfile,'.rho')>0) then
        call mapToGrid_modelParam(P_grid,P)
        call writeGrid_modelParam(P_grid,cfile)
        call deall_modelParam(P_grid)
    else if(index(cfile,'.prm')>0) then
        call writeLayers_modelParam(P,cfile)
    else
        call warning('(write_modelParam) unknown file format requested; writing to file skipped')
    end if

  end subroutine write_modelParam



  ! **********************************************************************
  ! * Use this to output resistivity model on the grid
  subroutine writeGrid_modelParam(P,cfile)

    type (modelParam_t), intent(in)         :: P
    character(*), intent(in)                :: cfile
    ! local
    real(8)                           :: lon,lat,depth
    integer                           :: ios,i,j,k

    if(.not.P%allocated) then
       call warning('(writeGrid_modelParam) parametrization not allocated yet')
       return
    else if(trim(P%type) .ne. 'grid') then
       call warning('(writeGrid_modelParam) skipped writing this model parametrization type')
       return
    else if(.not. associated(P%rho0)) then
       call errStop('(writeGrid_modelParam) background resistivity not associated')
    end if

    if (cfile /= '') then
      open(ioMdl,file=cfile,status='unknown',form='formatted',iostat=ios)
                write(ioMdl,*) '# lon(i), lat(j), depth(k), rho(ijk) of cell node'
                write(ioMdl,'(3i3)') P%grid%nx,P%grid%ny,P%grid%nzEarth
        do k=P%grid%nzAir+1,P%grid%nz
          do i=1,P%grid%nx
            do j=1,P%grid%ny
                lon = P%grid%ph(i)*r2d
                lat = 90.0d0-P%grid%th(j)*r2d
                depth = EARTH_R - P%grid%r(k)
              write(ioMdl,'(3g15.7)',advance='no') lon,lat,depth
              write(ioMdl,'(g15.7)') P%rho%v(i,j,k)
            end do
          end do
        end do
      close(ioMdl)
    else
      do k=P%grid%nzAir,P%grid%nz
        i=P%grid%nx
        j=1
        write(*,'(a10,i3,i3,i3)',advance='no') 'i,j,k = ',i,j,k
        write(*,*) ' rho(ijk) = ',P%rho%v(i,j,k)
      end do
      do k=P%grid%nzAir,P%grid%nz
        write(*,*) 'k,radius(k) = ',k,P%grid%r(k)
      end do
    end if

  end subroutine writeGrid_modelParam  ! outputModel

  ! **********************************************************************
  ! * This is an essential subroutine that is used to output the model
  ! * solution in e.g. NLCG. This solution is then used as the prior/initial
  ! * model, or to compute the fields.
  ! * We output both the final "smooth" model m = Cm^{1/2} mhat + m0,
  ! * and the "rough" model mhat, which is the result of the inverse search.
  ! * The smooth model m is best used to save and to plot the model; while
  ! * the rough model mhat is best used in conjunction with the prior to
  ! * restart the inversion, or compute the responses / derivatives.
  ! * Therefore, to avoid confusion we DO NOT write the regularization
  ! * parameters for a model if the smoothing has already been applied.
  ! * Regularization parameters are only needed to smooth a rough model,
  ! * so we leave them out for m (and leave them in for mhat).
  ! * Does not at present support spherical harmonic source output.
  subroutine writeLayers_modelParam(P,cfile)

    implicit none
    type (modelParam_t), intent(in)         :: P
    character(*), intent(in)                :: cfile

    integer                                 :: lmax,i,j,k,istat
    character(6)                            :: if_log_char,if_var_char

    if(.not.P%allocated) then
       call warning('(writeLayers_modelParam) parametrization not allocated yet')
       return
    else if(trim(P%type) .ne. 'harmonic') then
       call warning('(writeLayers_modelParam) skipped writing this model parametrization type')
       return
    end if

    open(unit=ioPrm, file=cfile, status='unknown', iostat=istat)

    lmax = getDegree(P)

    write(ioPrm,'(a24,i2,a8,i2)') 'Format: harmonic layers ',P%nL,' degree ',lmax
    write(ioPrm,*)

    do j=1,P%nL
        if (P%L(j)%if_log) then
            if_log_char = 'log'
        elseif (P%L(j)%if_tan) then
            if_log_char = 'tan'
        else
            if_log_char = 'linear'
        end if
        lmax = getDegree(P,j)
        write(ioPrm,'(a6)',advance='no') if_log_char
        write(ioPrm,'(a8,i2)',advance='no') ' degree ',lmax
        write(ioPrm,'(a7,g10.5)',advance='no') ' layer ',P%L(j)%depth
        if (P%smoothed) then
          ! OUTPUT SMOOTH MODEL SO DO NOT WRITE THE REGULARISATION PARAMETERS...
          ! INSTEAD, START NEW LINE:
          write(ioPrm,*)
        else
          write(ioPrm,'(a5,2g10.5,a8,g12.6)') ' reg ',P%L(j)%alpha, P%L(j)%beta,' weight ',P%L(j)%gamma
        end if
        write(ioPrm,*) '  l   m   value       min       max'
        do i=1,P%nF
            if (.not.P%c(j,i)%exists) then
                cycle
            end if
            write(ioPrm,'(2i4,g15.7)',advance='no') P%F(i)%l,P%F(i)%m,P%c(j,i)%value
            if (P%c(j,i)%frozen) then
                if_var_char = 'const'
            else
                if_var_char = 'range'
            end if
            write(ioPrm,'(2g15.7,a6)') P%c(j,i)%min,P%c(j,i)%max,if_var_char
        end do
        write(ioPrm,*)
    end do
    write(ioPrm,*)

    close(ioPrm)

  end subroutine writeLayers_modelParam

  ! **********************************************************************
  ! * Use to output a vector of model parameters, normally for the full
  ! * sensitivity matrix (NOTE ioSens HARDCODED).
  ! * Assume that the vector consists of parameters of the same type
  ! * (size, etc); could be verified on input.
  ! * Assume also that the file is already open; will be closed in the
  ! * calling subroutine.
  subroutine writeVec_modelParam(np,P,header,cfile)

    implicit none
    integer, intent(in)                     :: np
    type (modelParam_t), intent(in)         :: P(np)
    character(*), intent(in)                :: header, cfile
    ! * EOP

    integer                                 :: lmax,i,j,k,istat
    logical                                 :: opened
    character(6)                            :: if_log_char,if_var_char

    inquire(file=cfile, opened=opened)
    if (.not. opened) then
        open(unit=ioSens, file=cfile, status='unknown', iostat=istat)
    endif

    lmax = getDegree(P(1))

    write(ioSens,'(a24,i2,a8,i2)') 'Format: harmonic layers ',P(1)%nL,' degree ',lmax
    write(ioSens,*)

    do j=1,P(1)%nL
        if (P(1)%L(j)%if_log) then
            if_log_char = 'log'
        elseif (P(1)%L(j)%if_tan) then
            if_log_char = 'tan'
        else
            if_log_char = 'linear'
        end if
        lmax = getDegree(P(1),j)
        write(ioSens,'(a6)',advance='no') if_log_char
        write(ioSens,'(a8,i2)',advance='no') ' degree ',lmax
        write(ioSens,'(a7,g10.5)',advance='no') ' layer ',P(1)%L(j)%depth
        if (P(1)%smoothed) then
          ! OUTPUT SMOOTH MODEL SO DO NOT WRITE THE REGULARISATION PARAMETERS...
          ! INSTEAD, START NEW LINE:
          write(ioPrm,*)
        else
          write(ioPrm,'(a5,2g10.5,a8,g12.6)') ' reg ',P(1)%L(j)%alpha, P(1)%L(j)%beta,' weight ',P(1)%L(j)%gamma
        end if
        write(ioSens,*) '  l   m   value       min       max'
        do i=1,P(1)%nF
            if (.not.P(1)%c(j,i)%exists) then
                cycle
            end if
            write(ioSens,'(2i4)',advance='no') P(1)%F(i)%l,P(1)%F(i)%m
            do k=1,size(P)
                write(ioSens,'(g15.7)',advance='no') P(k)%c(j,i)%value
            end do
            if (P(1)%c(j,i)%frozen) then
                if_var_char = 'const'
            else
                if_var_char = 'range'
            end if
            write(ioSens,'(2g15.7,a6)') P(1)%c(j,i)%min,P(1)%c(j,i)%max,if_var_char
        end do
        write(ioSens,*)
    end do
    write(ioSens,*)

    !... do not close!
    !close(ioSens)

  end subroutine writeVec_modelParam

  ! **********************************************************************
  ! * BOP
  subroutine print_modelParam(P,verbose,comment)

    implicit none
    type (modelParam_t), intent(in)         :: P
    integer, intent(in)                     :: verbose
    character(*),intent(in),optional        :: comment
    ! * EOP

    type (modelCoeff_t)                     :: coeff
    integer                                 :: i,j

    if(.not.P%allocated) then
       call warning('(print_modelParam) parametrization not allocated yet')
       return
    else if(trim(P%type) .ne. 'harmonic') then
       call warning('(print_modelParam) skipped printing this model parametrization type')
       return
    end if

    write(0,*)
    if(present(comment)) then
        write(0,*) node_info,comment
    end if

    if (verbose>0) then
        write(0,'(a12,a50,i4)') node_info,'Number of layers in script: ',P%nL
        do j=1,P%nL
            write(0,'(a12,a46,i2,a2,i4)') node_info,'Number of coefficients in layer ',j,': ',count(.not.P%c(j,:)%frozen)
        end do
        write(0,'(a12,a50,i4)') node_info,'Number of variable parameters in script: ',count(.not.P%c%frozen)
        write(0,*)
    end if

    do j=1,P%nL
        write(0,'(a12,a46,i2,a2,g15.7)') node_info,'Degree and order zero coefficient in layer ',j,': ',P%c(j,1)%value
        end do
    write(0,*)

    if (verbose>0) then
        do j=1,P%nL
            write(0,'(a12,a46,i2,a2,g10.5)') node_info,'Horizontal regularization in layer ',j,': ',P%L(j)%alpha
        end do
        write(0,*)
    end if

    if (verbose>0) then
        do j=1,P%nL
            write(0,'(a12,a46,i2,a2,g10.5)') node_info,'Vertical regularization in layer ',j,': ',P%L(j)%beta
        end do
        write(0,*)
    end if

    if (verbose>0) then
        do j=1,P%nL
            write(0,'(a12,a46,i2,a2,g10.5)') node_info,'Relative layer weight in layer ',j,': ',P%L(j)%gamma
        end do
        write(0,*)
    end if

    if (verbose>3) then
        do i=1,P%nc
            coeff=getCoeff_modelParam(P,i)
            write (0,'(a12,2i8,g17.9)') node_info,coeff%L%num,coeff%F%num,coeff%value
        end do
    end if

  end subroutine print_modelParam

  ! **********************************************************************
  ! * BOP
  subroutine setCrust_modelParam(crust,P)

    type (modelShell_t), intent(in)          :: crust
    type (modelParam_t), intent(inout)       :: P
    ! * EOP
    integer                                  :: nx,ny,status

    if(.not.P%allocated) then
       call errStop('(setCrust_modelParam) parametrization not allocated yet')
    end if

    if(.not.crust%allocated) then
       ! crust not allocated intentionally; do nothing
       P%crust%allocated = .false.
       return
    end if

    if(.not. associated(crust%cond)) then
       call errStop('(setCrust_modelParam) crust not allocated yet')
    end if

    if(associated(P%crust%cond)) then
       deallocate(P%crust%cond, STAT=status)
    end if

    nx = size(crust%cond,1)
    ny = size(crust%cond,2)
    allocate(P%crust%cond(nx,ny), STAT=status)

    P%crust%cond = crust%cond
    P%crust%avg  = sum(crust%cond)/(nx*ny)
    P%crust%variable = crust%variable
    P%crust%allocated = .true.

  end subroutine setCrust_modelParam

  ! ***************************************************************************
  ! * setGrid_modelParam is the routine that sets the pointer to grid
  subroutine setGrid_modelParam(grid,P)

    type (grid_t), intent(in), target          :: grid
    type (modelParam_t), intent(inout)         :: P

    if(.not.P%allocated) then
       call warning('(setGrid_modelParam) parametrization not allocated yet')
       return
    else if(.not. grid%allocated) then
       call errStop('(setGrid_modelParam) input grid not allocated')
    end if

    P%grid => grid

  end subroutine setGrid_modelParam

  ! ***************************************************************************
  ! * setBackground_modelParam is the routine that sets the pointer to rho0
  subroutine setBackground_modelParam(rho0,P)

    type (rscalar), intent(in), target         :: rho0
    type (modelParam_t), intent(inout)         :: P

    if(.not.P%allocated) then
       call warning('(setBackground_modelParam) parametrization not allocated yet')
       return
    else if(.not. rho0%allocated) then
       call errStop('(setBackground_modelParam) input background resistivity not allocated')
    end if

    P%rho0 => rho0

  end subroutine setBackground_modelParam  ! initModel

end module ModelSpace

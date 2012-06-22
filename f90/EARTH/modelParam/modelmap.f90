! *****************************************************************************
module modelmap
  ! Module to initialize the model on the grid once the grid and parametrization
  ! information is available

  use griddef
  use sg_scalar
  use modeldef
  use paramfunc
  implicit none

  INTERFACE mapToGrid
     MODULE PROCEDURE mapToGrid_modelParam
  END INTERFACE

  ! ***************************************************************************
  ! * targetRho0: Target background resistivity for model mappings
  type (rscalar), target, save                  :: targetRho0

Contains


  ! ***************************************************************************
  ! * mapToGrid_modelParam is the routine that generates the 3-D resistivity
  ! * map on the grid using the information stored in the grid, the parametrization
  ! * and the background electrical resistivity on the grid (rho0)
  subroutine mapToGrid_modelParam(P1,P)

    type (modelParam_t), intent(in)         :: P
    type (modelParam_t), intent(inout)      :: P1
    ! local
    type (grid_t), pointer                  :: grid
    type (rscalar)                          :: resist
    !type (rscalar), intent(in), optional            :: resist0
    !real(8), dimension(:,:,:), intent(inout)       :: resist !(nx,ny,nz)
    integer                                         :: i,j,k,l,istat
    integer                                         :: iL,ip
    real(8)                                         :: value,coeff
    type (modelLayer_t), pointer                        :: this_layer
    type (modelFunc_t)                              :: func
    type (modelPoint_t)                             :: point
    real(8)                                         :: crust_depth
    logical                                         :: background

    if(.not.P%allocated) then
       call warning('(mapToGrid_modelParam) parametrization not allocated yet')
       return
    end if

    P1 = P

    ! First initialize resistivity in air and possibly crust, if given

    write(0,*) node_info,'Mapping to grid from model parameter of type: ',trim(P%type)

    if(associated(P%rho0)) then
       background = .true.
    else
       write(0,*) node_info,'Background resistivity not associated and will not be used for this mapping'
       background = .false.
    end if
    grid => P%grid

    crust_depth = grid%z(grid%nzAir+1) - grid%z(grid%nzAir+grid%nzCrust+1) !KM2M*(EARTH_R-CRUST_R)
    write(0,'(a12,a21,f12.3,a34)') node_info,'Using crustal depth ',crust_depth,' to define thinsheet resistivities'

    if (trim(P%type) .eq. 'grid') then

        ! set output to the resistivity on the grid
        resist = P%rho

        ! ... and insert thinsheet if allocated (if no crust, nzCrust == 0)
        do k=grid%nzAir+1,grid%nzAir+grid%nzCrust !
          do i=1,grid%nx
            do j=1,grid%ny
              if(P%crust%allocated) then
                  if(P%crust%cond(i,j) <= R_ZERO) then
                    write(0, *) 'Error: (mapToGrid_modelParam) negative or zero conductance at',i,j,k
                    stop
                  end if
                  resist%v(i,j,k) = crust_depth/P%crust%cond(i,j)
              else
                  ! if not allocated leave it alone
                  !resist%v(i,j,k) = crust_depth/P%crust%avg
              end if
            end do
          end do
        end do

    else


        call create_rscalar(grid,resist,CENTER)

        resist%v = R_ZERO

        forall (i=1:grid%nx, j=1:grid%ny, k=1:grid%nzAir)
            resist%v(i,j,k) = 1/SIGMA_AIR
        end forall

        do k=grid%nzAir+1,grid%nz

          ! Find current layer by locating the upper boundary of a cell
          do l=1,P%nL
            if (in_layer(grid%r(k),P%L(l))) then
              this_layer => P%L(l)
              exit
            end if
          end do

          iL = this_layer%num

          do i=1,grid%nx
            do j=1,grid%ny

              point%phi   = (grid%ph(i) + grid%ph(i+1))/2
              point%theta = (grid%th(j) + grid%th(j+1))/2
              point%r     = (grid%r(k) + grid%r(k+1))/2

              ! Sum up the coeffs * F_at_point in the given layer
              value = 0.0d0
              do ip=1,P%nF

                coeff = P%c(iL,ip)%value
                func = P%F(ip)
                if(coeff /= 0.0d0) then
                  value = value + coeff * F_at_point(func,point)
                end if

              end do

              if (background) then

                  ! use background resistivity AND parametrization to compute resistivity
                  ! on the grid; parameterize deviations only
                  if (this_layer%if_log) then
                    resist%v(i,j,k) = P%rho0%v(i,j,k) * exp(log(10.0d0)*value)
                  elseif (this_layer%if_tan) then
                    resist%v(i,j,k) = P%rho0%v(i,j,k) * ((2.0d0/PI)*atan(value)+1.0d0)
                  else
                    resist%v(i,j,k) = P%rho0%v(i,j,k) + value
                  end if

              else

                  ! compute background resistivity; generally reasonable to use log10
                  ! for background even if the variable parametrization is arctan
                  if (this_layer%if_log) then
                    resist%v(i,j,k) = exp(log(10.0d0)*value)
                  elseif (this_layer%if_tan) then
                    resist%v(i,j,k) = ((2.0d0/PI)*atan(value)+1.0d0)
                  else
                    resist%v(i,j,k) = value
                  end if

                  ! overwrite thinsheet (if no crust, nzCrust == 0)
                  if (k <= grid%nzAir+grid%nzCrust) then
                      if(P%crust%allocated) then
                          if(P%crust%cond(i,j) <= R_ZERO) then
                            write(0, *) 'Error: (mapToGrid_modelParam) negative or zero conductance at',i,j,k
                            stop
                          end if
                          resist%v(i,j,k) = crust_depth/P%crust%cond(i,j)
                      else
                          ! if not allocated use the average
                          !resist%v(i,j,k) = crust_depth/P%crust%avg
                      end if
                  end if

              end if

              if (resist%v(i,j,k) <= 0.0d0) then
                write(0, *) 'Error: (initModel) negative or zero resistivity at',i,j,k,'=',resist%v(i,j,k)
                if (background) then
                    write(0, *) 'Background and logicals: ',P%rho0%v(i,j,k),this_layer%if_log,this_layer%if_tan
                end if
                stop
              end if

            end do
          end do
        end do
     end if

     P1%rho = resist
     P1%type = 'grid'
     call deall_rscalar(resist)
     if (background) then
        write(0,'(a12,a50)') node_info,'Mapping to grid complete. Used a background model.'
     else
        write(0,'(a12,a51)') node_info,'Mapping to grid complete. No background model used.'
     end if

  end subroutine mapToGrid_modelParam  ! initModel


  ! ***************************************************************************
  ! * initResist initializes resistivity structure with the known values above
  ! * mantle; if crust distribution is known, this is also used; if not, this
  ! * case is handled by setting grid%nzCrust=0, so that no resistivity
  ! * values below air layers are initialized here.

  subroutine insertShell(grid,resist,crust)

    type (grid_t), intent(in)					:: grid
    type (rscalar), intent(inout)                     :: resist
	type (modelShell_t), intent(in)		            :: crust
	!real(8), dimension(:,:,:), intent(out)			:: resist !(nx,ny,nz)
	real(8)											:: crust_depth
	integer											:: i,j,k

	if(.not. resist%allocated) then
	   write(0, *) 'Error: (insertShell) resistivity not allocated'
	   stop
	end if

	crust_depth = grid%z(grid%nzAir+1) - grid%z(grid%nzAir+grid%nzCrust+1) !KM2M*(EARTH_R-CRUST_R)
	write(0,'(a12,a21,f12.3,a34)') node_info,'Using crustal depth ',crust_depth,' to define thinsheet resistivities'

	do i=1,grid%nx
	  do j=1,grid%ny
		do k=grid%nzAir+1,grid%nzAir+grid%nzCrust ! if no crust, nzCrust == 0
		  if(crust%allocated) then
			  if(crust%cond(i,j) <= R_ZERO) then
				write(0, *) 'Error: (insertShell) negative or infinite resistivity at',i,j,k
				stop
			  end if
			  resist%v(i,j,k) = crust_depth/crust%cond(i,j)
		  else
		      ! if not allocated use the average
              resist%v(i,j,k) = crust_depth/crust%avg
		  end if
		end do
	  end do
	end do


  end subroutine insertShell


  ! ***************************************************************************
  ! * initModel is the routine that generates the 3-D resistivity map on the
  ! * grid using the information stored in the grid and the parametrization

  subroutine initModel(param,resist,grid,resist0)

    type (grid_t), intent(in)						:: grid
    type (modelParam_t), intent(in)					:: param
	type (rscalar), intent(out)					    :: resist
    type (rscalar), intent(in), optional            :: resist0
	!real(8), dimension(:,:,:), intent(inout)		:: resist !(nx,ny,nz)
	integer											:: i,j,k,l,istat
	integer											:: iL,ip
	real(8)											:: value,coeff
	type (modelLayer_t), pointer						:: this_layer
	type (modelFunc_t)								:: func
	type (modelPoint_t)								:: point
    real(8)                                         :: crust_depth

	! First initialize resistivity in air and possibly crust, if given

	write(0,*) node_info,'Mapping to grid from model parameter of type: ',trim(param%type)

    crust_depth = grid%z(grid%nzAir+1) - grid%z(grid%nzAir+grid%nzCrust+1) !KM2M*(EARTH_R-CRUST_R)
    write(0,'(a12,a21,f12.3,a34)') node_info,'Using crustal depth ',crust_depth,' to define thinsheet resistivities'

    if (trim(param%type) .eq. 'grid') then

        ! set output to the resistivity on the grid
        resist = param%rho

        ! ... and insert thinsheet if allocated (if no crust, nzCrust == 0)
        do k=grid%nzAir+1,grid%nzAir+grid%nzCrust !
          do i=1,grid%nx
            do j=1,grid%ny
              if(param%crust%allocated) then
                  if(param%crust%cond(i,j) <= R_ZERO) then
                    write(0, *) 'Error: (initModel) negative or zero conductance at',i,j,k
                    stop
                  end if
                  resist%v(i,j,k) = crust_depth/param%crust%cond(i,j)
              else
                  ! if not allocated leave it alone
                  !resist%v(i,j,k) = crust_depth/param%crust%avg
              end if
            end do
          end do
        end do

        !call insertShell(grid,resist,param%crust)

    else


        call create_rscalar(grid,resist,CENTER)

        resist%v = R_ZERO

        forall (i=1:grid%nx, j=1:grid%ny, k=1:grid%nzAir)
            resist%v(i,j,k) = 1/SIGMA_AIR
        end forall

        do k=grid%nzAir+1,grid%nz

          ! Find current layer by locating the upper boundary of a cell
          do l=1,param%nL
            if (in_layer(grid%r(k),param%L(l))) then
              this_layer => param%L(l)
              exit
            end if
          end do

          iL = this_layer%num

          do i=1,grid%nx
            do j=1,grid%ny

              point%phi   = (grid%ph(i) + grid%ph(i+1))/2
              point%theta = (grid%th(j) + grid%th(j+1))/2
              point%r     = (grid%r(k) + grid%r(k+1))/2

              ! Sum up the coeffs * F_at_point in the given layer
              value = 0.0d0
              do ip=1,param%nF

                coeff = param%c(iL,ip)%value
                func = param%F(ip)
                if(coeff /= 0.0d0) then
                  value = value + coeff * F_at_point(func,point)
                end if

              end do

              if (present(resist0)) then

                  ! use background resistivity AND parametrization to compute resistivity
                  ! on the grid; parameterize deviations only
                  if (this_layer%if_log) then
                    resist%v(i,j,k) = resist0%v(i,j,k) * exp(log(10.0d0)*value)
                  elseif (this_layer%if_tan) then
                    resist%v(i,j,k) = resist0%v(i,j,k) * ((2.0d0/PI)*atan(value)+1.0d0)
                  else
                    resist%v(i,j,k) = resist0%v(i,j,k) + value
                  end if

              else

                  ! compute background resistivity; generally reasonable to use log10
                  ! for background even if the variable parametrization is arctan
                  if (this_layer%if_log) then
                    resist%v(i,j,k) = exp(log(10.0d0)*value)
                  elseif (this_layer%if_tan) then
                    resist%v(i,j,k) = ((2.0d0/PI)*atan(value)+1.0d0)
                  else
                    resist%v(i,j,k) = value
                  end if

                  ! overwrite thinsheet (if no crust, nzCrust == 0)
                  if (k <= grid%nzAir+grid%nzCrust) then
                      if(param%crust%allocated) then
                          if(param%crust%cond(i,j) <= R_ZERO) then
                            write(0, *) 'Error: (initModel) negative or zero conductance at',i,j,k
                            stop
                          end if
                          resist%v(i,j,k) = crust_depth/param%crust%cond(i,j)
                      else
                          ! if not allocated use the average
                          !resist%v(i,j,k) = crust_depth/param%crust%avg
                      end if
                  end if

              end if

              if (resist%v(i,j,k) <= 0.0d0) then
                write(0, *) 'Error: (initModel) negative or zero resistivity at',i,j,k,'=',resist%v(i,j,k)
                if (present(resist0)) then
                    write(0, *) 'Background and logicals: ',resist0%v(i,j,k),this_layer%if_log,this_layer%if_tan
                end if
                stop
              end if

            end do
          end do
        end do

!        ! if creating a background resistivity, overwrite with thinsheet (if no crust, nzCrust == 0)
!        if (.not. present(resist0)) then
!            forall (i=1:grid%nx, j=1:grid%ny, k=grid%nzAir+1:grid%nzAir+grid%nzCrust)
!
!              if(param%crust%allocated) then
!                  if(param%crust%cond(i,j) <= R_ZERO) then
!                    write(0, *) 'Error: (initModel) negative or zero conductance at',i,j,k
!                    stop
!                  end if
!                  resist%v(i,j,k) = crust_depth/param%crust%cond(i,j)
!              else
!                  ! if not allocated use the average
!                  resist%v(i,j,k) = crust_depth/param%crust%avg
!              end if
!
!            end forall
!        end if

        !        if (present(resist0)) then
!            forall (i=1:grid%nx, j=1:grid%ny, k=grid%nzAir+1:grid%nzAir+grid%nzCrust)
!                resist%v(i,j,k) = resist0%v(i,j,k)
!            end forall
!        else
!            call insertShell(grid,resist,param%crust)
!        end if

	end if


  end subroutine initModel	! initModel


end module modelmap

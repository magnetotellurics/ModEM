! *****************************************************************************
module modelmap
  ! Module to initialize the model on the grid once the grid and parametrization
  ! information is available

  use griddef
  use sg_scalar
  use modeldef
  use paramfunc
  implicit none


Contains


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
	write(0,'(a12,a20,g12.3,a34)') node_info,'Using crustal depth ',crust_depth,' to define thinsheet resistivities'

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
		      !Leave alone if not allocated ...
              !resist%v(i,j,k) = R_ZERO
		  end if
		end do
	  end do
	end do


  end subroutine insertShell


  ! ***************************************************************************
  ! * initModel is the routine that generates the 3-D resistivity map on the
  ! * grid using the information stored in the grid and the parametrization

  subroutine initModel(grid,param,resist)

    type (grid_t), intent(in)						:: grid
    type (modelParam_t), intent(in)					:: param
	type (rscalar), intent(out)					    :: resist
	!real(8), dimension(:,:,:), intent(inout)		:: resist !(nx,ny,nz)
	integer											:: i,j,k,l,istat
	integer											:: iL,ip
	real(8)											:: value,coeff
	type (modelLayer_t), pointer						:: this_layer
	type (modelFunc_t)								:: func
	type (modelPoint_t)								:: point

	! First initialize resistivity in air and possibly crust, if given

	write(0,*) node_info,'Mapping to grid from model parameter of type: ',trim(param%type)

    if (trim(param%type) .eq. 'grid') then

        resist = param%rho

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

              if (this_layer%if_log) then
                resist%v(i,j,k) = exp(log(10.0d0)*value)
              else
                resist%v(i,j,k) = value
              end if

              if (resist%v(i,j,k) <= 0.0d0) then
                write(0, *) 'Error: (initModel) negative or zero resistivity at',i,j,k
                stop
              end if

            end do
          end do
        end do

	end if

    ! Finally, insert thinsheet if allocated
    call insertShell(grid,resist,param%crust)


  end subroutine initModel	! initModel


end module modelmap

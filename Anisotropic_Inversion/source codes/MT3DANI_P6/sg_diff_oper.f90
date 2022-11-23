! ****************************************************************************
! Generic differential operators like Div, Grad, and Curl on stagggered grid.
! Not all routines are used at preseent, but are useful for debugging and
! ancilliary calculations. Belongs to SG_Basics class: staggered cartesian grid, 
! data types defined on this grid, and operations defined on these data types. Not
! specific to EM problem, no dependency on outside (from other classes) modules. 

module sg_diff_oper

  use math_constants
  use sg_vector
  use sg_scalar
  implicit none

  ! Div computes the divergence of a vector
  ! Grad computes the gradient of a scalar
  ! Curl computes the curl of a vector

  public                                :: Div, Grad

Contains

  ! ***************************************************************************
  ! * Div computes the divergence for a complex vector
  subroutine Div(inV, outSc)

    implicit none
    type (cvector), intent(in)                      :: inV    
    type (cscalar), intent(inout)                   :: outSc 
    integer                                         :: ix, iy, iz 
    
    IF(.not.inV%allocated) THEN
 	WRITE(0,*) 'inV not allocated in Div'
 	STOP
    ENDIF 
    
    IF(.not.outSc%allocated) THEN
 	WRITE(0,*) 'outSc not allocated in Div'
 	STOP
    ENDIF    
   
    ! Check whether all the inputs/ outputs involved 
    ! are even of the same size
    if ((inV%nx == outSc%nx).and.&
         (inV%ny == outSc%ny).and.&
         (inV%nz == outSc%nz)) then

       if ((inV%gridType == EDGE).and.(outSc%gridType == CORNER)) then

          ! computation done only for internal nodes
          do ix = 2, outSc%nx
             do iy = 2, outSc%ny
                do iz = 2, outSc%nz

                   outSc%v(ix, iy, iz) = (inV%x(ix, iy, iz) - &
                        inV%x(ix - 1, iy, iz))/ inV%grid%delX(ix) + &
                        (inV%y(ix, iy, iz) - inV%y(ix, iy - 1, iz))/&
                        inV%grid%delY(iy) + &
                        (inV%z(ix, iy, iz) - inV%z(ix, iy, iz - 1))/&
                        inV%grid%delZ(iz) 

                enddo   ! iz
             enddo      ! iy
          enddo         ! ix

       else if ((inV%gridType == FACE).and.(outSc%gridType == CENTER)) then

          ! computation done only for internal nodes
          do ix = 1, outSc%nx
             do iy = 1, outSc%ny
                do iz = 1, outSc%nz

                   outSc%v(ix, iy, iz) = (inV%x(ix+1, iy, iz) - &
                        inV%x(ix, iy, iz))/ inV%grid%dx(ix) + &
                        (inV%y(ix, iy+1, iz) - inV%y(ix, iy, iz))/&
                        inV%grid%dy(iy) + &
                        (inV%z(ix, iy, iz+1) - inV%z(ix, iy, iz))/&
                        inV%grid%dz(iz) 

                enddo   ! iz
             enddo      ! iy
          enddo         ! ix

       else
          write (0, *) 'Div: not compatible usage for existing data types'
       end if

    else
       write(11, *) 'Error-all input/ output in Div are not same size'
    endif

  end subroutine Div  ! Div


  ! ***************************************************************************
  ! * Grad computes the gradient for a complex scalar
  subroutine Grad(inSc, outV, job)

    implicit none
    type (rscalar), intent(in)            :: inSc   
    type (cvector), intent(inout)         :: outV    
    integer                               :: ix, iy, iz,job   
	real(kind=prec)                       :: HGrad
    complex(kind=prec)                    :: Grad1
       
    IF(.not.inSc%allocated) THEN
 	WRITE(0,*) 'inSc not allocated in Grad'
 	STOP
    ENDIF 
    
    IF(.not.outV%allocated) THEN
 	WRITE(0,*) 'outV not allocated in Grad'
 	STOP
    ENDIF 

    if ((inSc%nx == outV%nx).and.&
         (inSc%ny == outV%ny).and.&
         (inSc%nz == outV%nz)) then

       if ((outV%gridType == EDGE).and.(inSc%gridType == CORNER)) then

          ! the conversion in Grad is only done for interior nodes
          do ix = 1, outV%nx 
             do iy = 2, outV%ny
                do iz = 2, outV%nz

                   HGrad = (inSc%v(ix+1, iy, iz) - &
                        inSc%v(ix, iy, iz))/ inSc%grid%dx(ix)				    				 
                   if(job.eq.0) then
				     Grad1 = cmplx(HGrad,R_ZERO)
				   else
                     Grad1 = cmplx(R_ZERO,HGrad)
				   end if
                   outV%x(ix, iy, iz) = outV%x(ix, iy, iz) + Grad1

                enddo
             enddo
          enddo

          do ix = 2, outV%nx 
             do iy = 1, outV%ny
                do iz = 2, outV%nz

                   HGrad = (inSc%v(ix, iy+1, iz) - &
                        inSc%v(ix, iy, iz))/ inSc%grid%dy(iy)
                   if(job.eq.0) then
				     Grad1 = cmplx(HGrad,R_ZERO)
				   else
                     Grad1 = cmplx(R_ZERO,HGrad)
				   end if
                   outV%y(ix, iy, iz) = outV%y(ix, iy, iz) + Grad1

                enddo
             enddo
          enddo

          do ix = 2, outV%nx 
             do iy = 2, outV%ny
                do iz = 1, outV%nz  

                   HGrad = (inSc%v(ix, iy, iz+1) - &
                        inSc%v(ix, iy, iz))/ inSc%grid%dz(iz)
                   if(job.eq.0) then
				     Grad1 = cmplx(HGrad,R_ZERO)
				   else
                     Grad1 = cmplx(R_ZERO,HGrad)
				   end if
                   outV%z(ix, iy, iz) = outV%z(ix, iy, iz) + Grad1

                enddo
             enddo
          enddo

       else if ((outV%gridType == FACE).and.(inSc%gridType == CENTER)) then

          ! the conversion in Grad is only done for interior nodes
          do ix = 2, outV%nx 
             do iy = 1, outV%ny
                do iz = 1, outV%nz

                   HGrad = (inSc%v(ix, iy, iz) - &
                        inSc%v(ix-1, iy, iz))/ inSc%grid%delX(ix)
                   if(job.eq.0) then
				     Grad1 = cmplx(HGrad,R_ZERO)
				   else
                     Grad1 = cmplx(R_ZERO,HGrad)
				   end if
                   outV%x(ix, iy, iz) = outV%x(ix, iy, iz) + Grad1

                enddo
             enddo
          enddo

          do ix = 1, outV%nx 
             do iy = 2, outV%ny
                do iz = 1, outV%nz

                   HGrad = (inSc%v(ix, iy, iz) - &
                        inSc%v(ix, iy-1, iz))/ inSc%grid%delY(iy)
                   if(job.eq.0) then
				     Grad1 = cmplx(HGrad,R_ZERO)
				   else
                     Grad1 = cmplx(R_ZERO,HGrad)
				   end if
                   outV%y(ix, iy, iz) = outV%y(ix, iy, iz) + Grad1

                enddo
             enddo
          enddo

          do ix = 1, outV%nx 
             do iy = 1, outV%ny
                do iz = 2, outV%nz  

                   HGrad = (inSc%v(ix, iy, iz) - &
                        inSc%v(ix, iy, iz-1))/ inSc%grid%delZ(iz)
                   if(job.eq.0) then
				     Grad1 = cmplx(HGrad,R_ZERO)
				   else
                     Grad1 = cmplx(R_ZERO,HGrad)
				   end if
                   outV%z(ix, iy, iz) = outV%z(ix, iy, iz) + Grad1

                enddo
             enddo
          enddo

       else
          write (0, *) 'Grad: not compatible usage for existing data types'
       end if

    else
       write(11, *) 'Error-all input/ output in Grad are not same size'
    endif

  end subroutine Grad  ! Grad


end module sg_diff_oper

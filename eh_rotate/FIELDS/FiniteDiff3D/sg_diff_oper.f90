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
  public                                :: Curl

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
  subroutine Grad(inSc, outV)

    implicit none
    type (cscalar), intent(in)            :: inSc   
    type (cvector), intent(inout)         :: outV    
    integer                               :: ix, iy, iz   
    
       
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

                   outV%x(ix, iy, iz) = (inSc%v(ix+1, iy, iz) - &
                        inSc%v(ix, iy, iz))/ inSc%grid%dx(ix)

                enddo
             enddo
          enddo

          do ix = 2, outV%nx 
             do iy = 1, outV%ny
                do iz = 2, outV%nz

                   outV%y(ix, iy, iz) = (inSc%v(ix, iy+1, iz) - &
                        inSc%v(ix, iy, iz))/ inSc%grid%dy(iy)

                enddo
             enddo
          enddo

          do ix = 2, outV%nx 
             do iy = 2, outV%ny
                do iz = 1, outV%nz  

                   outV%z(ix, iy, iz) = (inSc%v(ix, iy, iz+1) - &
                        inSc%v(ix, iy, iz))/ inSc%grid%dz(iz)

                enddo
             enddo
          enddo

       else if ((outV%gridType == FACE).and.(inSc%gridType == CENTER)) then

          ! the conversion in Grad is only done for interior nodes
          do ix = 2, outV%nx 
             do iy = 1, outV%ny
                do iz = 1, outV%nz

                   outV%x(ix, iy, iz) = (inSc%v(ix, iy, iz) - &
                        inSc%v(ix-1, iy, iz))/ inSc%grid%delX(ix)

                enddo
             enddo
          enddo

          do ix = 1, outV%nx 
             do iy = 2, outV%ny
                do iz = 1, outV%nz

                   outV%y(ix, iy, iz) = (inSc%v(ix, iy, iz) - &
                        inSc%v(ix, iy-1, iz))/ inSc%grid%delY(iy)

                enddo
             enddo
          enddo

          do ix = 1, outV%nx 
             do iy = 1, outV%ny
                do iz = 2, outV%nz  

                   outV%z(ix, iy, iz) = (inSc%v(ix, iy, iz) - &
                        inSc%v(ix, iy, iz-1))/ inSc%grid%delZ(iz)

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


  ! ***************************************************************************
  ! * computes the curl of the vector V1 with results in V2 for complex vectors
  subroutine Curl(V1, V2)

    implicit none
    type(cvector), intent(in)               :: V1
    type(cvector), intent(inout)            :: V2
    integer                                 :: nx,ny,nz,ix,iy,iz
    
    IF(.not.V1%allocated) THEN
 	WRITE(0,*) 'V1 not allocated in Curl'
 	STOP
    ENDIF 
    
    IF(.not.V2%allocated) THEN
 	WRITE(0,*) 'V2 not allocated in Curl'
 	STOP
    ENDIF 
    

    Nx = V1%nx
    Ny = V1%ny
    Nz = V1%nz

    if((V1%gridType .eq. EDGE).and.(V2%gridType .eq. FACE)) then

       !   H = curl E (for E defined on edges)
       !  Hx
       do iy = 1,Ny
          do iz = 1,Nz
             V2%x(:,iy,iz) =  &
                  (V1%z(:,iy+1,iz)-V1%z(:,iy,iz))*V1%grid%dYinv(iy) &
                  -(V1%y(:,iy,iz+1)-V1%y(:,iy,iz))*V1%grid%dZinv(iz)
          enddo
       enddo

       !  Hy
       do iz = 1,Nz
          do ix = 1,Nx
             V2%y(ix,:,iz) = &
                  (V1%x(ix,:,iz+1)-V1%x(ix,:,iz))*V1%grid%dZinv(iz) &
                  -(V1%z(ix+1,:,iz)-V1%z(ix,:,iz))*V1%grid%dXinv(ix)
          enddo
       enddo

       !  Hz
       do ix = 1,Nx
          do iy = 1,Ny
             V2%z(ix,iy,:) = &
                  (V1%y(ix+1,iy,:)-V1%y(ix,iy,:))*V1%grid%dXinv(ix) &
                  -(V1%x(ix,iy+1,:)-V1%x(ix,iy,:))*V1%grid%dYinv(iy)
          enddo
       enddo

       elseif((V2%gridType .eq. EDGE).and.(V1%gridType .eq. FACE)) then

          !   E = curl H (for H defined on edges)
          !   NOTE: for this case boundary edge nodes are not computed
          !  Ex
       do iy = 2,Ny
          do iz = 2,Nz
             V2%x(:,iy,iz) =  &
                  (V1%z(:,iy,iz)-V1%z(:,iy-1,iz))*V1%grid%delYinv(iy) &
                  -(V1%y(:,iy,iz)-V1%y(:,iy,iz-1))*V1%grid%delZinv(iz)
          enddo
       enddo

       !  Ey
       do iz = 2,Nz
          do ix = 2,Nx
             V2%y(ix,:,iz) = &
                  (V1%x(ix,:,iz)-V1%x(ix,:,iz-1))*V1%grid%delZinv(iz) &
                  -(V1%z(ix,:,iz)-V1%z(ix-1,:,iz))*V1%grid%delXinv(ix)
          enddo
       enddo

       !  Ez
       do ix = 2,Nx
          do iy = 2,Ny
             V2%z(ix,iy,:) = &
                  (V1%y(ix,iy,:)-V1%y(ix-1,iy,:))*V1%grid%delXinv(ix) &
                  -(V1%x(ix,iy,:)-V1%x(ix,iy-1,:))*V1%grid%delYinv(iy)
          enddo
       enddo
    else

       !  incompatible grids ...
       write(0,*) 'Curl: vector types not compatible for curl operator'
    endif

  end subroutine Curl  ! Curl


end module sg_diff_oper

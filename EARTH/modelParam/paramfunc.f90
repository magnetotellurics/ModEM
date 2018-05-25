! *****************************************************************************
module paramfunc
	! This module contains subroutines that can be passed by reference to the
	! implementation of the parametrization from the original model parameters
	! to the resistivity on the grid

  use polpak, only: legendre_associated,spherical_harmonic,d_factorial
  use griddef
  implicit none

  ! Legendre polynomials are computed once and saved for efficiency;
  ! this becomes very important in an inversion when this module is called repeatedly
  ! The dimensions are (Nt+1,lmax+1,lmax+1).
  ! Once these are computed, legendre_allocated is set to .true.
  real(8), dimension(:,:,:), allocatable, public, save  :: center_P_lm
  real(8), dimension(:,:,:), allocatable, public, save  :: corner_P_lm
  logical, public, save                                 :: legendre_allocated_at_centers=.false.
  logical, public, save                                 :: legendre_allocated_at_corners=.false.
  integer, public, save                                 :: computed_at_centers_to_degree=0
  integer, public, save                                 :: computed_at_corners_to_degree=0

Contains

  ! ***************************************************************************
  ! * the following are the subroutines for the *unnormalized* expressions for
  ! * the spherical harmonics with *no radial dependence*; they are used in the
  ! * (externally defined) parametrizations for the inversion. The reason we
  ! * need them here is that for the operators P and Pt we need the mapping from
  ! * the original model parameters to the resistivity on the grid and backwards.
  ! * For the log10 spherical harmonic parametrization, we compute the resistivity
  ! * from the coefficient through the following expression:
  ! * $\log_10(\rho) = \sum_\i a_i \tau_i(\phi,\th,r)$,
  ! * where $\tau_i$ are the unnormalized spherical harmonics
  ! * sin(m*phi)*P_l^m(cos(theta)) and cos(m*phi)*P_l^m(cos(theta))
  ! * The coefficients and the number of harmonics can be different layer to layer.

subroutine legendre_deallocate_at_centers()
! For efficiency in an inversion, we save the Legendre polynomials once these
! have been computed. We reuse these whenever this module is invoked again.
! Use this subroutine to explicitly deallocate these arrays.

    integer     :: istat

    if (legendre_allocated_at_centers) then
        deallocate(center_P_lm, STAT=istat)
        legendre_allocated_at_centers = .false.
    end if

end subroutine


subroutine legendre_deallocate_at_corners()
! For efficiency in an inversion, we save the Legendre polynomials once these
! have been computed. We reuse these whenever this module is invoked again.
! Use this subroutine to explicitly deallocate these arrays.

    integer     :: istat

    if (legendre_allocated_at_corners) then
        deallocate(corner_P_lm, STAT=istat)
        legendre_allocated_at_corners = .false.
    end if

end subroutine


subroutine legendre_initialize(lmax,grid,gridType)
! usage: compute a version of Schmidt semi-normalised Legendre polynomials
! for all degrees and orders at either centers or corners of the grid.
! however, we have taken the computation of Pnm out of the subroutine to improve
! on efficiency, and call it only when cost is updated.

    integer, intent(in)                                         :: lmax
    type(grid_t), intent(in)                                    :: grid
    character(*), intent(in)                                    :: gridType
    ! local
    type(timer_t) :: timer
    real(8)     :: rone
    integer     :: j,l,m,Nt,istat
    real(8), dimension(grid%ny+1)  :: theta
    real(8), dimension(grid%ny+1,lmax+1,lmax+1)  :: P_lm
    real(8), dimension(lmax+1,lmax+1)  :: factor

    rone = 1.0d0
    factor(:,:) = rone
    P_lm(:,:,:) = 0.0d0
    Nt = grid%ny

    if (trim(gridType) == CENTER)  then
        if (legendre_allocated_at_centers) then
            if (computed_at_centers_to_degree >= lmax) then
                write(0,'(a12,a58,i3,a24)') node_info,'Legendre polynomials at centers pre-computed to degree ',&
                            computed_at_centers_to_degree,'; initialization skipped'
                return
            else
                call legendre_deallocate_at_centers()
            end if
        end if
        write(0,'(a12,a57,i3)') node_info,'Allocating Legendre polynomials at cell centers to degree ',lmax
        theta(1:Nt) = grid%th(1:Nt) + grid%dt/2.
        allocate(center_P_lm(Nt+1,lmax+1,lmax+1),STAT=istat)
        center_P_lm = 0.0d0
    else if (trim(gridType) == CORNER) then
        if (legendre_allocated_at_corners) then
            if (computed_at_corners_to_degree >= lmax) then
                write(0,'(a12,a58,i3,a24)') node_info,'Legendre polynomials at corners pre-computed to degree ',&
                            computed_at_corners_to_degree,'; initialization skipped'
                return
            else
                call legendre_deallocate_at_corners()
            end if
        end if
        write(0,'(a12,a57,i3)') node_info,'Allocating Legendre polynomials at cell corners to degree ',lmax
        theta = grid%th(1:Nt+1)
        allocate(corner_P_lm(Nt+1,lmax+1,lmax+1),STAT=istat)
        corner_P_lm = 0.0d0
    end if

    call reset_time(timer)

    do l = 1,lmax
        do m = 1,l
            factor(l+1,m+1) = ( (-1)**m ) * &
                sqrt ( 2.*rone * ( d_factorial(l - m) ) / ( d_factorial(l + m) ) )
        end do
    end do

    !compute Schmidt Seminormalized Associated Legendre Functions, defined as:
    !Q_n^m(x) = P_n(x) for m=0
    !Q_n^m(x) = (-1)^m sqrt ( 2 * (n-m)! / (n+m)! ) P_n^m for m>0
    do j = 1,Nt
        do m = 0,lmax
            call legendre_associated(lmax,m,cos(theta(j)),P_lm(j,:,m+1))
            P_lm(j,:,m+1) = P_lm(j,:,m+1) * factor(:,m+1)
        end do
    end do

    !save the Legendre polynomials in the module
    if (trim(gridType) == CENTER)  then
        center_P_lm = P_lm
        computed_at_centers_to_degree = lmax
        legendre_allocated_at_centers = .true.
    else if (trim(gridType) == CORNER) then
        corner_P_lm = P_lm
        computed_at_corners_to_degree = lmax
        legendre_allocated_at_corners = .true.
    end if

    if (output_level > 4) then
        write(*,*) 'LEGENDRE 1 0:',P_lm(:,2,1)
        write(*,*) 'LEGENDRE 1 1:',P_lm(:,2,2)
    end if

    if (output_level > 0) then
        write(*,*) node_info,'Done computing Legendre polynomials degree ',lmax,': ',elapsed_time(timer),' secs'
    end if

end subroutine

  ! ***************************************************************************
  ! * computes the value for the (l,m) spherical harmonic term in the expansion
  ! * THIS IS THE INEFFICIENT VERSION OF THE ROUTINE legendre_initialize
  ! * IT IS NO LONGER USED
  real(8) function SphHarm(l,m,phi,theta) result (Y)

	integer, intent(in)				:: l,m
	real(8), intent(in)				:: phi, theta
	real(8), dimension(0:l)				:: legendre
        real(8)                                         :: factor
        integer                                         :: m_abs

        m_abs = iabs(m)

	! I am not making efficient use of this subroutine, that outputs the first
	! l+1 values. However, this usage simplifies the parametrization by allowing
	! to immediately discard the coefficients which would otherwise be multiplied
	! by a zero-valued spherical harmonic term (with l<m)
	call legendre_associated(l,m_abs,cos(theta),legendre)

 !    We are using Schmidt Seminormalized Associated Legendre Functions, defined as:
 !    Q_n^m(x) = P_n(x) for m=0
 !    Q_n^m(x) = (-1)^m sqrt ( 2 * (n-m)! / (n+m)! ) P_n^m for m>0

        factor = ( (-1)**m_abs ) * &
             sqrt ( 2.0d0 * ( d_factorial(l - m_abs) ) / ( d_factorial(l + m_abs) ) )

	Y = legendre(l)

	if (m>0) then

	  Y = Y * factor * cos(m*phi)

	else if(m<0) then

	  Y = Y * factor * sin(iabs(m)*phi)

	end if

  end function SphHarm	! unnormalized spherical harmonics


  ! ***************************************************************************
  ! * to test whether xval belongs to the interval [xmin,xmax)
  integer function I1(xmin,xmax,xval)

	real(8),intent(in)	  :: xmin,xmax,xval

	if((xval>=xmin).and.(xval<xmax)) then

	  I1=1

	else

	  I1=0

	end if

  end function I1 ! identity function

  ! ***************************************************************************
  ! * to test whether (xval,yval) belongs to the rectangle area specified 
  integer function I2(xmin,xmax,xval,ymin,ymax,yval)

	real(8),intent(in)	  :: xmin,xmax,xval
	real(8),intent(in)	  :: ymin,ymax,yval

!	if((xval>=xmin).and.(xval<xmax).and.(yval>=ymin).and.(yval<ymax)) then
	I2 = I1(xmin,xmax,xval) * I1(ymin,ymax,yval)


  end function I2 ! identity function


end module paramfunc

! *****************************************************************************
module maxwells
  ! Central module in the project earth3d, containing the subroutine to solve
  ! numerically the integral form of the quasi-static approximation to Maxwell's
  ! equations on the given grid

  use iotypes
  use field_vectors ! contains global vectors for divergence correction
  use math_constants
  use coreFwd
  use wrapper
  implicit none

  public
      ! store the grid dimensions in the original format
  integer                                       :: nx, ny, nz, nzEarth, nzAir
      ! storing the (spherical) grid in Randie Mackie's format
      ! when allocated, dimensions will be x(nx), y(ny+1), z(nz+1)
  real(8), allocatable, dimension(:)            :: x,y,z


  public        :: SolveMaxwells

Contains

  ! ***************************************************************************
  ! * operator $A_{\rho,\omega}$ : E_i -> E_i stands for the numerical solution
  ! * of Maxwell's equations; A <h> = <b>, where <b> is a known vector of boundary
  ! * values, while <h> is unknown. Both <b> and <h> are complex edge vectors.
  ! * Both <h> and <b> are pre-multiplied by edge lengths, such that <h>=<lH> and
  ! * <b>=<lB>. This is the typical use of this operator. However, the subroutine
  ! * given here is the general version so that for any two vectors <x> and <y>,
  ! * where <y> is known, the output is <x>.
  ! *
  ! * SolveMaxwells is the routine to solve numerically the system of equations.
  !	* Input includes rho in cell centres, scalar angular frequency and the forcing
  !	* vector <y>. Output includes the vector <x> of line integrals lH, defined along
  !	* the edges of the cells, and the integer exit flag (zero if successful).
  ! * Whichever forcing vector <y> is used, the output <x> will by construction
  ! * represent the physical fields pre-multiplied by edge lengths. This allows
  ! * us to compute the divergence of H at each step of the relaxation process.
  ! * However, the divergence of the source term cannot be computed from <y> in
  ! * general, which makes it necessary to pass the source term s to the divergence
  ! * correction routine. This has been done by defining and pre-initializing
  ! * global vector sx,sy,sz externally.
  ! *
  ! * It is therefore most important to note that this subroutine uses the global
  ! * vectors hx,hy,hz,sx,sy,sz for the divergence correction. Care should be taken
  ! * to set them to reasonable initial values before calling SolveMaxwells.
  ! *
  ! * Last modified: Nov 7, 2005 by Anna Kelbert

  subroutine SolveMaxwells(xvec,yvec,omega,resist,ctrls,exitflg)

  ! np  = nhx+nhy+nhz, true dimension of most vectors here
    integer, save                            :: np
    integer                                  :: nhx,nhy,nhz
  ! vectors <x> and <y>, as in A <x> = <y>
	complex(8),dimension(np2),intent(inout)	 :: xvec
	complex(8),dimension(np2),intent(in)	 :: yvec
	real(8), intent(in)						 :: omega
	real(8), dimension(:,:,:), intent(in)	 :: resist
	type (fwdCtrl_t), intent(in)			 :: ctrls
	integer, intent(out)					 :: exitflg
  !---------------------------------------------------------
  ! Local variables
  !---------------------------------------------------------
  ! A info
    real(8),   dimension(:),allocatable      :: sp1	!np5
    integer,   dimension(:),allocatable      :: ija1  !np5
    complex(8),dimension(:),allocatable      :: dp1	!np5
    integer                                  :: rowlen ! used in atimes
  ! LU info
    complex(8),dimension(:),allocatable      :: sp	!np4
    integer,   dimension(:),allocatable      :: ija	!np4
  ! delta <h>
    complex(8),dimension(:),allocatable      :: dxvec	!np
  ! RHS during relaxations, as in A <x> = d
    complex(8),dimension(:),allocatable      :: dvec  !np2
  ! vector of residuals
    complex(8),dimension(:),allocatable      :: rvec,rvec0  !np2
  ! utility vectors, same purpose as described in BiCGSTAB
    complex(8),dimension(:),allocatable      :: pvec,qvec  !np2
    complex(8),dimension(:),allocatable      :: tvec,svec	!np2
  !---------------------------------------------------------
  ! utility variables used in BiCGSTAB algorithm
    complex(8)                               :: alpha ! <r0,r>/<r0,q>
    complex(8)                               :: beta  ! ~ <r0,r>/<r0,r(i-1)>
    complex(8)                               :: lambda	! <t,s>/<t,t>
  ! subsidiary variables (H stands for the Hermitian operator)
    complex(8)                               :: gamma ! <r0,r>
    complex(8)                               :: gamma0	! <r0,r(i-1)>
    complex(8)                               :: yHy,rHq,tHs,tHt
  !---------------------------------------------------------------------
  ! convergence criteria
    real(8)                                  :: errend,errend1,errend2
  ! estimating error in <x> = l*<h>
    complex(8)                               :: herr
  ! estimating error in RHS = A <x>
    real(8)                                  :: erre,errlast,errchk  !toh
  !---------------------------------------------------------------------
  ! variables needed for the divH=0 correction
    integer                                  :: npp
    integer,dimension(:),allocatable		 :: ijapot	!np6
    real(8),dimension(:),allocatable		 :: sppot,sppot1 !np6
    real(8),dimension(:),allocatable		 :: pot  !np1
    real(8),dimension(:),allocatable		 :: rpot,spot,ppot,dpot,zpot !np1
    real(8),dimension(:),allocatable		 :: qpot,wk1 !np1
    real(8)                                  :: errdiv
    real(8), save                            :: errdiv0
    integer                                  :: job
  !---------------------------------------------------------------------
  ! variables needed for the divH=0 correction
	integer									 :: IPOTLOOP,ICTPOT,ICOUNT
	integer									 :: IPOTMAX,IPOTINT,IPOTMAXMAX
	integer									 :: IPOTLOOPMAX
	integer									 :: NRELMAX,NRELDIVH
  !---------------------------------------------------------------------
  ! variables used for testing of initpot3 and set_op3
    integer,   dimension(100)                :: itest
    integer                                  :: i,j,iflag,idiv,itest_n

    complex(8)                               :: cone, czero  ! constant
    real(8)									 :: eps
!------------------------------------------------------------------
! The following variable determines if this routine has been called
! previously (e.g. for a previous frequency) - which if so, changes
! the initialisation of the values used in the relaxation
!------------------------------------------------------------------
    integer, save                            :: call_flag

    data call_flag/-1/

!------------------------------------------------------------------
! About to increment call_flag by one
!------------------------------------------------------------------
    call_flag = call_flag+1

    errlast   =  1.D99
	exitflg	  = -1


!toh: put the stopping criterion of the first frequency order of magnitude
!     more strict than others. hope this may work for a better accuracy of a
!     frequency loop in ASCENDING order. [22/FEB/1999]
!ak: don't know if greater accuracy is required for dh; errend = solution update tolerance
!ak(08/2011): don't change accuracy based on call_flag (deal with such issues at a higher level)

!	if ( nfreq > 1 .and. ifreq == 1 ) then
!	if ( call_flag == 0 ) then
!
!	  errend = ctrls%errend/10.d0
!
!	else

	  errend = ctrls%errend

!	end if



!-------------------------------------------------------------
!	  Initialize real 1 and 0 in the complex plane
!-------------------------------------------------------------
	cone  = dcmplx( 1.0d0,0.0d0)
	czero = dcmplx( 0.0d0,0.0d0)

	eps = 1.0e-30


!-------------------------------------------------------------
!     set up the operator needed for the divH=0 correction
!-------------------------------------------------------------
!toh:29/FEB/2000
	npp=(nz-1)*(2+(ny-1)*nx)
!           ----> from k=3 to N, pot should be specified....

	allocate(sppot(np6),sppot1(np6),ijapot(np6),pot(np1))
	allocate(rpot(np1),spot(np1),ppot(np1),dpot(np1),zpot(np1))
	allocate(qpot(np1),wk1(np1))

	call initpot3uBC3(nx,ny,nz,npp,sppot,sppot1,ijapot,x,y,z)

	! (AK) print *,'for initpot3......'
	itest_n=0
	do i=1,npp
	   idiv=ijapot(i+1)-ijapot(i)
	   iflag=0
	   do j=1,itest_n
		  if (itest(j) == idiv) iflag=1
	   end do
	   if (iflag == 0) then
		  itest_n=itest_n+1
		  itest(itest_n)=idiv
	   end if
	end do
	! (AK) do i=1,itest_n
	! (AK)   print *,itest(i)
	! (AK) end do

	call icdecompr(npp,ijapot,sppot,nx)

	!-------------------------------------------------------------
	! Dimensions used here are
	! np1=nx*(ny+1)*(nz+1); np6=4*np1; (for hx,hy,hz & potentials)
	! np2=3*np1; (for vectors)
	! np4=9*np1; np5=21*np1; (for matrices)
	!-------------------------------------------------------------
	! nz=number of layers in model including air layers
	! nz+1=number of vertical node positions at which hx and hy
	!		  are specified
	! nz-1=number of vertical node positions of variable hx,hy,hz
	!		(boundary values are fixed)
	!-------------------------------------------------------------
	! (nhx,nhy,nhz) <----- number of variables for Hx, Hy and Hz
	!           np  <----- total number of variables
	! For spherical problem,
	! note that Hz at the pole is irrespective of longitude.
	!-------------------------------------------------------------
	! Currently, to go along with the previous version of the code,
	! all vectors that should be length np are length np2, where
	!	np=nhx+nhy+nhz; np2=3*nx*(ny+1)*(nz+1) (ak, 12/Apr/2005)
	!
	! Lemma: sufficient condition for np<=np2 to hold is nx>=1.
	! Proof:
	! nhx+nhy<2*nx(ny+1)(nz+1); nx*(ny-1)+2<=nx(ny+1) iff nx>=1.
	!-------------------------------------------------------------

	nhx = nx*(ny-1)*(nz-1)
	nhy = nx*ny*(nz-1)
	nhz = ( nx*(ny-1)+2 )*(nz-1)
	np  = nhx+nhy+nhz


!toh: Dynamic allocation of sp1 & dp1 here.

	allocate(dp1(np))
	allocate(sp1(np+1:np5))

!ak:
	allocate(sp(np4),ija(np4),ija1(np5))
	allocate(dvec(np2),dxvec(np),rvec(np2),rvec0(np2))
	allocate(pvec(np2),qvec(np2),tvec(np2),svec(np2))


	! (AK) print*
	! (AK) print*,'rTr relaxation'
	! (AK) print*,'number of unknowns [np]',np
	! (AK) print 17,'relaxing',l,'x',m,'x',n,'model'
	! (AK) print 18,'freq=',omega/(2.0d0*pi)
17  format (a8,1x,i2,a1,i2,a1,i2,1x,a5)
18  format (a7,1x,e15.7)

!-------------------------------------------------------------
!  toh: set up 3d em operator
!-------------------------------------------------------------

	! compute SP1 of SP1 <h>=<p> and its diagonal portion SP
	call set_op3u(nx,ny,nz,sp1,ija1,sp,ija,resist,x,y,z, &
				 omega,dp1,np,rowlen)

	! (AK) print *,'for set_op3......'
	itest_n=0
	do i=1,np
	   idiv=ija(i+1)-ija(i)
	   iflag=0
	   do j=1,itest_n
		  if (itest(j) == idiv) iflag=1
	   end do
	   if (iflag == 0) then
		  itest_n=itest_n+1
		  itest(itest_n)=idiv
	   end if
	end do
	! (AK) do i=1,itest_n
	! (AK)   print *,itest(i)
	! (AK) end do


!----------------------------------------------------------------
! Incomplete Cholesky Factorization (ICF) of A
! A ~ LU = U^t*U, where U is an upper triangular matrix
! Note that U^t*U preserves A only in the first row, first column
! and the diagonal elements.
!----------------------------------------------------------------
	call icdecompc(ija,sp,np)

	errend1 = errend
!toh: introduce a new criterion for residual norms
!        errend2 = epsilon(erre)*epsilon(erre)
!toh: on 8/May/99
	errend2 = 1.0d-30
!toh: This was introduced on 8/MAR/2000 to stop iteration re erre
!toh: at a level of 10^(-30)

!-------------------------------------------------------------
!	if called for the first time, initialize divH variables
!-------------------------------------------------------------
	if (call_flag == 0) then

	   call divide_vec_by_l(nx,ny,nz,xvec,x,y,z)
	   call copyd1_d3_b(nx,ny,nz,hx,hy,hz,xvec,x,y,z)
	   call mult_vec_by_l(nx,ny,nz,xvec,x,y,z)
	 ! compute S{ Div[<H>] }dV
	   call calc_vdivBC3(nx,ny,nz,hx,hy,hz,divr,divi,x,y,z)
	   call calc_vdivBC3(nx,ny,nz,sx,sy,sz,divsr,divsi,x,y,z)
	   divr = divr - divsi/(omega*MU_0)	! -omega, ie complex conjugate
	   divi = divi + divsr/(omega*MU_0)	! of the differential operator
	 ! compute Sum( |S{ Div[<H>] }dV| )/Volume
	   call sum_divBC3(nx,ny,nz,divr,divi,z,errdiv0)

	end if

	if (output_level>2) then
	  print*
	  print*,trim(node_info),'initial divH =',errdiv0
	end if

!-------------------------------------------------------------
! Initialize some variables using the user-defined information
! stored in the global variable fwdCtrls. All initial values
! for utility variables specified here for clarity of the code.
!-------------------------------------------------------------

	IPOTINT=ctrls%ipotint
	IPOTMAXMAX=ctrls%ipot_max
	IPOTLOOPMAX=ctrls%ipotloopmax
	NRELMAX=ctrls%nrelmax
	NRELDIVH=ctrls%n_reldivh

	IPOTLOOP=1
	IPOTMAX=ctrls%ipot0


	! compute Norm2(<y>) = <y,y>
    yHy = dot_product(yvec(1:np),yvec(1:np))


	Main: do

	  ICTPOT=0
      IPOTMAX=IPOTMAX+IPOTINT
      if(IPOTMAX > IPOTMAXMAX) IPOTMAX=IPOTMAXMAX

!-------------------------------------------------------------
! If the forcing is zero, set solution to zero and exit
!-------------------------------------------------------------
	  if (cdabs(yHy) < eps) then
		if (output_level>2) then
			print*,trim(node_info),'initial forcing norm =',yHy
			print*,trim(node_info),'exiting with zero solution...'
		end if
		xvec = czero
		exitflg=0
		exit Main
	  end if

!-------------------------------------------------------------
! Calculate initial residual, and begin relaxation iterations
!-------------------------------------------------------------

	! compute <d> = A <x(0)>
      call atimes(xvec,dvec,sp1,ija1,dp1,np,rowlen)

	! compute <r> = <y> - <d>
      rvec(1:np) = cone*(yvec(1:np)-dvec(1:np))

!-------------------------------------------------------------
! toh: Hiro's Bi-CGSTAB... [6/MAR/2000]
!-------------------------------------------------------------

	  ! update vector <r(0)> = (LU)^{-1} (<b> - A<h(0)>)
      call asolvec(rvec,rvec0,sp,ija,np)
	  ! initialize vector <r> = <r(0)>
      rvec(1:np) = rvec0(1:np)
      pvec(1:np) = rvec0(1:np)
	  ! initialize gamma0 = <r(0),r>
	  gamma0 = dot_product(rvec0(1:np),rvec0(1:np))

	  BiCGSTAB: do ICOUNT=1,NRELMAX

	  ! compute <d> = A<p(i-1)>
		call atimes(pvec,dvec,sp1,ija1,dp1,np,rowlen)
	  ! solve (LU)<q(i)> = A<p(i-1)> for <q>
		call asolvec(dvec,qvec,sp,ija,np)
		rHq = dot_product(rvec0(1:np),qvec(1:np))

		alpha = gamma0/rHq
	  ! update vector <s(i)>
		svec(1:np) = rvec(1:np) - alpha*qvec(1:np)
	  ! compute <d> = A<s(i)>
		call atimes(svec,dvec,sp1,ija1,dp1,np,rowlen)
	  ! solve (LU)<t(i)> = A<s(i)> for <t>
		call asolvec(dvec,tvec,sp,ija,np)

		tHs = dot_product(tvec(1:np),svec(1:np))
		tHt = dot_product(tvec(1:np),tvec(1:np))
		lambda = tHs/tHt ! STAB-P => Anyway, SLOW !

		dxvec(1:np) = alpha*pvec(1:np) + lambda*svec(1:np)

	  ! update vector <h(i)>
		xvec(1:np) = xvec(1:np) + dxvec(1:np)

	  ! update vector <r(i)>
		rvec(1:np) = svec(1:np) - lambda*tvec(1:np)
		gamma = dot_product(rvec0(1:np),rvec(1:np))
		beta = (gamma/gamma0)*(alpha/lambda)

	  ! update vector <p(i)>
		pvec(1:np) = rvec(1:np) &
		  + beta*( pvec(1:np) - lambda*qvec(1:np) )

		gamma0 = gamma

	  ! convergence criterion (is the <x> update significant?)
		herr = cdsqrt( dot_product(dxvec(1:np),dxvec(1:np)) &
			/ dot_product(xvec(1:np),xvec(1:np)) )

	  ! estimate error <y> - A <x(i)>
		erre = cdabs( dot_product(rvec(1:np),rvec(1:np)) )/yHy

		if (output_level>2) then
		  print *, trim(node_info), ICOUNT, erre, dreal(herr)
		end if

		ICTPOT=ICTPOT+1
		if((ICOUNT > IPOTMAX).and.(IPOTLOOP <= IPOTLOOPMAX)) then
		! If i'th error > (i-1)'th error, update hx,hy,hz and start again;
		! otherwise either exit or update error and continue BiCGSTAB.
		  errchk=dlog10(errlast/erre)
		  ! If the error has increased in the last iteration, do divH correction
		  if(errchk < 0.0d0) then !Gary's version
		  ! if(errchk > 0.0d0) then !Adam's version
		  !----------------------------------------------------------------
		  ! Static Divergence Correction (SDC)
		  ! Finite values of divH sometimes emerge if the inductive effect
		  ! is weak, such as when either omega or sigma are very small, the
		  ! latter certainly true in the air layers. Correct by introducing
		  ! magnetic gauge grad(phi_m). This is seen to prevent convergence
		  ! to local minima rather than global solutions.
		  !----------------------------------------------------------------
			call divide_vec_by_l(nx,ny,nz,xvec,x,y,z)
			call copyd1_d3_b(nx,ny,nz,hx,hy,hz,xvec,x,y,z)

			call calc_vdivBC3(nx,ny,nz,hx,hy,hz,divr,divi,x,y,z)
			call calc_vdivBC3(nx,ny,nz,sx,sy,sz,divsr,divsi,x,y,z)
			divr = divr - divsi/(omega*MU_0)	! -omega, ie complex conjugate
			divi = divi + divsr/(omega*MU_0)	! of the differential operator

			call sum_divBC3(nx,ny,nz,divr,divi,z,errdiv)
			if (output_level>2) then
			  print*,trim(node_info),'divH before update =',errdiv,IPOTLOOP
			end if

			call relaxpot_chol3BC3(nx,ny,nz,pot,divr,rpot,zpot,qpot,ppot,spot, &
				  dpot,sppot,sppot1,ijapot,wk1,ctrls%n_reldivh)
			job=0
			call updateh3BC3(nx,ny,nz,pot,hx,hy,hz,x,y,z,job)
			call relaxpot_chol3BC3(nx,ny,nz,pot,divi,rpot,zpot,qpot,ppot,spot, &
				  dpot,sppot,sppot1,ijapot,wk1,ctrls%n_reldivh)
			job=1
			call updateh3BC3(nx,ny,nz,pot,hx,hy,hz,x,y,z,job)

			call calc_vdivBC3(nx,ny,nz,hx,hy,hz,divr,divi,x,y,z)
			call calc_vdivBC3(nx,ny,nz,sx,sy,sz,divsr,divsi,x,y,z)
			divr = divr - divsi/(omega*MU_0)	! -omega, ie complex conjugate
			divi = divi + divsr/(omega*MU_0)	! of the differential operator

			call sum_divBC3(nx,ny,nz,divr,divi,z,errdiv)
			if (output_level>2) then
			  print*,trim(node_info),'divH after  update =',errdiv,IPOTLOOP
			end if

			call copyd3_d1_b(nx,ny,nz,hx,hy,hz,xvec,x,y,z)
			call mult_vec_by_l(nx,ny,nz,xvec,x,y,z)

			IPOTLOOP=IPOTLOOP+1
			cycle Main
		  end if

		end if

		!--------------------------------------------------------------------
		! Check stopping criteria: if done, deallocate and exit SolveMaxwells
		!--------------------------------------------------------------------

		if(ICOUNT < NRELMAX) then

		!toh: modified such that erre is protected against numerical saturation
		!toh: by one digit error of the machine used [8/May/1999]

		  if( dreal(herr) <= errend1 .or. erre <= errend2 ) then
			inquire(ioERR,OPENED=opened)
			if (opened) then
			  write(ioERR,'(i7,2e15.7)')ICOUNT,erre,cdabs(herr)
			end if
			exitflg=0
			exit Main
		  end if

		else if (ICOUNT == NRELMAX) then

		  inquire(ioERR,OPENED=opened)
		  if (opened) then
			write(ioERR,'(i7,2e15.7)')ICOUNT,erre,cdabs(herr)
		  end if
		  exitflg=1
		  exit Main

		end if

		errlast = erre

	  end do BiCGSTAB

	end do Main


!toh: Deallocate dp1 & sp1 explicitly. 25/FEB/99

	deallocate(dp1)
	deallocate(sp1)

!ak: Deallocate all vectors except divr,divi,hx,hy,hz (which are global)

	deallocate(sp,ija,ija1)
    deallocate(dvec,dxvec,rvec,rvec0)
	deallocate(pvec,qvec,tvec,svec)

	deallocate(sppot,sppot1,ijapot,pot)
	deallocate(rpot,spot,ppot,dpot,zpot,qpot,wk1)

	return

  end subroutine SolveMaxwells	! SolveMaxwells



end module maxwells	! maxwells


!--------------------------------------------------------------------
!     Based on subroutine sph3d (A.Schultz 17 February 1999)
!
!        -----> no DivH check at each step to enhance the speed
!
!     3d H-field solution by relaxation with assigned
!     3-component h at top & tangential h at bottom nodes
!
!     geometry of diff equation is edge h and surface e(normals)
!     min rTr for symmetric matrices using preconditioning of m(-1)
!     preconditioning by incomplete cholesky decomposition of matrix
!     a correction is made to ensure divH=0 - this is done by incomplete
!          cholesky decomposition
!     real space values are used - no renormalization is done.
!
!     several air layers are added on top of the earth model. The
!          H fields are set to be uniform on the top of these layers.
!
!     modified to read boundary value file that may have TE mode
!     values for the Hx polarization.
!
!     15-12-97 last modified by M.Uyeshima
!     27-12-97 computing surface EM-field
!     30-12-97 air resistivity can be given
!     10-01-98 modified for divH correction (timing and number)
!     13-01-98 modified for diVH correction (position)
!              ---> pot=0 at the second-top shell
!     21-01-98 modified for Martinec-Period......
!      4-08-98 modified by A. Schultz for f90
!     12-12-98 modified by H. Toh for real-complex separation of A
!      4-02-99 modified by H. Toh into a BiCG algorithm
!     10-02-99 modified by H. Toh for dynamic allocation (sp1 & dp1)
!     17-02-99 modified by A. Schultz hx,hy,hz made static & relaxation changed
!     27-02-99 modified by H. Toh introduction of a weight function for a
!                                 stopping criterion
!     17-03-99 modified by A. Schultz hx,hy,hz were saved instead of being made
!                                 static
!     17-05-99 modified by H. Toh construct a scalar REFERENCE code.
!                                 1. with sorting subroutine [sortcdat]
!                                 2. Adam's 'atimes'
!                                 3. every dynamic allocation disabled
!                                 4. D-response is added with its average
!                                    transferred to the main routine
!                                 5. with 1D analytical solutions
!                                 6. no weight function on stopping criteria
!     27-05-99 modified by H. Toh for a return of dynamic allocations for
!                                 sp1, dp1 & h since dynamic allocations
!                                 turned out compatible with MPI. Also,
!                                 hyt & hzt in set_initial_p10 was tied up
!                                 with mm.
!     28-05-99 modified by H. Toh for Hx, Hy & Hz outputs in files [response]
!     10-06-99 modified by H. Toh for Jx, Jy & Jz outputs in a file [allj_out]
!     26-06-99 bug fix by H. Toh for the wrong boundary value of Hz at the top.
!              Now all the three components of the magnetic field on both
!              boundaries are correctly assigned. However, this boundary value
!              problem should be solved by assigning the tangential components
!              alone, so this still needs further modification...
!	  28-03-05 modified by Anna Kelbert to clean up. Bug fix: no longer updating
!			   the source and b.c. vector b every iteration of the cycle. Using
!			   omega for frequency internally in all calculations. Also modified
!			   to use in the module structure, so that the grid and the model are
!			   now global variables. Dimensions are now variables initialized from
!			   user-defined files, rather than parameters. All vectors modified
!			   for dynamic memory allocation. Notation and order of variable
!			   updates in BiCGSTAB redone to match the original BiCGSTAB algorithm
!			   (see Table 1 of Toh, Schultz and Uyeshima). Two new vectors
!			   introduced for clarity of notation, obsolete defns removed.
!	  09-11-05 modified by Anna Kelbert. SolveMaxwells now works for any input
!			   vector <y> to produce an output vector <x>. We are using the fact
!			   that by construction of operator A, <x> = lH, to perform the
!			   divergence correction. We have introduced the global vectors
!			   hx,hy,hz,sx,sy,sz (magnetic fields and the source term components).
!			   This is how we pass the information about the starting solution
!			   and about the source term forcing (which in practice has been used
!			   in generating <y>) to the divergence correction routine.
!
!	  Profiling shows 'atimes' is the longest by far routine. This subroutine
!	  is called up to IPOTMAX*(2*(NRELMAX+1)) times.
!--------------------------------------------------------------------

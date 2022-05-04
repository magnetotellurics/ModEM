! *************
! 
! Derived class to define a MT Source with boundary data computed by 1D solutions
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
!**
module SourceCSEM_Dipole1D
    !
	use dipole1d
	!
    use Constants
    use cVector3D_SG
	use rVector3D_SG
    use Grid3D_SG
    use Source
    use ModelOperator
	use ModelParameterCell_SG
    !
    type, extends( Source_t ) :: SourceCSEM_Dipole1D_t
        !
        real( kind=prec )  :: azimuth, dip, moment
        real( kind=prec )  :: location(3)
        !
		type( rVector3D_SG_t ) :: CondAnomaly_h
        contains
            !
            final :: SourceCSEM_Dipole1D_dtor
            !
            procedure, public :: setRHS => setRHSMT_1D
            procedure, public :: setE   => setESourceCSEM_Dipole1D
			procedure, public :: create_Ep_from_Dipole1D
			procedure, public :: set1DModel
            !
    end type SourceCSEM_Dipole1D_T
    !
    interface SourceCSEM_Dipole1D_t
        module procedure SourceCSEM_Dipole1D_ctor
    end interface SourceCSEM_Dipole1D_t
    !
contains
    !
    ! SourceCSEM_Dipole1D constructor
    function SourceCSEM_Dipole1D_ctor( model_operator, model_parameter, period, location, dip, azimuth, moment, E ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in )  :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: model_parameter
		!
        real( kind=prec ), intent( in )  :: period, azimuth, dip, moment, location(3)
        !
        class( cVector_t ), intent( in ), optional      :: E
        !
        type( SourceCSEM_Dipole1D_t ) :: self
        !
        !write( *, * ) "Constructor SourceCSEM_Dipole1D_t"
        !
        call self%init()
        !
        self%model_operator  => model_operator
        self%model_parameter => model_parameter
        !
		self%period = period
		self%location = location
		self%dip = dip
		self%azimuth = azimuth
		self%moment = moment
		!
        if ( present( E ) ) then
             !
             allocate( self%E, source = E )
             !
             call self%setRHS()
             !
        endif
        !
    end function SourceCSEM_Dipole1D_ctor
    !
    ! SourceCSEM_Dipole1D destructor
    subroutine SourceCSEM_Dipole1D_dtor( self )
        implicit none
        !
        type( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor SourceCSEM_Dipole1D_t"
        !
        call self%dealloc()
        !
    end subroutine SourceCSEM_Dipole1D_dtor
    !
    ! Set self%E from forward modelling 1D
    subroutine setESourceCSEM_Dipole1D( self, polarization )
		implicit none
		!
		class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
		integer, intent( in )                   :: polarization
		!
		complex(kind=prec)	    :: i_omega_mu 
		real(kind=prec)	    :: omega 
		!
		! Get the Transmitter setting:
		xTx1D = self%location(1)
		yTx1D = self%location(2)
		zTx1D = self%location(3)
		write (*,*) xTx1D,yTx1D,zTx1D,self%period,self%azimuth,self%moment,self%dip
		ftx1D = 1.0d0/self%period
		sdm1D = self%moment           ! (Am), dipole moment. Normalize to unit source moment
		azimuthTx1D = self%azimuth ! (degrees) 
		dipTx1D     = self%dip
		!
		HTmethod1D      = "kk_ht_201"    ! Use 201 point HT digital filters.
		outputdomain1D  = "spatial"      ! Assume spatial domain comps
		lbcomp          = .false.        ! This is changed to true if magnetics in data file
		lUseSpline1D    = .true.         ! Use spline interpolation for faster 1D computations
		linversion      = .false.        ! Compute derivatives with respect to self%model_parameter(layers)

		phaseConvention = "lag"          ! The usual default is lag, where phase becomes larger 
		!    positive values with increasing range.
		lenTx1D         = 00.d0        ! (m) Dipole length 0 = point dipole
		numIntegPts     = 0             ! Number of points to use for Gauss quadrature integration for finite dipole

		call self%set1DModel( xTx1D, yTx1D )
		write( *, * ) "FINISHES set1DModel"
		! Put the privte sig1D and zlay1D into public sig1D and zlay1D required in Dipole1D
		!nlay1D=nlay1D
		! if(allocated(zlay1D)) then
		! Deallocate(zlay1D, sig1D)
		! end if
		! allocate(zlay1D(nlay1D),sig1D(nlay1D))
		! zlay1D=zlay1D
		! sig1D=sig1D

		!call setAnomConductivity(self%model_parameter)
		call initilize_1d_vectors(self%model_parameter%grid)     ! Initilaize the 1D vectors where to compupte the E field
		write( *, * ) "FINISHES initilize_1d_vectors"
		call comp_dipole1D                  ! Calculate E-Field by Key"s code
		write( *, * ) "FINISHES comp_dipole1D"
		call self%create_Ep_from_Dipole1D(self%model_parameter%grid)
		write( *, * ) "FINISHES create_Ep_from_Dipole1D:", allocated( self%E ), self%E%is_allocated, self%CondAnomaly_h%is_allocated
		!
        omega = 2.0 * PI / self%period
        write( *, * ) "FINISHES omega"
		!
		i_omega_mu = cmplx( 0., real( -1.0d0*ISIGN*MU_0*omega ), kind=prec)
		write( *, * ) "FINISHES i_omega_mu"
		!
		!
		if( allocated( self%rhs ) ) deallocate( self%rhs )
		!
            select type( grid => self%model_parameter%grid )
                class is( Grid3D_SG_t )
                    !
                    allocate( self%rhs, source = cVector3D_SG_t( grid, EDGE ) )
					!
					self%rhs = self%E * self%CondAnomaly_h
					write( *, * ) "FINISHES self%E * self%CondAnomaly_h"
					!
					!call diagMult(CondAnomaly_h,self%E,source)
					self%rhs = self%rhs * i_omega_mu
					write( *, * ) "FINISHES self%rhs * i_omega_mu"
					!
					!call scMult(i_omega_mu,source,source)   


					!clean
					!call deall_rvector(CondAnomaly_h)
					write( *, * ) "FINISHES setESourceCSEM_Dipole1D"
					!
                    !
            end select
            !
		!

    end subroutine setESourceCSEM_Dipole1D
    !
!#############################################
subroutine initilize_1d_vectors(grid)
 class(Grid_t), intent(in)        :: grid 
!Local
integer counter,ix,iy,iz


	  n1D = (grid%Nx)*(grid%Ny+1)*(grid%Nz+1)
	  n1D = n1D + (grid%Nx+1)*(grid%Ny)*(grid%Nz+1)
	  n1D = n1D + (grid%Nx+1)*(grid%Ny+1)*(grid%Nz)

	  if (allocated (x1D)) then  
		 Deallocate(x1D, y1D, z1D)
		 Deallocate(ex1D,ey1D,jz1D)
		 Deallocate(bx1D,by1D,bz1D)
	  end if
	  
	  allocate (x1D(n1D), y1D(n1D), z1D(n1D))
	  allocate (ex1D(n1D),ey1D(n1D),jz1D(n1D))
	  allocate (bx1D(n1D),by1D(n1D),bz1D(n1D))
	
	 
	  
	
	  !====================================================================
	  ! Create position vector that the primary field has to be calculated
	  !====================================================================
	  counter = 1
	  ! E-field corresponing to these nodes is Ex
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny+1 !Edge Y
			  Do ix = 1,grid%Nx !Center X
				x1D(counter) = grid%xCenter(ix)
				y1D(counter) = grid%yEdge(iy)
				z1D(counter) = grid%zEdge(iz)
				counter = counter + 1
				End Do
		  End Do
	  End Do
	  
	  ! E-field corresponing to these nodes is Ey
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny !Center y
			  Do ix = 1,grid%Nx+1 !Edge x
				  x1D(counter) = grid%xEdge(ix)
				  y1D(counter) = grid%yCenter(iy)
				  z1D(counter) = grid%zEdge(iz)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
	  
	  ! E-field corresponing to these nodes is Ez
	  Do iz = 1,grid%Nz !Center Z
		  Do iy = 1,grid%Ny+1 !Edge y
			  Do ix = 1,grid%Nx+1 !Edge x
				  x1D(counter) = grid%xEdge(ix)
				  y1D(counter) = grid%yEdge(iy)
				  z1D(counter) = grid%zCenter(iz)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
end subroutine 	initilize_1d_vectors 

subroutine create_Ep_from_Dipole1D( self, grid )

!
		class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
		!
class(Grid_t), intent(in)        :: grid  
!Local 
integer ix,iy,iz,counter
!
            select type( grid  )
                class is( Grid3D_SG_t )
                    !
                    if( allocated( self%E ) ) deallocate( self%E )
      allocate( self%E, source = cVector3D_SG_t( grid, EDGE ) )
	  !
	select type( E => self%E )
		class is( cVector3D_SG_t )
			!
			!
 	  counter = 1
	  ! E-field corresponing to these nodes is Ex
	   Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny+1 !Edge Y
			  Do ix = 1,grid%Nx !Center X	  		  	  
				  E%x(ix,iy,iz) = ex1D(counter)				  
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
		  
	  ! E-field corresponing to these nodes is Ey
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny !Center y
			  Do ix = 1,grid%Nx+1 !Edge x	  				  
				  E%y(ix,iy,iz) = ey1D(counter)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
	  
	  ! E-field corresponing to these nodes is Ez
	  Do iz = 1,grid%Nz !Center Z
		  Do iy = 1,grid%Ny+1 !Edge y
			  Do ix = 1,grid%Nx+1 !Edge x
				  E%z(ix,iy,iz) = jz1D(counter)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
	  
		 Deallocate(x1D, y1D, z1D)
		 Deallocate(ex1D,ey1D,jz1D)
		 Deallocate(bx1D,by1D,bz1D)
		 
			!
	end select
	

                    !
            end select
            !
    end subroutine create_Ep_from_Dipole1D

    ! Set RHS from self%E
    subroutine setRHSMT_1D( self )
        implicit none
        !
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
        !
        select type( E => self%E )
            class is( cVector3D_SG_t )
                !
                if( allocated( self%rhs ) ) deallocate( self%rhs )
                allocate( self%rhs, source = cVector3D_SG_t( E%grid, EDGE ) )
                !
                call self%model_operator%MultAib( self%E%Boundary(), self%rhs )
                !
                self%rhs = self%rhs * C_MinusOne
                !
        end select
        !
    end subroutine setRHSMT_1D
	!
	subroutine set1DModel( self, xTx1D, yTx1D )
		!
        class( SourceCSEM_Dipole1D_t ), intent( inout ) :: self
        !
		real(kind=prec),intent(in)                  :: xTx1D, yTx1D 

		!   local variables ... this is an easy, but not necessarily most efficient
		!   way to get an average background layered conductivity ...
		!    could add routines to modelParameter module to do this more directly
		type( rScalar3D_SG_t )   ::	 sigmaCell 
		character(len=80)      :: paramtype
		type( rScalar3D_SG_t ) :: model
		type(ModelParameterCell_SG_t) :: aModel, Anomalous_model
		type( rVector3D_SG_t )        ::  cond, condNomaly_h


		integer	:: nzEarth,Nz,nzAir,i,j,k,ixTx,iyTx,izTx,counter
		real(kind=prec)	:: wt,vAir,asigma,temp_sigma_value
		character(len=256)   ::       PrimaryFile
		!character(len=20)    ::       get_1D_from



		!   first define conductivity on cells  
		!   (extract into variable which is public)
		!call modelParamToCell(model_parameter, sigmaCell, paramtype)
		!
		select type( model_parameter => self%model_parameter )
			class is( ModelParameterCell_SG_t )
				!
				!
		sigmaCell = model_parameter%cellCond
		nlay1D = sigmaCell%nz+sigmaCell%grid%nzAir
		nzEarth = sigmaCell%grid%nzEarth
		nzAir = sigmaCell%grid%nzAir

		ixTx= minNode(xTx1D, sigmaCell%grid%xEdge)  
		iyTx= minNode(yTx1D, sigmaCell%grid%yEdge)
		izTx= minNode(zTx1D, sigmaCell%grid%zEdge)  		 
		!get_1D_from= Trim(solverControl%get_1D_from)
		!   for layer boundaries use z-edges of 3D grid
		if(allocated(zlay1D)) then
		Deallocate(zlay1D, sig1D)
		end if

		allocate(zlay1D(nlay1D))
		allocate(sig1D(nlay1D))
		do k=1,nlay1D
		zlay1D(k) = sigmaCell%grid%zEdge(k)
		end do

		! For create sig1D, we divide this process into two parts (1) for air layers and 
		!    (2) for earth layers
		! For air layer, sig1D equal to air layer conductivity
		! For earth layer, The Geometric mean is be used to create sig1D

		sig1D(1:nzAir) = SIGMA_AIR !sigmaCell%v(1,1,1:nzAir)

		if (trim(get_1D_from) =="Geometric_mean") then
		   do k = nzAir+1,nlay1D
			wt = R_ZERO
			temp_sigma_value=R_ZERO
				do i = 1,sigmaCell%grid%Nx
					do j = 1,sigmaCell%grid%Ny
						  
							wt = wt + sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
							temp_sigma_value = temp_sigma_value + (sigmaCell%v(i,j,k-nzAir))* &
							sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
						
							
					end do
				end do
				   sig1D(k) = exp(temp_sigma_value/wt)
			   write(220,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
		   end do
		elseif (trim(get_1D_from) =="At_Tx_Position") then
		   do k = nzAir+1,nlay1D
			   sig1D(k)=sigmaCell%v(ixTx,iyTx,k)
			   !write(230,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
		   end do						
		elseif (trim(get_1d_from)=="Geometric_mean_around_Tx") then
			do k = nzAir+1,nlay1D
			  wt = R_ZERO
				do i = ixTx-5,ixTx+5
					do j = iyTx-5,iyTx+5
						if (log(sigmaCell%v(i,j,k)) .gt. -20.0) then
							wt = wt + sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
							sig1D(k) = sig1D(k) + log(sigmaCell%v(i,j,k))* &
							sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
						 end if
					end do
				end do           
				sig1D(k) = exp(sig1D(k)/wt)	
			  !write(240,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
			end do
		elseif (trim(get_1d_from)=="Full_Geometric_mean") then
			wt = R_ZERO
			temp_sigma_value=R_ZERO
			counter=0
			do k = nzAir+1,nlay1D
				do i = 1,sigmaCell%grid%Nx
					do j = 1,sigmaCell%grid%Ny
						if (log(sigmaCell%v(i,j,k)) .gt. -20.0  ) then
							counter=counter+1
							wt = wt + sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)*sigmaCell%grid%dz(k)
							temp_sigma_value = temp_sigma_value + log(sigmaCell%v(i,j,k))
						end if
					end do
				end do           
			end do
			do k = nzAir+1,nlay1D
			  sig1D(k) = exp(temp_sigma_value/counter)	
			  !write(250,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
			end do	
		elseif (trim(get_1d_from)=="Fixed_Value") then
			temp_sigma_value=sigmaCell%v(ixTx,iyTx,nzAir+1) !the value exactly below the Tx
			 do k = nzAir+1,nlay1D
			  sig1D(k) = temp_sigma_value
			  !write(260,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
			end do					
		end if
		!
		
do k = 1,nlay1D
write (*,*) k, sig1D(k),log(sig1D(k)),sigmaCell%v(1,1,k),zlay1D(k),sigmaCell%nz
end do
		model = sigmaCell
		!call getValue_modelParam(model_parameter,paramType,model,vAir)
		!
		! Put the background (Primary) "condNomaly" conductivities in ModEM model format
		model%v=R_ZERO
		do k = nzAir+1,nlay1D
			asigma = sig1D(k)
			
			if( trim(model_parameter%ParamType) == LOGE) asigma = log(asigma)
			!write (*,*) k, nzAir,k-nzAir,nzEarth,nlay1D, trim(model_parameter%ParamType), asigma!,log(asigma)
			do i = 1,sigmaCell%grid%Nx
				do j = 1,sigmaCell%grid%Ny	    
					model%v(i,j,k-nzAir) = (asigma)
				end do
			end do
		end do   

		!call copy_modelParam(amodel,model_parameter)  
		amodel = model_parameter
		!call setType_modelParam(amodel,paramType)
		!call setValue_modelParam(amodel,paramType,model,vAir)
		amodel%cellCond=model
		!call ModelParamToEdge(amodel,condNomaly_h)
		!condNomaly_h = amodel%PDEmapping()
		condNomaly_h = model_parameter%dPDEmapping( amodel )
		cond = model_parameter%PDEmapping()
		self%CondAnomaly_h=cond-condNomaly_h
				!
		end select
		

	end subroutine set1DModel
	!

  ! **************************************************************************
  function minNode(x, xNode) result(ix)
    !  This is a utility routine, used by several data functional
    !  set up routines, and for other interpolation functions
    !  Returns index ix such that  xNode(ix) <= x < xNode(ix+1)
    !  If x is out of range:
    !  x < xNode(1) returns 0; if x> xNode(nx) returns nx
    !  Assumes xNode is strictly increasing; does not check this
    !  NOTE: as presently coded, when xNode is called with center
    !  (face) node positions, this routine will return zero for
    !  the coordinates in the outer half cell nearest the boundary
    !  If evaluation over the complete model domain is to be allowed
    !  a more general interpolation rule will be required.
    !  A.K.: modified to allow input of any size, nx = size(xNode).

    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec), dimension(:), intent(in)     :: xNode

    integer                                     :: ix
    integer                                     :: i

    ix = size(xNode)
    do i = 1,size(xNode)
       if(clean(xNode(i)) .gt. clean(x)) then
          ix = i-1
          exit
       endif
    enddo

  end function minNode


  ! **************************************************************************
  function maxNode(x, xNode) result(ix)
    !  This is a utility routine, used by several data functional
    !  set up routines, and for other interpolation functions
    !  Returns index ix such that  xNode(ix) <= x < xNode(ix+1)
    !  If x is out of range:
    !  x > xNode(1) returns 0; if x< xNode(nx) returns nx
    !  Assumes xNode is strictly decreasing; does not check this
    !  NOTE: as presently coded, when xNode is called with center
    !  (face) node positions, this routine will return zero for
    !  the coordinates in the outer half cell nearest the boundary
    !  If evaluation over the complete model domain is to be allowed
    !  a more general interpolation rule will be required.
    !  A.K.: modified to allow input of any size, nx = size(xNode).

    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec), dimension(:), intent(in)     :: xNode

    integer                                     :: ix
    integer                                     :: i

    ix = size(xNode)
    do i = 1,size(xNode)
       if(clean(xNode(i)) .lt. clean(x)) then
          ix = i-1
          exit
       endif
    enddo

  end function maxNode
   ! **************************************************************************
  function clean(x)
    ! This is a utility routine that provides an expression used to battle
	! against machine error problems. It returns the same real or real(8)
	! as the input, but without the extra digits at the end that are often
	! a cause of wrong comparisons in the if statements. ALWAYS use clean(x)
	! instead of x in an inequality!!!
	! LARGE_REAL is defined in the module math_constants
	! A.K.
    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec)						       :: clean

	clean = dnint(x*LARGE_REAL)/LARGE_REAL

  end function clean
  
end module SourceCSEM_Dipole1D

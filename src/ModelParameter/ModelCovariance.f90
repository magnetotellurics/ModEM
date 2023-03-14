!
!> 3D_MT model covariance by Gary Egbert and Anna Kelbert: procedures.
!>
!> Recursive autoregression
!>
!
module ModelCovarianceRec
    !
    use Constants
    use ModelParameterCell_SG
    !
    type :: ModelCovarianceRec_t
        !
        !> Dimensions of the grid
        integer :: Nx, Ny, NzEarth
        !
        !> Number of times the smoothing operator should be applied
        integer :: N
        !
        !> General rules for smoothing in the X and Y-directions, dimension(NzEarth)
        real( kind=prec ), pointer, dimension(:) :: Sx, Sy
        !
        !> General rule for vertical smoothing
        real( kind=prec ) :: Sz
        !
        !> Special rules for smoothing across a surface stored as a sparse array
        !type( sparsevecc ) :: S
        !
        !> Integer array that defines regions for smoothing and scaling purposes
        !type( iscalar ) :: mask
        !
        !> True when all arrays are is_allocated, initialized and ready to use
        logical :: is_allocated
        !
        contains
            !
            final :: ModelCovarianceRec_dtor
            !
            procedure, public :: multBy_Cm
            procedure, public :: multBy_CmSqrt
            procedure, public :: multBy_CmSqrtInv
            procedure, public :: read_CmSqrt
            procedure, public :: RecursiveAR
            procedure, public :: RecursiveARInv
            procedure, public :: SmoothX
            procedure, public :: SmoothY
            procedure, public :: SmoothZ
            procedure, public :: Scaling
            !
    end type
    !
    interface ModelCovarianceRec_t
         module procedure ModelCovarianceRec_ctor
    end interface ModelCovarianceRec_t
    !
contains
    !
    !> Initializes self variable stored in RecursiveAR.hd. If cfile
    !> is specified, gets this information from file.
    function ModelCovarianceRec_ctor( m, cfile ) result( self )
        implicit none
        !
        class( ModelParameter_t ), allocatable, intent( in ) :: m
        character(*), intent( in ), optional :: cfile
        !
        type( ModelCovarianceRec_t ) :: self
        !
        integer :: istat
        logical :: exists
        !
        !> Initializing self
        self%Nx = m%metric%grid%Nx
        self%Ny = m%metric%grid%Ny
        self%NzEarth = m%metric%grid%NzEarth
        allocate( self%Sx(self%NzEarth), STAT=istat )
        allocate( self%Sy(self%NzEarth), STAT=istat )
        self%Sx = 0.3
        self%Sy = 0.3
        self%Sz = 0.3
        !call create_iscalar( m%grid, self%mask, CELL_EARTH )
        !self%mask%v = FREE
        self%N = 1
        self%is_allocated = .TRUE.
        !
        if( present( cfile ) ) then
            !
            !> Attempt to read self from cfile
            inquire( FILE=cfile, EXIST=exists )
            !
            if( exists ) then
                call self%read_CmSqrt( cfile )
            else
                !
                write( *, * ) "Error: Unable to find the input covariance file [", trim( cfile ), "]!"
                stop
                !
            endif
            !
            if( ( self%Nx /= m%metric%grid%Nx ) .OR. ( self%Ny /= m%metric%grid%Ny ) .OR. ( self%NzEarth /= m%metric%grid%NzEarth ) ) then
                !
                write( *, * ) "Error: Grid dimensions do not match in input covariance file [", trim( cfile ), "]!"
                stop
                !
            endif
            !
        endif
        !
    end function ModelCovarianceRec_ctor
    !
    !> No subroutine briefing
    subroutine ModelCovarianceRec_dtor( self )
        implicit none
        !
        type( ModelCovarianceRec_t ), intent( inout ) :: self
        !
        integer :: istat
        !
        deallocate( self%Sx, self%Sy, STAT=istat )
        !call deall_sparsevecc( self%S )
        !call deall_iscalar( self%mask )
        self%is_allocated = .FALSE.
        !
    end subroutine ModelCovarianceRec_dtor
    !
    !> Multiplies by the full model covariance,
    !> which is viewed as a smoothing operator. Intended
    !> to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
    !> For efficiency, self is a saved, private variable inside
    !> the modelParam module. Before this routine can be called,
    !> it has to be initialized by calling create_CmSqrt(m).
    !
    subroutine multBy_Cm( self, target_model )
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        class( ModelParameter_t ), allocatable, intent( inout ) :: target_model
        !
        class( ModelParameter_t ), allocatable :: temp_model
        !
        allocate( temp_model, source = target_model )
        !
        select type( temp_model )
            !
            class is( ModelParameterCell_SG_t )
                !
                select type( target_model )
                    !
                    class is( ModelParameterCell_SG_t )
                        !
                        call self%RecursiveAR( temp_model%cell_cond%v, target_model%cell_cond%v, 2 )
                        !
                    class default
                        stop "Error: multBy_Cm > Unclassified target_model"
                    !
                end select
                !
            class default
                stop "Error: multBy_Cm > Unclassified temp_model"
            !
        end select
        !
        deallocate( temp_model )
        !
    end subroutine multBy_Cm
    !
    !> Multiplies by the square root of the model covariance,
    !> which is viewed as a smoothing operator. Intended
    !> to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
    !> For efficiency, self is a saved, private variable inside
    !> the modelParam module. Before this routine can be called,
    !> it has to be initialized by calling create_CmSqrt(m).
    !
    subroutine multBy_CmSqrt( self, mhat, dsigma )
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        class( ModelParameter_t ), allocatable, intent( in ) :: mhat
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        if( allocated( dsigma ) ) deallocate( dsigma )
        allocate( dsigma, source = mhat )
        !
        select type( mhat )
            !
            class is( ModelParameterCell_SG_t )
                !
                select type( dsigma )
                    !
                    class is( ModelParameterCell_SG_t )
                        !
                        call self%RecursiveAR( mhat%cell_cond%v, dsigma%cell_cond%v, self%N )
                        !
                    class default
                        stop "Error: multBy_CmSqrt > Unclassified dsigma"
                    !
                end select
                !
            class default
                stop "Error: multBy_CmSqrt > Unclassified mhat"
            !
        end select
        !
    end subroutine multBy_CmSqrt
    !
    !> Multiplies by the inverse square root of the model covariance,
    !> which is viewed as a roughening operator. Intended
    !> to be used to compute \tilde{m} = C_m^{-1/2} ( m - m_0 ).
    !> For efficiency, self is a saved, private variable inside
    !> the modelParam module. Before this routine can be called,
    !> it has to be initialized by calling create_CmSqrt(m).
    !
    function multBy_CmSqrtInv( self, dm ) result ( mhat )
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        class( ModelParameter_t ), allocatable, intent( in ) :: dm
        class( ModelParameter_t ), allocatable :: mhat
        !
        mhat = dm
        !
        select type( mhat )
            !
            class is( ModelParameterCell_SG_t )
                !
                select type( dm )
                    !
                    class is( ModelParameterCell_SG_t )
                        !
                        call self%RecursiveARInv( dm%cell_cond%v, mhat%cell_cond%v, self%N )
                        !
                    class default
                        stop "Error: multBy_CmSqrtInv > Unclassified dm"
                    !
                end select
                !
            class default
                stop "Error: multBy_CmSqrtInv > Unclassified mhat"
            !
        end select
        !
    end function multBy_CmSqrtInv
    !
    !> ...  Copyright (C) 2008 Anna Kelbert. All rights reserved.
    !> 
    !> +-----------------------------------------------------------------------------+
    !> | This file defines model covariance for a recursive autoregression scheme.   |
    !> | The model space may be divided into distinct areas using integer masks.     |
    !> | Mask 0 is reserved for air; mask 1 is reserved for ocean. Smoothing between |
    !> | air, ocean and the rest of the model is turned off automatically. You can   |
    !> | also define exceptions to override smoothing between any two model areas.   |
    !> | To turn off smoothing set it to zero. This header is 16 lines long.         |
    !> | 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |
    !> | 2. Smoothing in the X direction (NzEarth real values)                       |
    !> | 3. Smoothing in the Y direction (NzEarth real values)                       |
    !> | 4. Vertical smoothing (1 real value)                                        |
    !> | 5. Number of times the smoothing should be applied (1 integer >= 0)         |
    !> | 6. Number of exceptions (1 integer >= 0)                                    |
    !> | 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 3 & 4) |
    !> | 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|
    !> +-----------------------------------------------------------------------------+
    !
    subroutine read_CmSqrt( self, cfile )
        implicit none
        !
        !> The minimal covariance information includes the AR parameters
        !> alpha(k), beta(k) for smoothing in x, y directions and gamma for
        !> the vertical smoothing. Both alpha and beta could depend on the
        !> vertical layer. The scaling is the identity when not specified.
        !> This information is read from a file. Also, we read an integer mask
        !> array that subdivides the model grid into different regions
        !> (AIR, OCEAN, EARTH) and a set of rules that overrides the default
        !> smoothing parameters across a particular surface between two
        !> distinct regions. We use this to set up the covariance self.
        !
        !> Strictly speaking, to define the smoothing across surfaces in
        !> full generality while maintaining efficiency, it has to be a sparse
        !> real vector defined on FACES (sparsevecr). We only have a complex
        !> sparse vector implemented (sparsevecc). We could either use that,
        !> or imitate the structure.
        class( ModelCovarianceRec_t ), intent( inout ) :: self
        character(*), intent( in ) :: cfile
        !
        !> Exception rules
        integer, pointer, dimension(:) :: mask1, mask2, ii, jj, kk, xyz
        real( kind=prec ), pointer, dimension(:) :: smoothing, S
        !
        integer :: Nx, Ny, NzEarth, nrules, nS, i, j, k, n, istat
        integer :: fid = 30
        !
        if( .NOT. self%is_allocated ) then
            stop "Error: Model covariance must be is_allocated before reading from file in read_CmSqrt"
        endif
        !
        open( unit=fid, file=cfile, form="formatted", status="old" )
        !
        !> skip the 16 lines header
        do j = 1,16
            read( fid, * )
        enddo
        !
        !> read grid dimensions
        read( fid, * ) Nx,Ny,NzEarth
        self%Nx = Nx
        self%Ny = Ny
        self%NzEarth = NzEarth
        !
        !> read smoothing parameters
        read( fid, * ) self%Sx
        read( fid, * ) self%Sy
        read( fid, * ) self%Sz
        !
        !> read number of times to apply the smoothing
        read( fid, * ) self%N
        !
        !> read exception rules for smoothing across surfaces
        read( fid, * ) nrules
        allocate( mask1( nrules ), mask2( nrules ), smoothing( nrules ),STAT=istat)
        !
        !>
        do n = 1, nrules
            read( fid, * ) mask1(n), mask2(n), smoothing(n)
        enddo
        !
        !> create and read the mask array
        !call read_iscalar( fid, self%mask )
        !
        close( fid )
        !
        !> create a huge sparse vector to make sure we accommodate all smoothing exceptions
        !call create_sparsevecc(Nx*Ny*NzEarth, self%S, FACE)
        ! !
        ! !> now, parse the exceptions
        ! nS = 0
        ! do k = 2,NzEarth
            ! do j = 2,Ny
                ! do i = 2,Nx
                    ! do n = 1,nrules
                        ! !
                        ! !> look back in the X-direction
                        ! if(((self%mask%v(i-1,j,k) == mask1(n)) .AND. (self%mask%v(i,j,k) == mask2(n))) &
                        ! .OR. ((self%mask%v(i-1,j,k) == mask2(n)) .AND. (self%mask%v(i,j,k) == mask1(n)))) &
                        ! then
                            ! nS = nS+1
                            ! self%S%i(nS) = i-1
                            ! self%S%j(nS) = j
                            ! self%S%k(nS) = k
                            ! self%S%xyz(nS) = 1
                            ! self%S%c(nS) = smoothing(n)
                        ! endif
                        ! !
                        ! !> look back in the Y-direction
                        ! if(((self%mask%v(i,j-1,k) == mask1(n)) .AND. (self%mask%v(i,j,k) == mask2(n))) &
                        ! .OR. ((self%mask%v(i,j-1,k) == mask2(n)) .AND. (self%mask%v(i,j,k) == mask1(n)))) &
                        ! then
                            ! nS = nS+1
                            ! self%S%i(nS) = i
                            ! self%S%j(nS) = j-1
                            ! self%S%k(nS) = k
                            ! self%S%xyz(nS) = 2
                            ! self%S%c(nS) = smoothing(n)
                        ! endif
                        ! !> look back in the Z-direction
                        ! if(((self%mask%v(i,j,k-1) == mask1(n)) .AND. (self%mask%v(i,j,k) == mask2(n))) &
                        ! .OR. ((self%mask%v(i,j,k-1) == mask2(n)) .AND. (self%mask%v(i,j,k) == mask1(n)))) &
                        ! then
                            ! nS = nS+1
                            ! self%S%i(nS) = i
                            ! self%S%j(nS) = j
                            ! self%S%k(nS) = k-1
                            ! self%S%xyz(nS) = 3
                            ! self%S%c(nS) = smoothing(n)
                        ! endif
                    ! enddo
                ! enddo
            ! enddo
        ! enddo
        ! !
        deallocate( mask1, mask2, smoothing )
        !
        !> now, truncate the smoothing vector to the correct number of components
        !call reall_sparsevecc(nS, self%S)
        !
    end subroutine read_CmSqrt
    !
    !> Implements the recursive autoregression algorithm for a 3D real array.
    !> In our case, the assumed-shape array would be e.g. conductivity
    !> in each cell of the Nx x Ny x NzEarth grid.
    !
    subroutine RecursiveAR( self, w, v, n )
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        real( kind=prec ), intent( in ) :: w(:,:,:)
        real( kind=prec ), intent( out ) :: v(:,:,:)
        integer, intent( in ) :: n
        integer :: Nx, Ny, NzEarth, i, j, k, iSmooth
        !
        Nx = size( w, 1 )
        Ny = size( w, 2 )
        NzEarth = size( w, 3 )
        !
        if( maxval( abs(shape(w) - shape(v) ) ) > 0 ) then
            !
            stop "Error: The input arrays should be of the same shapes in RecursiveAR!"
            !
        endif
        !
        v = w
        !
        do iSmooth = 1, n
            !
            !> smooth in the X-direction (Sx)
            do k = 1,NzEarth
                do j = 1,Ny
                    !v(1,j,k) = v(1,j,k)
                    do i = 2,Nx
                        v(i,j,k) = self%SmoothX(i-1,j,k) * v(i-1,j,k) + v(i,j,k)
                    enddo
                enddo
            enddo
            !
            !> smooth in the Y-direction (Sy)
            do k = 1,NzEarth
                do i = 1,Nx
                    !> v(i,1,k) = v(i,1,k)
                    do j = 2,Ny
                        v(i,j,k) = self%SmoothY(i,j-1,k) * v(i,j-1,k) + v(i,j,k)
                    enddo
                enddo
            enddo
            !
            !> smooth in the Z-direction (Sz)
            do j = 1,Ny
                do i = 1,Nx
                    !> v(i,j,1) = v(i,j,1)
                    do k = 2,NzEarth
                        v(i,j,k) = self%SmoothZ(i,j,k-1) * v(i,j,k-1) + v(i,j,k)
                    enddo
                enddo
            enddo
            !
            !> smooth in the Z-direction (Sz^T)
            do j = Ny,1,-1
                do i = Nx,1,-1
                    !> v(i,j,NzEarth) = v(i,j,NzEarth)
                    do k = NzEarth,2,-1
                        v(i,j,k-1) = v(i,j,k-1) + self%SmoothZ(i,j,k-1) * v(i,j,k)
                    enddo
                enddo
            enddo
            !
            !> smooth in the Y-direction (Sy^T)
            do k = NzEarth,1,-1
                do i = Nx,1,-1
                    !> v(i,Ny,k) = v(i,Ny,k)
                    do j = Ny,2,-1
                        v(i,j-1,k) = v(i,j-1,k) + self%SmoothY(i,j-1,k) * v(i,j,k)
                    enddo
                enddo
            enddo
            !
            !> smooth in the X-direction (Sx^T)
            do k = NzEarth,1,-1
                do j = Ny,1,-1
                    !> v(Nx,j,k) = v(Nx,j,k)
                    do i = Nx,2,-1
                        v(i-1,j,k) = v(i-1,j,k) + self%SmoothX(i-1,j,k) * v(i,j,k)
                    enddo
                enddo
            enddo
            !
        enddo
        !
        !> apply the scaling operator C
        do k = 1,NzEarth
            do j = 1,Ny
                do i = 1,Nx
                    v(i,j,k) = (self%Scaling(i,j,k)**n) * v(i,j,k)
                enddo
            enddo
        enddo
        !
    end subroutine RecursiveAR
    !
    !> ... and the inverse "roughening" operator useful for starting
    !> the inversion with an arbitrary model: \tilde{m} = C_m^{-1/2} (m - m_0).
    !> In our case, the assumed-shape array would be e.g. conductivity
    !> in each cell of the Nx x Ny x NzEarth grid.
    !> NOTE: the inverse covariance operator is poorly conditioned!!!
    !> e.g., at alpha=0.3 n=4 white noise completely overwhelmes the
    !> inverse model. Be extra careful when you use this function and
    !> always look at the result before using it to start the inversion.
    !> In the future, may want to stabilize this.
    subroutine RecursiveARInv( self, w, v, n)
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        real( kind=prec ), intent( in ) :: w(:,:,:)
        real( kind=prec ), intent( out ) :: v(:,:,:)
        integer, intent( in ) :: n
        !
        integer :: Nx, Ny, NzEarth, i, j, k, iSmooth, istat
        real( kind=prec ), allocatable :: u(:,:,:)
        !
        Nx = size( w, 1 )
        Ny = size( w, 2 )
        NzEarth = size( w, 3 )
        !
        if( maxval( abs( shape(w) - shape(v) ) ) >0 ) then
            stop "Error: The input arrays should be of the same shapes in RecursiveARInv"
        endif
        !
        allocate( u( Nx, Ny, NzEarth ), stat=istat )
        !
        v = w
        !
        do iSmooth = 1,n
            !
            u = v
            !
            !> invert smoothing in the X-direction (Sx^T)
            do k = NzEarth,1,-1
                do j = Ny,1,-1
                    v(Nx,j,k) = u(Nx,j,k)
                    do i = Nx,2,-1
                        v(i-1,j,k) = u(i-1,j,k) - self%SmoothX(i-1,j,k) * u(i,j,k)
                    enddo
                enddo
            enddo
            !
            u = v
            !
            !> invert smoothing in the Y-direction (Sy^T)
            do k = NzEarth,1,-1
                do i = Nx,1,-1
                    v(i,Ny,k) = u(i,Ny,k)
                    do j = Ny,2,-1
                        v(i,j-1,k) = u(i,j-1,k) - self%SmoothY(i,j-1,k) * u(i,j,k)
                    enddo
                enddo
            enddo
            !
            u = v
            !
            !> invert smoothing in the Z-direction (Sz^T)
            do j = Ny,1,-1
                do i = Nx,1,-1
                    v(i,j,NzEarth) = u(i,j,NzEarth)
                    do k = NzEarth,2,-1
                        v(i,j,k-1) = u(i,j,k-1) - self%SmoothZ(i,j,k-1) * u(i,j,k)
                    enddo
                enddo
            enddo
            !
            u = v
            !
            !> invert smoothing in the Z-direction (Sz)
            do j = 1,Ny
                do i = 1,Nx
                    v(i,j,1) = u(i,j,1)
                    do k = 2,NzEarth
                        v(i,j,k) = - self%SmoothZ(i,j,k-1) * u(i,j,k-1) +u(i,j,k)
                    enddo
                enddo
            enddo
            !
            u = v
            !
            !> invert smoothing in the Y-direction (Sy)
            do k = 1,NzEarth
                do i = 1,Nx
                    v(i,1,k) = u(i,1,k)
                    do j = 2,Ny
                        v(i,j,k) = - self%SmoothY(i,j-1,k) * u(i,j-1,k) + u(i,j,k)
                    enddo
                enddo
            enddo
            !
            u = v
            !
            !> invert smoothing in the X-direction (Sx)
            do k = 1,NzEarth
                do j = 1,Ny
                    v(1,j,k) = u(1,j,k)
                    do i = 2,Nx
                        v(i,j,k) = - self%SmoothX(i-1,j,k) * u(i-1,j,k) + u(i,j,k)
                    enddo
                enddo
            enddo
            !
        enddo
        !
        !> apply the inverse of the scaling operator C
        do k = 1,NzEarth
            do j = 1,Ny
                do i = 1,Nx
                    if( abs( self%Scaling( i, j, k ) ) < R_TINY ) then
                        v( i, j, k ) = 0
                    else
                        v( i, j, k ) =  v( i, j, k ) / ( self%Scaling( i, j, k )**n )
                    endif
                enddo
            enddo
        enddo
        !
        deallocate( u, stat=istat )
        !
    end subroutine RecursiveARInv
    !
    !> computes the smoothing coefficient in the x-direction based on self
    function SmoothX( self, i, j, k ) result( alpha )
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        integer, intent( in ) :: i, j, k
        real( kind=prec ) :: alpha
        !
        integer :: n
        !
        if( .NOT. associated( self%Sx ) ) then
            stop "Error: self%Sx has to be is_allocated before calling SmoothX"
        endif
        !
        alpha = self%Sx(k)
        !
        if((i < 1) .OR. (i > self%Nx)) then
            stop "Error: index i out of bounds in SmoothX(i,j,k)"
        else if((j < 1) .OR. (j > self%Ny)) then
            stop "Error: index j out of bounds in SmoothX(i,j,k)"
        else if((k < 1) .OR. (k > self%NzEarth)) then
            stop "Error: index k out of bounds in SmoothX(i,j,k)"
        endif
        ! !
        ! if(self%S%is_allocated .AND. (self%S%nCoeff > 0)) then
            ! !> scan through the special rules and possibly update the result
            ! do n = 1,self%S%nCoeff
                ! if(self%S%xyz(n) == 1) then
                    ! if((self%S%i(n) == i) .AND. (self%S%j(n) == j) .AND. (self%S%k(n) == k)) then
                        ! alpha = self%S%c(n)
                        ! exit
                    ! endif
                ! endif
            ! enddo
        ! endif
        ! !
    end function SmoothX
    !
    !> Computes the smoothing coefficient in the y-direction based on self
    !
    function SmoothY( self, i, j, k ) result( beta )
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        integer, intent( in ) :: i, j, k
        real( kind=prec ) :: beta
        !
        integer :: n
        !
        if( .NOT. associated( self%Sy ) ) then
            stop "Error: self%Sy has to be is_allocated before calling SmoothY"
        endif
        !
        beta = self%Sy(k)
        !
        if( ( i < 1 ) .OR. ( i > self%Nx ) ) then
            stop "Error: index i out of bounds in SmoothY(i,j,k)"
        else if((j < 1) .OR. (j > self%Ny)) then
            stop "Error: index j out of bounds in SmoothY(i,j,k)"
        else if((k < 1) .OR. (k > self%NzEarth)) then
            stop "Error: index k out of bounds in SmoothY(i,j,k)"
        endif
        ! !
        ! if( self%S%is_allocated .AND. ( self%S%nCoeff > 0 ) ) then
            ! !> scan through the special rules and possibly update the result
            ! do n = 1,self%S%nCoeff
                ! if( self%S%xyz(n) == 2) then
                    ! if( ( self%S%i(n) == i ) .AND. ( self%S%j(n) == j ) .AND. ( self%S%k(n) == k ) ) then
                        ! beta = self%S%c(n)
                        ! exit
                    ! endif
                ! endif
            ! enddo
        ! endif
        ! !
    end function SmoothY
    !
    !> computes the smoothing coefficient in the z-direction based on self
    !
    function SmoothZ( self, i, j, k ) result( gamma )
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        integer, intent( in ) :: i, j, k
        real( kind=prec ) :: gamma
        !
        integer :: n
        !
        gamma = self%Sz
        !
        if((i < 1) .OR. (i > self%Nx)) then
            stop "Error: index i out of bounds in SmoothZ(i,j,k)"
        else if((j < 1) .OR. (j > self%Ny)) then
            stop "Error: index j out of bounds in SmoothZ(i,j,k)"
        else if((k < 1) .OR. (k > self%NzEarth)) then
            stop "Error: index k out of bounds in SmoothZ(i,j,k)"
        endif
        ! !
        ! if( self%S%is_allocated .AND. ( self%S%nCoeff > 0 ) ) then
            ! !> scan through the special rules and possibly update the result
            ! do n = 1,self%S%nCoeff
                ! if( self%S%xyz(n) == 3 ) then
                    ! if( ( self%S%i(n) == i ) .AND. ( self%S%j(n) == j ) .AND. ( self%S%k(n) == k ) ) then
                        ! gamma = self%S%c(n)
                        ! exit
                    ! endif
                ! endif
            ! enddo
        ! endif
        ! !
    end function SmoothZ
    !
    !> computes the scaling coefficient based on self
    !
    function Scaling( self, i, j, k ) result( c )
        implicit none
        !
        class( ModelCovarianceRec_t ), intent( in ) :: self
        integer, intent( in ) :: i,j,k
        real( kind=prec ) :: c, alpha, beta, gamma
        !
        if((i < 1) .OR. (i > self%Nx)) then
            stop "Error: index i out of bounds in Scaling(i,j,k)"
        else if((j < 1) .OR. (j > self%Ny)) then
            stop "Error: index j out of bounds in Scaling(i,j,k)"
        else if((k < 1) .OR. (k > self%NzEarth)) then
            stop "Error: index k out of bounds in Scaling(i,j,k)"
        endif
        !
        alpha = self%SmoothX( i, j, k )
        beta = self%SmoothY( i, j, k )
        gamma = self%SmoothZ( i, j, k )
        !
        !if( self%mask%v( i, j, k ) == AIR ) then
            !c = 0.0
        !else if( self%mask%v(i,j,k) == OCEAN ) then
            !c = 0.0
        !else
            c = (1 - alpha)**2 * (1 - beta)**2 * (1 - gamma)**2
        !endif
        !
    end function Scaling
    !
end Module ModelCovarianceRec

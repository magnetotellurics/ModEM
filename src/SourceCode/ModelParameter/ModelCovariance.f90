!
!> 3D_MT model covariance by Gary Egbert and Anna Kelbert: procedures.
!>
!> Recursive autoregression
!>
!
module ModelCovariance
    !
    use Constants
    use ModelParameterCell_SG
    use cVectorSparse3D_SG
    use rScalar3D_SG
    use iScalar3D_SG
    !
    !> define the mask for air and ocean here. By default, we do not switch off
    !> the smoothing between the air & ocean and the rest of the model. However,
    !> we do set the scaling to zero for both of these regions.
    integer, parameter, private :: AIR   = 0
    integer, parameter, private :: OCEAN = 9
    integer, parameter, private :: FREE  = 1 ! anything 1-8
    !
    type :: ModelCovariance_t
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
        type( cVectorSparse3D_SG_t ) :: S
        !
        !> Integer array that defines regions for smoothing and scaling purposes
        type( iScalar3D_SG_t ) :: mask
        !
        !> True when all arrays are is_allocated, initialized and ready to use
        logical :: is_allocated
        !
        contains
            !
            final :: ModelCovariance_dtor
            !
            procedure, public :: multBy_Cm
            procedure, public :: multBy_CmSqrt
            procedure, public :: multBy_CmSqrtInv
            procedure, public :: read_CmSqrt
            procedure, public :: RecursiveAR
            procedure, public :: RecursiveARInv
            procedure, public :: SmoothX, SmoothY, SmoothZ
            procedure, public :: Scaling
            !
    end type
    !
    !> Public Global ModelCovariance object
    class( ModelCovariance_t ), allocatable :: model_cov
    !
    interface ModelCovariance_t
         module procedure ModelCovariance_ctor
    end interface ModelCovariance_t
    !
contains
    !
    !> Initializes self variable stored in RecursiveAR.hd. If cfile
    !> is specified, gets this information from file.
    function ModelCovariance_ctor( m, cfile ) result( self )
        implicit none
        !
        class( ModelParameter_t ), allocatable, intent( in ) :: m
        character(*), intent( in ), optional :: cfile
        !
        type( ModelCovariance_t ) :: self
        !
        integer :: istat
        logical :: exists
        !
        !> Initializing self
        self%Nx = m%metric%grid%Nx
        self%Ny = m%metric%grid%Ny
        self%NzEarth = m%metric%grid%NzEarth
        !
        allocate( self%Sx(self%NzEarth), STAT=istat )
        allocate( self%Sy(self%NzEarth), STAT=istat )
        !
        self%Sx = 0.3
        self%Sy = 0.3
        self%Sz = 0.3
        !
        self%S%is_allocated = .FALSE.
        !
        self%mask = iScalar3D_SG_t( m%metric%grid, CELL_EARTH )
        self%mask%v = FREE
        !
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
                call errStop( "ModelCovariance_ctor > Unable to find the input covariance file ["//trim( cfile )//"]!" )
                !
            endif
            !
            if( ( self%Nx /= m%metric%grid%Nx ) .OR. ( self%Ny /= m%metric%grid%Ny ) .OR. ( self%NzEarth /= m%metric%grid%NzEarth ) ) then
                !
                call errStop( "ModelCovariance_ctor > Grid dimensions do not match in input covariance file ["//trim( cfile )//"]!" )
                !
            endif
            !
        endif
        !
    end function ModelCovariance_ctor
    !
    !> No subroutine briefing
    subroutine ModelCovariance_dtor( self )
        implicit none
        !
        type( ModelCovariance_t ), intent( inout ) :: self
        !
        integer :: istat
        !
        deallocate( self%Sx, self%Sy, STAT=istat )
        !call deall_sparsevecc( self%S )
        !call deall_iscalar( self%mask )
        self%is_allocated = .FALSE.
        !
    end subroutine ModelCovariance_dtor
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
        class( ModelCovariance_t ), intent( in ) :: self
        class( ModelParameter_t ), allocatable, intent( inout ) :: target_model
        !
        integer :: i
        type( rScalar3D_SG_t ) :: target_cond, temp_cond
        !
        if( .NOT. target_model%is_allocated ) then
            call errStop( "multBy_Cm > target_model not allocated!" )
        endif
        !
        do i = 1, size( target_model%getCond() )
            !
            target_cond = target_model%getCond(i)
            !
            temp_cond = target_cond
            !
            call self%RecursiveAR( temp_cond%v, target_cond%v, 2 )
            !
            call target_model%setCond( target_cond, i )
            !
        enddo
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
        class( ModelCovariance_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in ) :: mhat
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        type( rScalar3D_SG_t ) :: mhat_cond, dsigma_cond
        integer :: i
        !
        if( .NOT. mhat%is_allocated ) then
            call errStop( "multBy_CmSqrt > mhat not allocated!" )
        endif
        !
        if( allocated( dsigma ) ) deallocate( dsigma )
        allocate( dsigma, source = mhat )
        dsigma = mhat
        !
        do i = 1, size( mhat%getCond() )
            !
            mhat_cond = mhat%getCond(i)
            !
            dsigma_cond = dsigma%getCond(i)
            !
            call self%RecursiveAR( mhat_cond%v, dsigma_cond%v, self%N )
            !
            call dsigma%setCond( dsigma_cond, i )
            !
        enddo
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
    function multBy_CmSqrtInv( self, dsigma ) result( mhat )
        implicit none
        !
        class( ModelCovariance_t ), intent( in ) :: self
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        class( ModelParameter_t ), allocatable :: mhat
        !
        type( rScalar3D_SG_t ) :: dsigma_cond, mhat_cond
        integer :: i
        !
        if( allocated( mhat ) ) deallocate( mhat )
        allocate( mhat, source = dsigma )
        mhat = dsigma
        !
        do i = 1, size( dsigma%getCond() )
            !
            dsigma_cond = dsigma%getCond(i)
            !
            mhat_cond = mhat%getCond(i)
            !
            call self%RecursiveARInv( dsigma_cond%v, mhat_cond%v, self%N )
            !
            call mhat%setCond( mhat_cond, i )
            !
        enddo
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
    !
    subroutine read_CmSqrt( self, cfile )
        implicit none
        !
        class( ModelCovariance_t ), intent( inout ) :: self
        character(*), intent( in ) :: cfile
        !
        !> Exception rules
        integer, pointer, dimension(:) :: mask1, mask2, ii, jj, kk, xyz
        real( kind=prec ), pointer, dimension(:) :: smoothing, S, temp
        !
        integer :: k1, k2, Nx, Ny, Nz, NzEarth, nrules, nS, i, j, k, n, ios, istat
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "read_CmSqrt > Model covariance must be is_allocated before reading from file!" )
        endif
        !
        open( unit=ioCovariance, file=cfile, form="formatted", status="old", iostat = ios )
        !
        if( ios == 0 ) then
            !
            !> skip the 16 lines header
            do j = 1,16
                read( ioCovariance, * )
            enddo
            !
            !> read grid dimensions
            read( ioCovariance, * ) Nx, Ny, NzEarth
            self%Nx = Nx
            self%Ny = Ny
            self%NzEarth = NzEarth
            !
            !> read smoothing parameters
            read( ioCovariance, * ) self%Sx
            read( ioCovariance, * ) self%Sy
            read( ioCovariance, * ) self%Sz
            !
            !> read number of times to apply the smoothing
            read( ioCovariance, * ) self%N
            !
            !> read exception rules for smoothing across surfaces
            read( ioCovariance, * ) nrules
            allocate( mask1( nrules ), mask2( nrules ), smoothing( nrules ),STAT=istat)
            !
            !>
            do n = 1, nrules
                read( ioCovariance, * ) mask1(n), mask2(n), smoothing(n)
            enddo
            !
            !> create and read the mask array
            !call self%mask%read( ioCovariance )
            !
            Nx = size(self%mask%v, 1)
            Ny = size(self%mask%v, 2)
            Nz = size(self%mask%v, 3)
            !
            allocate(temp(Ny), STAT = istat)
            !
            i = 1
            do
                 read( ioCovariance, *, iostat = istat ) k1, k2
                 !
                 if( istat /= 0) exit
                 !
                 if( ( k1 < 0 ) .OR. ( k2 > Nz ) ) then
                        write( *, * ) "Error: read_CmSqrt > While reading the ", i, "th block!"
                        stop
                 elseif( k1 > k2) then
                        write( *, * ) "     "//achar(27)//"[91m# Warning:"//achar(27)//"[0m read_CmSqrt > Block ", i, " will be ignored."
                 endif
                 !
                 do j = Nx, 1, -1
                        !
                        read(ioCovariance, *, iostat = istat) temp
                        !
                        if( istat /= 0) then
                             write( *, * ) "Error: read_CmSqrt > While reading the ", j, "th row in ", i,"th block."
                             stop
                        endif
                        !
                        do k = k1, k2
                             self%mask%v(j, :, k) = temp
                        enddo
                 enddo
                 !
                 if( k == Nz) exit
                 !
                 i = i + 1
                 !
            enddo
            !
            deallocate( temp )
            !
            close( ioCovariance )
            !
        else
            write( *, * ) "Error: read_CmSqrt > Cant open file [", cfile, "]!"
            stop
        endif
        !
        !> create a huge sparse vector to make sure we accommodate all smoothing exceptions
        self%S = cVectorSparse3D_SG_t( Nx * Ny * NzEarth, FACE )
        !
        !> now, parse the exceptions
        nS = 0
        do k = 2, NzEarth
            do j = 2, Ny
                do i = 2, Nx
                    do n = 1, nrules
                        !
                        !> look back in the X-direction
                        if(((self%mask%v(i-1,j,k) == mask1(n)) .AND. (self%mask%v(i,j,k) == mask2(n))) &
                        .OR. ((self%mask%v(i-1,j,k) == mask2(n)) .AND. (self%mask%v(i,j,k) == mask1(n)))) &
                        then
                            nS = nS+1
                            self%S%i(nS) = i-1
                            self%S%j(nS) = j
                            self%S%k(nS) = k
                            self%S%xyz(nS) = 1
                            self%S%c(nS) = smoothing(n)
                        endif
                        !
                        !> look back in the Y-direction
                        if(((self%mask%v(i,j-1,k) == mask1(n)) .AND. (self%mask%v(i,j,k) == mask2(n))) &
                        .OR. ((self%mask%v(i,j-1,k) == mask2(n)) .AND. (self%mask%v(i,j,k) == mask1(n)))) &
                        then
                            nS = nS+1
                            self%S%i(nS) = i
                            self%S%j(nS) = j-1
                            self%S%k(nS) = k
                            self%S%xyz(nS) = 2
                            self%S%c(nS) = smoothing(n)
                        endif
                        !> look back in the Z-direction
                        if(((self%mask%v(i,j,k-1) == mask1(n)) .AND. (self%mask%v(i,j,k) == mask2(n))) &
                        .OR. ((self%mask%v(i,j,k-1) == mask2(n)) .AND. (self%mask%v(i,j,k) == mask1(n)))) &
                        then
                            nS = nS+1
                            self%S%i(nS) = i
                            self%S%j(nS) = j
                            self%S%k(nS) = k-1
                            self%S%xyz(nS) = 3
                            self%S%c(nS) = smoothing(n)
                        endif
                    enddo
                enddo
            enddo
        enddo
        !
        deallocate( mask1, mask2, smoothing )
        !
        !> now, truncate the smoothing vector to the correct number of components
        call self%S%reall( nS )
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
        class( ModelCovariance_t ), intent( in ) :: self
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
            call errStop( "RecursiveAR > The input arrays should be of the same shapes!" )
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
    !
    subroutine RecursiveARInv( self, w, v, n)
        implicit none
        !
        class( ModelCovariance_t ), intent( in ) :: self
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
        if( maxval( abs( shape(w) - shape(v) ) ) > 0 ) then
            !
            call errStop( "RecursiveARInv > The input arrays should be of the same shapes in RecursiveARInv" )
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
        class( ModelCovariance_t ), intent( in ) :: self
        integer, intent( in ) :: i, j, k
        real( kind=prec ) :: alpha
        !
        integer :: n
        !
        if( .NOT. associated( self%Sx ) ) then
            call errStop( "SmoothX > self%Sx has to be is_allocated before calling this" )
        endif
        !
        alpha = self%Sx(k)
        !
        if((i < 1) .OR. (i > self%Nx)) then
            call errStop( "SmoothX > index i out of bounds in SmoothX(i,j,k)" )
        elseif((j < 1) .OR. (j > self%Ny)) then
            call errStop( "SmoothX > index j out of bounds in SmoothX(i,j,k)" )
        elseif((k < 1) .OR. (k > self%NzEarth)) then
            call errStop( "SmoothX > index k out of bounds in SmoothX(i,j,k)" )
        endif
        !
        if( self%S%is_allocated .AND. (self%S%nCoeff > 0) ) then
            !> scan through the special rules and possibly update the result
            do n = 1,self%S%nCoeff
                if(self%S%xyz(n) == 1) then
                    if((self%S%i(n) == i) .AND. (self%S%j(n) == j) .AND. (self%S%k(n) == k)) then
                        alpha = self%S%c(n)
                        exit
                    endif
                endif
            enddo
        endif
        !
    end function SmoothX
    !
    !> Computes the smoothing coefficient in the y-direction based on self
    !
    function SmoothY( self, i, j, k ) result( beta )
        implicit none
        !
        class( ModelCovariance_t ), intent( in ) :: self
        integer, intent( in ) :: i, j, k
        real( kind=prec ) :: beta
        !
        integer :: n
        !
        if( .NOT. associated( self%Sy ) ) then
            call errStop( "SmoothY > self%Sy has to be is_allocated before calling this" )
        endif
        !
        beta = self%Sy(k)
        !
        if( ( i < 1 ) .OR. ( i > self%Nx ) ) then
            call errStop( "SmoothY > index i out of bounds in Y(i,j,k)" )
        elseif( (j < 1) .OR. ( j > self%Ny) ) then
            call errStop( "SmoothY > index j out of bounds in Y(i,j,k)" )
        elseif( (k < 1) .OR. ( k > self%NzEarth) ) then
            call errStop( "SmoothY > index k out of bounds in Y(i,j,k)" )
        endif
        !
        if( self%S%is_allocated .AND. ( self%S%nCoeff > 0 ) ) then
            !> scan through the special rules and possibly update the result
            do n = 1,self%S%nCoeff
                if( self%S%xyz(n) == 2) then
                    if( ( self%S%i(n) == i ) .AND. ( self%S%j(n) == j ) .AND. ( self%S%k(n) == k ) ) then
                        beta = self%S%c(n)
                        exit
                    endif
                endif
            enddo
        endif
        !
    end function SmoothY
    !
    !> computes the smoothing coefficient in the z-direction based on self
    !
    function SmoothZ( self, i, j, k ) result( gamma )
        implicit none
        !
        class( ModelCovariance_t ), intent( in ) :: self
        integer, intent( in ) :: i, j, k
        real( kind=prec ) :: gamma
        !
        integer :: n
        !
        gamma = self%Sz
        !
        if((i < 1) .OR. (i > self%Nx)) then
            call errStop( "SmoothZ > index i out of bounds in Z(i,j,k)" )
        elseif((j < 1) .OR. (j > self%Ny)) then
            call errStop( "SmoothZ > index j out of bounds in Z(i,j,k)" )
        elseif((k < 1) .OR. (k > self%NzEarth)) then
            call errStop( "SmoothZ > index k out of bounds in Z(i,j,k)" )
        endif
        !
        if( self%S%is_allocated .AND. ( self%S%nCoeff > 0 ) ) then
            !> scan through the special rules and possibly update the result
            do n = 1,self%S%nCoeff
                if( self%S%xyz(n) == 3 ) then
                    if( ( self%S%i(n) == i ) .AND. ( self%S%j(n) == j ) .AND. ( self%S%k(n) == k ) ) then
                        gamma = self%S%c(n)
                        exit
                    endif
                endif
            enddo
        endif
        !
    end function SmoothZ
    !
    !> computes the scaling coefficient based on self
    !
    function Scaling( self, i, j, k ) result( c )
        implicit none
        !
        class( ModelCovariance_t ), intent( in ) :: self
        integer, intent( in ) :: i,j,k
        real( kind=prec ) :: c, alpha, beta, gamma
        !
        if((i < 1) .OR. (i > self%Nx)) then
            call errStop( "Scaling > index i out of bounds in Scaling(i,j,k)" )
        elseif((j < 1) .OR. (j > self%Ny)) then
            call errStop( "Scaling > index j out of bounds in Scaling(i,j,k)" )
        elseif((k < 1) .OR. (k > self%NzEarth)) then
            call errStop( "Scaling > index k out of bounds in Scaling(i,j,k)" )
        endif
        !
        alpha = self%SmoothX( i, j, k )
        beta = self%SmoothY( i, j, k )
        gamma = self%SmoothZ( i, j, k )
        !
        if( self%mask%v( i, j, k ) == AIR ) then
            c = 0.0
        elseif( self%mask%v(i,j,k) == OCEAN ) then
            c = 0.0
        else
            c = ( 1 - alpha ) ** 2 * ( 1 - beta ) ** 2 * ( 1 - gamma ) ** 2
        endif
        !
    end function Scaling
    !
end Module ModelCovariance

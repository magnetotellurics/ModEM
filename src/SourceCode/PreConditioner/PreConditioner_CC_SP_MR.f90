!
!> Derived class to define a Curl-Curl PreConditioner
!> Using Sparse Matrices System
!
module PreConditioner_CC_SP_MR
    !
    use PreConditioner
    use ModelOperator_SP_V1
    use ModelOperator_SP_V2
    !
    type, extends( PreConditioner_t ) :: PreConditioner_CC_SP_MR_t
        !
        !> upper and lower triangular
        type( spMatCSR_Cmplx ) :: L, U
        !
        type( spMatCSR_Cmplx ) :: LH, UH
        !
        contains
            !
            final :: PreConditioner_CC_SP_MR_dtor
            !
            procedure, public :: setPreConditioner => setPreConditioner_CC_SP_MR !> This needs to be called by Solver    object
            !
            procedure, public :: LTSolve => LTSolvePreConditioner_CC_SP_MR !> These are left(M1) and right(M2)
            procedure, public :: UTSolve => UTSolvePreConditioner_CC_SP_MR !> preconditioning matrices for curl-curl equation.
            procedure, public :: LUSolve => LUSolvePreConditioner_CC_SP_MR !> preconditoner for symmetric divCGrad operator
            !
    end type PreConditioner_CC_SP_MR_t
    !
    interface PreConditioner_CC_SP_MR_t
        module procedure PreConditioner_CC_SP_MR_ctor
    end interface PreConditioner_CC_SP_MR_t
    !
contains
    !
    !> No subroutine briefing
    !
    function PreConditioner_CC_SP_MR_ctor( model_operator ) result( self ) 
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( PreConditioner_CC_SP_MR_t ) :: self
        !
        !write( *, * ) "Constructor PreConditioner_CC_SP_MR_t"
        !
        self%omega = R_ZERO
        !
        self%model_operator => model_operator
        !
    end function PreConditioner_CC_SP_MR_ctor
    !
    !> PreConditioner_DC_SP destructor
    !
    subroutine PreConditioner_CC_SP_MR_dtor( self )
        implicit none
        !
        type( PreConditioner_CC_SP_MR_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor PreConditioner_CC_SP_MR_t"
        !
        call deall_spMatCSR( self%L )
        call deall_spMatCSR( self%U )
        call deall_spMatCSR( self%LH )
        call deall_spMatCSR( self%UH )
        !
    end subroutine PreConditioner_CC_SP_MR_dtor
    !
    !> Block DILU preconditioner for CC operator, should be
    !> comparable to what is implemented for matrix-free version
    !
    subroutine setPreConditioner_CC_SP_MR( self, omega )
        implicit none
        !
        class( PreConditioner_CC_SP_MR_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega
        !
        integer, allocatable, dimension(:) :: ix, iy, iz
        real( kind=prec ), allocatable, dimension(:) :: d
        integer :: nx, ny, na, nz
        integer :: nEdge, nEdgeT, n, i, j
        type( spMatCSR_real ) :: CCxx
        type( spMatCSR_Cmplx ) :: Axx
        type( spMatCSR_Cmplx ), pointer  :: Lblk(:), Ublk(:)
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_SP_V1_t )
                !
                call errStop( "setPreConditioner_CC_SP_MR > ModelOperator_SP_V1_t NEED CHECK" )
                !
            class is( ModelOperator_SP_V2_t )
                !
                !>  block DILU preconditioner for CC operator, should be
                !> comparable to what is implemented for matrix-free verson
                !if( model_operator%metric%grid%is_initialized ) then
                    !> This generates indicies(in list of interior edges)
                    !> for x, y, z edges.
                    !
                    call model_operator%metric%grid%numberOfEdges( nx, ny, nz )
                    !
                    !> ix
                    nEdgeT = 0
                    nEdge = 0
                    do i = 1, size( model_operator%metric%grid%EDGEi )
                        if( model_operator%metric%grid%EDGEi(i) <= nx ) nEdge = nEdge + 1
                    enddo
                    !
                    allocate( ix(nEdge) )
                    !
                    ix =(/(j, j = nEdgeT + 1, nEdgeT + nEdge) /)
                    !
                    !> iy
                    nEdgeT = nEdgeT + nEdge
                    nEdge = 0
                    do i = 1, size( model_operator%metric%grid%EDGEi )
                        !
                        if( model_operator%metric%grid%EDGEi(i) > nx .AND. &
                            model_operator%metric%grid%EDGEi(i) <= nx + ny ) then
                            nEdge = nEdge + 1
                        endif
                        !
                    enddo
                    !
                    allocate( iy(nEdge) )
                    !
                    iy =(/(j,j = nEdgeT + 1, nEdgeT + nEdge) /)
                    !
                    !> iz
                    nEdgeT = nEdgeT+nEdge
                    nEdge = 0
                    do i = 1, size( model_operator%metric%grid%EDGEi )
                        !
                        if( model_operator%metric%grid%EDGEi(i) > nx + ny ) then
                            nEdge = nEdge + 1
                        end if
                        !
                    end do
                    !
                    allocate( iz(nEdge) )
                    !
                    iz =(/(j, j = nEdgeT + 1, nEdgeT + nEdge) /)
                    ! !
                ! else
                    ! !
                    ! nEdgeT = 0
                    ! call model_operator%metric%grid%setLimits( XEDGE, nx, ny, nz )
                    ! nEdge = nx*(ny-2)*(nz-2)
                    ! allocate(ix(nEdge))
                    ! ix =(/(j,j=nEdgeT+1,nEdgeT+nEdge) /)
                    ! !
                    ! nEdgeT = nEdgeT+nEdge
                    ! call model_operator%metric%grid%setLimits( YEDGE, nx, ny, nz )
                    ! nEdge =(nx-2)*ny*(nz-2)
                    ! allocate(iy(nEdge))
                    ! iy =(/(j,j=nEdgeT+1,nEdgeT+nEdge) /)
                    ! !
                    ! nEdgeT = nEdgeT+nEdge
                    ! call model_operator%metric%grid%setLimits( ZEDGE, nx, ny, nz )
                    ! nEdge =(nx-2)*(ny-2)*nz
                    ! allocate(iz(nEdge))
                    ! iz =(/(j,j=nEdgeT+1,nEdgeT+nEdge) /)
                    ! !
                ! end if
                ! !
                !> construct submatrices for x, y, z components
                allocate( Lblk(3) )
                allocate( Ublk(3) )
                !
                call SubMatrix_Real( model_operator%AAii, ix, ix, CCxx )
                !
                n = size(ix)
                allocate(d(n))
                d = model_operator%VomegaMuSig(ix)
                !
                call CSR_R2Cdiag( CCxx, d, Axx )
                call ilu0_Cmplx( Axx, Lblk(1), Ublk(1) )
                !
                deallocate(d)
                !
                call SubMatrix_Real( model_operator%AAii, iy, iy, CCxx )
                n = size(iy)
                allocate(d(n))
                d = model_operator%VomegaMuSig(iy)
                !
                call CSR_R2Cdiag( CCxx, d, Axx )
                call ilu0_cmplx( Axx, Lblk(2), Ublk(2) )
                !
                deallocate(d)
                !
                call SubMatrix_Real( model_operator%AAii, iz, iz, CCxx )
                n = size(iz)
                allocate(d(n))
                d = model_operator%VomegaMuSig(iz)
                !
            class default
                call errStop( "setPreConditioner_CC_SP_MR > Unclassified target" )
            !
        end select
        !
        call CSR_R2Cdiag( CCxx, d, Axx )
        call ilu0_Cmplx( Axx, Lblk(3), Ublk(3) )
        !
        deallocate(d)
        !
        !> could merge into a single LT and UT matrix, or solve systems individually
        call BlkDiag_Cmplx( Lblk, self%L )
        call BlkDiag_Cmplx( Ublk, self%U )
        call CMATtrans( self%L, self%LH )
        call CMATtrans( self%U, self%UH )
        deallocate(ix)
        deallocate(iy)
        deallocate(iz)
        call deall_spMatCSR(CCxx)
        call deall_spMatCSR(Axx)
        !
        do j = 1, 3
            call deall_spMatCSR(Lblk(j))
            call deall_spMatCSR(Ublk(j))
        enddo
        !
        deallocate(Lblk)
        deallocate(Ublk)
        !
    end subroutine setPreConditioner_CC_SP_MR
    !
    !> Implement the sparse matrix solve for curl-curl operator
    !
    subroutine LTSolvePreConditioner_CC_SP_MR( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_CC_SP_MR_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ) :: adjoint
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_e_v, out_e_v_int
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "LTSolvePreConditioner_CC_SP_MR > in_e not allocated yet" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "LTSolvePreConditioner_CC_SP_MR > out_e not allocated" )
        endif
        !
        in_e_v = in_e%getArray()
        !
        out_e_v = out_e%getArray()
        out_e_v_int = out_e_v( out_e%indInterior() )
        !
        if( adjoint ) then
            !
            call UTsolve_Cmplx( self%LH, in_e_v( in_e%indInterior() ), out_e_v_int )
            !
        else
            !
            out_e_v_int = C_ZERO
            !
            call LTsolve_Cmplx( self%L, in_e_v( in_e%indInterior() ), out_e_v_int )
            !
        endif
        !
        out_e_v( out_e%indInterior() ) = out_e_v_int
        !
        call out_e%setArray( out_e_v )
        !
    end subroutine LTSolvePreConditioner_CC_SP_MR
    !
    !> Procedure UTSolvePreConditioner_CC_SP_MR
    !> Purpose: to solve the upper triangular system(or it"s adjoint);
    !> for the d-ilu pre-condtioner
    !
    subroutine UTSolvePreConditioner_CC_SP_MR( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_CC_SP_MR_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ) :: adjoint
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_e_v
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v_int, out_e_v_int
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "UTSolvePreConditioner_CC_SP_MR > in_e not allocated yet" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "UTSolvePreConditioner_CC_SP_MR > out_e not allocated" )
        endif
        !
        in_e_v = in_e%getArray()
        in_e_v_int = in_e_v( in_e%indInterior() )
        !
        out_e_v = out_e%getArray()
        out_e_v_int = out_e_v( out_e%indInterior() )
        !
        if( adjoint ) then
            !
            call LTsolve_Cmplx( self%UH, in_e_v_int, out_e_v_int )
            !
        else
            !
            call UTsolve_Cmplx( self%U, in_e_v_int, out_e_v_int )
            !
        endif
        !
        out_e_v( out_e%indInterior() ) = out_e_v_int
        !
        call out_e%setArray( out_e_v )
        !
    end subroutine UTSolvePreConditioner_CC_SP_MR
    !
    !> Procedure LUSolvePreConditioner_CC_SP_MR
    !> this is dummy routine required by abstract preconditioner class
    !
    subroutine LUSolvePreConditioner_CC_SP_MR( self, in_phi, out_phi )
        implicit none
        !
        class( PreConditioner_CC_SP_MR_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: in_phi
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        call errStop( "LUSolvePreConditioner_CC_SP_MR not implemented" )
        !
    end subroutine LUSolvePreConditioner_CC_SP_MR
    !
end module PreConditioner_CC_SP_MR


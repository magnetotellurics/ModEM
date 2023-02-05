!
!> Derived class to define a PreConditioner_CC_SP
!
module PreConditioner_CC_SP
    !
    use PreConditioner
    use ModelOperator_SP
    !
    type, extends( PreConditioner_t ) :: PreConditioner_CC_SP_t
        !
        !> upper and lower triangular
        type( spMatCSR_Cmplx ) :: L, U
        !
        type( spMatCSR_Cmplx ) :: LH, UH
        !
        contains
            !
            final :: PreConditioner_CC_SP_dtor
            !
            procedure, public :: setPreConditioner => setPreConditioner_CC_SP !> This needs to be called by Solver    object
            !
            procedure, public :: LTSolve => LTSolvePreConditioner_CC_SP !> These are left (M1) and right (M2)
            procedure, public :: UTSolve => UTSolvePreConditioner_CC_SP !> preconditioning matrices for curl-curl equation.
            procedure, public :: LUSolve => LUSolvePreConditioner_CC_SP !> preconditoner for symmetric divCgrad operator
            !
    end type PreConditioner_CC_SP_t
    !
    interface PreConditioner_CC_SP_t
        module procedure PreConditioner_CC_SP_ctor
    end interface PreConditioner_CC_SP_t
    !
contains
    !
    !> No subroutine briefing
    !
    function PreConditioner_CC_SP_ctor( model_operator ) result( self ) 
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( PreConditioner_CC_SP_t ) :: self
        !
        !write( *, * ) "Constructor PreConditioner_CC_SP_t"
        !
        self%omega = R_ZERO
        !
        self%model_operator => model_operator
        !
    end function PreConditioner_CC_SP_ctor
    !
    !> PreConditioner_DC_SP destructor
    subroutine PreConditioner_CC_SP_dtor( self )
        implicit none
        !
        type( PreConditioner_CC_SP_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor PreConditioner_CC_SP_t"
        !
        call deall_spMatCSR( self%L )
        call deall_spMatCSR( self%U )
        call deall_spMatCSR( self%LH )
        call deall_spMatCSR( self%UH )
        !
    end subroutine PreConditioner_CC_SP_dtor
    !
    !> Block DILU preconditioner for CC operator, should be
    !> comparable to what is implemented for matrix-free version
    !
    subroutine setPreConditioner_CC_SP( self, omega )
        implicit none
        !
        class( PreConditioner_CC_SP_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega
        !
        integer, allocatable, dimension(:) :: ix, iy, iz
        real(kind=prec), allocatable, dimension(:) :: d
        integer :: nx, ny, na, nz
        integer :: nEdge, nEdgeT, n, j
        type( spMatCSR_real ) :: CCxx
        type( spMatCSR_Cmplx ) :: Axx
        type( spMatCSR_Cmplx ), pointer  :: Lblk(:), Ublk(:)
        !
        ! find indicies of x, y, z elements
        !
        ! this generates indicies (in list of interior edges)
        ! for x, y, z edges
        nEdgeT = 0
        !
        call setLimitsSP( XEDGE, self%model_operator%metric%grid, nx, ny, nz )
        !
        nEdge = nx * ( ny - 2 ) * ( nz - 2 )
        !
        allocate( ix( nEdge ) )
        !
        ix = (/ (j, j=nEdgeT+1, nEdgeT+nEdge) /)
        !
        nEdgeT = nEdgeT + nEdge
        !
        call setLimitsSP( YEDGE, self%model_operator%metric%grid, nx, ny, nz )
        !
        nEdge = ( nx - 2 ) * ny * ( nz - 2 )
        !
        allocate( iy( nEdge ) )
        !
        iy = (/ (j, j=nEdgeT+1, nEdgeT+nEdge) /)
        !
        nEdgeT = nEdgeT+nEdge
        !
        call setLimitsSP( ZEDGE, self%model_operator%metric%grid, nx, ny, nz )
        !
        nEdge = ( nx - 2 ) * ( ny - 2 ) * nz
        !
        allocate( iz( nEdge ) )
        !
        iz = (/ (j, j=nEdgeT+1, nEdgeT+nEdge) /)
        !
        ! Construct submatrices for x, y, z components
        allocate( Lblk(3) )
        allocate( Ublk(3) )
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_SP_t )
                !
                call SubMatrix_Real( model_operator%CCii, ix, ix, CCxx )
                !
                n = size(ix)
                !
                allocate( d(n) )
                !
                d = model_operator%VomegaMuSig(ix)
                !
                call CSR_R2Cdiag( CCxx, d, Axx )
                !
                call dilu_Cmplx( Axx, Lblk(1), Ublk(1) )
                !
                deallocate(d)
                !
                call SubMatrix_Real( model_operator%CCii, iy, iy, CCxx )
                !
                n = size(iy)
                !
                allocate( d(n) )
                !
                d = model_operator%VomegaMuSig(iy)
                !
                call CSR_R2Cdiag( CCxx, d, Axx )
                !
                call dilu_Cmplx( Axx, Lblk(2), Ublk(2) )
                !
                deallocate(d)
                !
                call SubMatrix_Real( model_operator%CCii, iz, iz, CCxx )
                !
                n = size(iz)
                !
                allocate( d(n) )
                !
                d = model_operator%VomegaMuSig(iz)
                !
                call CSR_R2Cdiag( CCxx, d, Axx )
                !
                call dilu_Cmplx( Axx, Lblk(3), Ublk(3) )
                !
                deallocate(d)
                !
                ! Could merge into a single LT and UT matrix, or solve systems individually
                call BlkDiag_Cmplx( Lblk, self%L )
                !
                call BlkDiag_Cmplx( Ublk, self%U )
                !
                call CMATtrans( self%L, self%LH )
                call CMATtrans( self%U, self%UH )
                !
                deallocate( ix, iy, iz )
                !
                call deall_spMatCSR( CCxx )
                !
                call deall_spMatCSR( Axx )
                !
                do j = 1, 3
                    call deall_spMatCSR( Lblk(j) )
                    call deall_spMatCSR( Ublk(j) )
                enddo
                !
                deallocate( Lblk, Ublk )
                !
            class default
                stop "LTSolvePreConditioner_CC_SP: Unclassified ModelOperator"
            !
        end select
        !
    end subroutine setPreConditioner_CC_SP
    !
    !> Implement the sparse matrix solve for curl-curl operator
    !
    subroutine LTSolvePreConditioner_CC_SP( self, inE, outE, adjoint )
        implicit none
        !
        class( PreConditioner_CC_SP_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: inE
        class( Vector_t ), intent( inout ) :: outE
        logical, intent( in ) :: adjoint
        !
        complex( kind=prec ), allocatable, dimension(:) :: temp_array_inE, temp_array_outE
        !
        temp_array_inE = inE%getArray()
        !
        if( adjoint ) then
            !
            call UTsolve_Cmplx( self%LH, temp_array_inE, temp_array_outE )
            !
        else
            !
            call LTsolve_Cmplx( self%L, temp_array_inE, temp_array_outE )
            !
        endif
        !
        deallocate( temp_array_inE )
        !
        call outE%setArray( temp_array_outE )
        !
        deallocate( temp_array_outE )
        !
    end subroutine LTSolvePreConditioner_CC_SP
    !
    !> Procedure UTSolvePreConditioner_CC_SP
    !> Purpose: to solve the upper triangular system (or it"s adjoint);
    !> for the d-ilu pre-condtioner
    !
    subroutine UTSolvePreConditioner_CC_SP( self, inE, outE, adjoint )
        implicit none
        !
        class( PreConditioner_CC_SP_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: inE
        class( Vector_t ), intent( inout ) :: outE
        logical, intent( in ) :: adjoint
        !
        complex( kind=prec ), allocatable, dimension(:) :: temp_array_inE, temp_array_outE
        !
        temp_array_inE = inE%getArray()
        !
        if( adjoint ) then
            !
            call LTsolve_Cmplx( self%UH, temp_array_inE, temp_array_outE )
            !
        else
            !
            call LTsolve_Cmplx( self%U, temp_array_inE, temp_array_outE )
            !
        endif
        !
        deallocate( temp_array_inE )
        !
        call outE%setArray( temp_array_outE )
        !
        deallocate( temp_array_outE )
        !
    end subroutine UTSolvePreConditioner_CC_SP
    !
    !> Procedure LUSolvePreConditioner_CC_SP
    !> this is dummy routine required by abstract preconditioner class
    subroutine LUSolvePreConditioner_CC_SP( self, inPhi, outPhi )
        implicit none
        !
        class( PreConditioner_CC_SP_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: inPhi
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        STOP "Error: LUsolve is not coded for this pre-conditioner class"
        !
    end subroutine LUSolvePreConditioner_CC_SP
    !
end module PreConditioner_CC_SP


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
        integer :: sum_x, sum_y, sum_z
        real( kind=prec ), allocatable, dimension(:) :: d
        integer :: nx, ny, na, nz
        integer :: nEdge, nEdgeT, n, i, j
        type( spMatCSR_real ) :: Mii, CCxx
        type( spMatCSR_Cmplx ) :: Axx
        type( spMatCSR_Cmplx ), pointer, dimension(:)  :: Lblk, Ublk
        !
        !> Generate the index arrays ix, iy and iz from grid%iXYZinterior
        select type( grid => self%model_operator%metric%grid )
            !
            class is( Grid3D_MR_t )
                !
                !> Find and sum all x=1, y=2, z=3; for ix, iy, iz allocation sizes
                sum_x = 0; sum_y = 0; sum_z = 0
                !
                do i = 1, size( grid%iXYZinterior )
                    !
                    if( grid%iXYZinterior(i) == 1 ) sum_x = sum_x + 1
                    !
                    if( grid%iXYZinterior(i) == 2 ) sum_y = sum_y + 1
                    !
                    if( grid%iXYZinterior(i) == 3 ) sum_z = sum_z + 1
                    !
                enddo
                !
                allocate( ix( sum_x ) )
                allocate( iy( sum_y ) )
                allocate( iz( sum_z ) )
                !
                !> Populate ix, iy, iz
                sum_x = 1; sum_y = 1; sum_z = 1
                !
                do i = 1, size( grid%iXYZinterior )
                    !
                    if( grid%iXYZinterior(i) == 1 ) then
                        ix(sum_x) = i
                        sum_x = sum_x + 1
                    endif
                    !
                    if( grid%iXYZinterior(i) == 2 ) then
                        iy(sum_y) = i
                        sum_y = sum_y + 1
                    endif
                    !
                    if( grid%iXYZinterior(i) == 3 ) then
                        iz(sum_z) = i
                        sum_z = sum_z + 1
                    endif
                    !
                enddo
                !
            class default
                call errStop( "setPreConditioner_CC_SP_MR > Grid must be MR" )
            !
        end select
        !
        !> Construct submatrices for x, y, z components
        allocate( Lblk(3) )
        allocate( Ublk(3) )
        !
        !> Ensure that the model_operator is SP
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_SP_t )
                !
                !> Differentiation between SP1 (CCii) and SP2 (AAii)
                select type( model_operator )
                    !
                    class is( ModelOperator_SP_V1_t )
                        !
                        Mii = model_operator%CCii
                        !
                        write(6666,*) "##### setPreConditioner_CC_SP_MR_V1"
                        !
                    class is( ModelOperator_SP_V2_t )
                        !
                        Mii = model_operator%AAii
                        !
                        write(6666,*) "##### setPreConditioner_CC_SP_MR_V2"
                        !
                    class default
                        call errStop( "setPreConditioner_CC_SP_MR > model_operator must be SP V1 or V2" )
                    !
                end select
                !
                call SubMatrix_Real( Mii, ix, ix, CCxx )
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
                call SubMatrix_Real( Mii, iy, iy, CCxx )
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
                call SubMatrix_Real( Mii, iz, iz, CCxx )
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
                call BlkDiag_Cmplx( Ublk, self%U )
                !
                call CMATtrans( self%L, self%LH )
                call CMATtrans( self%U, self%UH )
                !
                deallocate( ix, iy, iz )
                !
                call deall_spMatCSR( CCxx )
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
                call errStop( "setPreConditioner_CC_SP_MR > Unclassified model_operator" )
            !
        end select
        !
    end subroutine setPreConditioner_CC_SP_MR
    !
    !> Implement the sparse matrix LT solve for curl-curl operator
    !
    subroutine LTSolvePreConditioner_CC_SP_MR( self, in_e, out_e, adjoint )
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
            call errStop( "LTSolvePreConditioner_CC_SP_MR > in_e not allocated yet" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "LTSolvePreConditioner_CC_SP_MR > out_e not allocated" )
        endif
        !
        in_e_v = in_e%getArray()
        in_e_v_int = in_e_v( in_e%indInterior() )
        !
        !
        out_e_v = out_e%getArray()
        out_e_v_int = out_e_v( out_e%indInterior() )
        !
        if( adjoint ) then
            !
            call UTsolve_Cmplx( self%LH, in_e_v_int, out_e_v_int )
            !
        else
            !
            out_e_v_int = C_ZERO
            !
            call LTsolve_Cmplx( self%L, in_e_v_int, out_e_v_int )
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
    !> Purpose: to solve the upper triangular system(or its adjoint);
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


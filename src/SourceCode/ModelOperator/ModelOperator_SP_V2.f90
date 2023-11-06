!
!> Derived class to define a ModelOperator
!> with basic operations for Sparse Matrices
!
module ModelOperator_SP_V2
    !
    use ModelOperator_SP
    !
    type, extends( ModelOperator_SP_t ) :: ModelOperator_SP_V2_t
        !
        !> sparse matrix representation
        !> of the modified system equation
        type( spMatCSR_Real ) :: AAii, AAit
        !
        !> save only the Earth part of GD
        type( spMatCSR_Real ) :: GDii
        !
        contains
            !
            final :: ModelOperator_SP_V2_dtor
            !
            !> Setup
            procedure, public :: setCond => setCond_ModelOperator_SP_V2
            !
            procedure, private :: GradDivSetup2, airNIndex
            !
            !> Operations
            procedure, public :: amult => amult_ModelOperator_SP_V2
            !
            !> Miscellaneous
            procedure, public :: print => print_ModelOperator_SP_V2
            !
    end type ModelOperator_SP_V2_t
    !
    interface ModelOperator_SP_V2_t
        module procedure ModelOperator_SP_V2_ctor
    end interface ModelOperator_SP_V2_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelOperator_SP_V2_ctor( grid ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        !
        type( ModelOperator_SP_V2_t ) :: self
        !
        call self%baseInit
        !
        call self%create( grid )
        !
    end function ModelOperator_SP_V2_ctor
    !
    !> ModelOperator_SP_V2 destructor
    !
    subroutine ModelOperator_SP_V2_dtor( self )
        implicit none
        !
        type( ModelOperator_SP_V2_t ), intent( inout ) :: self
        !
        call self%baseDealloc
        !
        call self%dealloc
        !
        call deall_spMatCSR( self%AAii )
        call deall_spMatCSR( self%AAit )
        call deall_spMatCSR( self%GDii )
        !
    end subroutine ModelOperator_SP_V2_dtor
    !
    !> No subroutine briefing
    !
    subroutine setCond_ModelOperator_SP_V2( self, sigma, omega_in )
        implicit none
        !
        class( ModelOperator_SP_V2_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: omega_in
        !
        class( Vector_t ), allocatable :: sigma_edge
        class( Scalar_t ), allocatable :: sigma_node
        real( kind=prec ), allocatable, dimension(:) :: sigma_edge_array, sigma_edge_b0_array
        real( kind=prec ), allocatable, dimension(:) :: v_edge_array, sigma_node_array
        !
        !> SigEdge
        call self%metric%createVector( real_t, EDGE, sigma_edge )
        !
        call sigma%PDEmapping( sigma_edge )
        !
        sigma_edge_array = sigma_edge%getArray()
        !
        deallocate( sigma_edge )
        !
        sigma_edge_b0_array = sigma_edge_array
        !
        !> SigNode
        call self%metric%createScalar( real_t, NODE, sigma_node )
        !
        call sigma%nodeCond( sigma_node )
        !
        sigma_node_array = sigma_node%getArray()
        !
        deallocate( sigma_node )
        !
        !> force the boundary to be zeros
        sigma_edge_b0_array( self%metric%grid%EDGEb ) = R_ZERO
        sigma_node_array( self%metric%grid%NODEb ) = R_ZERO
        !
        !> modify the system equation here,
        !> as the GD should be updated whenever the omega or the conductivity
        !> is updated
        !> could find a better place for this
        call self%GradDivSetup2( sigma_edge_b0_array, sigma_node_array )
        !
        v_edge_array = self%metric%v_edge%getArray()
        !
        self%VomegaMuSig = MU_0 * omega_in * sigma_edge_array( self%metric%grid%EDGEi ) * v_edge_array( self%metric%grid%EDGEi )
        !
        self%omega = omega_in
        !
    end subroutine setCond_ModelOperator_SP_V2
    !
    !> Implement the sparse matrix multiply for curl-curl operator
    !> for interior elements
    !> assume output y is already allocated
    !
    subroutine amult_ModelOperator_SP_V2( self, in_e, out_e, omega, adjoint )
        implicit none
        !
        class( ModelOperator_SP_V2_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        real( kind=prec ), intent( in ) :: omega
        logical, intent( in ) :: adjoint
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_e_v
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v_int, out_e_v_int
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "amult_ModelOperator_SP_V2 > in_e not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "amult_ModelOperator_SP_V2 > out_e not allocated" )
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
            call RMATxCVEC( self%AAit, in_e_v_int, out_e_v_int )
            !
            out_e_v_int = out_e_v_int - ONE_I * ISIGN * self%VomegaMuSig * in_e_v( in_e%indInterior() )
            !
        else
            !
            call RMATxCVEC( self%AAii, in_e_v_int, out_e_v_int )
            !
            out_e_v_int = out_e_v_int + ONE_I * ISIGN * self%VomegaMuSig * in_e_v( in_e%indInterior() )
            !
        endif
        !
        out_e_v( in_e%indInterior() ) = out_e_v_int
        !
        call out_e%setArray( out_e_v )
        !
    end subroutine amult_ModelOperator_SP_V2
    !
    !> a subroutine to add divergence enforcement into system Matrix A to avoid
    !> using divergence correction with iterative solvers.
    !> essentially this builds Grad(Div(sigma*E)) = 0 and add it to the CC
    !> operator
    !
    !> The curl-Curl operator in:
    !> curl(curl(E))+ i*omega*sigma*E = 0
    !> can be rewritten as:
    !>    A = grad(div(E))-laplacian(E)
    !> as the divergence in the Air or the Earth requires:
    !>    GD = grad(div(sigma*E)) == 0
    !> the modified A can be expressed as:
    !>    A = A - GD
    !> in practice, VGD is added instead of GD
    !> for the AIR:
    !>    VGaDa
    !> for the EARTH:
    !>    VGeDe
    !
    subroutine GradDivSetup2( self, SigEdge, SigNode )
        implicit none
        !
        class( ModelOperator_SP_V2_t ), intent( inout ) :: self
        real( kind=prec ), dimension(:), intent( in ) :: SigEdge, SigNode
        !
        type( spMatCSR_Real ) :: Dt, GDa, GDe, GD
        real( kind=prec ), allocatable, dimension(:) :: M1air, M2air
        real( kind=prec ), allocatable, dimension(:) :: M0earth, M1earth, M2earth
        real( kind=prec ), allocatable, dimension(:) :: M3earth, M4earth
        real( kind=prec ), allocatable, dimension(:) :: Nair, Nearth
        real( kind=prec ), allocatable, dimension(:) :: dual_face_area_v, v_node_v, edge_length_v
        real( kind=prec ) :: tol
        integer :: Ne, Nei, Nn, fid, i
        !
        Nn = size( self%metric%grid%NODEi ) + size( self%metric%grid%NODEb )
        Nei = size( self%metric%grid%EDGEi )
        Ne = Nei + size( self%metric%grid%EDGEb )
        !
        tol = 1E-8
        !
        allocate( M1air(Ne) )
        allocate( M2air(Nn) )
        allocate( M0earth(Ne) )
        allocate( M1earth(Ne) )
        allocate( M2earth(Nn) )
        allocate( M3earth(Ne) )
        allocate( M4earth(Ne) )
        allocate( Nearth(Nn) )
        allocate( Nair(Nn) )
        !
        call self%airNIndex( SigEdge, SIGMA_AIR, Nair, Nearth )
        !
        dual_face_area_v = self%metric%dual_face_area%getArray()
        !
        v_node_v = self%metric%v_node%getArray()
        !
        edge_length_v = self%metric%edge_length%getArray()
        !
        !> (for Air sigma is not necessary as it is constant everywhere)
        M1air = dual_face_area_v
        !
        M2air = Nair / v_node_v
        !
        ! rescale earth part with lambda = 1./SigNode
        ! we don't have to distinguish between air and earth edges with node based
        ! scaling factor (as node is either air or earth)
        M0earth = dual_face_area_v
        M1earth = SigEdge * dual_face_area_v
        !
        !> WORKAROUND FOR LINE 295 DIV ????
        !SigNode( self%metric%grid%NODEb ) = R_ONE
        !
        M2earth = ( R_ONE / SigNode ) * ( Nearth / v_node_v )
        !
        M3earth = R_ONE / edge_length_v
        !
        M4earth = dual_face_area_v
        !
        call RMATtrans( self%topology%G, Dt )
        !> build the Air part...
        !
        call DIAGxRMAT( M1air, self%topology%G, GD )
        call RMATxDIAG( GD, M2air, GDa )
        call RMATxRMAT( GDa, Dt, GD )
        call RMATxDIAG( GD, M1air, GDa )
        !
        ! build the Earth part...
        call DIAGxRMAT( M0earth, self%topology%G, GD )
        !
        call RMATxDIAG( GD, M2earth, GDe )
        !
        call RMATxRMAT( GDe, Dt, GD )
        !
        call RMATxDIAG( GD, M1earth, GDe )
        !
        ! now assemble the GradDiv matrix...
        call RMATplusRMAT( GDa, GDe, GD )
        !
        call deall_spMatCSR( GDa )
        !
        call subMatrix_Real( GD, self%metric%grid%EDGEi, self%metric%grid%EDGEi, self%GDii )
        !
        call RMATplusRMAT( self%CCii, self%GDii, self%AAii )
        !
        call RMATtrans( self%AAii, self%AAit )
        !
        ! build the GradDiv matrix for additional terms in RHS...
        call DIAGxRMAT( M3earth, self%topology%G, GDe )
        !
        call RMATxDIAG( GDe, M2earth, GD )
        !
        call RMATxRMAT( GD, Dt, GDe )
        !
        call RMATxDIAG( GDe, M4earth, GD )
        !
        call subMatrix_Real( GD, self%metric%grid%EDGEi, self%metric%grid%EDGEi, self%GDii )
        !
        !> 
        call deall_spMatCSR( GDe )
        call deall_spMatCSR( GD )
        call deall_spMatCSR( Dt )
        !
        !> no need to keep this if GDii is symmetric
        !> call RMATtrans(AAii,ATii)
        deallocate( M1air, M2air )
        deallocate( M0earth, M1earth, M2earth, M3earth, M4earth )
        deallocate( Nearth )
        deallocate( Nair )
        !
    end subroutine GradDivSetup2
    !
    !> this generate the air/earth domain index for Nodes only,
    !> these indexes will be used to construct distinct air/earth grad and div
    !> operators, with "node-based" scaling factors
    !
    subroutine airNIndex( self, SigEdge, AirCond, Nair, Nearth )
        implicit none
        !
        class( ModelOperator_SP_V2_t ), intent( in ) :: self
        real( kind=prec ), intent(in) :: SigEdge(:), AirCond
        real( kind=prec ), intent(inout), dimension(:) :: Nair, Nearth
        real( kind=prec ), dimension(:), allocatable :: Eair, Eearth
        !
        type( spMatCSR_Real ) :: Dt
        integer :: Ne,Nn,i
        !
        Ne = size( SigEdge )
        Nn = size( self%metric%grid%NODEi ) + size( self%metric%grid%NODEb )
        !
        call RMATtrans( self%topology%G, Dt )
        !
        Dt%val = abs( Dt%val )
        !
        !> set edge indexes
        allocate( Eair( Ne ) )
        allocate( Eearth( Ne ) )
        !
        Eair = R_ZERO
        Eearth = R_ZERO
        !
        do i = 1, Ne
            !
            if( SigEdge(i) .LT. 1.1 * AirCond ) then
                Eair(i) = R_ONE
            else
                Eearth(i) = R_ONE
            endif
            !
            if( SigEdge(i) .EQ. R_ZERO ) then ! boundary should always be zero
                Eair(i) = R_ZERO
                Eearth(i) = R_ZERO
            endif
            !
        enddo
        !
        !> set node indexes
        call RMATxRVEC( Dt, Eearth, Nair )
        !
        Nearth = R_ZERO
        ! any node that connects to at least one earth edge is an earth node
        do i = 1, Nn
            !
            if( Nair(i) .GT. R_ZERO ) then
                Nearth(i) = R_ONE
                Nair(i) = R_ZERO
            else
                Nair(i) = R_ONE
            endif
            !
        enddo
        !
        Nearth( self%metric%grid%NODEb ) = R_ZERO
        Nair( self%metric%grid%NODEb ) = R_ZERO
        !
        call deall_spMatCSR_Real( Dt )
        !
        deallocate( Eearth )
        deallocate( Eair )
        !
    end subroutine airNIndex
    !
    !> No subroutine briefing
    !
    subroutine print_ModelOperator_SP_V2( self )
        implicit none
        !
        class( ModelOperator_SP_V2_t ), intent( in ) :: self
        !
        call errStop( "print_ModelOperator_SP_V2 not implemented" )
        !
    end subroutine print_ModelOperator_SP_V2
    !
end module ModelOperator_SP_V2

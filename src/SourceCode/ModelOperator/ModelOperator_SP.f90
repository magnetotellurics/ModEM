!
!> Derived class to define a ModelOperator_SP
!
module ModelOperator_SP
    !
    use ModelOperator
    use SpOpTopology_SG
    use MetricElements_CSG
    use ModelParameterCell_SG
    !
    type, extends( ModelOperator_t ) :: ModelOperator_SP_t
        !
        type( SpOpTopology_SG_t ) :: topology_sg
        !
        integer, allocatable, dimension(:) :: EDGEi, EDGEb
        integer, allocatable, dimension(:) :: NODEi, NODEb
        !
        type( spMatCSR_Real ) :: CC, CCii, CCib
        !
        real( kind=prec ) :: omega
        !
        real( kind=prec ), allocatable, dimension(:) :: VomegaMuSig
        !
        type( spMatCSR_Real ) :: VDiv ! div : edges->nodes (interior only)
        type( spMatCSR_Real ) :: VDsG ! operator for div correction
        type( spMatCSR_Real ) :: VDs  ! divergence of current operator
        !
        type( spMatCSR_Real ) :: VDsG_L, VDsG_U
        !
        contains
            !
            final :: ModelOperator_SP_dtor
            !
            !> Setup
            procedure, public :: setEquations => setEquations_ModelOperator_SP
            procedure, public :: setCond => setCond_ModelOperator_SP
            !
            procedure, public :: divCorInit => divCorInit_ModelOperator_SP
            procedure, public :: divCorSetUp => divCorSetUp_ModelOperator_SP
            !
            !> Operations
            procedure, public :: amult => amult_ModelOperator_SP
            procedure, public :: multAib => multAib_ModelOperator_SP
            !
            procedure, public :: div => div_ModelOperator_SP
            procedure, public :: divC => divC_ModelOperator_SP
            procedure, public :: divCGrad => divCGrad_ModelOperator_SP
            !
            procedure, public :: grad => grad_ModelOperator_SP
            !
            !> Alloc/Dealloc
            procedure :: create => create_ModelOperator_SP 
            procedure :: dealloc => deallocate_ModelOperator_SP
            !
            !> Miscellaneous
            procedure, public :: print => print_ModelOperator_SP
            !
    end type ModelOperator_SP_t
    !
    interface ModelOperator_SP_t
        module procedure ModelOperator_SP_ctor
    end interface ModelOperator_SP_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelOperator_SP_ctor( grid ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        !
        type( ModelOperator_SP_t ) :: self
        !
        write( *, * ) "Constructor ModelOperator_SP"
        !
        call self%baseInit
        !
        !> Instantiation of the specific object MetricElements
        allocate( self%metric, source = MetricElements_CSG_t( grid ) )
        !
        call self%create
        !
    end function ModelOperator_SP_ctor
    !
    !> ModelOperator_SP destructor
    subroutine ModelOperator_SP_dtor( self )
        implicit none
        !
        type( ModelOperator_SP_t ), intent( inout ) :: self
        !
        write( *, * ) "Destructor ModelOperator_SP_t"
        !
        call self%baseDealloc
        !
        call self%dealloc
        !
    end subroutine ModelOperator_SP_dtor
    !
    !> No subroutine briefing
    !> using existing curl operator, create sparse matrix CC
    !> Note: this is the symmetric form, multiplied by edge volume elements
    !
    subroutine setEquations_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        type( spMatCSR_Real ) :: Temp, Ttrans
        integer :: m, n, nz
        real( kind=prec ), allocatable, dimension(:) :: Dtemp
        integer :: fid
        !
        write(*,*) "setEquations_ModelOperator_SP"
        !
        m = T%nRow
        n = T%nCol
        nz = T%row( T%nRow + 1 ) - 1
        !
        call create_spMatCSR( m, n, nz, Temp )
        call create_spMatCSR( n, m, nz, Ttrans )
        call create_spMatCSR( m, n, nz, self%CC )
        !
        call RMATxDIAG( T, real( self%metric%edge_length%getSV(), kind=prec ), Temp )
        !
        Dtemp = ( self%metric%dual_edge_length%getSV() / self%metric%face_area%getSV() )
        !
        call DIAGxRMAT( Dtemp, Temp, self%CC )
        !
        call RMATtrans( T, Ttrans )
        !
        call RMATxRMAT( Ttrans, self%CC, Temp )
        !
        call DIAGxRMAT( real( self%metric%edge_length%getSV(), kind=prec ), Temp, self%CC )
        !
        call subMatrix_Real( self%CC, self%EDGEi, self%EDGEi, self%CCii )
        !
        call subMatrix_Real( self%CC, self%EDGEi, self%EDGEb, self%CCib )
        !
        call deall_spMatCSR( Temp )
        call deall_spMatCSR( Ttrans )
        !
        self%eqset = .TRUE.
        !
        call self%divCorInit
        !
    end subroutine setEquations_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine setCond_ModelOperator_SP( self, sigma, omega )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( inout ) :: sigma
        real( kind=prec ), intent( in ), optional :: omega
        !
        integer :: i
        class( Scalar_t ), allocatable :: temp_cell_cond
        class( Vector_t ), allocatable :: sig_temp
        complex( kind=prec ), allocatable, dimension(:) :: v_edge_sv, sig_temp_sv
        !
		if( present( omega ) ) then
			write(*,*) "setCond_ModelOperator_SP :", omega
		else
			write(*,*) "setCond_ModelOperator_SP (no omega )"
		endif
        !
        if( present( omega ) ) then
            self%omega = omega
        else
            self%omega = 1.0
        endif
        !
        !> Sigma
        call self%metric%createVector( real_t, EDGE, sig_temp )
        !
        call sigma%PDEmapping( sig_temp )
        !
        sig_temp_sv = sig_temp%getSV()
        !
        !> vEdge
        v_edge_sv = self%metric%v_edge%getSV()
        !
        self%VomegaMuSig = MU_0 * self%omega * sig_temp_sv( self%EDGEi ) * v_edge_sv( self%EDGEi )
        !
        !call self%divCorSetUp
        !
    end subroutine setCond_ModelOperator_SP
    !
    !> To complete setup conductivity is required
    !> DivCorInit has to be called before this routine
    !
    subroutine divCorInit_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        type( spMatCSR_Real ) :: Temp, Temp2
        real( kind=prec ), allocatable, dimension(:) :: d, aux_vec
        integer, allocatable, dimension(:) :: allNodes
        integer :: i, m
        !
        write(*,*) "divCorInit_ModelOperator_SP"
        !
        !> set indexes for interior and boundary nodes
        call boundaryIndexSP( NODE, self%metric, self%NODEb, self%NODEi )
        !
        !> (1) first construct VDiv operator transpose of topology
        call RMATtrans( G, Temp )
        !
        !> pre-multiply by dual face area
        aux_vec = real( self%metric%dual_face_area%getSV(), kind=prec )
        !
        call RMATxDIAG( Temp, aux_vec, Temp2 )
        !
        call deall_spMatCSR( Temp )
        !
        !> select out interior nodes, edges
        call subMatrix_Real( Temp2, self%NODEi, self%EDGEi, self%VDiv )
        !
        call deall_spMatCSR( Temp2 )
        !
        !> (2) next turn G into actual gradient (not just topology,
        !> all nodes-> interior edges (not clear this is what we want!)
        allocate( d( G%nRow ) )
        !
        aux_vec = real( self%metric%edge_length%getSV(), kind=prec )
        !
        do i = 1, G%nRow
            d(i) = 1. / aux_vec(i)
        enddo
        !
        allocate( allNodes( G%nCol ) )
        !
        do i=1,G%nCol
            allNodes(i) = i
        enddo
        !
        call DIAGxRMAT( d, G, Temp )
        call subMatrix_Real( Temp, self%EDGEi, allNodes, G )
        call deall_spMatCSR( Temp )
        !
        deallocate( allNodes, d )
        !
    end subroutine divCorInit_ModelOperator_SP
    !
    !> To complete setup conductivity is required
    !> DivCorInit has to be called before this routine
    !
    subroutine divCorSetUp_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        type( spMatCSR_Real ) :: temp_matrix
        complex( kind=prec ), allocatable, dimension(:) :: v_edge_sv
        real( kind=prec ), allocatable, dimension(:) :: d
        integer, allocatable, dimension(:) :: allNodes
        integer :: n, i
        !
        write(*,*) "divCorSetUp_ModelOperator_SP"
        !
        !> Construct VDs .. multiply VDiv by Conductivity on edges; can use VomegaMuSig
        n = self%VDiv%nCol
        !
        d = ( self%VomegaMuSig / ( mu_0 * self%omega ) )
        !
        v_edge_sv = self%metric%v_edge%getSV()
        !
        d = ( d / v_edge_sv( self%EDGEi ) )
        !
        call RMATxDIAG( self%VDiv, d, self%VDs )
        !
        !>Construct VDsG: symmetric operator for divergence correction solver
        allocate( allNodes( G%nRow ) )
        !
        do i = 1, G%nRow
            allNodes( i ) = i
        enddo
        !
        call subMatrix_Real( G, allNodes, self%NODEi, temp_matrix )
        !
        call RMATxRMAT( self%VDs, temp_matrix, self%VDsG )
        !
        ! Setup preconditioner
        call Dilu_Real( self%VDsG, self%VDsG_L, self%VDsG_U )
        !
        !call CholInc_real(VDsG,VDsG_L)
        !call RMATtrans(VDsG_L,VDsG_U)
        !
        call deall_spMatCSR( temp_matrix )
        !
        deallocate( allNodes )
        !
    end subroutine divCorSetUp_ModelOperator_SP
    !
    !> Implement the sparse matrix multiply for curl-curl operator
    !> for interior elements
    !> assume output y is already allocated
    !
    subroutine amult_ModelOperator_SP( self, omega, in_e, out_e, p_adjoint )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        real( kind=prec ), intent( in ), optional :: omega
        class( Vector_t ), intent( inout ) :: in_e, out_e
        logical, intent( in ), optional :: p_adjoint
        !
        logical :: adjoint
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_e_v
        complex( kind=prec ), allocatable, dimension(:) :: out_e_v_int
        !
        in_e_v = in_e%getSV()
        !
        out_e_v = out_e%getSV()
        out_e_v_int = out_e_v( out_e%ind_interior )
        !
        call RMATxCVEC( self%CCii, in_e_v( in_e%ind_interior ), out_e_v_int )
        !
        if( present( p_adjoint ) ) then
            adjoint = p_adjoint
        else
            adjoint = .FALSE.
        endif
        !
        if( adjoint ) then
            out_e_v_int = out_e_v_int - ONE_I * ISIGN * self%VomegaMuSig * in_e_v( in_e%ind_interior )
        else
            out_e_v_int = out_e_v_int + ONE_I * ISIGN * self%VomegaMuSig * in_e_v( in_e%ind_interior )
        endif
        !
        write(*,*) "amult_ModelOperator_SP: ", self%CCii%nCol, size( in_e_v( in_e%ind_interior ) ), size( out_e_v_int ), adjoint
        !
        out_e_v( in_e%ind_interior ) = out_e_v_int
        !
        call out_e%setSV( out_e_v )
        !
    end subroutine amult_ModelOperator_SP
    !
    !> Implement the sparse matrix multiply for curl-curl operator
    !> for interior/boundary elements
    !> assume output y is already allocated
    !
    subroutine multAib_ModelOperator_SP( self, in_e, out_e )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e, out_e
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_e_v
        complex( kind=prec ), allocatable, dimension(:) :: out_e_v_int
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "multAib_ModelOperator_SP > in_e not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "multAib_ModelOperator_SP > out_e not allocated" )
        endif
        !
        in_e_v = in_e%getSV()
        !
        out_e_v = out_e%getSV()
        out_e_v_int = out_e_v( out_e%ind_interior )
        !
        write(*,*) "multAib_ModelOperator_SP: ", self%CCib%nCol, size( in_e_v( in_e%ind_boundaries ) ), size( out_e_v_int )
        !
        call RMATxCVEC( self%CCib, in_e_v( in_e%ind_boundaries ), out_e_v_int )
        !
        out_e_v( in_e%ind_interior ) = out_e_v_int
        !
        call out_e%setSV( out_e_v )
        !
    end subroutine multAib_ModelOperator_SP
    !
    !> No subroutine briefing
    subroutine div_ModelOperator_SP( self, in_e, out_phi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        type( spMatCSR_Real ) :: D
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_phi_v
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "div_ModelOperator_SP > in_e not allocated" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "div_ModelOperator_SP > out_phi not allocated" )
        endif
        !
        call RMATtrans( G, D )
        !
        in_e_v = in_e%getSV()
        !
        out_phi_v = out_phi%getSV()
        !
        write(*,*) "div_ModelOperator_SP: ", D%nCol, size( in_e_v( in_e%ind_interior ) ), size( out_phi_v )
        !
        call RMATxCVEC( D, in_e_v( in_e%ind_interior ), out_phi_v )
        !
        call deall_spMatCSR( D )
        !
        call out_phi%setSV( out_phi_v )
        !
    end subroutine div_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine divC_ModelOperator_SP( self, in_e, out_phi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_phi_v
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "divC_ModelOperator_SP > in_e not allocated" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "divC_ModelOperator_SP > out_phi not allocated" )
        endif
        !
        in_e_v = in_e%getSV()
        !
        out_phi_v = out_phi%getSV()
        !
        write(*,*) "divC_ModelOperator_SP: ", self%VDs%nCol, size( in_e_v( in_e%ind_interior ) ), size( out_phi_v )
        !
        call RMATxCVEC( self%VDs, in_e_v( in_e%ind_interior ), out_phi_v )
        !
        call out_phi%setSV( out_phi_v )
        !
    end subroutine divC_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine divCGrad_ModelOperator_SP( self, in_phi, out_phi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi, out_phi
        !
        complex( kind=prec ), allocatable, dimension(:) :: out_phi_v
        !
        if( .NOT. in_phi%is_allocated ) then
            call errStop( "divCGrad_ModelOperator_SP > in_phi not allocated" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "divCGrad_ModelOperator_SP > out_phi not allocated" )
        endif
        !
        out_phi_v = out_phi%getSV()
        !
        write(*,*) "divCGrad_ModelOperator_SP: ", self%VDsG%nCol, size( in_phi%getSV() ), size( out_phi_v )
        !
        call RMATxCVEC( self%VDsG, in_phi%getSV(), out_phi_v )
        !
        call out_phi%setSV( out_phi_v )
        !
    end subroutine divCGrad_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine grad_ModelOperator_SP( self, in_phi, out_e )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi
        class( Vector_t ), intent( inout ) :: out_e
        !
        complex( kind=prec ), allocatable, dimension(:) :: out_e_v, out_e_v_int
        !
        if( .NOT. in_phi%is_allocated ) then
            call errStop( "grad_ModelOperator_SP > in_phi not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "grad_ModelOperator_SP > out_e not allocated" )
        endif
        !
        out_e_v = out_e%getSV()
		out_e_v = C_ZERO
		out_e_v_int = out_e_v( out_e%ind_interior )
        !
        write(*,*) "grad_ModelOperator_SP: ", G%nCol, size( in_phi%getSV() ), size( out_e_v_int )
        !
        call RMATxCVEC( G, in_phi%getSV(), out_e_v_int )
        !
		out_e_v( out_e%ind_interior ) = out_e_v_int
		!
        call out_e%setSV( out_e_v )
        !
    end subroutine grad_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    !> Set sparse matrices for curl (T) and grad (G)
    !> operator topologies; these sparse matrices are stored
    !> in module spOpTopology
    !
    subroutine create_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        integer :: nInterior
        !
        self%is_allocated = .FALSE.
        !
        self%topology_sg = SpOpTopology_SG_t( self%metric%grid )
        !
        call self%topology_sg%curl( T )
        !
        call self%topology_sg%grad( G )
        !
        call boundaryIndexSP( EDGE, self%metric, self%EDGEb, self%EDGEi )
        !
        nInterior = size( self%EDGEi )
        !
        !> Find indexes (in vector of all) of boundary and interior edges
        !> allocate for diagonal part of curl-curl operator
        !> (maybe this should just be for interior edges)
        !> here for all edges
        allocate( self%VomegaMuSig( nInterior ) )
        !
        !> set a default omega
        self%omega = 0.0
        !
        self%is_allocated = .TRUE.
        !
    end subroutine create_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine deallocate_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        !> interior and edge indexes
        deallocate( self%EDGEi, self%EDGEb )
        deallocate( self%NODEi, self%NODEb )
        !
        call deall_spMatCSR( self%CC )
        call deall_spMatCSR( self%CCii )
        call deall_spMatCSR( self%CCib )
        !
        !> and the edge conductivities
        if( allocated( self%VomegaMuSig ) ) deallocate( self%VomegaMuSig )
        !
        !> and the curl and grad topology matrices
        call deall_spMatCSR( T )
        call deall_spMatCSR( G )
        !
        call deall_spMatCSR( self%VDiv )
        call deall_spMatCSR( self%VDsG )
        call deall_spMatCSR( self%VDs )
        call deall_spMatCSR( self%VDsG_L )
        call deall_spMatCSR( self%VDsG_U )
        !
        self%is_allocated = .FALSE.
        !
    end subroutine deallocate_ModelOperator_SP
    !
    !> No subroutine briefing
    subroutine print_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        !
        stop "Subroutine print not implemented for ModelOperator_SP"
        !
    end subroutine print_ModelOperator_SP
    !
end module ModelOperator_SP

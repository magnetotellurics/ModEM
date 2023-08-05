!
!> Derived class to define a ModelOperator_SP
!
module ModelOperator_SP
    !
    use ModelOperator
    use SpOpTopology_SG
    use MetricElements_CSG
    use ModelParameterCell
    !
    type, extends( ModelOperator_t ) :: ModelOperator_SP_t
        !
        type( SpOpTopology_SG_t ) :: topology_sg
        !
        integer, allocatable, dimension(:) :: EDGEi, EDGEb
        integer, allocatable, dimension(:) :: NODEi, NODEb
        !
        type( spMatCSR_Real ) :: D, Gd, CCii, CCib
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
        !write( *, * ) "Constructor ModelOperator_SP"
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
        !write( *, * ) "Destructor ModelOperator_SP_t"
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
        type( spMatCSR_Real ) :: matrix_1, Ttrans, CC
        integer :: m, n, nz
        real( kind=prec ), allocatable, dimension(:) :: Dtemp
        integer :: fid
        !
        m = T%nRow
        n = T%nCol
        nz = T%row( T%nRow + 1 ) - 1
        !
        call create_spMatCSR( m, n, nz, matrix_1 )
        call create_spMatCSR( n, m, nz, Ttrans )
        call create_spMatCSR( m, n, nz, CC )
        !
        call RMATxDIAG( T, real( self%metric%edge_length%getArray(), kind=prec ), matrix_1 )
        !
        Dtemp = ( self%metric%dual_edge_length%getArray() / self%metric%face_area%getArray() )
        !
        call DIAGxRMAT( Dtemp, matrix_1, CC )
        !
        call RMATtrans( T, Ttrans )
        !
        call RMATxRMAT( Ttrans, CC, matrix_1 )
        !
        Dtemp = self%metric%edge_length%getArray()
        !
        call DIAGxRMAT( Dtemp, matrix_1, CC )
        !
        call subMatrix_Real( CC, self%EDGEi, self%EDGEi, self%CCii )
        !
        call subMatrix_Real( CC, self%EDGEi, self%EDGEb, self%CCib )
        !
        call deall_spMatCSR( matrix_1 )
        call deall_spMatCSR( Ttrans )
        call deall_spMatCSR( CC )
        !
        self%eqset = .TRUE.
        !
        call self%divCorInit
        !
    end subroutine setEquations_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine setCond_ModelOperator_SP( self, sigma, omega_in )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: omega_in
        !
        class( Vector_t ), allocatable:: sig_temp
        complex( kind=prec ), allocatable, dimension(:) :: sig_vec_v, v_edge_v
        !
        call self%metric%createVector( real_t, EDGE, sig_temp )
        !
        call sigma%PDEmapping( sig_temp )
        !
        sig_vec_v = sig_temp%getArray()
        !
        deallocate( sig_temp )
        !
        v_edge_v = self%metric%v_edge%getArray()
        !
        self%VomegaMuSig = MU_0 * omega_in * sig_vec_v( self%EDGEi ) * v_edge_v( self%EDGEi )
        !
        self%omega = omega_in
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
        type( spMatCSR_Real ) :: matrix_1, matrix_2
        real( kind=prec ), allocatable, dimension(:) :: d, aux_vec, aux_vec_int
        integer, allocatable, dimension(:) :: allNodes
        integer :: i, m
        !
        !> #Part 1. Construction of VDiv (pre-Vds matrix) and D (div operator)
        !
        !> set indexes for interior and boundary nodes
        call boundaryIndexSP( NODE, self%metric, self%NODEb, self%NODEi )
        !
        !> matrix_1 -> transpose of topology G
        call RMATtrans( G, matrix_1 )
        !
        !> matrix_2 -> Multiply matrix_1 by dual face area
        aux_vec = self%metric%dual_face_area%getArray()
        !
        call RMATxDIAG( matrix_1, aux_vec, matrix_2 )
        !
        call deall_spMatCSR( matrix_1 )
        !
        !> Select matrix_2 interior nodes, edges to create VDiv
        call subMatrix_Real( matrix_2, self%NODEi, self%EDGEi, self%VDiv )
        !
        call deall_spMatCSR( matrix_2 )
        !
        !> v_node
        aux_vec = self%metric%v_node%getArray()
        !
        !> self%D -> Divide self%VDiv by v_node interior
        aux_vec_int = ( 1. / aux_vec( self%metric%v_node%ind_interior ) )
        !
        call DIAGxRMAT( aux_vec_int, self%VDiv, self%D )
        !
        !> #Part 2. Construction of self%Gd (grad operator)
        !
        !> turn G into actual gradient (not just topology,
        !> all nodes-> interior edges (not clear this is what we want!)
        allocate( d( G%nRow ) )
        !
        !> edge_length
        aux_vec = self%metric%edge_length%getArray()
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
        call DIAGxRMAT( d, G, matrix_1 )
        call subMatrix_Real( matrix_1, self%EDGEi, allNodes, self%Gd )
        call deall_spMatCSR( matrix_1 )
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
        type( spMatCSR_Real ) :: matrix
        complex( kind=prec ), allocatable, dimension(:) :: v_edge_v
        real( kind=prec ), allocatable, dimension(:) :: d
        integer, allocatable, dimension(:) :: allNodes
        integer :: n, i
        !
        !> Construct VDs .. multiply VDiv by Conductivity on edges; can use VomegaMuSig
        n = self%VDiv%nCol
        !
        d = ( self%VomegaMuSig / ( mu_0 * self%omega ) )
        !
        v_edge_v = self%metric%v_edge%getArray()
        !
        d = ( d / v_edge_v( self%EDGEi ) )
        !
        call RMATxDIAG( self%VDiv, d, self%VDs )
        !
        !> Construct VDsG: symmetric operator for divergence correction solver
        allocate( allNodes( self%Gd%nRow ) )
        !
        do i = 1, self%Gd%nRow
            allNodes( i ) = i
        enddo
        !
        call subMatrix_Real( self%Gd, allNodes, self%NODEi, matrix )
        !
        call RMATxRMAT( self%VDs, matrix, self%VDsG )
        !
        ! Setup preconditioner
        call Dilu_Real( self%VDsG, self%VDsG_L, self%VDsG_U )
        !
        call deall_spMatCSR( matrix )
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
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ), optional :: p_adjoint
        !
        logical :: adjoint
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_e_v
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v_int, out_e_v_int
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "amult_ModelOperator_SP > in_e not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "amult_ModelOperator_SP > out_e not allocated" )
        endif
        !
        in_e_v = in_e%getArray()
        in_e_v_int = in_e_v( in_e%ind_interior )
        !
        out_e_v = out_e%getArray()
        out_e_v_int = out_e_v( out_e%ind_interior )
        !
        !write(*,*) "amult_ModelOperator_SP: ", self%CCii%nCol, self%CCii%nRow, size( in_e_v_int ), size( out_e_v_int ), adjoint
        !
        call RMATxCVEC( self%CCii, in_e_v_int, out_e_v_int )
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
        out_e_v( in_e%ind_interior ) = out_e_v_int
        !
        call out_e%setArray( out_e_v )
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
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_e_v
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v_bry, out_e_v_int
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "multAib_ModelOperator_SP > in_e not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "multAib_ModelOperator_SP > out_e not allocated" )
        endif
        !
        in_e_v = in_e%getArray()
        in_e_v_bry = in_e_v( in_e%ind_boundary )
        !
        out_e_v = out_e%getArray()
        out_e_v_int = out_e_v( out_e%ind_interior )
        !
        !write(*,*) "multAib_ModelOperator_SP: ", self%CCib%nCol, self%CCib%nRow, size( in_e_v_bry ), size( out_e_v_int )
        !
        call RMATxCVEC( self%CCib, in_e_v_bry, out_e_v_int )
        !
        out_e_v( out_e%ind_interior ) = out_e_v_int
        !
        call out_e%setArray( out_e_v )
        !
    end subroutine multAib_ModelOperator_SP
    !
    !> No subroutine briefing
    subroutine div_ModelOperator_SP( self, in_e, out_phi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_phi_v
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v_int, out_phi_v_int
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "div_ModelOperator_SP > in_e not allocated" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "div_ModelOperator_SP > out_phi not allocated" )
        endif
        !
        in_e_v = in_e%getArray()
        in_e_v_int = in_e_v( in_e%ind_interior )
        !
        out_phi_v = out_phi%getArray()
        out_phi_v_int = out_phi_v( out_phi%ind_interior )
        !
        !write(*,*) "div_ModelOperator_SP: ", self%D%nCol, self%D%nRow, size( in_e_v_int ), size( out_phi_v_int )
        !
        call RMATxCVEC( self%D, in_e_v_int, out_phi_v_int )
        !
        out_phi_v( out_phi%ind_interior ) = out_phi_v_int
        !
        call out_phi%setArray( out_phi_v_int )
        !
    end subroutine div_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine divC_ModelOperator_SP( self, in_e, out_phi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v, out_phi_v
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v_int, out_phi_v_int
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "divC_ModelOperator_SP > in_e not allocated" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "divC_ModelOperator_SP > out_phi not allocated" )
        endif
        !
        in_e_v = in_e%getArray()
        in_e_v_int = in_e_v( in_e%ind_interior )
        !
        out_phi_v = out_phi%getArray()
        out_phi_v_int = out_phi_v( out_phi%ind_interior )
        !
        !write(*,*) "divC_ModelOperator_SP: ", self%VDs%nCol, self%VDs%nRow, size( in_e_v_int ), size( out_phi_v_int )
        !
        call RMATxCVEC( self%VDs, in_e_v_int, out_phi_v_int )
        !
        out_phi_v( out_phi%ind_interior ) = out_phi_v_int
        !
        call out_phi%setArray( out_phi_v )
        !
    end subroutine divC_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine divCGrad_ModelOperator_SP( self, in_phi, out_phi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( in ) :: in_phi
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_phi_v, out_phi_v
        complex( kind=prec ), allocatable, dimension(:) :: in_phi_v_int, out_phi_v_int
        !
        if( .NOT. in_phi%is_allocated ) then
            call errStop( "divCGrad_ModelOperator_SP > in_phi not allocated" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "divCGrad_ModelOperator_SP > out_phi not allocated" )
        endif
        !
        in_phi_v = in_phi%getArray()
        in_phi_v_int = in_phi_v( in_phi%ind_interior )
        !
        out_phi_v = out_phi%getArray()
        out_phi_v_int = out_phi_v( out_phi%ind_interior )
        !
        !write(*,*) "divCGrad_ModelOperator_SP: ", self%VDsG%nCol, self%VDsG%nRow, size( in_phi_v_int ), size( out_phi_v_int )
        !
        call RMATxCVEC( self%VDsG, in_phi_v_int, out_phi_v_int )
        !
        out_phi_v( out_phi%ind_interior ) = out_phi_v_int
        !
        call out_phi%setArray( out_phi_v )
        !
    end subroutine divCGrad_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine grad_ModelOperator_SP( self, in_phi, out_e )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( in ) :: in_phi
        class( Vector_t ), intent( inout ) :: out_e
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_phi_v, out_e_v
        complex( kind=prec ), allocatable, dimension(:) :: in_phi_v_int, out_e_v_int
        !
        if( .NOT. in_phi%is_allocated ) then
            call errStop( "grad_ModelOperator_SP > in_phi not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "grad_ModelOperator_SP > out_e not allocated" )
        endif
        !
        in_phi_v = in_phi%getArray()
        !
        out_e_v = out_e%getArray()
        out_e_v_int = out_e_v( out_e%ind_interior )
        !
        !write(*,*) "grad_ModelOperator_SP: ", self%Gd%nCol, self%Gd%nRow, size( in_phi_v ), size( out_e_v_int )
        !
        call RMATxCVEC( self%Gd, in_phi_v, out_e_v_int )
        !
        out_e_v( out_e%ind_interior ) = out_e_v_int
        !
        call out_e%setArray( out_e_v )
        !
    end subroutine grad_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    !> Set sparse matrices for curl (T) and grad (G)
    !> operator topologies; these sparse matrices are stored
    !> in module spOpTopology
    !
    !> Find indexes (in vector of all) of boundary and interior edges
    !> allocate for diagonal part of curl-curl operator
    !> (maybe this should just be for interior edges)
    !> here for all edges
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
        call deall_spMatCSR( self%Gd )
        call deall_spMatCSR( self%D )
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

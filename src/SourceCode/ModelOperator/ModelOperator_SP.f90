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
        type( spMatCSR_Real ) :: VDiv        ! div : edges->nodes (interior only)
        type( spMatCSR_Real ) :: VDsG        ! operator for div correction
        type( spMatCSR_Real ) :: VDs         ! divergence of current operator
        !
        type( spMatCSR_Real ) :: VDsG_L, VDsG_U
        !
        contains
            !
            final :: ModelOperator_SP_dtor
            !
            procedure, public :: setEquations => setEquations_ModelOperator_SP
            procedure, public :: setCond => setCond_ModelOperator_SP
            procedure, public :: amult => amultModelOperator_SP
            procedure, public :: multAib => multAib_ModelOperator_SP
            procedure, public :: multCurlT => multCurlT_ModelOperator_SP
            procedure, public :: divCorSetUp => divCorSetUp_ModelOperator_SP
            !
            procedure, public :: divCorInit => divCorInitModelOperator_SP
            !
            procedure :: divCGrad => divCGrad_ModelOperator_SP
            procedure :: divC => divC_ModelOperator_SP
            procedure :: grad => grad_ModelOperator_SP
            procedure :: div => div_ModelOperator_SP
            !
            procedure :: create => create_ModelOperator_SP 
            procedure :: dealloc => deallocate_ModelOperator_SP
            !
            procedure, public :: print => print_ModelOperator_SP
            !
            procedure, private :: updateOmegaMuSig
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
    !> No subroutine briefing
    subroutine create_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        integer :: nInterior
        !
        self%is_allocated = .FALSE.
        !
        !> Set sparse matrices for curl (T) and grad (G)
        !> operator topologies; these sparse matrices are stored
        !> in module spOpTopology
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
    !> using existing curl operator, create sparse matrix CC
    !> Note: this is the symmetric form, multiplied by edge volume elements
    !
    subroutine setEquations_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        type( spMatCSR_Real ) :: Temp, Ttrans
        integer :: m, n, nz
        real( kind=prec ), allocatable, dimension(:) :: Dtemp
        integer :: fid
        !
        m = T%nRow
        n = T%nCol
        nz = T%row( T%nRow + 1 ) - 1
        !
        !allocate( Dtemp( m ) )
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
        call self%divCorSetUp
        !
    end subroutine setEquations_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine setCond_ModelOperator_SP( self, sigma )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( inout ) :: sigma
        !
        integer :: i
        class( Scalar_t ), allocatable :: temp_cell_cond
        class( ModelParameter_t ), allocatable :: model
        class( Vector_t ), allocatable :: sig_temp
        complex( kind=prec ), allocatable, dimension(:) :: v_edge_sv, sig_temp_sv, cVomegaMuSig
        !
        self%omega = 1.0
        !
        !> vEdge
        v_edge_sv = self%metric%v_edge%getSV()
        !
        !> Sigma
        call self%metric%createVector( real_t, EDGE, sig_temp )
        !
        call sigma%PDEmapping( sig_temp )
        !
        sig_temp_sv = sig_temp%getSV()
        !
        self%VomegaMuSig = v_edge_sv( self%EDGEi ) * sig_temp_sv( self%EDGEi ) * mu_0 * self%omega
        !
        cVomegaMuSig = cmplx( self%VomegaMuSig, 0.0, kind=prec )
        !
        allocate( model, source = sigma )
        model = sigma
        !
        do i = 1, model%anisotropic_level
            !
            call self%metric%createScalar( real_t, CELL_EARTH, temp_cell_cond )
            !
            call temp_cell_cond%setSV( cVomegaMuSig )
            !
            call model%setCond( temp_cell_cond, i )
            !
        enddo
        !
        deallocate( model )
        !
    end subroutine setCond_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine updateOmegaMuSig( self, in_omega, model_param )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        real( kind=prec ), intent ( in ) :: in_omega
        class( ModelParameter_t ), intent( inout ), optional :: model_param
        !
        class( Vector_t ), allocatable :: sig_temp
        complex( kind=prec ), allocatable, dimension(:) :: sig_temp_sv, v_edge_sv
        !
        if( present( model_param ) ) then
            !
            !> Sigma
            call self%metric%createVector( real_t, EDGE, sig_temp )
            !
            call model_param%PDEmapping( sig_temp )
            !
            sig_temp_sv = sig_temp%getSV()
            !
            !> vEdge
            v_edge_sv = self%metric%v_edge%getSV()
            !
            self%VomegaMuSig = MU_0 * in_omega * sig_temp_sv( self%EDGEi ) * v_edge_sv( self%EDGEi )
            !
            self%omega = in_omega
            !
        else
            if( self%omega .gt. 0 ) then
                self%VomegaMuSig = ( self%VomegaMuSig / self%omega )
            endif
            !
            self%VomegaMuSig = ( self%VomegaMuSig * in_omega )
            !
            self%omega = in_omega
            !
        endif
        !
    end subroutine updateOmegaMuSig
    !
    !> To complete setup conductivity is required
    !> DivCorInit has to be called before this routine
    !
    subroutine divCorInitModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        type( spMatCSR_Real ) :: Temp, Temp2
        real( kind=prec ), allocatable, dimension(:) :: d, aux_vec
        integer, allocatable, dimension(:) :: allNodes
        integer :: i, m
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
    end subroutine divCorInitModelOperator_SP
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
    subroutine amultModelOperator_SP( self, omega, in_e, out_e, p_adjoint )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        real( kind=prec ), intent( in ), optional :: omega
        class( Vector_t ), intent( inout ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ), optional :: p_adjoint
        !
        logical :: adjoint
        complex( kind=prec ), allocatable, dimension(:) :: array_inE, array_outE
        complex( kind=prec ), allocatable, dimension(:) :: array_inE_int, array_result
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "amultModelOperator_SP > in_e not allocated" )
        endif
        !
        array_inE = in_e%getSV()
        array_outE = out_e%getSV()
        !
        array_inE_int = array_inE( in_e%ind_interior )
        !
        array_result = array_inE_int
        array_result = C_ZERO
        !
        call RMATxCVEC( self%CCii, array_inE_int, array_result )
        !
        array_outE( in_e%ind_interior ) = array_result
        !
        if( present( p_adjoint ) ) then
            adjoint = p_adjoint
        else
            adjoint = .FALSE.
        endif
        !
        if( adjoint ) then
            array_outE = array_outE - ONE_I * ISIGN * self%VomegaMuSig * array_inE
        else
            array_outE = array_outE + ONE_I * ISIGN * self%VomegaMuSig * array_inE
        endif
        !
        call out_e%setSV( array_outE )
        !
    end subroutine amultModelOperator_SP
    !
    !> Implement the sparse matrix multiply for curl-curl operator
    !> for interior/boundary elements
    !> assume output y is already allocated
    !
    subroutine multAib_ModelOperator_SP( self, in_e, out_e )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        !
        complex( kind=prec ), allocatable, dimension(:) :: array_inE, array_outE
        complex( kind=prec ), allocatable, dimension(:) :: array_inE_bdry, array_outE_int
        !
        if( .NOT. in_e%is_allocated ) then
            stop "Error: amultModelOperator_SP > in_e not allocated"
        endif
        !
        array_inE = in_e%getSV()
        array_inE_bdry = array_inE( in_e%ind_boundaries )
        !
        array_outE = out_e%getSV()
        array_outE_int = array_outE( in_e%ind_interior )
        array_outE_int = C_ZERO
        !
        call RMATxCVEC( self%CCib, array_inE_bdry, array_outE_int )
        !
        array_outE( in_e%ind_interior ) = array_outE_int
        !
        call out_e%setSV( array_outE )
        !
    end subroutine multAib_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine multCurlT_ModelOperator_SP( self, in_e, out_e )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Vector_t ), allocatable, intent( inout ) :: out_e
        !
        integer :: ix, iy, iz
        complex( kind=prec ), allocatable, dimension(:, :, :) :: in_e_x, in_e_y, in_e_z
        complex( kind=prec ), allocatable, dimension(:, :, :) :: out_e_x, out_e_y, out_e_z
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "multCurlT_ModelOperator_SP > output vector not allocated" )
        endif
        !
        out_e_x = out_e%getX()
        out_e_y = out_e%getY()
        out_e_z = out_e%getZ()
        !
        call in_e%div( self%Metric%face_area )
        !
        in_e_x = in_e%getX()
        in_e_y = in_e%getY()
        in_e_z = in_e%getZ()
        !
        !> Ex
        do iy = 2, in_e%Ny
            do iz = 2, in_e%Nz
                out_e_x(:, iy, iz) = (in_e_z(:, iy, iz) - &
                in_e_z(:, iy - 1, iz)) - &
                (in_e_y(:, iy, iz) - in_e_y(:, iy, iz - 1))
            enddo
        enddo
        !
        !> Ey
        do iz = 2, in_e%Nz
            do ix = 2, in_e%Nx
                out_e_y(ix, :, iz) = (in_e_x(ix, :, iz) - &
                in_e_x(ix, :, iz - 1)) - &
                (in_e_z(ix, :, iz) - in_e_z(ix - 1, :, iz))
            enddo
        enddo
        !
        !> Ez
        do ix = 2, in_e%Nx
            do iy = 2, in_e%Ny
                out_e_z(ix,iy,:) = (in_e_y(ix, iy, :) - &
                in_e_y(ix - 1, iy, :)) - &
                (in_e_x(ix, iy, :) - in_e_x(ix, iy - 1, :))
            enddo
        enddo
        !
        call out_e%setX( out_e_x )
        call out_e%setY( out_e_y )
        call out_e%setZ( out_e_z )
        !
        call out_e%mult( self%metric%edge_length )
        !
        call out_e%switchStoreState( singleton )
        !
    end subroutine multCurlT_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine divCGrad_ModelOperator_SP( self, in_phi, out_phi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi, out_phi
        !
        complex( kind=prec ), allocatable, dimension(:) :: array_inPhi, array_outPhi
        !
        array_inPhi = in_phi%getSV()
        array_outPhi = out_phi%getSV()
        array_outPhi = C_ZERO
        !
        call RMATxCVEC( self%VDsG, array_inPhi, array_outPhi )
        !
        call out_phi%setSV( array_outPhi )
        !
    end subroutine divCGrad_ModelOperator_SP
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
        complex( kind=prec ), allocatable, dimension(:) :: array_inE, array_outPhi
        complex( kind=prec ), allocatable, dimension(:) :: array_inE_int, array_result
        !
        array_inE = in_e%getSV()
        array_inE_int = array_inE( in_e%ind_interior )
        !
        array_outPhi = out_phi%getSV()
        array_outPhi = C_ZERO
        array_result = array_inE_int
        array_result = C_ZERO
        !
        call RMATxCVEC( self%VDs, array_inE_int, array_result )
        !
        array_outPhi( out_phi%ind_interior ) = array_result
        !
        call out_phi%setSV( array_outPhi )
        !
    end subroutine divC_ModelOperator_SP
    !
    !> No subroutine briefing
    subroutine grad_ModelOperator_SP( self, in_phi, out_e )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi
        class( Vector_t ), intent( inout ) :: out_e
        !
        complex( kind=prec ), allocatable, dimension(:) :: array_inPhi, array_outE
        complex( kind=prec ), allocatable, dimension(:) :: array_outE_int
        !
        array_inPhi = in_phi%getSV()
        !
        array_outE = out_e%getSV()
        array_outE_int = array_outE( out_e%ind_interior )
        !
        call RMATxCVEC( G, array_inPhi, array_outE_int )
        !
        array_outE( out_e%ind_interior ) = array_outE_int
        !
        call out_e%setSV( array_outE )
        !
    end subroutine grad_ModelOperator_SP
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
        complex( kind=prec ), allocatable, dimension(:) :: array_inE, array_outPhi
        !
        call RMATtrans( G, D )
        !
        array_inE = in_e%getSV()
        !
        array_outPhi = array_inE
        array_outPhi = C_ZERO
        !
        call RMATxCVEC( D, array_inE, array_outPhi )
        !
        call deall_spMatCSR( D )
        !
        call out_phi%setSV( array_outPhi )
        !
    end subroutine div_ModelOperator_SP
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

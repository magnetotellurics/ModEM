!
!> Derived class to define a ModelOperator
!> with basic operations for Sparse Matrices
!
module ModelOperator_SP
    !
    use ModelOperator
    use SpOptopology
    use SpOpTopology_MR
    use MetricElements_SG
    use MetricElements_MR
    !
    type, abstract, extends( ModelOperator_t ) :: ModelOperator_SP_t
        !
        class( SpOpTopology_t ), allocatable :: topology
        !
        type( spMatCSR_Real ) :: D, Gd, TCC, CCii, CCib
        !
        type( spMatCSR_Real ) :: VDiv ! div : edges->nodes (interior only)
        type( spMatCSR_Real ) :: VDsG ! operator for div correction
        type( spMatCSR_Real ) :: Ds ! divergence of current operator
        !
        type( spMatCSR_Real ) :: VDsG_L, VDsG_U
        !
        real( kind=prec ) :: omega
        !
        real( kind=prec ), allocatable, dimension(:) :: VomegaMuSig
        !
        contains
            !
            !> Setup
            procedure, public :: create => create_ModelOperator_SP
            !
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
            procedure, public :: multCurlT => multCurlT_ModelOperator_SP
            !
            procedure, public :: div => div_ModelOperator_SP
            procedure, public :: divC => divC_ModelOperator_SP
            procedure, public :: divCGrad => divCGrad_ModelOperator_SP
            !
            procedure, public :: grad => grad_ModelOperator_SP
            !
            !> Alloc/Dealloc
            procedure, public :: dealloc => deallocate_ModelOperator_SP
            !
            !> Miscellaneous
            procedure, public :: print => print_ModelOperator_SP
            !
    end type ModelOperator_SP_t
    !
contains
    !
    !> No subroutine briefing
    !
    !> Set sparse matrices for curl (topology%T) and grad (topology%G)
    !> operator topologies; these sparse matrices are stored
    !> in module SpOpTopology
    !
    !> Find indexes (in vector of all) of boundary and interior edges
    !> allocate for diagonal part of curl-curl operator
    !> (maybe this should just be for interior edges)
    !> here for all edges
    !
    subroutine create_ModelOperator_SP( self, grid )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        class( Grid_t ), target, intent( in ) :: grid
        !
        self%is_allocated = .FALSE.
        !
        !> Instantiation of the specific object MetricElements
        select type( grid )
            !
            class is ( Grid3D_SG_t )
                !
                allocate( self%metric, source = MetricElements_SG_t( grid ) )
                !
                allocate( self%topology, source = SpOpTopology_SG_t( grid ) )
                !
            class is ( Grid3D_MR_t )
                !
                allocate( self%metric, source = MetricElements_MR_t( grid ) )
                !
                allocate( self%topology, source = SpOpTopology_MR_t( grid ) )
                !
            class default
                !
                call errStop( "create_ModelOperator_SP > grid_format not provided" )
                !
        end select
        !
        call self%topology%curl( self%topology%T )
        !
        !call writeIJS_Matrix( self%topology%T, 6666 )
        !
        call self%topology%grad( self%topology%G )
        !
        !call writeIJS_Matrix( self%topology%G, 6667 )
        !
        allocate( self%VomegaMuSig( size( self%metric%grid%EDGEi ) ) )
        !
        !> set a default omega
        self%omega = 0.0
        !
        self%is_allocated = .TRUE.
        !
    end subroutine create_ModelOperator_SP
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
        type( spMatCSR_Real ) :: temp_matrix, Ttrans, CC
        integer :: m, n, nz
        real( kind=prec ), allocatable, dimension(:) :: temp_array
        integer :: fid
        !
        m = self%topology%T%nRow
        n = self%topology%T%nCol
        nz = self%topology%T%row( self%topology%T%nRow + 1 ) - 1
        !
        call create_spMatCSR( m, n, nz, temp_matrix )
        call create_spMatCSR( n, m, nz, Ttrans )
        call create_spMatCSR( m, n, nz, CC )
        !!
        call RMATxDIAG( self%topology%T, real( self%metric%edge_length%getArray(), kind=prec ), temp_matrix )
        !
        !> Create TCC for for multCurlT
        temp_array = 1. / self%metric%face_area%getArray()
        !
        call DIAGxRMAT( temp_array, temp_matrix, CC )
        !
        call RMATtrans( CC, self%TCC )
        !
        !> Create CCii for amult and CCib for multAib
        temp_array = self%metric%dual_edge_length%getArray() / self%metric%face_area%getArray()
        !
        call DIAGxRMAT( temp_array, temp_matrix, CC )
        !
        call RMATtrans( self%topology%T, Ttrans )
        !
        call RMATxRMAT( Ttrans, CC, temp_matrix )
        !
        temp_array = self%metric%edge_length%getArray()
        !
        call DIAGxRMAT( temp_array, temp_matrix, CC )
        !
        call subMatrix_Real( CC, self%metric%grid%EDGEi, self%metric%grid%EDGEi, self%CCii )
        !
        call subMatrix_Real( CC, self%metric%grid%EDGEi, self%metric%grid%EDGEb, self%CCib )
        !
        call deall_spMatCSR( temp_matrix )
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
        class( Grid_t ), allocatable :: temp_grid_mr
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
        self%VomegaMuSig = MU_0 * omega_in * sig_vec_v( self%metric%grid%EDGEi ) * v_edge_v( self%metric%grid%EDGEi )
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
        integer, allocatable, dimension(:) :: all_nodes
        integer :: i, m
        !
        !> #Part 1. Construction of VDiv (pre-Vds matrix) and D (div operator)
        !
        !> matrix_1 -> transpose of self%topology self%topology%G
        call RMATtrans( self%topology%G, matrix_1 )
        !
        !> matrix_2 -> Multiply matrix_1 by dual face area
        aux_vec = self%metric%dual_face_area%getArray()
        !
        call RMATxDIAG( matrix_1, aux_vec, matrix_2 )
        !
        call deall_spMatCSR( matrix_1 )
        !
        !> Select matrix_2 interior nodes, edges to create VDiv
        call subMatrix_Real( matrix_2, self%metric%grid%NODEi, self%metric%grid%EDGEi, self%VDiv )
        !
        call deall_spMatCSR( matrix_2 )
        !
        !> v_node
        aux_vec = self%metric%v_node%getArray()
        !
        !> self%D -> Divide self%VDiv by v_node interior
        aux_vec_int = ( 1. / aux_vec( self%metric%v_node%indInterior() ) )
        !
        call DIAGxRMAT( aux_vec_int, self%VDiv, self%D )
        !
        !> #Part 2. Construction of self%Gd (grad operator)
        !
        !> turn self%topology%G into actual gradient (not just self%topology,
        !> all nodes-> interior edges (not clear this is what we want!)
        allocate( d( self%topology%G%nRow ) )
        !
        !> edge_length
        aux_vec = self%metric%edge_length%getArray()
        !
        do i = 1, self%topology%G%nRow
            d(i) = 1. / aux_vec(i)
        enddo
        !
        allocate( all_nodes( self%topology%G%nCol ) )
        !
        do i=1,self%topology%G%nCol
            all_nodes(i) = i
        enddo
        !
        call DIAGxRMAT( d, self%topology%G, matrix_1 )
        call subMatrix_Real( matrix_1, self%metric%grid%EDGEi, all_nodes, self%Gd )
        call deall_spMatCSR( matrix_1 )
        !
        deallocate( all_nodes, d )
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
        type( spMatCSR_Real ) :: matrix, VDs
        complex( kind=prec ), allocatable, dimension(:) :: v_edge_v
        real( kind=prec ), allocatable, dimension(:) :: d
        integer, allocatable, dimension(:) :: all_nodes
        integer :: i
        !
        d = ( self%VomegaMuSig / ( mu_0 * self%omega ) )
        !
        v_edge_v = self%metric%v_edge%getArray()
        !
        d = ( d / v_edge_v( self%metric%grid%EDGEi ) )
        !
        call RMATxDIAG( self%D, d, self%Ds )
        !
        !> VDs -> multiply self%VDiv by Conductivity on edges using VomegaMuSig
        call RMATxDIAG( self%VDiv, d, VDs )
        !
        !> Construct VDsG: symmetric operator for divergence correction solver
        allocate( all_nodes( self%Gd%nRow ) )
        !
        do i = 1, self%Gd%nRow
            all_nodes( i ) = i
        enddo
        !
        call subMatrix_Real( self%Gd, all_nodes, self%metric%grid%NODEi, matrix )
        !
        call RMATxRMAT( VDs, matrix, self%VDsG )
        !
        call deall_spMatCSR( VDs )
        call deall_spMatCSR( matrix )
        !
        ! Setup preconditioner matrices: self%VDsG_L and self%VDsG_U
        call Dilu_Real( self%VDsG, self%VDsG_L, self%VDsG_U )
        !
        deallocate( all_nodes )
        !
    end subroutine divCorSetUp_ModelOperator_SP
    !
    !> Implement the sparse matrix multiply for curl-curl operator
    !> for interior elements
    !> assume output y is already allocated
    !
    subroutine amult_ModelOperator_SP( self, in_e, out_e, omega, adjoint )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        real( kind=prec ), intent( in ) :: omega
        logical, intent( in ) :: adjoint
        !
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
        in_e_v_int = in_e_v( in_e%indInterior() )
        !
        out_e_v = out_e%getArray()
        out_e_v_int = out_e_v( out_e%indInterior() )
        !
        !write(*,*) "amult_ModelOperator_SP: ", self%CCii%nCol, self%CCii%nRow, size( in_e_v_int ), size( out_e_v_int ), adjoint
        !
        call RMATxCVEC( self%CCii, in_e_v_int, out_e_v_int )
        !
        if( adjoint ) then
            out_e_v_int = out_e_v_int - ONE_I * ISIGN * self%VomegaMuSig * in_e_v( in_e%indInterior() )
        else
            out_e_v_int = out_e_v_int + ONE_I * ISIGN * self%VomegaMuSig * in_e_v( in_e%indInterior() )
        endif
        !
        out_e_v( in_e%indInterior() ) = out_e_v_int
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
        in_e_v_bry = in_e_v( in_e%indBoundary() )
        !
        out_e_v = out_e%getArray()
        out_e_v_int = out_e_v( out_e%indInterior() )
        !
        !write(*,*) "multAib_ModelOperator_SP: ", self%CCib%nCol, self%CCib%nRow, size( in_e_v_bry ), size( out_e_v_int )
        !
        call RMATxCVEC( self%CCib, in_e_v_bry, out_e_v_int )
        !
        out_e_v( out_e%indInterior() ) = out_e_v_int
        !
        call out_e%setArray( out_e_v )
        !
    end subroutine multAib_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine multCurlT_ModelOperator_SP( self, in_b, out_e )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_b
        class( Vector_t ), allocatable, intent( out ) :: out_e
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_b_v, out_e_v
        !
        if( .NOT. in_b%is_allocated ) then
            call errStop( "multAib_ModelOperator_SP > in_b not allocated" )
        endif
        !
        call self%metric%createVector( complex_t, EDGE, out_e )
        call out_e%zeros
        !
        in_b_v = in_b%getArray()
        !
        out_e_v = out_e%getArray()
        !
        !write(*,*) "multCurlT_ModelOperator_SP: ", self%TCC%nCol, self%TCC%nRow, size( in_b_v ), size( out_e_v )
        !
        call RMATxCVEC( self%TCC, in_b_v, out_e_v )
        !
        call out_e%setArray( out_e_v )
        !
    end subroutine multCurlT_ModelOperator_SP
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
        complex( kind=prec ), allocatable, dimension(:) :: in_e_v_int, out_phi_int
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
        in_e_v_int = in_e_v( in_e%indInterior() )
        !
        out_phi_v = out_phi%getArray()
        out_phi_int = out_phi_v( out_phi%indInterior() )
        !
        !write(*,*) "div_ModelOperator_SP: ", self%D%nCol, self%D%nRow, size( in_e_v_int ), size( out_phi_int )
        !
        call RMATxCVEC( self%D, in_e_v_int, out_phi_int )
        !
        out_phi_v( out_phi%indInterior() ) = out_phi_int
        !
        call out_phi%setArray( out_phi_v )
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
        in_e_v_int = in_e_v( in_e%indInterior() )
        !
        out_phi_v = out_phi%getArray()
        out_phi_v_int = out_phi_v( out_phi%indInterior() )
        !
        !write(*,*) "divC_ModelOperator_SP: ", self%Ds%nCol, self%Ds%nRow, size( in_e_v_int ), size( out_phi_v_int )
        !
        call RMATxCVEC( self%Ds, in_e_v_int, out_phi_v_int )
        !
        out_phi_v( out_phi%indInterior() ) = out_phi_v_int
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
        in_phi_v_int = in_phi_v( in_phi%indInterior() )
        !
        out_phi_v = out_phi%getArray()
        out_phi_v_int = out_phi_v( out_phi%indInterior() )
        !
        !write(*,*) "divCGrad_ModelOperator_SP: ", self%VDsG%nCol, self%VDsG%nRow, size( in_phi_v_int ), size( out_phi_v_int )
        !
        call RMATxCVEC( self%VDsG, in_phi_v_int, out_phi_v_int )
        !
        out_phi_v( out_phi%indInterior() ) = out_phi_v_int
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
        out_e_v_int = out_e_v( out_e%indInterior() )
        !
        !write(*,*) "grad_ModelOperator_SP: ", self%Gd%nCol, self%Gd%nRow, size( in_phi_v ), size( out_e_v_int )
        !
        call RMATxCVEC( self%Gd, in_phi_v, out_e_v_int )
        !
        out_e_v( out_e%indInterior() ) = out_e_v_int
        !
        call out_e%setArray( out_e_v )
        !
    end subroutine grad_ModelOperator_SP
    !
    !> No subroutine briefing
    !
    subroutine deallocate_ModelOperator_SP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        call self%baseDealloc
        !
        call deall_spMatCSR( self%CCii )
        call deall_spMatCSR( self%CCib )
        !
        !> and the edge conductivities
        if( allocated( self%VomegaMuSig ) ) deallocate( self%VomegaMuSig )
        !
        !> and the curl and grad self%topology matrices
        call deall_spMatCSR( self%topology%T )
        call deall_spMatCSR( self%topology%G )
        !
        if( allocated( self%topology ) ) deallocate( self%topology )
        !
        call deall_spMatCSR( self%TCC )
        call deall_spMatCSR( self%Gd )
        call deall_spMatCSR( self%D )
        call deall_spMatCSR( self%VDiv )
        call deall_spMatCSR( self%VDsG )
        call deall_spMatCSR( self%Ds )
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
!
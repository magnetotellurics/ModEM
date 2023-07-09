!
!> Derived class to define a ModelOperator_SP
!
module ModelOperator_SP
    !
    use ModelOperator
    use spOpTopology_SG
    use MetricElements_CSG
    use ModelParameterCell_SG
    !
    type, extends( ModelOperator_t ) :: ModelOperator_SP_t
        !
        type( spOpTopology_SG_t ) :: topology_sg
        !
        class( Vector_t ), allocatable :: sigma_C
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
            procedure, public :: setEquations => setEquationsModelOperatorSP
            procedure, public :: setCond => setCondModelOperatorSP
            procedure, public :: amult => amultModelOperatorSP
            procedure, public :: multAib => multAibModelOperatorSP
            procedure, public :: multCurlT => multCurlTModelOperatorSP
            procedure, public :: divCorSetUp => divCorSetUpModelOperatorSP
            !
            procedure, public :: divCorInit => divCorInitModelOperatorSP
            !
            procedure, private :: updateOmegaMuSig
            !
            procedure :: divCgrad => divCgradModelOperatorSP
            procedure :: divC => divCModelOperatorSP
            procedure :: grad => gradModelOperatorSP
            procedure :: div => divModelOperatorSP
            !
            procedure :: create => createModelOperatorSP 
            procedure :: deallocate => deallocateModelOperatorSP
            !
            procedure, public :: print => printModelOperatorSP
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
        class( Grid3D_SG_t ), target, intent( in ) :: grid
        !
        type( ModelOperator_SP_t ) :: self
        !
        !write( *, * ) "Constructor ModelOperator_SP"
        !
        call self%init
        !
        !> Instantiation of the specific object MetricElements
        allocate( self%metric, source = MetricElements_CSG_t( grid ) )
        !
        call self%create
        !
    end function ModelOperator_SP_ctor
    !
    !> No subroutine briefing
    subroutine createModelOperatorSP( self )
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
        self%topology_sg = spOpTopology_SG_t( self%metric%grid )
        !
        call self%topology_sg%curl( T )
        !
        call self%topology_sg%grad( G )
        !
        call boundaryIndexSP( EDGE, self%metric%grid, self%EDGEb, self%EDGEi )
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
    end subroutine createModelOperatorSP
    !
    !> ModelOperator_SP destructor
    subroutine ModelOperator_SP_dtor( self )
        implicit none
        !
        type( ModelOperator_SP_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelOperator_SP_t"
        !
        call self%dealloc
        !
        call self%deallocate()
        !
    end subroutine ModelOperator_SP_dtor
    !
    !> No subroutine briefing
    !
    subroutine deallocateModelOperatorSP( self )
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
    end subroutine deallocateModelOperatorSP
    !
    !> No subroutine briefing
    !> using existing curl operator, create sparse matrix CC
    !> Note: this is the symmetric form, multiplied by edge volume elements
    !
    subroutine setEquationsModelOperatorSP( self )
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
        call RMATxDIAG( T, real( self%metric%Edgelength%getArray(), kind=prec ), Temp )
        !
        Dtemp = ( self%metric%DualEdgelength%getArray() / self%metric%FaceArea%getArray() )
        !
        call DIAGxRMAT( Dtemp, Temp, self%CC )
        !
        call RMATtrans( T, Ttrans )
        !
        call RMATxRMAT( Ttrans, self%CC, Temp )
        !
        call DIAGxRMAT( real( self%metric%Edgelength%getArray(), kind=prec ), Temp, self%CC )
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
    end subroutine setEquationsModelOperatorSP
    !
    !> No subroutine briefing
    !
    subroutine setCondModelOperatorSP( self, sigma )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: sigma
        !
        integer :: i
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: cell_cond
        class( ModelParameter_t ), allocatable :: model
        type( rVector3D_SG_t ) :: temp_vec, sig_temp_vec
        !
        sig_temp_vec = rVector3D_SG_t( self%metric%grid, EDGE )
        !
        !> ON -> call ModelParamToEdge( sigma, sig_temp_vec )
        call sigma%PDEmapping( sig_temp_vec )
        !
        call sig_temp_vec%switchStoreState()
        !
        self%omega = 1.0
        !
        temp_vec = self%metric%VEdge
        !
        call temp_vec%switchStoreState()
        !
        self%VomegaMuSig = temp_vec%sv( self%EDGEi ) * sig_temp_vec%sv( self%EDGEi ) * mu_0 * self%omega
        !
        !> TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
        !> set static array for cell conductivities
        !> this stores conductivity values in a module structure
        !> that is readily accesible to boundary condition routines
        !> rvector sigma_C is created if it is not yet allocated
        !
        allocate( model, source = sigma )
        model = sigma
        !
        allocate( cell_cond( model%anisotropic_level ) )
        !
        do i = 1, model%anisotropic_level
            !
            cell_cond(i) = rScalar3D_SG_t( self%metric%grid, CELL_EARTH )
            !
            call cell_cond(i)%setArray( cmplx( self%VomegaMuSig, 0.0, kind=prec ) )
            !
        enddo
        !
        call model%setCond( cell_cond )
        !
        call sigma%dPDEmapping( model, self%sigma_C )
        !
        deallocate( model )
        !
    end subroutine setCondModelOperatorSP
    !
    !> No subroutine briefing
    !
    subroutine updateOmegaMuSig( self, inOmega, model_param )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        real( kind=prec ), intent ( in ) :: inOmega
        class( ModelParameter_t ), intent( in ), optional :: model_param
        !
        type( rVector3D_SG_t ) :: sigTemp
        complex( kind=prec ), allocatable, dimension(:) :: sigVec, temp_vec
        !
        if( present( model_param ) ) then
            !
            sigTemp = rVector3D_SG_t( self%metric%grid, EDGE )
            !
            call model_param%PDEmapping( sigTemp )
            !
            sigVec = sigTemp%getArray()
            !
            temp_vec = self%metric%VEdge%getArray()
            !
            self%VomegaMuSig = MU_0 * inOmega * sigVec( self%EDGEi ) * temp_vec( self%EDGEi )
            !
            self%omega = inOmega
            !
        else
            if( self%omega .gt. 0 ) then
                self%VomegaMuSig = ( self%VomegaMuSig / self%omega )
            endif
            !
            self%VomegaMuSig = ( self%VomegaMuSig * inOmega )
            !
            self%omega = inOmega
            !
        endif
        !
    end subroutine updateOmegaMuSig
    !
    !> To complete setup conductivity is required
    !> DivCorInit has to be called before this routine
    !
    subroutine divCorInitModelOperatorSP( self )
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
        call boundaryIndexSP( NODE, self%metric%grid, self%NODEb, self%NODEi )
        !
        !> (1) first construct VDiv operator transpose of topology
        call RMATtrans( G, Temp )
        !
        !> pre-multiply by dual face area
        aux_vec = real( self%metric%DualFaceArea%getArray(), kind=prec )
        !
        call RMATxDIAG( Temp, aux_vec, Temp2 )
        !
        !write( *, * ) "size( aux_vec ), Temp2%nRow, Temp2%nCol: ", size( aux_vec ), Temp2%nRow, Temp2%nCol
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
        aux_vec = real( self%metric%Edgelength%getArray(), kind=prec )
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
    end subroutine divCorInitModelOperatorSP
    !
    !> To complete setup conductivity is required
    !> DivCorInit has to be called before this routine
    !
    subroutine divCorSetUpModelOperatorSP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        !
        type( spMatCSR_Real ) :: temp_matrix
        type( rVector3D_SG_t ) :: temp_vec
        real( kind=prec ), allocatable, dimension(:) :: d
        integer, allocatable, dimension(:) :: allNodes
        integer :: n, i
        !
        !> Construct VDs .. multiply VDiv by Conductivity on edges; can use VomegaMuSig
        n = self%VDiv%nCol
        !
        d = ( self%VomegaMuSig / ( mu_0 * self%omega ) )
        !
        temp_vec = self%metric%VEdge
        !
        call temp_vec%switchStoreState
        !
        d = ( d / temp_vec%sv( self%EDGEi ) )
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
    end subroutine divCorSetUpModelOperatorSP
    !
    !> Implement the sparse matrix multiply for curl-curl operator
    !> for interior elements
    !> assume output y is already allocated
    !
    subroutine amultModelOperatorSP( self, omega, inE, outE, p_adjoint )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        real( kind=prec ), intent( in ), optional :: omega
        class( Field_t ), intent( in ) :: inE
        class( Vector_t ), intent( inout ) :: outE
        logical, intent( in ), optional :: p_adjoint
        !
        logical :: adjoint
        complex( kind=prec ), allocatable, dimension(:) :: array_inE, array_outE
        complex( kind=prec ), allocatable, dimension(:) :: array_inE_int, array_result
        !
        if( .NOT. inE%is_allocated ) then
            stop "Error: amultModelOperatorSP > inE not allocated"
        endif
        !
        array_inE = inE%getArray()
        array_outE = outE%getArray()
        !
        array_inE_int = array_inE( inE%ind_interior )
        !
        array_result = array_inE_int
        array_result = C_ZERO
        !
        call RMATxCVEC( self%CCii, array_inE_int, array_result )
        !
        if( present( p_adjoint ) ) then
            adjoint = p_adjoint
        else
            adjoint = .FALSE.
        endif
        !
        if( adjoint ) then
            array_result = array_result - ONE_I * ISIGN * self%VomegaMuSig * array_inE_int
        else
            array_result = array_result + ONE_I * ISIGN * self%VomegaMuSig * array_inE_int
        endif
        !
        !array_inE = C_ZERO
        array_outE( inE%ind_interior ) = array_result
        !
        call outE%setArray( array_outE )
        !
    end subroutine amultModelOperatorSP
    !
    !> Implement the sparse matrix multiply for curl-curl operator
    !> for interior/boundary elements
    !> assume output y is already allocated
    !
    subroutine multAibModelOperatorSP( self, inE, outE )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: inE
        class( Vector_t ), intent( inout ) :: outE
        !
        complex( kind=prec ), allocatable, dimension(:) :: array_inE, array_outE
        complex( kind=prec ), allocatable, dimension(:) :: array_inE_bdry, array_outE_int
        !
        if( .NOT. inE%is_allocated ) then
            stop "Error: amultModelOperatorSP > inE not allocated"
        endif
        !
        array_inE = inE%getArray()
        array_inE_bdry = array_inE( inE%ind_boundaries )
        !
        array_outE = outE%getArray()
        array_outE_int = array_outE( inE%ind_interior )
        array_outE_int = C_ZERO
        !
        ! ON CCib DIFFERENT length ???? 
        call RMATxCVEC( self%CCib, array_inE_bdry, array_outE_int )
        !
        array_outE( inE%ind_interior ) = array_outE_int
        !
        call outE%setArray( array_outE )
        !
    end subroutine multAibModelOperatorSP
    !
    !> No subroutine briefing
    !
    subroutine multCurlTModelOperatorSP( self, inH, outE )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Field_t ), intent( inout ) :: inH
        class( Vector_t ), allocatable, intent( inout ) :: outE
        !
        integer :: ix, iy, iz
        !
        select type( inH )
            !
            class is( cVector3D_SG_t )
                !
                call inH%div( self%Metric%FaceArea )
                !
                if( .NOT. outE%is_allocated ) then
                     write( *, * ) "Error:  multCurlTModelOperatorSP > outE not allocated"
                endif
                !
                select type( outE )
                    !
                    class is( cVector3D_SG_t )
                        !
                        !> Ex
                        do iy = 2, inH%Ny
                            do iz = 2, inH%Nz
                                outE%x(:, iy, iz) =   (inH%z(:, iy, iz) - &
                                inH%z(:, iy - 1, iz)) - &
                               (inH%y(:, iy, iz) - inH%y(:, iy, iz - 1))
                            enddo
                        enddo
                        !
                        !> Ey
                        do iz = 2, inH%Nz
                            do ix = 2, inH%Nx
                                outE%y(ix, :, iz) =(inH%x(ix, :, iz) - &
                                inH%x(ix, :, iz - 1)) - &
                               (inH%z(ix, :, iz) - inH%z(ix - 1, :, iz))
                            enddo
                        enddo
                        !
                        !> Ez
                        do ix = 2, inH%Nx
                            do iy = 2, inH%Ny
                                outE%z(ix, iy,:) =(inH%y(ix, iy, :) - &
                                inH%y(ix - 1, iy, :)) - &
                               (inH%x(ix, iy, :) - inH%x(ix, iy - 1, :))
                            enddo
                        enddo
                        !
                    class default
                        stop "Error: multCurlTModelOperatorSP > Unclassified outE"
                end select
                    !> 
            class default
                stop "Error: multCurlTModelOperatorSP > Unclassified input inH"
                !
        end select
        !
        call outE%mult( self%metric%Edgelength )
        !
    end subroutine multCurlTModelOperatorSP
    !
    !> No subroutine briefing
    !
    subroutine divCgradModelOperatorSP( self, inPhi, outPhi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( in ) :: inPhi
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        complex( kind=prec ), allocatable, dimension(:) :: array_inPhi, array_outPhi
        complex( kind=prec ), allocatable, dimension(:) :: array_inPhi_int, array_result
        !
        array_inPhi = inPhi%getArray()
        array_inPhi_int = array_inPhi( inPhi%ind_interior )
        !
        array_outPhi = outPhi%getArray()
        array_result = array_inPhi_int
        array_result = C_ZERO
        !
        call RMATxCVEC( self%VDsG, array_inPhi_int, array_result )
        !
        array_outPhi( inPhi%ind_interior ) = array_result
        !
        call outPhi%setArray( array_outPhi )
        !
    end subroutine divCgradModelOperatorSP
    !
    !> No subroutine briefing
    !
    subroutine divCModelOperatorSP( self, inE, outPhi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: inE
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        complex( kind=prec ), allocatable, dimension(:) :: array_inE, array_outPhi
        complex( kind=prec ), allocatable, dimension(:) :: array_inE_int, array_outPhi_int
        !
        array_inE = inE%getArray()
        array_inE_int = array_inE( inE%ind_interior )
        !
        array_outPhi = outPhi%getArray()
        array_outPhi_int = array_outPhi( outPhi%ind_interior )
        array_outPhi_int = C_ZERO
        !
        call RMATxCVEC( self%VDs, array_inE_int, array_outPhi_int )
        !
        array_outPhi( outPhi%ind_interior ) = array_outPhi_int
        !
        call outPhi%setArray( array_outPhi )
        !
    end subroutine divCModelOperatorSP
    !
    !> No subroutine briefing
    subroutine gradModelOperatorSP( self, inPhi, outE )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( in ) :: inPhi
        class( Vector_t ), intent( inout ) :: outE
        !
        complex( kind=prec ), allocatable, dimension(:) :: array_inPhi, array_outE
        complex( kind=prec ), allocatable, dimension(:) :: array_outE_int
        !
        array_inPhi = inPhi%getArray()
        !
        array_outE = outE%getArray()
        array_outE_int = array_outE( outE%ind_interior )
        !
        call RMATxCVEC( G, array_inPhi, array_outE_int )
        !
        array_outE( outE%ind_interior ) = array_outE_int
        !
        call outE%setArray( array_outE )
        !
    end subroutine gradModelOperatorSP
    !
    !> No subroutine briefing
    subroutine divModelOperatorSP( self, inE, outPhi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: inE
        class( Field_t ), intent( inout ) :: outPhi
        !
        type( spMatCSR_Real ) :: D
        !
        complex( kind=prec ), allocatable, dimension(:) :: array_inE, array_outPhi
        !
        call RMATtrans( G, D )
        !
        array_inE = inE%getArray()
        !
        array_outPhi = array_inE
        array_outPhi = C_ZERO
        !
        call RMATxCVEC( D, array_inE, array_outPhi )
        !
        call deall_spMatCSR( D )
        !
        call outPhi%setArray( array_outPhi )
        !
    end subroutine divModelOperatorSP
    !
    !> No subroutine briefing
    subroutine printModelOperatorSP( self )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        !
        stop "Subroutine print not implemented for ModelOperator_SP"
        !
    end subroutine printModelOperatorSP
    !
end module ModelOperator_SP

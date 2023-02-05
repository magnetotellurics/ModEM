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
        class( Scalar_t ), allocatable :: sigma_C
        integer, allocatable, dimension(:) :: EDGEi, EDGEb
        integer, allocatable, dimension(:) :: NODEi, NODEb
        !
        type( spMatCSR_Real ) :: CCii, CCib
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
        call self%create( grid )
        !
    end function ModelOperator_SP_ctor
    !
    !> No subroutine briefing
    subroutine createModelOperatorSP( self, grid )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( inout ) :: self
        class( Grid_t ), target, intent( in ) :: grid
        !
        integer :: nz, nInterior
        !
        self%is_allocated = .FALSE.
        !
        self%metric%grid => grid
        !
        !> Set sparse matrices for curl (T) and grad (G)
        !> operator topologies; these sparse matrices are stored
        !> in module spOpTopology
        call setCurlTopology( grid )
        !
        call setGradTopology( grid )
        !
        call boundaryIndexSP( EDGE, grid, self%EDGEb, self%EDGEi )
        !
        nInterior = size( self%EDGEi )
        !
        !> Find indexes (in vector of all) of boundary and interior edges
        !> allocate for diagonal part of curl-curl operator
        !> (maybe this should just be for interior edges)
        !> here for all edges
        allocate( self%VomegaMuSig( nInterior ) )
        !
        allocate( self%metric, source = MetricElements_CSG_t( grid ) )
        !
        !> set a default omega
        self%omega = 0.0
        !
        call self%setEquations
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
        call self%dealloc
        !
        !> and the curl and grad topology matrices
        call deall_spMatCSR(T)
        call deall_spMatCSR(G)
        !
        call deall_spMatCSR( self%VDiv )
        call deall_spMatCSR( self%VDs )
        call deall_spMatCSR( self%VDsG )
        call deall_spMatCSR( self%VDsG_L )
        call deall_spMatCSR( self%VDsG_U )
        !
        !> interior and edge indicies
        deallocate( self%EDGEi )
        deallocate( self%EDGEb )
        deallocate( self%NODEi )
        deallocate( self%NODEb )
        !
        !> and the edge conductivities
        if( allocated( self%VomegaMuSig ) ) then
            deallocate( self%VomegaMuSig )
        endif
        !
        !> and the cell conductivities
        !> note that sigma_C is only needed to set up boundary conditions
        deallocate( self%sigma_C )
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
        type( spMatCSR_Real )       :: Temp, CC, Ttrans
        integer :: m, n, nz
        real( kind=prec ), allocatable, dimension(:) :: Dtemp
        integer :: fid
        !
        m = T%nRow
        n = T%nCol
        nz = T%row( T%nRow + 1 ) - 1
        allocate( Dtemp( m ) )
        call create_spMatCSR( m, n, nz, Temp )
        call create_spMatCSR( n, m, nz, Ttrans )
        call create_spMatCSR( m, n, nz, CC )
        call RMATxDIAG( T, real( self%metric%EdgeLength%getArray(), kind=prec ), Temp )
        Dtemp = self%metric%DualEdgeLength%getArray() / self%metric%FaceArea%getArray()
        call DIAGxRMAT( Dtemp, Temp, CC )
        call RMATtrans( T, Ttrans )
        call RMATxRMAT( Ttrans, CC, Temp )
        call DIAGxRMAT( real( self%metric%EdgeLength%getArray(), kind=prec ), Temp, CC )
        call subMatrix_Real( CC, self%EDGEi, self%EDGEi, self%CCii )
        call subMatrix_Real( CC, self%EDGEi, self%EDGEb, self%CCib )
        call deall_spMatCSR( Temp )
        call deall_spMatCSR( Ttrans )
        call deall_spMatCSR( CC )
        deallocate( Dtemp )
        !
        self%eqset = .TRUE.
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
        type( rvector3D_SG_t ) :: temp_vec, sig_temp_vec
        type( rScalar3D_SG_t ) :: cell_cond
        !
        sig_temp_vec = rvector3D_SG_t( self%metric%grid, EDGE )
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
        !> ON -> call ModelParamToCell( sigma, sigma_C )
        cell_cond = rScalar3D_SG_t( self%metric%grid, EDGE )
        caLL cell_cond%setArray( cmplx( self%VomegaMuSig, 0.0, kind=prec ) )
        call sigma%dPDEmapping( ModelParameterCell_SG_ctor( self%metric%grid, cell_cond ), self%sigma_C )
        !
    end subroutine setCondModelOperatorSP
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
        type( rvector3D_SG_t ) :: temp_vec
        real( kind=prec ), allocatable, dimension(:) :: d
        integer, allocatable, dimension(:) ::  allNodes
        integer :: n, i
        !
        !> Construct VDs .. multiply VDiv by Conductivity on edges; can use VomegaMuSig
        n = self%VDiv%nCol
        !
        allocate( d( n ) )
        !
        d = self%VomegaMuSig / ( mu_0 * self%omega )
        !
        temp_vec = self%metric%VEdge
        !
        call temp_vec%switchStoreState
        !
        d = d / temp_vec%sv( self%EDGEi )
        !
        call RMATxDIAG( self%VDiv, d, self%VDs )
        !
        !>Construct VDsG: symmetric operator for divergence correction solver
        allocate( allNodes( G%nRow ) )
        do i = 1, G%nRow
            allNodes( i ) = i
        enddo
        !
        call subMatrix_Real( G, allNodes, self%NODEi, temp_matrix )
        !
        call RMATxRMAT( self%VDs, temp_matrix, self%VDsG )
        ! Setup preconditioner
        call Dilu_Real( self%VDsG, self%VDsG_L, self%VDsG_U )
        !
        !call CholInc_real(VDsG,VDsG_L)
        !call RMATtrans(VDsG_L,VDsG_U)
        !
        call deall_spMatCSR( temp_matrix )
        !
        deallocate( d, allNodes )
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
        complex( kind=prec ), allocatable, dimension(:) :: temp_array_inE, temp_array_outE
        logical :: adjoint
        !
        if( .NOT. inE%is_allocated ) then
            stop "Error: amultModelOperatorSP > inE not allocated"
        endif
        !
        temp_array_inE = inE%getArray()
        !
        call RMATxCVEC( self%CCii, temp_array_inE, temp_array_outE )
        !
        if( present( p_adjoint ) ) then
            adjoint = p_adjoint
        else
            adjoint = .FALSE.
        endif
        !
        if( adjoint ) then
            temp_array_outE = temp_array_outE - ONE_I * ISIGN * self%VomegaMuSig * temp_array_inE
        else
            temp_array_outE = temp_array_outE + ONE_I * ISIGN * self%VomegaMuSig * temp_array_inE
        endif
        !
        deallocate( temp_array_inE )
        !
        outE = cVector3D_SG_t( self%metric%grid, EDGE )
        !
        call outE%setArray( temp_array_outE )
        !
        deallocate( temp_array_outE )
        !
    end subroutine amultModelOperatorSP
    !
    !> Implement the sparse matrix multiply for curl-curl operator
    !> for interior/boundary elements
    !> assume output y is already allocated
    !
    subroutine multAibModelOperatorSP( self, bdry, outE )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: bdry
        class( Vector_t ), intent( inout ) :: outE
        !
        stop "Error: multAibModelOperatorSP not implemented > Code ready but commented out!"
        ! !
        ! logical :: adjoint
        ! type( spMatCSR_Real ) :: CCibt
        ! complex( kind=prec ), allocatable, dimension(:) :: temp_array, temp_array_inE, temp_array_outE
        ! !
        ! if( .NOT. bdry%is_allocated ) then
            ! stop "Error: amultModelOperatorSP > bdry not allocated"
        ! endif
        ! !
        ! temp_array_inE = bdry%getArray()
        ! !
        ! call RMATxCVEC( self%CCii, temp_array_inE, temp_array_outE )
        ! !
        ! if( present( p_adjoint ) ) then
            ! adjoint = p_adjoint
        ! else
            ! adjoint = .FALSE.
        ! endif
        ! !
        ! !> ON ORIGINAL OMPLEMENTATION
        ! if( adjoint ) then
            ! !
            ! call RMATtrans( self%CCib, CCibt )
            ! !
            ! allocate( temp_array( size( self%EDGEb ) ) )
            ! !
            ! call RMATxCVEC( CCibt, temp_array_inE, temp_array )
            ! !
            ! temp_array_outE( self%EDGEb ) = temp_array;
            ! !
            ! call deall_spMATcsr( CCibt )
            ! !
            ! deallocate( temp_array )
            ! !
        ! else
            ! call RMATxCVEC( self%CCib, temp_array_inE, temp_array_outE )
        ! endif
        ! !
        ! deallocate( temp_array_inE )
        ! !
        ! outE = rVector3D_SG_t( self%metric%grid, EDGE )
        ! !
        ! call outE%setArray( temp_array_outE )
        ! !
        ! deallocate( temp_array_outE )
        ! !
    end subroutine multAibModelOperatorSP
    !
    !> No subroutine briefing
    subroutine multCurlTModelOperatorSP( self, inH, outE )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: inH
        class( Vector_t ), allocatable, intent( inout ) :: outE
        !
        integer :: ix, iy, iz
        !
        select type( inH )
            class is( cVector3D_SG_t )
                !
                call inH%div( self%Metric%FaceArea )
                !
                if(.NOT.outE%is_allocated) then
                     write( *, * ) "Error:  multCurlTModelOperatorSP > output vector not allocated"
                endif
                !
                select type( outE )
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
                        stop "Error: multCurlTModelOperatorSP > Incompatible input [outE]"
                end select
                    !> 
            class default
                stop "Error: multCurlTModelOperatorSP > Incompatible input [inH]"
                !
        end select
        !
        call outE%mult( self%metric%EdgeLength )
        !
    end subroutine multCurlTModelOperatorSP
    !
    !> No subroutine briefing
    subroutine divCgradModelOperatorSP( self, inPhi, outPhi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Scalar_t ), intent( in ) :: inPhi
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        complex( kind=prec ), allocatable, dimension(:) :: temp_array_inPhi, temp_array_outPhi
        !
        temp_array_inPhi = inPhi%getArray()
        !
        call RMATxCVEC( self%VDsG, inPhi%getArray(), temp_array_outPhi )
        !
        deallocate( temp_array_inPhi )
        !
        outPhi = rScalar3D_SG_t( self%metric%grid, EDGE )
        !
        call outPhi%setArray( temp_array_outPhi )
        !
        deallocate( temp_array_outPhi )
        !
    end subroutine divCgradModelOperatorSP
    !
    !> No subroutine briefing
    subroutine divCModelOperatorSP( self, inE, outPhi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: inE
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        complex( kind=prec ), allocatable, dimension(:) :: temp_array_inE, temp_array_outPhi
        !
        temp_array_inE = inE%getArray()
        !
        call RMATxCVEC( self%VDs, temp_array_inE, temp_array_outPhi )
        !
        deallocate( temp_array_inE )
        !
        outPhi = rVector3D_SG_t( inE%grid, inE%grid_type )
        !
        call outPhi%setArray( temp_array_outPhi )
        !
        deallocate( temp_array_outPhi )
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
        complex( kind=prec ), allocatable, dimension(:) :: temp_array_inPhi, temp_array_outE
        !
        temp_array_inPhi = inPhi%getArray()
        !
        call RMATxCVEC( G, temp_array_inPhi, temp_array_outE )
        !
        deallocate( temp_array_inPhi )
        !
        outE = rVector3D_SG_t( inPhi%grid, inPhi%grid_type )
        !
        call outE%setArray( temp_array_outE )
        !
        deallocate( temp_array_outE )
        !
    end subroutine gradModelOperatorSP
    !
    !> No subroutine briefing
    subroutine divModelOperatorSP( self, inE, outPhi )
        implicit none
        !
        class( ModelOperator_SP_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: inE
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        type( spMatCSR_Real ) :: D
        !
        complex( kind=prec ), allocatable, dimension(:) :: temp_array_inE, temp_array_outPhi
        !
        call RMATtrans( G, D )
        !
        temp_array_inE = inE%getArray()
        !
        call RMATxCVEC( D, temp_array_inE, temp_array_outPhi )
        !
        call deall_spMatCSR( D )
        !
        deallocate( temp_array_inE )
        !
        outPhi = rScalar3D_SG_t( inE%grid, inE%grid_type )
        !
        call outPhi%setArray( temp_array_outPhi )
        !
        deallocate( temp_array_outPhi )
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

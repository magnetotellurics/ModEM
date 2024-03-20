!
!> SUMMARY
!> 
!> Defines differential operators of 3D rectangular grid
!> Grad, Curl Fields are discretized on the edges.
!> This is the variant of the "classic" symmetric multi-res
!> operators, for which fine grid interface edges are redundant, and
!> coarse grid edges are active;
!> 
!> This code is based on the ModEMM version:
!>  file TOperatorTopology_cmrb.m
!
module SpOpTopology_MR
    !
    use SpOpTopology_SG
    !
    use rVector3D_MR
    use rScalar3D_MR
    !
    type, extends( SpOpTopology_t ) :: SpOpTopology_MR_t
        !
        type( Grid3D_MR_t ), pointer :: grid
        !
        contains
            !
            procedure, public :: curl => curl_SpOpTopology_MR
            procedure, public :: grad => grad_SpOpTopology_MR
            !
            procedure, private :: getT1, getG1
            !
    end type SpOpTopology_MR_t

    interface SpOpTopology_MR_t
        module procedure SpOpTopology_MR_ctor
    end interface SpOpTopology_MR_t
    !
    type :: Comp
        !
        character :: xyz
        integer :: xmin, xstep, xmax
        integer :: ymin, ystep, ymax
        integer :: zmin, zstep, zmax
        !
    end type Comp
    !
contains
    !
    function SpOpTopology_MR_ctor( grid ) result( self )
        implicit none
        !
        type( Grid3D_MR_t ), target, intent( in ) :: grid
        !
        type( SpOpTopology_MR_t ) :: self
        !
        self%grid => grid
        !
    end function SpOpTopology_MR_ctor
    !
    !> The multi-resolution version can be  viewed  as a product of
    !>   three sparse matrices:   curl = T3*T2*T1  where
    !>   --> T1 maps from the active edges, to a full set of edges
    !>          for all subVectors; all active edges are just copied,
    !>          and redundant edges in the outputs are set as averages
    !>          of the pairs of corresponding fine grid inputs
    !> 
    !>   --> T2 is the block diagonal mapping defined on each
    !>           subgrid, mapping from all edges to all faces
    !> 
    !>   --> T3 just selects out the active faces; this does not
    !>          actually need to be implemented as a matrix; just reduce
    !>          the product T2*T1 to the submatrix of rows corresponding to
    !>          active faces
    !
    subroutine curl_SpOpTopology_MR( self, curl ) 
        implicit none
        !
        class( SpOpTopology_MR_t ), intent( in ) :: self
        type( spMatCSR_Real ), intent( inout ) :: curl
        !
        type( spMatCSR_Real ), pointer, dimension(:) :: T2_array
        type( spMatCSR_Real ) :: T1, T2, Ctmp
        type( SpOpTopology_SG_t ) :: TOp1
        type( rVector3D_MR_t ) :: vecC
        integer, allocatable, dimension(:) :: col
        integer :: i
        !
        call self%getT1(T1)
        !
        allocate( T2_array( self%grid%n_grids ) )
        !
        do i = 1, self%grid%n_grids
            TOp1 = SpOpTopology_SG_t( self%grid%sub_grid(i) )
            call Top1%curl( T2_array(i) )
        enddo
        !
        call BlkDiag_Real( T2_array, T2 )
        call RMATxRMAT( T2, T1, Ctmp )
        !
        !> Pick out rows associated with active faces ...
        vecC = rVector3D_MR_t( self%grid, FACE )
        allocate( col( Ctmp%nCol ) )
        col =(/(i, i = 1, Ctmp%nCol)/)
        call subMatrix_Real( Ctmp, self%grid%FACEa, col, curl )
        !
        !> Clean up
        do i = 1, self%grid%n_grids
            call deall_spMatCSR_Real( T2_array(i) )
        enddo
        deallocate( T2_array )
        !
        call deall_spMatCSR( T1 )
        call deall_spMatCSR( T2 )
        call deall_spMatCSR( Ctmp )
        !
    end subroutine curl_SpOpTopology_MR
    !
    !> The multi-resolution version can be  viewed  as a product of
    !> three sparse matrices:   grad = G3*G2*G1  where
    !>   --> G1 maps from the active nodes, to a full set of nodes
    !>       for all subVectors; all active nodes are just copied,
    !>       and redundant nodes on interfaces are copied from the
    !>       corresponding(non-redundent) edge on the fine-grid side
    !>       of the interface.
    !> 
    !>   --> G2 is the block diagonal mapping defined on each
    !>       subgrid, mapping from all nodes to all edges
    !> 
    !>   --> G3 just selects out the active edges; this does not
    !>       actually need to be implemented as a matrix; just reduce
    !>       the product T2*T1 to the submatrix of rows corresponding to
    !>       active edges.
    !
    subroutine grad_SpOpTopology_MR( self, grad )
        implicit none
        !
        class( SpOpTopology_MR_t ), intent( in ) :: self
        type( spMatCSR_Real ), intent( inout ) :: grad
        !
        type( spMatCSR_Real ), pointer, dimension(:) :: G2_array
        type( spMatCSR_Real ) :: G1, G2, Gtmp
        type( SpOpTopology_SG_t ) :: TOp1
        integer, allocatable, dimension(:) :: col
        integer :: i
        !
        call self%getG1( G1 )
        !
        allocate( G2_array( self%grid%n_grids ) )
        !
        do i = 1, self%grid%n_grids
            TOp1 = SpOpTopology_SG_t( self%grid%sub_grid(i) )
            call Top1%grad( G2_array(i) )
        enddo
        !
        call BlkDiag_Real( G2_array, G2 )
        call RMATxRMAT( G2, G1, Gtmp )
        !
        !> Pick out rows associated with active edges ...
        allocate( col( Gtmp%nCol ) )
        col =(/(i, i = 1, Gtmp%nCol)/)
        call subMatrix_Real( Gtmp, self%grid%EDGEa, col, grad )
        !
        do i = 1, self%grid%n_grids
            call deall_spMatCSR_Real( G2_array(i) )
        enddo
        deallocate( G2_array )
        !
        call deall_spMatCSR( G1 )
        call deall_spMatCSR( G2 )
        call deall_spMatCSR( Gtmp )
        !
    end subroutine grad_SpOpTopology_MR
    !
    !> No private subroutine briefing
    !
    subroutine getT1( self, T1 )
        implicit none
        !
        class( SpOpTopology_MR_t ), intent( in ) :: self
        type( spMatCSR_Real ), intent( out ) :: T1
        !
        type( spMatIJS_Real) :: T1_ijs
        type(rVector3D_MR_t) :: vecR
        !
        integer, allocatable, dimension(:) :: R, C
        real( kind=prec ), allocatable, dimension(:) :: S
        integer, allocatable, dimension(:) :: Rtmp, Ctmp
        real( kind=prec ), allocatable, dimension(:) :: Stmp
        integer, allocatable, dimension(:) :: indSet

        type( Comp ) :: compR1(4), compR2(4), &
        compC1(2), compC2(2), compC3(2)
        type( rVector3D_MR_t ) :: vecR1, vecR2, vecC1, vecC2, vecC3

        integer, allocatable, dimension(:) :: indXcoarse, indYcoarse
        integer, allocatable, dimension(:) :: indXcoarse1, indYcoarse1
        integer, allocatable, dimension(:) :: indXcoarse2, indYcoarse2

        real( kind=prec ) :: Cc(2), cR(4)
        integer :: i, k, n, n_grids, iComp
        integer :: kVecC, kVecR
        integer :: n_rows, n_cols
        logical :: XY
        integer :: nXedge, nYedge, nZedge
        !
        !> Total number of edges
        call self%grid%numberOfEdges(nXedge, nYedge, nZedge)
        n_rows = nXedge + nYedge + nZedge
        !
        !> This extra argument set to true results in redundant x and
        !> y edges distinguished by different values.
        XY = .TRUE.;
        !
        n_cols = size( self%grid%EDGEa )
        !
        allocate( R(n_cols), C(n_cols), S(n_cols) )
        !
        R = self%grid%EDGEa
        C =(/(i, i = 1, n_cols)/)
        S = 1
        !
        !> findValuefine grid interface edges which subdivide coarse grid
        !>     edges:  these all need to be set
        !
        !> Below,
        !>     =  0 means 'end' location in the array'(Matlab notation).
        !>     = -1 means 'end -1' locaton in the array'(Matlab notation)
        !>
        !> First 4 are fine grid edges aligned with  coarse edges.
        compR1(1)%xyz = 'x';
        compR1(1)%xmin = 1; compR1(1)%xstep = 2; compR1(1)%xmax = 0 
        compR1(1)%ymin = 1; compR1(1)%ystep = 2; compR1(1)%ymax = 0 
        !
        compR1(2)%xyz = 'x';
        compR1(2)%xmin = 2; compR1(2)%xstep = 2; compR1(2)%xmax = 0 
        compR1(2)%ymin = 1; compR1(2)%ystep = 2; compR1(2)%ymax = 0 
        !
        compR1(3)%xyz = 'y';
        compR1(3)%xmin = 1; compR1(3)%xstep = 2; compR1(3)%xmax = 0 
        compR1(3)%ymin = 1; compR1(3)%ystep = 2; compR1(3)%ymax = 0 
        !
        compR1(4)%xyz = 'y';
        compR1(4)%xmin = 1; compR1(4)%xstep = 2; compR1(4)%xmax = 0 
        compR1(4)%ymin = 2; compR1(4)%ystep = 2; compR1(4)%ymax = 0 
        !
        !> Next 4 are fine grid edges that subdivide coarse grid faces.
        compR2(1)%xyz = 'x';
        compR2(1)%xmin = 1; compR2(1)%xstep = 2; compR2(1)%xmax = 0 
        compR2(1)%ymin = 2; compR2(1)%ystep = 2; compR2(1)%ymax = 0 
        !
        compR2(2)%xyz = 'x';
        compR2(2)%xmin = 2; compR2(2)%xstep = 2; compR2(2)%xmax = 0 
        compR2(2)%ymin = 2; compR2(2)%ystep = 2; compR2(2)%ymax = 0 
        !
        compR2(3)%xyz = 'y';
        compR2(3)%xmin = 2; compR2(3)%xstep = 2; compR2(3)%xmax = 0 
        compR2(3)%ymin = 1; compR2(3)%ystep = 2; compR2(3)%ymax = 0 
        !
        compR2(4)%xyz = 'y';
        compR2(4)%xmin = 2; compR2(4)%xstep = 2; compR2(4)%xmax = 0 
        compR2(4)%ymin = 2; compR2(4)%ystep = 2; compR2(4)%ymax = 0 
        !
        cR =(/1, 2, 3, 4/)
        !
        !> These are the coarse grid edges;(1) and(4) are used for
        !> "copying"; others for averaging.
        !> Below,
        !>     =  0 means 'end' location in the array'(Matlab notation).
        !>     = -1 means 'end -1' locaton in the array'(Matlab notation)
        compC1(1)%xyz = 'x';
        compC1(1)%xmin = 1; compC1(1)%xstep = 1; compC1(1)%xmax = 0 
        compC1(1)%ymin = 1; compC1(1)%ystep = 1; compC1(1)%ymax = 0 
        !
        compC2(1)%xyz = 'x';
        compC2(1)%xmin = 1; compC2(1)%xstep = 1; compC2(1)%xmax = 0 
        compC2(1)%ymin = 1; compC2(1)%ystep = 1; compC2(1)%ymax = -1 
        !
        compC3(1)%xyz = 'x';
        compC3(1)%xmin = 1; compC3(1)%xstep = 1; compC3(1)%xmax = 0 
        compC3(1)%ymin = 2; compC3(1)%ystep = 1; compC3(1)%ymax = 0 
        !
        compC1(2)%xyz = 'y';
        compC1(2)%xmin = 1; compC1(2)%xstep = 1; compC1(2)%xmax = 0 
        compC1(2)%ymin = 1; compC1(2)%ystep = 1; compC1(2)%ymax = 0 
        !
        compC2(2)%xyz = 'y';
        compC2(2)%xmin = 1; compC2(2)%xstep = 1; compC2(2)%xmax = -1 
        compC2(2)%ymin = 1; compC2(2)%ystep = 1; compC2(2)%ymax = 0  
        !
        compC3(2)%xyz = 'y';
        compC3(2)%xmin = 2; compC3(2)%xstep = 1; compC3(2)%xmax = 0 
        compC3(2)%ymin = 1; compC3(2)%ystep = 1; compC3(2)%ymax = 0 
        !
        vecR1 = rVector3D_MR_t( self%grid, EDGE )
        vecR2 = rVector3D_MR_t( self%grid, EDGE )
        vecC1 = rVector3D_MR_t( self%grid, EDGE )
        vecC2 = rVector3D_MR_t( self%grid, EDGE )
        vecC3 = rVector3D_MR_t( self%grid, EDGE )
        !
        Cc =(/1, 2/)
        !
        !> Looping over interfaces, set zlev for coarse and fine grids
        !> vecC(for columns) is for coarse grid, vecR(rows) is for fine.
        !
        SubGrids: do k = 2, self%grid%n_grids
            !
            if( self%grid%coarseness( k-1, 1 ) < self%grid%coarseness( k, 1 ) ) then
                !
                !> Fine grid is on top -- z-level to average to is at
                !> bottom of grid k - 1.
                !
                do i = 1, 4
                    !
                    compR1(i)%zmin = 0
                    compR1(i)%zstep = 1
                    compR1(i)%zmax = 0 ! zlev = 'end'
                    !
                    compR2(i)%zmin = 0
                    compR2(i)%zstep = 1
                    compR2(i)%zmax = 0 ! zlev = 'end'
                    !
                enddo
                !
                kVecR = k - 1
                !
                do i = 1, 2
                    !
                    compC1(i)%zmin = 1
                    compC1(i)%zstep = 1
                    compC1(i)%zmax = 1 ! zlev = 1
                    !
                    compC2(i)%zmin = 1
                    compC2(i)%zstep = 1
                    compC2(i)%zmax = 1 ! zlev = 1
                    !
                    compC3(i)%zmin = 1
                    compC3(i)%zstep = 1
                    compC3(i)%zmax = 1 ! zlev = 1
                    !
                enddo
                !
                kVecC = k
                !
            else
            !
            do i = 1, 4
                !
                compR1(i)%zmin = 1
                compR1(i)%zstep = 1
                compR1(i)%zmax = 1 ! zlev = 1
                !
                compR2(i)%zmin = 1
                compR2(i)%zstep = 1
                compR2(i)%zmax = 1 ! zlev = 1
                !
            enddo
            !
            kVecR = k
            !
            do i = 1, 2
                !
                compC1(i)%zmin = 0
                compC1(i)%zstep = 1
                compC1(i)%zmax = 0 ! zlev = 'end'
                !
                compC2(i)%zmin = 0
                compC2(i)%zstep = 1
                compC2(i)%zmax = 0 ! zlev = 'end'
                !
                compC3(i)%zmin = 0
                compC3(i)%zstep = 1
                compC3(i)%zmax = 0 ! zlev = 'end'
                !
            enddo
            !
            kVecC = k - 1
            endif
            !
            do i = 1, 2
                !
                call vecC1%sub_vector(kVecC)%setVecComponents( compC1(i)%xyz, &
                compC1(i)%xmin, compC1(i)%xstep, compC1(i)%xmax, &
                compC1(i)%ymin, compC1(i)%ystep, compC1(i)%ymax, &
                compC1(i)%zmin, compC1(i)%zstep, compC1(i)%zmax, cmplx( Cc(i), 0.0, kind=prec ) )
                !
                call vecC2%sub_vector(kVecC)%setVecComponents(compC2(i)%xyz, &
                compC2(i)%xmin, compC2(i)%xstep, compC2(i)%xmax, &
                compC2(i)%ymin, compC2(i)%ystep, compC2(i)%ymax, &
                compC2(i)%zmin, compC2(i)%zstep, compC2(i)%zmax, cmplx( Cc(i), 0.0, kind=prec ) )
                !
                call vecC3%sub_vector(kVecC)%setVecComponents(compC3(i)%xyz, &
                compC3(i)%xmin, compC3(i)%xstep, compC3(i)%xmax, &
                compC3(i)%ymin, compC3(i)%ystep, compC3(i)%ymax, &
                compC3(i)%zmin, compC3(i)%zstep, compC3(i)%zmax, cmplx( Cc(i), 0.0, kind=prec ) )
            enddo
            !
            do i = 1, 4
                !
                call vecR1%sub_vector(kVecR)%setVecComponents(compR1(i)%xyz, &
                compR1(i)%xmin, compR1(i)%xstep, compR1(i)%xmax, &
                compR1(i)%ymin, compR1(i)%ystep, compR1(i)%ymax, &
                compR1(i)%zmin, compR1(i)%zstep, compR1(i)%zmax, cmplx( cR(i), 0.0, kind=prec ) )
                !
                call vecR2%sub_vector(kVecR)%setVecComponents( compR2(i)%xyz, &
                compR2(i)%xmin, compR2(i)%xstep, compR2(i)%xmax, &
                compR2(i)%ymin, compR2(i)%ystep, compR2(i)%ymax, &
                compR2(i)%zmin, compR2(i)%zstep, compR2(i)%zmax, cmplx( cR(i), 0.0, kind=prec ) )
                !
            enddo
            !
        enddo SubGrids
        !
        !> First do the edges that coincide with coarse(active) edges.
        !
        indXcoarse = vecC1%findValue( 1.0_prec )
        indYcoarse = vecC1%findValue( 2.0_prec )
        !
        n = size(C) + 2 * size( indXcoarse ) + 2 * size( indYcoarse )
        allocate(Ctmp(n))
        Ctmp =(/C, indXcoarse, indXcoarse, indYcoarse, indYcoarse /)
        deallocate(C)
        call move_alloc(Ctmp, C)
        !
        do iComp = 1, 4
            !
            indSet = vecR1%findFull(real(iComp, prec))
            !
            allocate(Rtmp(size(R) + size(indSet)))
            Rtmp =(/R, indSet/)
            deallocate(R)
            call move_alloc(Rtmp, R)
            !
            allocate(Stmp(size(S) + size(indSet)))
            Stmp(1:size(S)) = S; Stmp(size(S)+1:) = 0.5
            deallocate(S)
            call move_alloc(Stmp, S)
            !
            deallocate(indSet)
            !
        enddo
        !
        !> Next fill in fine grid edges that subdivide coarse face.
        !
        indXcoarse1 = vecC2%findValue( 1.0_prec )
        indXcoarse2 = vecC3%findValue( 1.0_prec )
        indYcoarse1 = vecC2%findValue( 2.0_prec )
        indYcoarse2 = vecC3%findValue( 2.0_prec )
        !
        n = size(C) + 2 * ( size( indXcoarse1 ) + size( indXcoarse2 ) ) + &
        2 * ( size( indYcoarse1 ) + size( indYcoarse2 ) )
        !
        allocate( Ctmp(n) )
        !
        Ctmp =(/C, indXcoarse1, indXcoarse2, indXcoarse1, indXcoarse2,&
        indYcoarse1, indYcoarse2, indYcoarse1, indYcoarse2/)
        deallocate(C)
        !
        call move_alloc( Ctmp, C )
        !
        do iComp = 1, 4
            !
            indSet = vecR2%findFull( real(iComp, prec) )
            !
            allocate( Rtmp( size(R) + 2 * size(indSet) ) )
            Rtmp = (/R, indSet, indSet/)
            deallocate(R)
            !
            call move_alloc( Rtmp, R )
            !
            allocate( Stmp( size(S) + 2 * size(indSet) ) )
            Stmp( 1:size(S) ) = S; Stmp( size(S) + 1:) = 0.25
            deallocate(S)
            !
            call move_alloc(Stmp, S)
            !
            deallocate(indSet)
            !
        enddo
        !
        !> Create output curl matrix
        !
        !> First temporary IJS format
        n = size(S)
        !
        call create_spMatIJS( n_rows, n_cols, n, T1_ijs )
        !
        do i = 1, n
            !
            T1_ijs%I(i) = R(i)
            T1_ijs%J(i) = C(i)
            T1_ijs%S(i) = S(i)
            !
        enddo
        !>
        !> Finally in CSR format
        call create_spMatCSR( n_rows, n_cols, n, T1 )
        call ijs2csr( T1_ijs, T1 )
        !
        call deall_spMatIJS( T1_ijs )
        !
    end subroutine getT1
    !
    !> No private subroutine briefing
    !
    subroutine getG1( self, G1 )
        implicit none
        !
        class(SpOpTopology_MR_t), intent( in ) :: self
        type(spMatCSR_Real), intent( out ) :: G1
        ! Local variables
        type(spMatIJS_Real) :: G1_ijs
        type(rScalar3D_MR_t) :: vecR, vecC(9)
        !
        integer, allocatable, dimension(:) :: R, C
        real(kind=prec), dimension(:), allocatable :: S
        integer, allocatable, dimension(:) :: Rtmp, Ctmp
        real(kind=prec), dimension(:), allocatable :: Stmp
        integer, allocatable, dimension(:) :: indR, indC_i
        type(Comp) :: compR(4), compC(9)
        real(kind=prec) :: Cc(9), cR(4)

        integer :: kVecC, kVecR
        integer :: i, k, n, n_grids
        integer :: n_rows, n_cols
        !
        n_rows = self%grid%numberOfNodes()
        !
        n_cols = size( self%grid%NODEa )
        allocate(R(n_cols), C(n_cols), S(n_cols))

        R = self%grid%NODEa 
        C =(/(i, i = 1, n_cols)/)
        S = 1

        !**
        ! Fine grid interface nodes to interpolate to(rows of sparse
        ! matrix).
        !   = 0 means 'end' of the array'
        compR(1)%xmin = 1; compR(1)%xstep = 2; compR(1)%xmax = 0 
        compR(1)%ymin = 1; compR(1)%ystep = 2; compR(1)%ymax = 0 
        !**
        ! 2 and 3 subdivide x and y edges respectively, 4 are
        !  centers of faces.
        !   = 0 means 'end' of the array'
        compR(2)%xmin = 2; compR(2)%xstep = 2; compR(2)%xmax = 0 
        compR(2)%ymin = 1; compR(2)%ystep = 2; compR(2)%ymax = 0

        compR(3)%xmin = 1; compR(3)%xstep = 2; compR(3)%xmax = 0 
        compR(3)%ymin = 2; compR(3)%ystep = 2; compR(3)%ymax = 0 

        compR(4)%xmin = 2; compR(4)%xstep = 2; compR(4)%xmax = 0 
        compR(4)%ymin = 2; compR(4)%ystep = 2; compR(4)%ymax = 0 

        cR =(/1, 2, 3, 4/)

        !**
        ! Coarse grid interface nodes to average.
        ! Copy adjacent nodes.
        !   = 0 means 'end' of the array'
        compC(1)%xmin = 1; compC(1)%xstep = 1; compC(1)%xmax = 0 
        compC(1)%ymin = 1; compC(1)%ystep = 1; compC(1)%ymax = 0
        ! 2:5 are for subdivided edges(x then y), last 4 are
        ! for center nodes.
        !   =  0 means 'end' element of the array(MATLAB notation)
        !   = -1 means 'end - 1' element of the array(MATLAB notation)
        compC(2)%xmin = 1; compC(2)%xstep = 1; compC(2)%xmax = -1 
        compC(2)%ymin = 1; compC(2)%ystep = 1; compC(2)%ymax = 0

        compC(3)%xmin = 2; compC(3)%xstep = 1; compC(3)%xmax = 0 
        compC(3)%ymin = 1; compC(3)%ystep = 1; compC(3)%ymax = 0

        compC(4)%xmin = 1; compC(4)%xstep = 1; compC(4)%xmax = 0 
        compC(4)%ymin = 1; compC(4)%ystep = 1; compC(4)%ymax = -1

        compC(5)%xmin = 1; compC(5)%xstep = 1; compC(5)%xmax = 0
        compC(5)%ymin = 2; compC(5)%ystep = 1; compC(5)%ymax = 0

        compC(6)%xmin = 1; compC(6)%xstep = 1; compC(6)%xmax = -1 
        compC(6)%ymin = 1; compC(6)%ystep = 1; compC(6)%ymax = -1

        compC(7)%xmin = 2; compC(7)%xstep = 1; compC(7)%xmax = 0
        compC(7)%ymin = 1; compC(7)%ystep = 1; compC(7)%ymax = -1

        compC(8)%xmin = 1; compC(8)%xstep = 1; compC(8)%xmax = -1 
        compC(8)%ymin = 2; compC(8)%ystep = 1; compC(8)%ymax = 0

        compC(9)%xmin = 2; compC(9)%xstep = 1; compC(9)%xmax = 0
        compC(9)%ymin = 2; compC(9)%ystep = 1; compC(9)%ymax = 0

        cC =(/(i, i = 1, 9)/)

        do i = 1, 9
            vecC(i) = rScalar3D_MR_t( self%grid, NODE )
        enddo

        !**
        ! Looping over interfaces, set zlev for coarse and fine grids
        ! vecC(for columns) is for coarse grid, vecR(rows) is for
        ! fine.
        !*
        
        ! Create a vector consisting of all nodes
        vecR = rScalar3D_MR_t( self%grid, NODE )


        SubGrids:do k = 2, self%grid%n_grids
        if(self%grid%coarseness(k - 1, 1) < &
        self%grid%coarseness(k, 1)) then
        ! Fine grid is on top -- z-level for averaging is at
        ! bottom of grid k - 1.
        do i = 1, 4
        ! zlev = 'end' ModEMM notation
        compR(i)%zmin = 0; compR(i)%zstep = 1; compR(i)%zmax = 0 
        enddo

        kVecR = k - 1

        do i = 1, 9
        ! zlev = '1' ModEMM notation
        compC(i)%zmin = 1; compC(i)%zstep = 1; compC(i)%zmax = 1
        enddo

        kVecC = k
        else
        do i = 1, 4
        ! zlev = '1' ModEMM notation
        compR(i)%zmin = 1; compR(i)%zstep = 1; compR(i)%zmax = 1 
        enddo

        kVecR = k

        do i = 1, 9
        ! zlev = 'end' ModEMM notation
        compC(i)%zmin = 0; compC(i)%zstep = 1; compC(i)%zmax = 0 
        enddo

        kVecC = k - 1
        endif

        do i = 1, 9
        call vecC(i)%sub_scalar(kVecC)%setVecComponents(&
        compC(i)%xyz, &
        compC(i)%xmin, compC(i)%xstep, compC(i)%xmax, &
        compC(i)%ymin, compC(i)%ystep, compC(i)%ymax, &
        compC(i)%zmin, compC(i)%zstep, compC(i)%zmax, cC(i))
        enddo

        do i = 1, 4
        call vecR%sub_scalar(kVecR)%setVecComponents(&
        compR(i)%xyz, &
        compR(i)%xmin, compR(i)%xstep, compR(i)%xmax, &
        compR(i)%ymin, compR(i)%ystep, compR(i)%ymax, &
        compR(i)%zmin, compR(i)%zstep, compR(i)%zmax, cR(i))
        enddo
        enddo SubGrids

        !**
        ! First copy
        !*
        indR = vecR%findFull(1.0_prec)
        allocate(Rtmp(size(R) + size(indR)))
        Rtmp =(/R, indR/)
        deallocate(R)
        call move_alloc(Rtmp, R)

        allocate(Stmp(size(S) + size(indR)))
        Stmp = 1
        Stmp(1:size(S)) = S
        deallocate(S)
        call move_alloc(Stmp, S)

        indC_i = vecC(1)%findValue(1.0_prec)
        allocate(Ctmp(size(C) + size(indC_i)))
        Ctmp =(/C, indC_i/)
        deallocate(C)
        call move_alloc(Ctmp, C)
        deallocate(indC_i)

        !***
        ! Nodes on x edges.
        !*

        deallocate(indR)
        indR = vecR%findFull(2.0_prec)
        allocate(Rtmp(size(R) + 2*size(indR)))
        Rtmp =(/R, indR, indR/)
        deallocate(R)
        call move_alloc(Rtmp, R)

        allocate(Stmp(size(S) + 2*size(indR)))
        Stmp = 0.5
        Stmp(1:size(S)) = S    
        deallocate(S)
        call move_alloc(Stmp, S)

        do i = 2, 3
        indC_i = vecC(i)%findValue(real(i, prec))
        allocate(Ctmp(size(C) + size(indC_i)))
        Ctmp =(/C, indC_i/)
        deallocate(C)
        call move_alloc(Ctmp, C)
        deallocate(indC_i)
        enddo

        !***
        ! Nodes on y edges.
        !*
        deallocate(indR)
        indR = vecR%findFull(3.0_prec)
        allocate(Rtmp(size(R) + 2*size(indR)))
        Rtmp =(/R, indR, indR/)
        deallocate(R)
        call move_alloc(Rtmp, R)

        allocate(Stmp(size(S) + 2*size(indR)))
        Stmp = 0.5
        Stmp(1:size(S)) = S
        deallocate(S)
        call move_alloc(Stmp, S)

        do i = 4, 5
        indC_i = vecC(i)%findValue(real(i, prec))
        allocate(Ctmp(size(C) + size(indC_i)))
        Ctmp =(/C, indC_i/)
        deallocate(C)
        call move_alloc(Ctmp, C)
        deallocate(indC_i)
        enddo

        !***
        ! Nodes in coarse grid face centers.
        !*
        deallocate(indR)
        indR = vecR%findFull(4.0_prec)

        allocate(Rtmp(size(R) + 4*size(indR)))
        Rtmp =(/R, indR, indR, indR, indR/)
        deallocate(R)
        call move_alloc(Rtmp, R)

        allocate(Stmp(size(S) + 4*size(indR)))
        Stmp = 0.25
        Stmp(1:size(S)) = S
        deallocate(S)
        call move_alloc(Stmp, S)

        do i = 6, 9
        indC_i = vecC(i)%findValue(real(i, prec))
        allocate(Ctmp(size(C) + size(indC_i)))
        Ctmp =(/C, indC_i/)
        deallocate(C)
        call move_alloc(Ctmp, C)
        deallocate(indC_i)
        enddo

        !**
        ! Create output gradient matrix
        !
        ! First temporary IJS format
        n = size(S)
        call create_spMatIJS(n_rows, n_cols, n, G1_ijs)
        do i = 1, n
        G1_ijs%I(i) = R(i)
        G1_ijs%J(i) = C(i)
        G1_ijs%S(i) = S(i)
        enddo

        ! Finally in CSR format
        call create_spMatCSR(n_rows, n_cols, n, G1)
        call ijs2csr(G1_ijs, G1)

    end subroutine getG1

end module SpOpTopology_MR

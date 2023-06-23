!
!> Abstract Base class to define a Source
!
module SourceCSEM
    !
    use Source
    use dipole1d
    use cVector3D_SG
    use rVector3D_SG
    use Grid3D_SG
    use EM1D
    use ModelParameterCell_SG
    use ModelParameterCell_SG_VTI
    !
    character(:), allocatable :: source_type_csem
    character( len=15 ), parameter :: SRC_CSEM_EM1D = "SourceCSEM_EM1D"
    character( len=19 ), parameter :: SRC_CSEM_DIPOLE1D = "SourceCSEM_Dipole1D"
    !
    character(:), allocatable :: get_1d_from
    character( len=5 ), parameter :: FROM_FIXED_VALUE = "Fixed"
    character( len=14 ), parameter :: FROM_GEO_MEAN = "Geometric_mean"
    character( len=14 ), parameter :: FROM_TX_GEO_MEAN = "Mean_around_Tx"
    character( len=11 ), parameter :: FROM_TX_LOCATION = "Tx_Position"
    !
    class( Vector_t ), allocatable :: E_p
    !
    type, abstract, extends( Source_t ) :: SourceCSEM_t
        !
        ! No derived properties
        !
        contains
            !
            procedure( interface_set_1d_model_source ), deferred, public :: set1DModel
            !
            procedure, public :: setCondAnomally
            !
    end type SourceCSEM_t
    !
    abstract interface
        !
        subroutine interface_set_1d_model_source( self )
            !
            import :: SourceCSEM_t, prec
            class( SourceCSEM_t ), intent( inout ) :: self
            !
        end subroutine interface_set_1d_model_source
        !
    end interface
    !
    contains
    !
    subroutine setCondAnomally( self, sigma_cell, cond_anomaly )
        implicit none
        !
        class( SourceCSEM_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: sigma_cell
        class( Vector_t ), allocatable, intent( inout ) :: cond_anomaly
        !
        class( ModelParameter_t ), allocatable :: aModel
        class( Scalar_t ), allocatable :: cond_cell
        complex( kind=prec ), allocatable :: v(:, :, :)
        real( kind=prec ) :: wt, sigma_1d
        integer :: nzAir, i, j, k
        class( Vector_t ), allocatable :: cond_edges, cond_nomaly
        !
        !>
        nzAir = sigma_cell%grid%nzAir
        !
        !> Set 1D conductivities to an auxiliary cond_cell
        allocate( cond_cell, source = sigma_cell )
        !
        v = cond_cell%getV()
        !
        v = R_ZERO
        !
        do k = nzAir+1, nlay1D
            !
            sigma_1d = sig1D(k)
            !
            if( trim( self%sigma%param_type ) == LOGE ) sigma_1d = log( sigma_1d )
            !
            do i = 1, cond_cell%grid%Nx
                do j = 1, cond_cell%grid%Ny
                    v( i, j, k-nzAir ) = sigma_1d
                enddo
            enddo
            !
        enddo
        !
        call cond_cell%setV( v )
        !
        !> Create ModelParam - aModel with these conductivities
        allocate( aModel, source = self%sigma )
        call aModel%setCond( cond_cell )
        call aModel%setType( self%sigma%param_type )
        !
        !> Map sigma to cond_edges vector
        allocate( cond_edges, source = rVector3D_SG_t( self%sigma%metric%grid, EDGE ) )
        call self%sigma%PDEmapping( cond_edges )
        !
        !> Map aModel to cond_nomaly vector
        allocate( cond_nomaly, source = cond_edges )
        call aModel%PDEmapping( cond_nomaly )
        !
        deallocate( aModel )
        !
        !> Subtract cond_edges from cond_nomaly
        allocate( cond_anomaly, source = cond_nomaly )
        !
        deallocate( cond_nomaly )
        !
        call cond_anomaly%sub( cond_edges )
        !
        deallocate( cond_edges )
        !
    end subroutine setCondAnomally
    !
end module SourceCSEM

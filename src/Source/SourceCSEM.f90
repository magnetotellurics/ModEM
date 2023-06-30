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
        real( kind=prec ) :: location(3)
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
    subroutine setCondAnomally( self, cond_anomaly, ani_level )
        implicit none
        !
        class( SourceCSEM_t ), intent( in) :: self
        type( rVector3D_SG_t ), intent( out ) :: cond_anomaly
        integer, intent( in ) :: ani_level
        !
        class( Scalar_t ), allocatable :: sigma_cell
        complex( kind=prec ), allocatable :: v(:, :, :)
        class( ModelParameter_t ), allocatable :: aModel
        real( kind=prec ) :: wt, temp_sigma_1d
        integer :: nzAir, nzEarth, i, j, k
        type( rVector3D_SG_t ) :: cond_nomaly
        !
        call self%sigma%getCond( sigma_cell, ani_level )
        !
        nzAir = sigma_cell%grid%nzAir
        !
        nlay1D = sigma_cell%nz + nzAir
        !
        if( allocated( zlay1D ) ) deallocate( zlay1D )
        allocate( zlay1D( nlay1D ) )
        !
        do i = 1, nlay1D
            zlay1D(i) = sigma_cell%grid%zEdge(i)
        enddo
        !
        !> For create sig1D, we divide this process into two parts (1) for air layers and 
        !>    (2) for earth layers
        !> For air layer, sig1D equal to air layer conductivity
        !> For earth layer, The Geometric mean is be used to create sig1D
        !
        if( allocated( sig1D ) ) deallocate( sig1D )
        allocate( sig1D( nlay1D ) )
        !
        sig1D(1:nzAir) = SIGMA_AIR
        !
        !> Verbose
        write( *, "( a39, a14 )" ) "- Get 1D according to: ", trim( get_1d_from )
        !
        if( trim( get_1D_from ) == "Geometric_mean" ) then
            !
            v = sigma_cell%getV()
            !
            do k = nzAir+1, nlay1D
                !
                wt = R_ZERO
                temp_sigma_1d = R_ZERO
                !
                do i = 1, sigma_cell%grid%Nx
                    do j = 1, sigma_cell%grid%Ny
                        !
                        wt = wt + sigma_cell%grid%dx(i) * sigma_cell%grid%dy(j)
                        !
                        temp_sigma_1d = temp_sigma_1d + v(i,j,k-nzAir) * &
                        sigma_cell%grid%dx(i) * sigma_cell%grid%dy(j)
                        !
                    enddo
                enddo
                !
                sig1D(k) = exp( temp_sigma_1d / wt )
                !
            enddo
            !
        elseif( trim( get_1D_from ) == "At_Tx_Position" ) then
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > At_Tx_Position not implemented yet"
            !
        elseif( trim(get_1d_from) == "Geometric_mean_around_Tx" ) then
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > Geometric_mean_around_Tx not implemented yet"
            !
        elseif( trim(get_1d_from) == "Full_Geometric_mean" ) then
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > Full_Geometric_mean not implemented yet"
            !
        elseif( trim( get_1d_from ) == "Fixed_Value" ) then
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > Fixed_Value not implemented yet"
            !
        else
            !
            stop "Error: setTemp_SourceCSEM_Dipole1D > Unknown get_1d_from"
            !
        endif
        !
        nzEarth = sigma_cell%grid%nzEarth
        !
        v = SIGMA_AIR
        !
        do k = 1, nzEarth
            !
            temp_sigma_1d = sig1D( k + nzAir )
            !
            if( trim( self%sigma%param_type ) == LOGE ) temp_sigma_1d = log( temp_sigma_1d )
            !
            do i = 1, sigma_cell%grid%Nx
                do j = 1, sigma_cell%grid%Ny
                    v( i, j, k ) = temp_sigma_1d
                enddo
            enddo
            !
        enddo
        !
        call sigma_cell%setV( v )
        !
        !> Create ModelParam from 1D: aModel
        !> with sigma_cell conductivity in the proper anisotropic direction
        allocate( aModel, source = self%sigma )
        call aModel%setCond( sigma_cell, ani_level )
        call aModel%setType( self%sigma%param_type )
        !
        deallocate( sigma_cell )
        !
        cond_nomaly = rVector3D_SG_t( self%sigma%metric%grid, EDGE )
        call aModel%PDEmapping( cond_nomaly )
        !
        !> Map sigma to cond_edges vector
        cond_anomaly = rVector3D_SG_t( self%sigma%metric%grid, EDGE )
        call self%sigma%PDEmapping( cond_anomaly )
        !
        !cond_anomaly = cond_anomaly - cond_nomaly
        call cond_anomaly%sub( cond_nomaly )
        !
    end subroutine setCondAnomally
    !
end module SourceCSEM

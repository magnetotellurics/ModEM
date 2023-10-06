!
!> Abstract Base class to define a Source
!
module SourceCSEM
    !
    use Source
    use dipole1d
    use ModelOperator_SP_V2
    !
    character(:), allocatable :: source_type_csem
    character( len=15 ), parameter :: SRC_CSEM_EM1D = "SourceCSEM_EM1D"
    character( len=19 ), parameter :: SRC_CSEM_DIPOLE1D = "SourceCSEM_Dipole1D"
    !
    character(:), allocatable :: get_1d_from
    character( len=5 ), parameter :: FROM_FIXED_VALUE = "Fixed"
    character( len=14 ), parameter :: FROM_GEOM_MEAN = "Geometric_mean"
    character( len=14 ), parameter :: FROM_TX_GEOM_MEAN = "Mean_around_Tx"
    character( len=11 ), parameter :: FROM_TX_LOCATION = "Tx_Position"
    !
    type, abstract, extends( Source_t ) :: SourceCSEM_t
        !
        real( kind=prec ) :: location(3)
        !
        contains
            !
            procedure( interface_set_1d_model_source ), deferred, public :: set1DModel
            !
            procedure, public :: setCondAnomally => setCondAnomally_SourceCSEM
            !
            procedure, public :: createRHS => createRHS_SourceCSEM
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
    subroutine setCondAnomally_SourceCSEM( self, cond_anomaly, ani_level )
        implicit none
        !
        class( SourceCSEM_t ), intent( in) :: self
        class( Vector_t ), allocatable, intent( out ) :: cond_anomaly
        integer, intent( in ) :: ani_level
        !
        type( rScalar3D_SG_t ) :: sigma_cell
        class( ModelParameter_t ), allocatable :: aModel
        real( kind=prec ) :: wt, temp_sigma_1d
        integer :: nzAir, nzEarth, i, j, k
        class( Vector_t ), allocatable :: cond_nomaly
        !
        !>
        sigma_cell = self%sigma%getCond( ani_level )
        !
        select type( grid => sigma_cell%grid )
            !
            class is( Grid3D_SG_t )
                !
                nzAir = grid%nzAir
                !
                nlay1D = sigma_cell%nz + nzAir
                !
                if( allocated( zlay1D ) ) deallocate( zlay1D )
                allocate( zlay1D( nlay1D ) )
                !
                do i = 1, nlay1D
                    zlay1D(i) = grid%z_edge(i)
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
                    do k = nzAir+1, nlay1D
                        !
                        wt = R_ZERO
                        temp_sigma_1d = R_ZERO
                        !
                        do i = 1, grid%Nx
                            do j = 1, grid%Ny
                                !
                                wt = wt + grid%dx(i) * grid%dy(j)
                                !
                                temp_sigma_1d = temp_sigma_1d + sigma_cell%v(i,j,k-nzAir) * &
                                grid%dx(i) * grid%dy(j)
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
                    call errStop( "setCondAnomally_SourceCSEM > At_Tx_Position not implemented yet" )
                    !
                elseif( trim( get_1d_from ) == "Geometric_mean_around_Tx" ) then
                    !
                    call errStop( "setCondAnomally_SourceCSEM > Geometric_mean_around_Tx not implemented yet" )
                    !
                elseif( trim( get_1d_from ) == "Full_Geometric_mean" ) then
                    !
                    call errStop( "setCondAnomally_SourceCSEM > Full_Geometric_mean not implemented yet" )
                    !
                elseif( trim( get_1d_from ) == "Fixed_Value" ) then
                    !
                    call errStop( "setCondAnomally_SourceCSEM > Fixed_Value not implemented yet" )
                    !
                else
                    !
                    call errStop( "setCondAnomally_SourceCSEM > Unknown get_1d_from" )
                    !
                endif
                !
                nzEarth = grid%nzEarth
                !
                sigma_cell%v = SIGMA_AIR
                !
                do k = 1, nzEarth
                    !
                    temp_sigma_1d = sig1D( k + nzAir )
                    !
                    if( trim( self%sigma%param_type ) == LOGE ) temp_sigma_1d = log( temp_sigma_1d )
                    !
                    do i = 1, grid%Nx
                        do j = 1, grid%Ny
                            sigma_cell%v( i, j, k ) = temp_sigma_1d
                        enddo
                    enddo
                    !
                enddo
                !
            class default
                call errStop( "setCondAnomally_SourceCSEM: undefined grid" )
            !
        end select
        !
        !> Create ModelParam from 1D: aModel
        !> with sigma_cell conductivity in the proper anisotropic direction
        allocate( aModel, source = self%sigma )
        aModel = self%sigma
        call aModel%setCond( sigma_cell, ani_level )
        call aModel%setType( self%sigma%param_type )
        !
        call self%sigma%metric%createVector( real_t, EDGE, cond_nomaly )
        call aModel%PDEmapping( cond_nomaly )
        !
        !> Map sigma to cond_edges vector
        call self%sigma%metric%createVector( real_t, EDGE, cond_anomaly )
        call self%sigma%PDEmapping( cond_anomaly )
        !
        !cond_anomaly = cond_anomaly - cond_nomaly
        call cond_anomaly%sub( cond_nomaly )
        !
    end subroutine setCondAnomally_SourceCSEM
    !
    !> Set RHS from self%E
    !
    subroutine createRHS_SourceCSEM( self )
        implicit none
        !
        class( SourceCSEM_t ), intent( inout ) :: self
        !
        type( cVector3D_MR_t ) :: temp_vec_mr
        complex( kind=prec ) :: iOmegaMuInv
        complex( kind=prec ), allocatable, dimension(:) :: v_edge_v, rhs_v, rhs_v_int, in_e, in_e_int, out_e
        !
        !if( allocated( self%rhs ) ) deallocate( self%rhs )
        allocate( self%rhs(1) )
        !
        !> Check if grid is MR 
        !> RHS calculated as MR vector
        select type( grid => self%model_operator%metric%grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( self%rhs(1)%v, source = self%E(1) )
                !
            class is( Grid3D_MR_t )
                !
                temp_vec_mr = cVector3D_MR_t( grid, self%E(1)%grid_type )
                !
                call temp_vec_mr%fromSG( self%E(1) )
                !
                allocate( self%rhs(1)%v, source = temp_vec_mr )
                !
            class default
                call errStop( "createRHS_SourceCSEM > model_operator must be SP V1 or V2" )
            !
        end select
        !
        !> Check if model_operator is SP2
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_SP_V2_t )
                !
                !RHS = RHS + (i omega mu)^-1 GDii* E(EDGEi)
                iOmegaMuInv = isign / ( cmplx( 0.0, real( 2.0 * ( PI / self%period ) * MU_0, kind=prec ), kind=prec ) )
                !
                in_e = self%E(1)%getArray()
                !
                in_e_int = in_e( self%E(1)%indInterior() )
                !
                out_e = in_e_int
                out_e = C_ZERO
                !
                call RMATxCVEC( model_operator%GDii, in_e_int, out_e )
                !
                v_edge_v = self%model_operator%metric%v_edge%getArray()
                !
                out_e = out_e * iOmegaMuInv
                !
                rhs_v = self%rhs(1)%v%getArray()
                !
                rhs_v_int = ( rhs_v( self%rhs(1)%v%indInterior() ) + out_e )
                !
                rhs_v = C_ZERO
                rhs_v( self%rhs(1)%v%indInterior() ) = rhs_v_int
                !
                call self%rhs(1)%v%setArray( rhs_v )
                !
        end select
        !
        call self%rhs(1)%v%mult( self%model_operator%metric%v_edge )
        !
    end subroutine createRHS_SourceCSEM
    !
end module SourceCSEM

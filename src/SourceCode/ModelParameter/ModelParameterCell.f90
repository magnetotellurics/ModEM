!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only.
!
module ModelParameterCell
    !
    use ModelParameter
    use rScalar3D_SG
    use FileUnits
    !
    type, abstract, extends( ModelParameter_t ) :: ModelParameterCell_t
        !
        class( Grid_t ), allocatable :: param_grid
        !
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: cell_cond
        !
        contains
            !
            !> Procedures
            procedure, public :: baseDealloc => baseDealloc_ModelParameter
            !
            procedure, public :: getOneCond => getOneCond_ModelParameterCell
            procedure, public :: getAllCond => getAllCond_ModelParameterCell
            !
            procedure, public :: setOneCond => setOneCond_ModelParameterCell
            procedure, public :: setAllCond => setAllCond_ModelParameterCell
            !
            procedure, public :: zeros => zeros_ModelParameterCell
            !
            procedure, public :: copyFrom => copyFrom_ModelParameterCell
            !
            procedure, public :: countModel => countModel_ModelParameterCell
            !
            procedure, public :: dotProd => dotProd_ModelParameterCell
            !
            procedure, public :: linComb => linComb_ModelParameterCell
            !
            procedure, public :: setType => setType_ModelParameterCell
            !
            procedure, public :: print => print_ModelParameterCell
            !
    end type ModelParameterCell_t
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine baseDealloc_ModelParameter( self )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        !
        if( allocated( self%param_grid ) ) deallocate( self%param_grid )
        !
    end subroutine baseDealloc_ModelParameter
    !
    !> No function briefing
    !
    function getOneCond_ModelParameterCell( self, i_cond ) result( cond )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        integer, intent( in ) :: i_cond
        !
        type( rScalar3D_SG_t ) :: cond
        !
        if( i_cond .GT. self%anisotropic_level ) then
            !
            call errStop( "getOneCond_ModelParameterCell > conductivity level too high" )
            !
        endif
        !
        cond = self%cell_cond( i_cond )
        !
    end function getOneCond_ModelParameterCell
    !
    !> No function briefing
    !
    function getAllCond_ModelParameterCell( self ) result( cond )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        !
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: cond
        !
        cond = self%cell_cond
        !
    end function getAllCond_ModelParameterCell
    !
    !> No interface subroutine briefing
    !
    subroutine setOneCond_ModelParameterCell( self, cond, i_cond )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: cond
        integer, intent( in ) :: i_cond
        !
        if( .NOT. cond%is_allocated ) then
            call errStop( "setOneCond_ModelParameterCell > cond not allocated" )
        endif
        !
        if( i_cond .LE. self%anisotropic_level ) then
            !
            self%cell_cond( i_cond ) = cond
            !
        else
            !
            call errStop( "setOneCond_ModelParameterCell > Unsupport general anisotropy yet" )
            !
        endif
        !
    end subroutine setOneCond_ModelParameterCell
    !
    !> No interface subroutine briefing
    !
    subroutine setAllCond_ModelParameterCell( self, cond )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        type( rScalar3D_SG_t ), dimension(:), intent( in ) :: cond
        !
        integer :: i
        !
        do i = 1, self%anisotropic_level
            !
            if( .NOT. cond(i)%is_allocated ) then
                !
                call errStop( "setAllCond_ModelParameterCell > cond not allocated" )
                !
            elseif( .NOT. self%cell_cond(i)%is_allocated )then
                !
                call errStop( "setAllCond_ModelParameterCell > self%cell_cond not allocated" )
                !
            endif
            !
            call self%setCond( cond(i), i )
            !
        enddo
        !
    end subroutine setAllCond_ModelParameterCell
    !
    !> No subroutine briefing
    !
    subroutine zeros_ModelParameterCell( self )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        !
        integer :: i
        !
        do i = 1, self%anisotropic_level
            !
            call self%cell_cond(i)%zeros
            !
        enddo
        !
    end subroutine zeros_ModelParameterCell
    !
    !> No subroutine briefing
    !
    subroutine copyFrom_ModelParameterCell( self, rhs )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "copyFrom_ModelParameterCell > rhs not allocated." )
        endif
        !
        select type( rhs )
            !
            class is( ModelParameterCell_t )
                !
                self%metric => rhs%metric
                !
                self%mKey = rhs%mKey
                !
                self%anisotropic_level = rhs%anisotropic_level
                !
                self%air_cond = rhs%air_cond
                !
                self%param_type = rhs%param_type
                !
                self%is_allocated = rhs%is_allocated
                !
                self%param_grid = rhs%param_grid
                !
                self%cell_cond = rhs%cell_cond
                !
                self%sigmap_ptr => rhs%sigmap_ptr
                !
            class default
                call errStop( "copyFrom_ModelParameterCell > Unclassified rhs." )
            !
        end select
        !
    end subroutine copyFrom_ModelParameterCell
    !
    !> No function briefing
    !
    function countModel_ModelParameterCell( self ) result( counter )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        !
        integer :: i, counter, nx, ny, nz, nzAir, nz_earth
        !
        counter = 0
        !
        do i = 1, self%anisotropic_level
            !
            if( .NOT. self%cell_cond(i)%is_allocated ) then
                write( *, * ) "Error: countModel_ModelParameterCell > cell_cond (", i, ") not allocated!"
                stop
            endif
            !
            call self%cell_cond(i)%grid%getDimensions( nx, ny, nz, nzAir )
            nz_earth = nz - nzAir
            !
            counter = counter + self%cell_cond(i)%Nx * self%cell_cond(i)%Ny * nz_earth
            !
        enddo
        !
    end function countModel_ModelParameterCell
    !
    !> No subroutine briefing
    !
    subroutine linComb_ModelParameterCell( self, a1, a2, rhs )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: a1, a2
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        select type( rhs )
            !
            class is( ModelParameterCell_t )
                !
                do i = 1, self%anisotropic_level
                    !
                    if( self%cell_cond(i)%isCompatible( rhs%cell_cond(i) ) ) then
                        !
                        self%cell_cond(i)%v = a1 * self%cell_cond(i)%v + a2 * rhs%cell_cond(i)%v
                        !
                    else
                        write( *, * ) "Error: linComb_ModelParameterCell > Incompatible rhs cell_cond (", i, ")!"
                        stop
                    endif
                    !
                enddo
                !
            class default
                call errStop( "linComb_ModelParameterCell > undefined rhs" )
            !
        end select
        !
        !> ????
        !self%air_cond = rhs%air_cond
        !
    end subroutine linComb_ModelParameterCell
    !
    !> No subroutine briefing
    !
    function dotProd_ModelParameterCell( self, rhs ) result( rvalue )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        real( kind=prec ) :: rvalue
        !
        integer :: i
        !
        rvalue = R_ZERO
        !
        select type( rhs )
            !
            class is( ModelParameterCell_t )
                !
                do i = 1, self%anisotropic_level
                    !
                    if( self%cell_cond(i)%isCompatible( rhs%cell_cond(i) ) ) then
                        !
                        rvalue = rvalue + sum( self%cell_cond(i)%v * rhs%cell_cond(i)%v )
                        !
                    else
                        write( *, * ) "Error: dotProd_ModelParameterCell > Incompatible rhs cell_cond (", i, ")!"
                        stop
                    endif
                    !
                enddo
                !
            class default
                call errStop( "dotProd_ModelParameterCell > Unclassified rhs" )
            !
        end select
        !
    end function dotProd_ModelParameterCell
    !
    !> No subroutine briefing
    !
    subroutine setType_ModelParameterCell( self, param_type )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        character(:), allocatable, intent( in ) :: param_type
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
                call errStop( "setType_ModelParameterCell > Self not allocated." )
        endif
        !
        do i = 1, self%anisotropic_level
            !
            if( trim( param_type ) .EQ. trim( self%param_type ) ) then
                ! Nothing to be done
            elseif( self%param_type == "" ) then
                self%param_type = trim( param_type )
            elseif( self%param_type == LINEAR ) then
                !
                if( param_type == LOGE ) then
                    !
                    self%cell_cond(i)%v = log( self%cell_cond(i)%v )
                    !
                elseif( param_type == LOG_10) then
                    !
                    self%cell_cond(i)%v = log10( self%cell_cond(i)%v )
                    !
                endif
                !
            elseif( param_type == LINEAR ) then
                !
                if( self%param_type == LOGE ) then
                    !
                    self%cell_cond(i)%v = exp( self%cell_cond(i)%v )
                    !
                elseif( self%param_type == LOG_10 ) then
                    !
                    self%cell_cond(i)%v = exp( self%cell_cond(i)%v * log(10.) )
                    !
                endif
                !
            elseif( ( self%param_type == LOGE ) .AND. ( param_type == LOG_10 ) ) then
                !
                self%cell_cond(i)%v = self%cell_cond(i)%v / log(10.)
                !
            elseif( ( self%param_type == LOG_10 ) .AND. ( param_type == LOGE ) ) then
                !
                self%cell_cond(i)%v = self%cell_cond(i)%v * log(10.)
                !
            else
                call errStop( "setType_ModelParameterCell > Unknown param_type." )
            endif
            !
            !cself%cell_cond(i)%v = cmplx( self%cell_cond(i)%v, 0.0, kind=prec )
            !
            !call self%cell_cond(i)%setArray( cself%cell_cond(i)%v )
            !
        enddo
        !
        self%param_type = param_type 
        !
    end subroutine setType_ModelParameterCell
    !
    !> No subroutine briefing
    !
    subroutine print_ModelParameterCell( self )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        !
        integer :: i
        !
        write( *, * ) "ModelParameterCell_t:", self%mKey, self%air_cond, self%param_type, &
        self%is_allocated, self%param_grid%nx, self%param_grid%ny, self%param_grid%nz, self%param_grid%nzAir
        !
        do i = 1, self%anisotropic_level
            !
            call self%cell_cond(i)%print
            !
        enddo
        !
    end subroutine print_ModelParameterCell
    !
end Module ModelParameterCell

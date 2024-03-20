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
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: cell_cond
        !
        class( Grid_t ), pointer :: param_grid
        !
        contains
            !
            !> Procedures
            procedure, public :: deallocCell => deallocCell_ModelParameterCell
            !
            procedure, public :: setMetric => setMetric_ModelParameterCell
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
            procedure, public :: write => write_ModelParameterCell
            !
            procedure, public :: print => print_ModelParameterCell
            !
    end type ModelParameterCell_t
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine deallocCell_ModelParameterCell( self )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        !
        if( allocated( self%cell_cond ) ) deallocate( self%cell_cond )
        !
    end subroutine deallocCell_ModelParameterCell
    !
    !> No subroutine briefing
    !
    subroutine setMetric_ModelParameterCell( self, metric )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        class( MetricElements_t ), target, intent( in ) :: metric
        !
        self%metric => metric
        !
    end subroutine setMetric_ModelParameterCell
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
        type( rScalar3D_SG_t ), intent( in ) :: cond
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
                self%param_grid => rhs%param_grid
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
        !> NEED THIS LINE????
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
                call errStop( "setType_ModelParameterCell > self not allocated." )
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
    !> opens cfile on unit ioModelParam, writes out object of
    !> type modelParam in Weerachai Siripunvaraporn"s format,
    !> closes file.
    !
    subroutine write_ModelParameterCell( self, file_name, comment )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        character(*), intent( in ) :: file_name
        character(*), intent( in ), optional :: comment
        !
        real( kind=prec ), allocatable, dimension(:,:,:) :: cond_v
        integer :: Nx, Ny, NzEarth, ii, i, j, k, ios
        !
        ! Verbose
        !write( *, * ) "     > Write Model to file: [", file_name, "]"
        !
        open( ioModelParam, file = file_name, action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            if( present( comment ) ) then
                write( ioModelParam, * ) "# ", trim( comment )
            else
                write( ioModelParam, * ) "# 3D MT model written by ModEM in WS format"
            endif
            !
            !> Write grid geometry definitions
            Nx = self%metric%grid%nx
            Ny = self%metric%grid%ny
            NzEarth = self%metric%grid%nz - self%metric%grid%nzAir
            !
            write( ioModelParam, "(4i5)", advance = "no" ) Nx, Ny, NzEarth, 0
            !
            write( ioModelParam, "(a10)", advance = "no" ) trim( self%param_type )
            !
            if( self%anisotropic_level == 2 ) then
                !
                write( ioModelParam, * ) " VTI"
                !
            else
                !
                write( ioModelParam, * )
                !
            endif
            !
            !> Write self%metric%grid spacings
            do j = 1, self%metric%grid%nx
                write( ioModelParam, "(f12.3)", advance = "no" ) self%metric%grid%dx(j)
            enddo
            !
            write( ioModelParam, * )
            !
            do j = 1, self%metric%grid%ny
                write( ioModelParam, "(f12.3)", advance = "no" ) self%metric%grid%dy(j)
            enddo
            !
            write( ioModelParam, * )
            !
            do j = self%metric%grid%nzAir + 1, self%metric%grid%nz
                write( ioModelParam, "(f12.3)", advance = "no" ) self%metric%grid%dz(j)
            enddo
            !
            write( ioModelParam, * )
            !
            do ii = 1, self%anisotropic_level
                !
                !> Convert (horizontal) conductivity to resistivity
                !
                cond_v = self%cell_cond(ii)%v
                !
                if( index( self%param_type, "LOGE" ) > 0 .OR. index( self%param_type, "LOG10" ) > 0 ) then
                    cond_v = -cond_v
                elseif( index(self%param_type, "LINEAR" ) > 0 ) then
                    cond_v = ONE / cond_v
                endif
                !
                !> Write the (horizontal) resistivity
                !
                write( ioModelParam, * )
                !
                do k = 1, nzEarth
                    do j = 1, Ny
                        do i = Nx, 1, -1
                            write( ioModelParam, "(es13.5)", iostat = ios, advance = "no" ) cond_v(i,j,k)
                        enddo
                        !
                        write( ioModelParam, * )
                        !
                    enddo
                    !
                    write( ioModelParam, * )
                    !
                enddo
                !
            enddo
            !
            !> Note that our standard subroutine doesn"t work with Weerachai"s
            !> real value format. It is still better than either Mackie"s or WS"s...
            !> call write_rscalar(ioModelParam,rho)
            !> Also write the self%metric%grid origin (in metres!) and rotation (in degrees)...
            !
            write( ioModelParam, "(3f16.3)", iostat = ios) self%metric%grid%ox, self%metric%grid%oy, self%metric%grid%oz
            write( ioModelParam, "(f9.3)", iostat = ios)  self%metric%grid%rotdeg
            !
            close( ioModelParam )
            !
        else
            call errStop( "write_ModelParameterCell > Error opening file ["//file_name//"]!" )
        endif
        !
    end subroutine write_ModelParameterCell
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
        !write( *, * ) "ModelParameterCell_t:", self%mKey, self%air_cond, self%param_type, &
        !self%is_allocated, self%param_grid%nx, self%param_grid%ny, self%param_grid%nz, self%param_grid%nzAir
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

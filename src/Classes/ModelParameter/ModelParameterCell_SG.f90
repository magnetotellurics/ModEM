!**
! Basic model parameter class --
! cell conductivities defined on numerical grid
! in Earth cells only -- as coded this is specific
! to usual CSG.
!*
module ModelParameterCell_SG
    !
    use Constants
    use Grid
    use Grid3D_SG
    use rScalar
    use rVector
    use rScalar3D_SG
    use rVector3D_SG 
    use ModelParameter
    use Grid1D
    use Grid2D
    use ModelParameter1D
    use ModelParameter2D
    !
    type, extends( ModelParameter_t ) :: ModelParameterCell_SG_t
         !
         class( Grid_t ), allocatable :: paramGrid
         !
         type( rScalar3D_SG_t ) :: cellCond
         !
        contains
              !
              final :: ModelParameterCell_SG_dtor
              !
              procedure, public :: zeros    => zerosModelParameterCell
              procedure, public :: copyFrom => copyFromModelParameterCell
              !
              ! Model mapping methods
              procedure, public :: PDEmapping   => PDEmappingModelParameterCell
              procedure, public :: dPDEmapping  => dPDEmappingModelParameterCell
              procedure, public :: dPDEmappingT => dPDEmappingTModelParameterCell
              !
              procedure, public :: slice1D => slice1DModelParameterCell
              procedure, public :: slice2D => slice2DModelParameterCell
              !
              procedure, public :: avgModel1D => avgModel1DModelParameterCell
              !
              procedure, public :: setType => setTypeModelParameterCell
              !
    end type ModelParameterCell_SG_t
    !
    interface ModelParameterCell_SG_t
         module procedure ModelParameterCell_SG_ctor
    end interface ModelParameterCell_SG_t
    !
contains
    
    !**
    !
    !*
    function ModelParameterCell_SG_ctor( grid, ccond, paramType ) result( self )
        implicit none
        !
        class( Grid3D_SG_t ), target, intent( in )        :: grid
        type( rScalar3D_SG_t ), intent( in )              :: ccond
        character(:), allocatable, optional, intent( in ) :: paramType
        !
        type( ModelParameterCell_SG_t ) :: self
        !
        integer :: nx, ny, nz, nzAir
        !
        write(*,*) "Constructor ModelParameterCell_SG_t"
        !
        call self%init()
        !
        if ( .NOT. present( paramType ) ) then
              self%paramType = LOGE
        else
              self%paramType = trim( paramType )
        end if
        !
        ! Point to the original grid
        self%grid => grid
        !
        nx = grid%nx
        ny = grid%ny
        nz = grid%nz - grid%nzAir
        !
        nzAir = 0
        !
        self%ParamGrid = Grid3D_SG_t( nx, ny, nzAir, nz, &
          grid%dx, grid%dy, &
          grid%dz( grid%nzAir+1:grid%nz ) )
        !
        self%cellCond = ccond
        !
        if ( present( paramType ) ) then
              call self%SetSigMap( paramType )
              ! We initially specify airCond as linear conductivity!
        !    self%AirCond = self%SigMap( self%airCond, "inverse" )
        !    GDE:  we should always keep AirCond as actual linear conductvity!
        !        Never apply SigMap to air layers!!!
        end if
        !
        self%is_allocated = .true.
        !
    end function ModelParameterCell_SG_ctor
    !
    ! ModelOperator_MF destructor
    subroutine ModelParameterCell_SG_dtor( self )
        implicit none
        !
        type( ModelParameterCell_SG_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor ModelParameterCell_SG"
        !
        deallocate( self%ParamGrid )
        !
    end subroutine ModelParameterCell_SG_dtor
    !
    function slice1DModelParameterCell( self, ix, iy ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        integer, intent( in ) :: ix, iy
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        real( kind=prec ), allocatable, dimension(:) :: CondSlice
        !
        !    create 1D model parameter
        model_param_1D = ModelParameter1D_t( self%grid%Slice1D() )
        !    comnductivity slice
        allocate( CondSlice( model_param_1D%grid%nz ) )
        !
        !  extract slice; convert to linear conductivity:
        CondSlice = self%SigMap( self%cellCond%v( ix, iy, : ) )
        !
        call model_param_1D%SetConductivity( CondSlice, self%AirCond, self%paramType, self%mKey )
        !
        deallocate( CondSlice )
        !
    end function slice1DModelParameterCell
    !
    function avgModel1DModelParameterCell( self ) result( model_param_1D )
        !    extracts slice corresponding to column j of model parameter
        implicit none
        ! Arguments
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        real( kind=prec ), allocatable, dimension(:) :: CondSlice
        real( kind=prec ) :: wt, temp_sigma_value
        integer :: i, j, k
        !
        model_param_1D = ModelParameter1D_t( self%grid%Slice1D() )
        !
        allocate( CondSlice( self%grid%nzEarth ) )
        !
        do k = 1, self%grid%nzEarth
            !
            wt = R_ZERO
            temp_sigma_value = R_ZERO
            do i = 1, self%grid%Nx
                do j = 1, self%grid%Ny
                    wt = wt + self%grid%dx(i) * self%grid%dy(j)
                    temp_sigma_value = temp_sigma_value + self%CellCond%v( i, j, k ) * &
                    self%grid%dx(i) * self%grid%dy(j)
                end do
            end do
            !
            CondSlice( k ) = self%SigMap( temp_sigma_value / wt )
            !
        end do
        !
        call model_param_1D%SetConductivity( CondSlice, self%AirCond, self%paramType, self%mKey )
        !
        deallocate( CondSlice )
        !
    end function avgModel1DModelParameterCell
    !
    function slice2DModelParameterCell( self, axis, j ) result( m2D )
        !    extracts slice corresponding to column j of model parameter
        implicit none
        ! Arguments
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        integer, intent( in )                          :: axis, j
        !
        type( ModelParameter2D_t ) :: m2D 
        !
        character(:), allocatable :: paramType
        real( kind=prec ), allocatable, dimension(:,:) :: CondSlice
        !
        paramType = LINEAR
        !
        !    create 2D model parameter
        m2D = ModelParameter2D_t( self%grid%Slice2D() )
        !    comnductivity slice
        allocate( CondSlice( self%grid%ny, self%grid%nzEarth ) )
        !
        if( axis == 1 ) then
            CondSlice = self%SigMap(Self%cellCond%v(j,:,:))
        else if( axis == 2 ) then
            CondSlice = self%SigMap(Self%cellCond%v(:,j,:))
        else if( axis == 3 ) then
            CondSlice = self%SigMap(Self%cellCond%v(:,:,j))
        else
            stop "ModelParameter:Slice2D: wrong axis"
        endif
        !
        call m2D%SetConductivity( CondSlice, self%AirCond, paramType, self%mKey )
        !
        deallocate( CondSlice )
        !
    end function slice2DModelParameterCell
    !**
    ! Zeros
    ! Zero model parameter
    !*
    subroutine zerosModelParameterCell( self )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent(inout) :: self
        !
        call self%cellCond%zeros()
        !
    end subroutine zerosModelParameterCell

    !**
    ! Copy rhs to self.
    !*
    subroutine copyFromModelParameterCell( self, rhs )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in )           :: rhs
        !
        select type( rhs )
            class is( ModelParameterCell_SG_t )
                self%ParamGrid = rhs%ParamGrid
                self%cellCond = rhs%cellCond
                self%airCond = rhs%airCond
                self%paramType = rhs%paramType
                self%mKey = rhs%mKey
                self%metric => rhs%metric
                self%grid => rhs%grid
            class default
                write(*, *) "ERROR:ModelParameterCell:CopyFrom"
                stop "              Incompatible input. Exiting."
        end select
        !
    end subroutine copyFromModelParameterCell
    !
    !    NOT SURE WE WANT THESE MAPPINGS TO BE FUNCTIONS ...
    function PDEmappingModelParameterCell( self ) result( eVec )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( rVector_t ), allocatable                :: eVec
        !
        type( rScalar3D_SG_t ) :: SigmaCell
        integer :: i, j, k, k0, k1, k2
        !
        type( rVector3D_SG_t ) :: length, area
        !
        select type( grid => self%grid )
            class is( Grid3D_SG_t )
                !
                allocate( eVec, source = rVector3D_SG_t( grid, EDGE ) )
                !
				call eVec%zeros()
				!
                SigmaCell = rScalar3D_SG_t( grid, CELL )
                !
                k0 = self%grid%nzAir
                k1 = k0 + 1
                k2 = self%grid%Nz
                SigmaCell%v(:, :, 1:k0) = self%airCond
                !
                ! Note: AirCond should always be in linear domain, but conductivity
                ! in cells is generally transformed -- SigMap converts to linear
                SigmaCell%v(:, :, k1:k2) = self%SigMap(self%cellCond%v)
                !
                ! Form Conductivity--cell volume product  -- now using Vcell from MetricElements
                call sigmaCell%mults( self%metric%Vcell )
                !
                ! Sum onto edges
                call eVec%SumCells( SigmaCell )
                !
                ! Divide by total volume -- sum of 4 cells
                ! surrounding edge -- just 4*V_E        
                call eVec%divs( self%metric%Vedge )
                !
                !  still need to divide by 4 ...
                !call eVec%mults( 0.25_prec )
                !
            class default
                write(*, *) "ERROR:ModelParameterCell_SG:PDEmapping:"
                stop "              Incompatible grid. Exiting."
                !
        end select
        !
    end function PDEmappingModelParameterCell
    
    !**
    ! PDE mapping linearized at background model
    ! parameter m0, applied to dm result is an edge-vector eVec.
    !*
    function dPDEmappingModelParameterCell( self, dm ) result( eVec )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in )        :: dm
        class( rVector_t ), allocatable                :: eVec
        !
        type( rScalar3D_SG_t ) :: SigmaCell
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: k0, k1, k2
        !
        select type( dm )
               class is( modelParameterCell_SG_t )
                    !
                    select type( grid => self%grid )
                    class is( Grid3D_SG_t )
                          !
                          allocate( eVec, source = rVector3D_SG_t( grid, EDGE ) )
                          !
                          SigmaCell = rScalar3D_SG_t( grid, CELL )
                        
                          ! Set Earth cells using m0, SigMap and dm
                          ! I am doing this explicitly -- could make SigmaCell on ParamGrid
                          ! then move Earth part to a Vector on ModelGrid (this is how we
                          ! would do this more generally, when model space was really different
                          ! from modeling grid.
                          k0 = self%ParamGrid%NzAir
                          k1 = k0 + 1
                          k2 = self%ParamGrid%Nz
                        
                          call SigmaCell%zeros()     !    need to zero to make sure values in air are zero
                          SigmaCell%v(:,:,k1:k2) = self%SigMap(self%cellCond%v, JOB)
                          SigmaCell%v(:,:,k1:k2) = SigmaCell%v(:,:,k1:k2)*dm%cellCond%v
                        
                          ! Average onto edges, as in PDEmapping ...
                          !
                          call sigmaCell%multS( self%metric%Vcell )
                          ! Sum onto edges
                          call eVec%SumCells( SigmaCell )
                          !
                          ! Divide by total volume -- sum of 4 cells
                          ! surrounding edge -- just 4*V_E
                          call eVec%divs( self%metric%Vedge )
                          !  still need to divide by 4 ...
                          call evec%mults( 0.25_prec )
                          !
                    class default
                        write(*, *) "ERROR:ModelParameterCell_SG:dPDEmapping:"
                        stop "              Incompatible grid. Exiting."
                        !
                    end select
                    !
               class default
                    write(*, *) "ERROR:ModelParameterCell_SG:dPDEmapping:"
                    stop "              Incompatible input [dm]. Exiting."
                    !
        end select
        !
    end function dPDEmappingModelParameterCell
    !**
    ! Transpose (adjoint) of dPDEmapping, applied to an edge-vector eVec
    ! result is a model parameter dm.
    !*
    function dPDEmappingTModelParameterCell( self, eVec ) result( dm )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( rVector_t ), intent( in )               :: eVec
        class( ModelParameter_t ), allocatable         :: dm
        !
        type( rScalar3D_SG_t ) :: sigmaCell
        type( rVector3D_SG_t ) :: vTemp
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: k0, k1, k2
        !
        !
        select type( param_grid => self%paramGrid )
            class is( Grid3D_SG_t )
                !
                dm = ModelParameterCell_SG_t( param_grid, self%cellCond, self%paramType )
                !
            class default
                   write(*, *) "ERROR:ModelParameterCell_SG:dPDEmappingT:"
                   stop "              Unknow grid"
        end select
        !
        select type( eVec )
              class is( rVector3D_SG_t )
                !
                select type( dm )
                      !
                      class is( ModelParameterCell_SG_t )
                           ! Create local temporary scalar and vector
                           vTemp = eVec%interior()
                           ! Divide by total volume -- sum of 4 cells
                           ! surrounding edge -- just 4*V_E
                           call vTemp%divs(self%metric%Vedge)
                           !  still need to divide by 4 ...
                           call vTemp%mults(0.25_prec)

                           sigmaCell = vTemp%SumEdges()
                           !
                           !deallocate( vTemp )
                           !
                           call sigmaCell%multS( self%metric%Vcell )

                           dm%cellCond%v = self%SigMap( self%cellCond%v,JOB ) 

                           k0 = self%ParamGrid%NzAir
                           k1 = k0 + 1
                           k2 = self%ParamGrid%Nz
                    
                           !     deleted select type for SigmaCell -- can"t see how we need this, since
                           !        this is local variable declared with explicit type!
                           dm%cellCond%v = dm%cellCond%v * SigmaCell%v(:,:,k1:k2)
                           !
                           !deallocate( SigmaCell )
                           !
                      class default
                           write(*, *) "ERROR:ModelParameterCell_SG:dPDEmappingT:"
                           stop "              Incompatible input [eVec]. Exiting."
                end select
        end select
        !
    end function dPDEmappingTModelParameterCell
    !**
    ! SetType
    !*
    subroutine setTypeModelParameterCell( self, paramType )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        character(:), allocatable, intent( in )           :: paramType
        !
        if (.NOT.(self%is_allocated)) then
              write(*, *) "ERROR:ModelPArameterCell_t:SetType:"
              stop "              Not allocated."
        end if
        !  NOTE: always keep AirCond linear (actual conductivity)
        !    parameter transformation is only needed for inversion,
        !    and we do not invert for AirCond!
        
        if (trim(paramType) .eq. trim(self%paramType)) then
              ! We are done
        else if (self%paramType == "") then
              self%paramType = trim(paramType)
        else if (self%paramType == LINEAR) then
              ! Convert to log
              if (paramType == LOGE) then
                self%cellCond%v = log(self%cellCond%v)
                !self%airCond = log(self%AirCond)
              else if (paramType == LOG_10) then
                self%cellCond%v = log10(self%cellCond%v)
                !self%airCond = log10(self%airCond)
              end if
        else if (paramType == LINEAR) then
              ! Convert from log to linear
              if (self%paramType == LOGE) then
                self%cellCond%v = exp(self%cellCond%v)
                !self%airCond = exp(self%airCond)
              else if(self%paramType == LOG_10) then
                self%cellCond%v = exp(self%cellCond%v * log(10.))
                !self%airCond = exp(self%AirCond * log(10.))
              end if
        else if ((self%paramType == LOGE) .and. (paramType == LOG_10)) then
              ! Convert from natural log to log10
              self%cellCond%v = self%cellCond%v / log(10.)
              !self%airCond = self%airCond / log(10.)
        else if ((self%paramType == LOG_10) .and. (paramType == LOGE)) then
              ! Convert from log10 to natural log
              self%cellCond%v = self%cellCond%v * log(10.)
              !self%airCond = self%airCond * log(10.)
        else
              write(*, *) "ERROR:ModelParameterCell_t:SetType:"
              stop "              Unknown paramType."
        end if
        !
        self%paramType = paramType 
        !
    end subroutine setTypeModelParameterCell
    
end Module ModelParameterCell_SG

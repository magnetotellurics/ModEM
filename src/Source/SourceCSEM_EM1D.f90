!
!> Derived class to define a CSEM Source with E_p computed using Dipole1D
!
module SourceCSEM_EM1D
    !
    use EM1D
    !
    use Constants
    use cVector3D_SG
    use rVector3D_SG
    use Grid3D_SG
    use Source
    use ModelOperator
    use ModelParameterCell_SG
    !
    type, extends( Source_t ) :: SourceCSEM_EM1D_t
        !
        integer :: i_tx
        !
        real( kind=prec ) :: location(3)
        !
        type( rVector3D_SG_t ) :: cond_anomaly_h, cond_anomaly_v
        !
        contains
            !
            final :: SourceCSEM_EM1D_dtor
            !
            procedure, public :: createE => createE_SourceCSEM_EM1D
            !
            procedure, public :: createRHS => createRHS_SourceCSEM_EM1D
            !
            procedure, private :: create_Ep_from_EM1D
            !
    end type SourceCSEM_EM1D_T
    !
    interface SourceCSEM_EM1D_t
        module procedure SourceCSEM_EM1D_ctor
    end interface SourceCSEM_EM1D_t
    !
contains
    !
    !> SourceCSEM_EM1D constructor
    function SourceCSEM_EM1D_ctor( model_operator, sigma, period, location, i_tx ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period, location(3)
        integer, intent( in ) :: i_tx
        !
        type( SourceCSEM_EM1D_t ) :: self
        !
        !write( *, * ) "Constructor SourceCSEM_EM1D_t"
        !
        call self%init
        !
        self%model_operator => model_operator
        !
        self%sigma => sigma
        !
        self%period = period
        !
        self%location = location
        !
        self%i_tx = i_tx
        !
        self%non_zero_source = .TRUE.
        !
        self%non_zero_bc = .FALSE.
        !
    end function SourceCSEM_EM1D_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    subroutine SourceCSEM_EM1D_dtor( self )
        implicit none
        !
        type( SourceCSEM_EM1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor SourceCSEM_EM1D_t"
        !
        call self%dealloc
        !
    end subroutine SourceCSEM_EM1D_dtor
    !
    !> Set self%E from forward modeling 1D
    subroutine createE_SourceCSEM_EM1D( self )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        !
        type( backgrounddata ) :: bgdat                    !model description, coordinates, data
        type( sorec ) :: src                               !source specification
        type( sorec ), dimension(:), pointer :: receivers  !receiver specification
        type( freqdata ) :: freqdat                        !frequency dependent specifications
        type( refl_struct ) :: refl_var                    !all variables that have to be remembered while computing 1D fields
        !
        integer :: ifreq, icur, comm, ix, iy, iz
        complex( kind=prec ) :: i_omega_mu 
        integer :: status
        !
        call set1DModel_VTI( self%sigma, self%location(1), self%location(2) )
        !
        call setAnomConductivity_VTI( self%sigma )
        !
        bgdat%omega = ( 2.0 * PI / self%period )
        bgdat%dowhat = 1
        call create_background_data( self%sigma%metric%grid, bgdat )
        call create_source_data( self%i_tx, src, freqdat )
        !
        ifreq=1
        icur=1
        comm=1
        refl_var%nzrecHxy=0
        refl_var%nzrecHz=0
        call reflectivity_unified( src, bgdat, refl_var, ifreq, icur, comm ) ! Output field will be saved in bgdat
        call self%create_Ep_from_EM1D( self%sigma%metric%grid, bgdat )            ! Put the 1D field into E_p
        !
        i_omega_mu = cmplx( 0., real( -1.0d0 * isign * mu_0 * ( 2.0 * PI / self%period ), kind=prec ), kind=prec )
        !
        allocate( self%E(1), source = E_p )
        !
        !> Fill e_vector (cVector3D_SG) from E2D (Esoln2DTM_t)
        select type( E => self%E(1) )
            !
            class is( cVector3D_SG_t )
                !
                !> Fill e_vector (cVector3D_SG) from E2D (Esoln2DTM_t)
                select type( E_P )
                    !
                    class is( cVector3D_SG_t )
                        !
                        E%x = self%cond_anomaly_h%x * E_P%x
                        E%y = self%cond_anomaly_h%y * E_P%y
                        E%z = self%cond_anomaly_v%z * E_P%z
                        !
                    class default
                        stop "createE_SourceCSEM_EM1D > Unclassified E_P"
                end select
                !
            class default
                stop "createE_SourceCSEM_EM1D > Unclassified E"
        end select
        !
        call self%E(1)%mult( i_omega_mu )
        !
        deallocate( bgdat%sigv, STAT=status )
        deallocate( bgdat%sigh, STAT=status )
        deallocate( bgdat%epsrv, STAT=status )
        deallocate( bgdat%epsrh, STAT=status )
        deallocate( bgdat%zbound, STAT=status )
        deallocate( bgdat%Expos, STAT=status )
        deallocate( bgdat%Eypos, STAT=status )
        deallocate( bgdat%Ezpos, STAT=status )
        deallocate( bgdat%Ex, STAT=status )
        deallocate( bgdat%Ey, STAT=status )
        deallocate( bgdat%Ez, STAT=status )
        !
        deallocate( src%nelem, STAT=status )
        deallocate( src%pos, STAT=status )
        deallocate( src%ljx, STAT=status )
        deallocate( src%ljy, STAT=status )
        deallocate( src%ljz, STAT=status )
        deallocate( src%akx, STAT=status )
        deallocate( src%aky, STAT=status )
        deallocate( src%akz, STAT=status )
        !
        call self%createRHS
        !
    end subroutine createE_SourceCSEM_EM1D
    !
    !> No subroutine briefing
    !
    subroutine create_Ep_from_EM1D( self, grid, bgdat )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        class( Grid_t ), intent( in ) :: grid
        type( backgrounddata ), intent( in ) :: bgdat
        !
        integer ix, iy, iz, counter
        !
        allocate( E_p, source = cVector3D_SG_t( grid, EDGE ) )
        !
        !> Fill e_vector (cVector3D_SG) from E2D (Esoln2DTM_t)
        select type( E_P )
            !
            class is( cVector3D_SG_t )
                !
                counter = 1
                ! E-field corresponing to these nodes is Ex
                do iz = 1,grid%Nz+1    !Edge Z
                    do iy = 1,grid%Ny+1     !Edge Y
                        do ix = 1,grid%Nx       !Center X
                            E_p%x(ix,iy,iz) = bgdat%Ex(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
                counter = 1
                ! E-field corresponing to these nodes is Ey
                do iz = 1, grid%Nz+1    !Edge Z
                    do iy = 1, grid%Ny      !Center y
                        do ix = 1, grid%Nx+1    !Edge x
                            E_p%y(ix,iy,iz) = bgdat%Ey(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
                counter = 1
                ! E-field corresponing to these nodes is Ez
                do iz = 1,grid%Nz !Center Z
                    do iy = 1,grid%Ny+1 !Edge y
                        do ix = 1,grid%Nx+1 !Edge x
                            E_p%z(ix,iy,iz) = bgdat%Ez(counter)
                            counter = counter + 1
                        enddo
                    enddo
                enddo
                !
            class default
                stop "create_Ep_from_EM1D > Unclassified E_P"
        end select
        !
    end subroutine create_Ep_from_EM1D
    !
    !> Set RHS from self%E
    !
    subroutine createRHS_SourceCSEM_EM1D( self )
        implicit none
        !
        class( SourceCSEM_EM1D_t ), intent( inout ) :: self
        !
        if( allocated( self%rhs ) ) deallocate( self%rhs )
        allocate( cVector3D_SG_t :: self%rhs(1) )
        !
        self%rhs(1) = self%E(1)
        !
        call self%rhs(1)%mult( self%model_operator%metric%Vedge )
        !
    end subroutine createRHS_SourceCSEM_EM1D
    !
    !>    This is a utility routine, used by several data functional
    !>    set up routines, and for other interpolation functions
    !>    Returns index ix such that    xNode(ix) <= x < xNode(ix+1)
    !>    If x is out of range:
    !>    x < xNode(1) returns 0; if x> xNode(nx) returns nx
    !>    Assumes xNode is strictly increasing; does not check this
    !>    NOTE: as presently coded, when xNode is called with center
    !>    (face) node positions, this routine will return zero for
    !>    the coordinates in the outer half cell nearest the boundary
    !>    If evaluation over the complete model domain is to be allowed
    !>    a more general interpolation rule will be required.
    !>    A.K.: modified to allow input of any size, nx = size(xNode).
    !
    function minNode( x, xNode ) result( ix )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( in ) :: xNode
        !
        integer :: ix, i
        !
        do i = 1, size( xNode )
            if( clean( xNode(i) ) .GT. clean(x) ) then
                ix = i-1
                exit
            endif
        enddo
        !
    end function minNode
    !
    !>    This is a utility routine, used by several data functional
    !>    set up routines, and for other interpolation functions
    !>    Returns index ix such that    xNode(ix) <= x < xNode(ix+1)
    !>    If x is out of range:
    !>    x > xNode(1) returns 0; if x< xNode(nx) returns nx
    !>    Assumes xNode is strictly decreasing; does not check this
    !>    NOTE: as presently coded, when xNode is called with center
    !>    (face) node positions, this routine will return zero for
    !>    the coordinates in the outer half cell nearest the boundary
    !>    If evaluation over the complete model domain is to be allowed
    !>    a more general interpolation rule will be required.
    !>    A.K.: modified to allow input of any size, nx = size(xNode).
    !
    function maxNode(x, xNode) result(ix)
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( in ) :: xNode
        !
        integer :: ix, i
        !
        do i = 1, size(xNode)
           if( clean( xNode(i)) .LT. clean(x) ) then
                ix = i-1
                exit
           endif
        enddo
        !
    end function maxNode
    !
    !> This is a utility routine that provides an expression used to battle
    !> against machine error problems. It returns the same real or real(8)
    !> as the input, but without the extra digits at the end that are often
    !> a cause of wrong comparisons in the if statements. ALWAYS use clean(x)
    !> instead of x in an inequality!!!
    !> R_LARGE is defined in the module math_constants
    !> A.K.
    !
    function clean(x)
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ) :: clean
        !
        clean = dnint(x*R_LARGE)/R_LARGE
        !
    end function clean
    !
end module SourceCSEM_EM1D

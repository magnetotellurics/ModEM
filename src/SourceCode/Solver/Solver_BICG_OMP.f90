!
!> module containing iterative equation solvers using OpenMP parallelization
!
module Solver_BICG_OMP
    !
    !use DeclarationMPI
    use FileUnits
    use Solver_CC
    use ModelOperator_MF_SG
    use ModelOperator_SP_V1
    use ModelOperator_SP_V2
    use PreConditioner_CC_MF_SG
    use PreConditioner_CC_SP_SG
    use PreConditioner_CC_SP_MR
    !
    use omp_lib
    !
    type, extends( Solver_CC_t ) :: Solver_BICG_OMP_t
        !
        integer :: OMP_integer_kind
        !integer( kind=OMP_integer_kind ), private :: Nthreads, Ithread
        integer :: Nthreads, Ithread
        complex( kind=prec ), allocatable, dimension(:) :: VomegaMuSig2
        !
    contains
        !
        procedure, public :: solve => solve_Solver_BICG_OMP
        !
        procedure, public :: BICG_OMP_ADJ
        procedure, public :: BICG_OMP_FWD
        procedure, public :: Solver_BICG_OMP_time
        !
    end type Solver_BICG_OMP_t
    !
    interface Solver_BICG_OMP_t
        module procedure Solver_BICG_OMP_ctor
    end interface Solver_BICG_OMP_t
    !
contains
    !
    !> No subroutine briefing
    !
    function Solver_BICG_OMP_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        !
        type( Solver_BICG_OMP_t ) :: self
        !
        !write( *, * ) "Constructor Solver_BICG_OMP_t"
        !
        call self%baseInit
        !
        self%OMP_integer_kind = 4
        !
        !> Instantiate the PreConditioner object according to the ModelOperator type
        select type( grid => model_operator%metric%grid )
            !
            class is( Grid3D_SG_t )
                !
                !> Instantiate the PreConditioner object according to the ModelOperator type
                select type( model_operator )
                    !
                    class is( ModelOperator_MF_SG_t )
                        !
                        allocate( self%preconditioner, source = PreConditioner_CC_MF_SG_t( model_operator ) )
                        !
                    class is( ModelOperator_SP_t )
                        !
                        allocate( self%preconditioner, source = PreConditioner_CC_SP_SG_t( model_operator ) )
                        !
                    class default
                        call errStop( "Solver_BICG_OMP_ctor > Unclassified SG model_operator" )
                    !
                end select
                !
            class is( Grid3D_MR_t )
                !
                !> Instantiate the PreConditioner object according to the ModelOperator type
                select type( model_operator )
                    !
                    class is( ModelOperator_MF_SG_t )
                        !
                        call errStop( "Solver_BICG_OMP_ctor > For MR use model_operator SP" )
                        !
                    class is( ModelOperator_SP_t )
                        !
                        allocate( self%preconditioner, source = PreConditioner_CC_SP_MR_t( model_operator ) )
                        !
                    class default
                        call errStop( "Solver_BICG_OMP_ctor > Unclassified MR model_operator" )
                    !
                end select
                !
            class default
                call errStop( "Solver_BICG_OMP_ctor > Unclassified grid" )
            !
        end select
        !
        call self%set( max_solver_iters, tolerance_solver )
        !
        call self%zeroDiagnostics
        !
    end function Solver_BICG_OMP_ctor
    !
    !> No subroutine briefing
    !
    subroutine Solver_BICG_OMP_time( it )
        implicit none
        !
        integer, intent( in ) :: it
        !
        real( kind=prec ), save :: ttime
        !
        ttime = OMP_get_wtime() ! start time
        !
        if( it == 1 ) return
        !
        ttime = OMP_get_wtime() - ttime ! end time
        !
        write( ioNodeInfo, "(a5,i3.3,a4)") "node[", taskid, "]:  "
        write( 6, 1 ) ioNodeInfo, ttime
        !
        1 format( a, "runtime(sec)=", f12.4 )
        !
    end subroutine Solver_BICG_OMP_time
    !
    !> No subroutine briefing
    !
    subroutine Solver_BICG_OMP_info( self, Nt_num )
        implicit none
        !
        class( Solver_BICG_OMP_t ), intent( inout ) :: self
        integer, intent( out ) :: Nt_num
        !
        character*10 :: cinfo
        !> Uncomment for method 1
        !call get_environment_variable("OMP_NUM_THREADS ",cinfo)
        !read(cinfo,fmt=*,iostat=ierr)Nt_num
        !if(ierr/=0)then
        !>  Nt_num=0
        !>  return
        !endif
        !self%Nthreads=OMP_get_num_threads()
        self%Nthreads = OMP_get_max_threads()
        Nt_num = int( self%Nthreads )
        !write( *, * )"Solver_BICG_OMP_info: num=",Nt_num
        !
    end subroutine Solver_BICG_OMP_info
    !
    !> Combined PC_Lsolve,PC_Usolve
    !> implement the sparse matrix solve for curl-curl operator
    !
    subroutine solve_Solver_BICG_OMP( self, b, x )
        implicit none
        !
        class( Solver_BICG_t ), intent(inout) :: self
        class( Vector_t ), intent(in) :: b
        class( Vector_t ), intent(inout) :: x
        !
        class( Vector_t ), allocatable :: R, RT, V, T
        class( Vector_t ), allocatable :: P, PT, PH, S, ST, SH, AX
        class( Vector_t ), allocatable :: xhalf, xmin
        real( kind=prec ) :: rnorm, bnorm, rnormin, btol
        complex( kind=prec ) :: RHO, ALPHA, BETA, OMEGA
        complex( kind=prec ) :: RTV, TT, RHO1
        integer :: iter, imin
        logical :: adjoint, ilu_adjt
        !
        if( .NOT. x%is_allocated ) then
            call errStop( "solve_Solver_BICG > x not allocated yet" )
        endif
        !
        if( .NOT. b%is_allocated ) then
            call errStop( "solve_Solver_BICG > b not allocated yet" )
        endif
        !
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, xhalf )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, xmin)
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, AX )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, R )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, RT )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, P )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, PT )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, PH )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, S )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, ST )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, SH )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, V )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, T )
        !
        self%iter = 1
        !
        !> Norm of rhs, residual
        bnorm = SQRT( b%dotProd( b ) )
        !
        !> this usually means an inadequate model, in which case Maxwell"s fails
        if( isnan( abs( bnorm ) ) ) then
            !
            call errStop( "solve_Solver_BICG > b contains NaNs" )
            !
        elseif( bnorm .EQ. 0.0 ) then ! zero rhs -> zero solution
            !
            call warning( "b in BICG has all zeros, returning zero solution" )
            !
            x = b
            !
            self%iter = 1
            self%n_iter = 1
            self%relErr(1) = 0.0
            !
            return
            !
        endif
        !
        !> now calculate the (original) residual
        adjoint = .FALSE.
        !
        !call self%preconditioned%model_operator%multA_N( x, R, adjoint )
        call self%preconditioner%model_operator%amult( x, R, self%omega, adjoint )
        !
        !> R= b - Ax, for initial guess x, that has been inputted to the routine
        rnorm = CDSQRT( R%dotProd( R ) )
        !
        call R%linComb( b, C_MinusOne, C_ONE )
        !
        !> Norm of residual
        rnorm = CDSQRT( R%dotProd( R ) )
        !
        btol = self%tolerance * bnorm
        !
        if( rnorm .LE. btol ) then ! the first guess is already good enough
            !
            self%n_iter = 1
            !
            call warning( "solve_Solver_BICG > the first guess is already good enough" )
            !
            self%relErr(1) = rnorm / bnorm
            !
            return
            !
        endif
        !
        !> ================= Now start configuring the iteration ===================
        !
        !> the adjoint (shadow) residual
        !
        rnormin = rnorm
        xmin = x
        self%relErr(1) = real( rnormin / bnorm )
        !
        self%converged = .FALSE.
        imin = 1
        RHO = C_ONE
        OMEGA = C_ONE
        RT = R ! use the overloaded =
        !
        !============================== looooops! ================================
        !
        do iter = 2, self%max_iters
            !
            self%iter = iter
            !
            RHO1 = RHO
            RHO = RT%dotProd( R )
            !
            if( RHO .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > RHO .EQ. 0.0" )
                !
            endif
            !
            if( self%iter .EQ. 2 ) then
                P = R
            else
                !
                BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
                !
                if( BETA .EQ. 0.0 ) then
                    !
                    call errStop( "solve_Solver_BICG > BETA .EQ. 0.0" )
                    !
                endif
                !
                !> P= R + BETA * (P - OMEGA * V);
                call P%linComb( V, C_One, -OMEGA )
                call P%linComb( R, BETA, C_One )
                !
            endif
            !
            !> L
            ilu_adjt = .FALSE.
            call PT%zeros
            call self%preconditioner%LTsolve( P, PT , ilu_adjt )
            !
            ! U
            ilu_adjt = .FALSE.
            call PH%zeros
            call self%preconditioner%UTsolve( PT, PH , ilu_adjt )
            !
            ! PH = P
            adjoint = .FALSE.
            !
            call V%zeros
            call self%preconditioner%model_operator%amult( PH, V, self%omega, adjoint )
            !
            RTV = RT%dotProd( V )
            !
            if( RTV .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > RTV .EQ. 0.0" )
                !
            endif
            !
            ALPHA = RHO / RTV
            !
            if( ALPHA .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > ALPHA .EQ. 0.0" )
                !
            endif
            !
            !> xhalf = x + ALPHA*PH ! the first half of iteration
            xhalf = x
            call xhalf%linComb( PH, C_One, ALPHA )
            !
            adjoint = .FALSE.
            !
            call AX%zeros
            call self%preconditioner%model_operator%amult( xhalf, AX, self%omega, adjoint )
            !
            call AX%linComb( b, C_MinusOne, C_One )
            rnorm = CDSQRT( AX%dotProd( AX ) )
            !
            if( rnorm .LT. btol ) then
                !
                x = xhalf
                self%n_iter = self%iter
                self%converged = .TRUE.
                self%relErr( self%iter ) = real( rnorm / bnorm )
                !
                exit
                !
            endif
            !
            if( rnorm .LT. rnormin) then
                !
                rnormin = rnorm
                xmin = xhalf
                imin = self%iter
                !
            endif
            !
            !> S = R - ALPHA*V  !residual for the 0.5 x
            S = R
            call S%linComb( V, C_One, -ALPHA )
            !
            !> L
            ilu_adjt = .FALSE.
            call self%preconditioner%LTsolve( S, ST, ilu_adjt )
            !
            !> U
            ilu_adjt = .FALSE.
            call self%preconditioner%UTsolve( ST, SH, ilu_adjt )
            !
            !> SH = S
            adjoint = .FALSE.
            call T%zeros
            call self%preconditioner%model_operator%amult( SH, T, self%omega, adjoint )
            !
            TT = T%dotProd( T )
            !
            if( TT .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > TT .eq. 0.0" )
                !
            endif
            !
            OMEGA = T%dotProd(S) / TT
            !
            if( OMEGA .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > OMEGA .eq. 0.0" )
                !
            endif
            !
            !> x = xhalf + OMEGA * SH  ! the second half (shadow) of iteration
            x = xhalf
            call x%linComb( SH, C_One, OMEGA )
            !
            adjoint = .FALSE.
            !
            call AX%zeros
            call self%preconditioner%model_operator%amult( x, AX, self%omega, adjoint )
            !
            call AX%linComb( b, C_MinusOne, C_One )
            !
            rnorm = CDSQRT( AX%dotProd( AX ) )
            !
            self%relErr( self%iter ) = real( rnorm / bnorm )
            !
            if( rnorm .LT. btol ) then
                !
                self%n_iter = self%iter
                self%converged = .TRUE.
                !
                exit
                !
            endif
            !
            if( rnorm .LT. rnormin) then
                !
                rnormin = rnorm
                xmin = x
                imin = self%iter
                !
            endif
            !
            !R = S - OMEGA * T  !residual for the 1.0 x
            R = S
            call R%linComb( T, C_One, -OMEGA )
            !
            !> Verbose
            !write( *, "( a36, i6, a3, es12.3 )" ) "BICG iter: ", self%iter, " : ", self%relErr( self%iter )
            !
        enddo
        !
        if( self%converged ) then
            write( *, "( a52, i6, a7, es12.3 )" ) "->Solver BICG converged within ", self%iter, ": err= ", self%relErr( self%iter )
        else
            write( *, "( a52, i6, a7, es12.3 )" ) "->Solver BICG not converged in ", self%max_iters, ": err= ", self%relErr( self%max_iters )
        endif
        !
        if( .NOT. self%converged ) then
            ! it should be noted that this is the way my matlab version works
            ! the bicg will return the 'best' (smallest residual) iteration
            x = xmin; !comment this line
            self%n_iter = self%max_iters
            self%relErr( self%max_iters ) = self%relErr( imin) ! and this line
            ! to use the last iteration result instead of the 'best'
        endif
        !
        deallocate( xhalf )
        deallocate( xmin)
        deallocate( AX )
        deallocate( R )
        deallocate( RT )
        deallocate( P )
        deallocate( PT )
        deallocate( PH )
        deallocate( S )
        deallocate( ST )
        deallocate( SH )
        deallocate( V )
        deallocate( T )
        !
        ! ################# MICHAEL IMPLEMENTATION #################
        !
        ! class( Vector_t ), intent( in ) :: b
        ! class( Vector_t ), intent( inout ) :: x
        ! !
        ! !logical, intent( in ):: adjt
        ! integer :: i, j, j2
        ! !> Adjoint call:
        ! !> 1) PC_Lsolve(b,adjt,x) -> call UTsolve_Cmplx(LH,b,x) = subroutine UTsolve_Cmplx(U,b,b)
        ! !>             P     PT                           P PT
        ! !> 2) PC_Usolve(b,adjt,x) -> call LTsolve_Cmplx(UH,b,x)
        ! !>            PT     PH                          PT PH
        ! if( adjt ) then ! 1) UTsolve_Cmplx(LH,b,x), 2) LTsolve_Cmplx(UH,b,x)
            ! !
            ! do i = LH%nRow,1,-1 ! UTsolve_Cmplx(LH,b,x)
                ! x = b
                ! j2 = LH%row(i)
                ! do j = j2+1,LH%row(i+1)-1
                ! x = x-LH%val(j)*x(LH%col(j))
                ! enddo
                ! x = x/LH%val(j2)
                ! enddo
                ! do i = 1,UH%nRow ! LTsolve_Cmplx(UH,b,b)
                ! j2=UH%row(i+1)-2
                ! do j = UH%row(i),j2
                ! x = x-UH%val(j)*x(UH%col(j))
                ! enddo
                ! x = x/UH%val(j2+1)
                ! !
            ! enddo
            ! !> Forward call:
            ! !> 1) PC_Lsolve(b,adjt,x) -> LTsolve_Cmplx(L,b,x) = subroutine LTsolve_Cmplx(L,b,b)
            ! !>             P     PT                           P PT
            ! !> 2) PC_Usolve(b,adjt,x) -> call UTsolve_Cmplx(U,b,x)
            ! !>            PT     PH                          PT PH
        ! else ! adjt=.FALSE.
            ! !
            ! do i = 1, L%nRow ! L: b=ST, b=S
                ! !
                ! x = b
                ! j2=L%row(i+1)-2
                ! !write(98,"(i7,a)")i," row"
                ! !
                ! do j = L%row(i), j2 !Lmat%row(i+1)-2
                    ! x = x-L%val(j)*x(L%col(j))
                    ! !write(98,"(i7,a)")Lmat%col(j)," col"
                ! enddo !j
                ! !
                ! x = x/L%val(j2+1)
                ! !
            ! enddo !i
            ! !
            ! do i = U%nRow, 1, -1 ! U: b=SH, b=ST
                ! !
                ! j2 = U%row(i)
                ! !
                ! do j = j2+1,U%row(i+1)-1
                    ! x = x-U%val(j)*x(U%col(j))
                ! enddo
                ! !
                ! x = x/U%val(j2)
                ! !
            ! enddo
            ! !
        ! endif
        ! !
    end subroutine solve_Solver_BICG_OMP
    !
    !> OpenMP solver versions
    !
    !> *****************************************************************************
    !> Solver contains subroutines for: 
    !> a) PCG- a quasi-generic pre-conditioned conjugate gradient, and 
    !> b) QMR - Quasi-Minimal Residual method (pre-conditioned, no look-ahead)
    !> c) BICG - bicgstab Stabilized Bi-conjugate Gradient method
    !
    !> Stabilized version of BiConjugate Gradient, set up for solving
    !> A x = b using routines in  mult_Aii.
    !> solves for the interior (edge) field
    !
    !> modified from my matlab version of BICGstab...
    !> so the naming might sound a little different from conventional ones
    !> also added the optional adjoint to solve adjoint system A^Tx = b
    !> 
    !> NOTE: BICG actually performs two sub (or half) line searches within a 
    !>       iteration, but here we only store the relerr for the second sub
    !>       just to be compatibility with QMR
    !
    !> 
    !> interface...........
    !> redefining some of the interfaces for our convenience (locally)
    !> generic routines for vector operations for edge/ face nodes
    !> in a staggered grid
    !
    !> NOTE: this has not been extensively tested - I believe it feels a 
    !> little unstable (despite the name)...
    !> if you have time reading this, test it!
    !> USE omp_lib
    !
    subroutine BICG_OMP_FWD( self, b, x )
        implicit none
        !
        class( Solver_BICG_OMP_t ), intent( inout ) :: self
        !> b is right hand side
        class( Vector_t ), intent( in ) :: b
        !
        !> solution vector is x ... on input is provided with the initial
        !> guess, on output is the iterate with smallest residual. 
        class( Vector_t ), intent( inout ) :: x
        !
        type( spMatCSR_real ) :: Mii
        complex( kind=prec ) :: R, RT, V, T, P, S, SH
        complex( kind=prec ) :: xhalf, xmin
        real( kind=prec ) :: rnorm, bnorm, rnormin, btol
        complex( kind=prec ) :: RHO, ALPHA, BETA, OMEGA, AX1
        complex( kind=prec ) :: RTV, TT, RHO1
        integer :: iter, imin, maxiter
        logical :: lmem(2)
        !
        !> CMC: For new OMP version
        integer :: i, j, j2
        complex( kind=prec ) :: cnorm
        !
        R = x%getArray(); R = C_ZERO
        RT = x%getArray(); RT = C_ZERO
        V = x%getArray(); V = C_ZERO
        T = x%getArray(); T = C_ZERO
        P = x%getArray(); P = C_ZERO
        S = x%getArray(); S = C_ZERO
        SH = x%getArray(); SH = C_ZERO
        xhalf = x%getArray(); xhalf = C_ZERO
        xmin = x%getArray(); xmin = C_ZERO
        !
        self%VomegaMuSig2 = x%getArray(); self%VomegaMuSig2 = C_ZERO
        !
        lmem = .FALSE.
        lmem(1) = .TRUE.
        !
        cnorm = C_ZERO; TT = C_ZERO
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%preconditioner%model_operator )
            !
            class is( ModelOperator_SP_t )
                !
                !> Square Matrix Mii: differentiated Differentiation between SP1 (CCii) and SP2 (AAii)
                select type( model_operator )
                    !
                    class is( ModelOperator_SP_V1_t )
                        !
                        Mii = model_operator%CCii
                        !
                    class is( ModelOperator_SP_V2_t )
                        !
                        Mii = model_operator%AAii
                        !
                    class default
                        call errStop( "BICG_OMP_FWD > model_operator must be SP V1 or V2" )
                    !
                end select
                !
                !$OMP PARALLEL PRIVATE(self%Ithread)
                !> self%Nthreads=OMP_get_num_threads(); self%Ithread=OMP_get_thread_num()
                !> if(self%Ithread==0)write( *, * )"OpenMP2 threads num,max=",self%Nthreads,OMP_get_max_threads()
                !
                !$OMP DO REDUCTION(+:cnorm,TT) PRIVATE(j)
                do i = 1, Mii%nRow
                    !
                    self%VomegaMuSig2(i) = ONE_I * ISIGN * model_operator%VomegaMuSig(i)
                    !
                    R(i) = self%VomegaMuSig2(i) * x(i)
                    !
                    do j = Mii%row(i), Mii%row(i+1)-1 
                        R = R + Mii%val(j) * x( Mii%col(j) )
                    enddo
                    !
                    R(i) = b(i) - R(i)
                    cnorm = cnorm + R * conjg( R )
                    TT = TT + b * b%conjugate() ! bnorm = SQRT(dot_product(b, b))
                    !
                enddo
                !$OMP enddo
                !
                !$OMP END PARALLEL
                !
                !> Norm of residual
                rnorm = DREAL( CDSQRT( cnorm ) )
                bnorm = SQRT( DREAL( TT ) )
                !
                write( *, * ) "Initial bnorm=", bnorm
                !
                if( bnorm < 1d-16 ) then ! zero rhs -> zero solution
                    !
                    call warning( "b in BICG_OMP_FWD has all zeros, returning zero solution" )
                    x = b
                    self%n_iter = 1
                    self%relErr = R_ZERO
                    goto 9
                    !
                endif
                !
                btol = self%tolerance * bnorm
                !
                if( rnorm .LE. btol ) then ! the first guess is already good enough
                    !
                    !> returning
                    self%n_iter = 1
                    self%relErr(1) = real( rnorm/bnorm, kind=prec )
                    goto 9
                    !
                endif
                !
                !> allocate the rest here
                !allocate( xhalf(xsize), xmin(xsize) )
                !allocate( RT(xsize), P(xsize), S(xsize), SH(xsize), V(xsize), T(xsize) )
                lmem(2) = .TRUE.
                !
                !================= Now start configuring the iteration ===================!
                !> the adjoint (shadow) residual
                rnormin = rnorm
                self%relErr(1) = real( rnormin/bnorm, kind=prec )
                !> write(6,*) "initial residual",  self%relErr(1)
                self%converged = .FALSE.
                maxiter = self%max_iters 
                imin = 0
                RHO = C_ONE; OMEGA = C_ONE
                RT = R
                xmin = x
                imin = 1
                !> subroutine Mult_Aii: File FWD_SP2/modelOperator3D.f90, L. 449
                !> call RMATxCVEC(Mii,x,y)
                !> subroutine RMATxCVEC: File SP_Topology/spOpTools.f90, L. 404
                !
                !============================== looooops! ================================!
                !write( *, * )Lmat%nRow,Umat%nRow,Mii%nRow,"size" ! all equal
                !
                do iter = 1, maxiter
                    !
                    RHO1 = RHO; RHO = C_ZERO
                    !
                    !$OMP PARALLEL
                    !$OMP DO REDUCTION(+:RHO)
                    do i = 1, Mii%nRow
                        RHO = RHO + conjg( RT ) * R
                    enddo
                    !
                    !$OMP enddo
                    if( iter == 1 )then
                        BETA = C_ZERO
                    else
                        BETA = (RHO/RHO1) * (ALPHA/OMEGA) 
                    endif
                    !
                    !$OMP DO
                    do i = 1, Mii%nRow
                        !
                        P = R + BETA * ( P - OMEGA * V )
                        V = P
                        !
                    enddo
                    !
                    !$OMP enddo
                    !$OMP END PARALLEL
                    !
                    !if(cdabs(RHO) < 1d-16) then
                    !>   BICGiter%failed = .TRUE.
                    !>   exit
                    !endif 
                    !call M1solve(P,ilu_adjt,PT) ! P->x  PT->y
                    !call M2solve(PT,ilu_adjt,PH)
                    !call solve_Solver_BICG_OMP( .FALSE., V, SH )
                    call self%solve( .FALSE., V, SH )
                    !
                    RTV = C_ZERO; cnorm = C_ZERO
                    !
                    !$OMP PARALLEL PRIVATE(self%Ithread)
                    self%Ithread = OMP_get_thread_num()
                    !
                    !$OMP DO REDUCTION(+:RTV) PRIVATE(j)
                    do i = 1, Mii%nRow
                        !
                        V = self%VomegaMuSig2(i) * SH
                        !
                        do j = Mii%row(i), Mii%row(i+1) - 1 
                            V = V + Mii%val(j) * SH !SH(Mii%col(j)) ????
                        enddo
                        !
                        RTV = RTV + RT%conjugate() * V
                        !
                    enddo
                    !$OMP enddo
                    !
                    if( self%Ithread == 0 ) ALPHA = RHO / RTV
                    !
                    !$OMP BARRIER
                    !
                    !$OMP DO
                    do i = 1, Mii%nRow
                        xhalf = x + ALPHA * SH ! the first half of iteration
                    enddo
                    !$OMP enddo
                    !
                    !$OMP DO REDUCTION(+:cnorm) PRIVATE(j,AX1)
                    do i = 1, Mii%nRow
                        !
                        call AX1%setArray( self%VomegaMuSig2(i) * xhalf*getArray()
                        !
                        do j = Mii%row(i), Mii%row(i+1)-1 
                            AX1 = AX1 + Mii%val(j) * xhalf !xhalf(Mii%col(j)) ????
                        enddo
                        !
                        AX1 = b - AX1
                        cnorm = cnorm + AX1 * AX1%conjugate()
                        S = R - ALPHA * V !residual for the 0.5 x
                        !
                    enddo
                    !$OMP enddo
                    !$OMP END PARALLEL
                    !
                    if( RTV .EQ. 0.0 .OR. ALPHA .EQ. 0.0 ) then
                        !
                        exit
                        !
                    endif
                    !
                    cnorm = CDSQRT( cnorm ); rnorm = DREAL( cnorm )
                    !
                    self%relErr(iter) = rnorm / bnorm
                    !write( *, * ) "iter # ",iter," xhalf residual: ", self%relErr(iter)
                    !
                    if( rnorm .LT. btol ) then
                        !
                        x = xhalf
                        self%n_iter = iter
                        self%converged = .TRUE.
                        exit
                        !
                    endif
                    !
                    if( rnorm .LT. rnormin ) then
                        !
                        rnormin = rnorm
                        xmin = xhalf
                        imin = iter
                        !
                    endif
                    !
                    !> ========================== M1,M2solve ==========================
                    !> call M1solve(S,ilu_adjt,ST) = call PC_Lsolve(x,adjt,y) 
                    !>                                call LTsolve_Cmplx(L,x,y)
                    !>                          subroutine LTsolve_Cmplx(L,b,x)
                    !> map variables:                    S->x->b, ST->y->x
                    !> matrix = Lmat => L (lower triang.)
                    !call LTsolve_Cmplx(L,x,y) S
                    !!$OMP PARALLEL PRIVATE(self%Ithread)
                    !>        self%Ithread =OMP_get_thread_num()
                    !>        if(self%Ithread==0)write(*,*)Lmat%nRow,"Lmat%nRow"
                    !!$OMP DO SCHEDULE(STATIC)
                    !>        do i = 1,Lmat%nRow
                    !>           if(iter==1)then
                    !>              if(self%Ithread==0)write(*,*)i,"THREAD0"
                    !>              if(self%Ithread==1)write(*,*)i,"THREAD1"
                    !>           endif
                    !>        enddo
                    !
                    !call solve_Solver_BICG_OMP( .FALSE., S, SH )
                    call self%solve( .FALSE., S, SH )
                    !
                    TT = C_ZERO; OMEGA = C_ZERO; cnorm = C_ZERO
                    !
                    !$OMP PARALLEL PRIVATE(self%Ithread)
                    self%Ithread = OMP_get_thread_num()
                    !
                    !$OMP DO REDUCTION(+:TT,OMEGA) PRIVATE(j)
                    do i = 1, Mii%nRow
                        !
                        T = self%VomegaMuSig2(i) * SH
                        do j = Mii%row(i), Mii%row(i+1) - 1 
                            T = T + Mii%val(j) * SH !SH(Mii%col(j)) ????
                        enddo
                        !
                        TT = TT + T * T%conjugate() ! complex dot-product
                        !
                        OMEGA = OMEGA + T%conjugate() * S ! complex dot-product
                        !
                    enddo
                    !$OMP enddo
                    !
                    if( self%Ithread == 0 ) OMEGA = OMEGA / TT
                    !if(iter==1)print*,"OMEGA",self%Ithread,OMEGA
                    !
                    !$OMP BARRIER
                    !$OMP DO
                    do i = 1, Mii%nRow
                        x = xhalf + OMEGA * SH  ! the second half of iteration
                    enddo
                    !
                    !$OMP enddo
                    !$OMP DO REDUCTION(+:cnorm) PRIVATE(j,AX1)
                    do i = 1, Mii%nRow
                        !
                        AX1 = self%VomegaMuSig2(i) * x
                        !
                        do j = Mii%row(i), Mii%row(i+1) - 1 
                            AX1 = AX1 + Mii%val(j) * x !x( Mii%col(j) ) ????
                        enddo
                        !
                        AX1 = b - AX1
                        cnorm = cnorm + AX1 * AX1%conjugate()
                        R = S - OMEGA * T  !residual for the 1.0 x
                        !
                    enddo
                    !$OMP enddo
                    !$OMP END PARALLEL
                    !
                    if( TT .EQ. 0.0 .OR. OMEGA .EQ. 0.0 ) then
                        !
                        exit
                        !
                    endif
                    !
                    rnorm = CDSQRT( cnorm )
                    self%relErr(iter) = real( rnorm / bnorm, kind=prec )
                    !
                    !> write(6,*) "iter # ",iter," x residual: ", self%relErr(2*iter)
                    if( rnorm .LT. btol ) then
                        !
                        self%n_iter = iter
                        self%converged = .TRUE.
                        exit
                        !
                    endif
                    !
                    if( rnorm .LT. rnormin ) then
                        !
                        rnormin = rnorm
                        xmin = x
                        imin = iter
                        !
                    endif
                    !
                enddo
                !
                if( .NOT. self%converged ) then 
                    !
                    !> it should be noted that this is the way my matlab version works
                    !> the bicg will return the "best" (smallest residual) iteration
                    !> x = xmin;  ! comment this line 
                    self%n_iter = self%max_iters
                    !> self%relErr(self%max_iters) = self%relErr(imin)  ! and this line
                    !> to use the last iteration result instead of the "best"
                endif
                !
                9 continue
                !
                if( lmem(1) ) deallocate( R, self%VomegaMuSig2 )
                !
                if( lmem(2) ) then
                    !
                    deallocate( xhalf, xmin )
                    deallocate( RT, P, S, SH, V, T )
                    !
                endif
                !
            class default
                call errStop( "BICG_OMP_FWD: Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine BICG_OMP_FWD ! BICG
    !
    !> Stabilized version of BiConjugate Gradient, set up for solving
    !> A x = b using routines in  mult_Aii.
    !> solves for the interior (edge) field
    !
    !> modified from my matlab version of BICGstab...
    !> so the naming might sound a little different from conventional ones
    !> also added the optional adjoint to solve adjoint system A^Tx = b
    !> 
    !> NOTE: BICG actually performs two sub (or half) line searches within a 
    !>       iteration, but here we only store the relerr for the second sub
    !>       just to be compatible with QMR
    !
    !> 
    !> interface...........
    !> redefining some of the interfaces for our convenience (locally)
    !> generic routines for vector operations for edge/ face nodes
    !> in a staggered grid
    !
    !> NOTE: this has not been extensively tested - I believe it feels a 
    !> little unstable (despite the name)...
    !> if you have time reading this, test it!
    !
    !>  b is right hand side
    !
    !>  solution vector is x ... on input is provided with the initial
    !>  guess, on output is the iterate with smallest residual.
    !
    subroutine BICG_OMP_ADJ( self, b, x )
        implicit none
        !
        class( Solver_BICG_OMP_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: b
        class( Vector_t ), intent( inout ) :: x
        !
        type( spMatCSR_real ) :: Mii
        real( kind=prec ) :: rnorm, bnorm, rnormin, btol
        class( Vector_t ), allocatable :: R, RT, V, T, P, S, SH
        class( Vector_t ), allocatable :: xhalf, xmin
        complex( kind=prec ) :: RHO, ALPHA, BETA, OMEGA, AX1
        complex( kind=prec ) :: RTV, TT, RHO1
        integer :: iter, xsize, imin, maxiter
        logical :: lmem(2)
        !
        !CMC: For new OMP version
        integer :: i, j, j2
        complex( kind=prec ) :: cnorm
        !
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, R )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, RT )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, V )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, T )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, P )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, S )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, SH )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, xhalf )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, xmin )
        !
        lmem(2) = .FALSE.
        !
        xsize = x%length()
        allocate( self%VomegaMuSig2(xsize) ); lmem(1) = .TRUE.
        cnorm = C_ZERO; TT = C_ZERO
        !
        !$OMP PARALLEL PRIVATE(self%Ithread)
        self%Nthreads = OMP_get_num_threads(); self%Ithread = OMP_get_thread_num()
        !
        if( self%Ithread == 0 ) write( *, * ) "OpenMP threads num,max=", self%Nthreads, OMP_get_max_threads()
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%preconditioner%model_operator )
            !
            class is( ModelOperator_SP_t )
                !
                !> Square Matrix Mii: differentiated Differentiation between SP1 (CCii) and SP2 (AAii)
                select type( model_operator )
                    !
                    class is( ModelOperator_SP_V1_t )
                        !
                        Mii = model_operator%CCii
                        !
                    class is( ModelOperator_SP_V2_t )
                        !
                        Mii = model_operator%AAii
                        !
                    class default
                        call errStop( "BICG_OMP_ADJ > model_operator must be SP V1 or V2" )
                    !
                end select
                !
                !$OMP DO REDUCTION(+:cnorm,TT) PRIVATE(j)
                do i = 1, Mii%nRow
                    !
                    self%VomegaMuSig2(i) = ONE_I * ISIGN * model_operator%VomegaMuSig(i)
                    !
                    R = self%VomegaMuSig2(i) * x
                    !
                    do j = Mii%row(i), Mii%row(i+1) - 1 
                        R = R + Mii%val(j) * x !x( Mii%col(j) ) ????
                    enddo
                    !
                    R = b - R
                    cnorm = cnorm + R * R%conjugate()
                    TT = TT + b * b%conjugate() ! bnorm = SQRT(dot_product(b, b))
                    !
                enddo
                !$OMP enddo
                !$OMP END PARALLEL
                !
                !> Norm of residual
                rnorm = DREAL( CDSQRT( cnorm ) )
                bnorm = SQRT( DREAL( TT ) )
                !
                write( *, * ) "Initial bnorm=",bnorm
                !
                if( bnorm < 1d-16 ) then ! zero rhs -> zero solution
                    !
                    call warning( "b in BICG_OMP_ADJ has all zeros, returning zero solution" )
                    !
                    x = b 
                    self%n_iter = 1
                    self%relErr = 0.0
                    goto 9
                    !
                endif
                !
                btol = self%tolerance * bnorm
                !
                if( rnorm .LE. btol ) then ! the first guess is already good enough
                    !
                    !> returning
                    self%n_iter = 1
                    self%relErr(1) = real( rnorm / bnorm, kind=prec )
                    goto 9
                    !
                endif
                !
                !> allocate the rest here
                !allocate(xhalf(xsize),xmin(xsize))
                !allocate(RT(xsize))
                !allocate(P(xsize))
                !allocate(S(xsize))
                !allocate(SH(xsize))
                !allocate(V(xsize))
                !allocate(T(xsize))
                !
                lmem(2) = .TRUE.
                !================= Now start configuring the iteration ===================!
                !> the adjoint (shadow) residual
                rnormin = rnorm
                self%relErr(1) = real( rnormin / bnorm, kind=prec )
                !> write(6,*) "initial residual",  self%relErr(1)
                self%converged = .FALSE.
                maxiter = self%max_iters 
                imin = 0
                RHO = C_ONE; OMEGA = C_ONE
                RT = R
                xmin = x
                imin = 1
                !> subroutine Mult_Aii: File FWD_SP2/modelOperator3D.f90, L. 449
                !>   call RMATxCVEC(Mii,x,y)
                !> subroutine RMATxCVEC: File SP_Topology/spOpTools.f90, L. 404
                !
                !============================== looooops! ================================!
                !write( *, * )Lmat%nRow,Umat%nRow,Mii%nRow,"size" ! all equal
                !
                do iter = 1, maxiter
                    !
                    RHO1 = RHO; RHO = C_ZERO
                    !
                    !$OMP PARALLEL
                    !$OMP DO REDUCTION(+:RHO)
                    do i=1,Mii%nRow
                        RHO = RHO + RT%conjugate() * R
                    enddo
                    !$OMP enddo
                    !
                    if( iter == 1 ) then
                        BETA = C_ZERO
                    else
                        BETA = (RHO/RHO1) * (ALPHA/OMEGA) 
                    endif
                    !
                    !$OMP DO
                    do i = 1, Mii%nRow
                        !
                        P = R + BETA * (P - OMEGA * V)
                        V = P
                        !
                    enddo
                    !
                    !$OMP enddo
                    !$OMP END PARALLEL
                    !if(cdabs(RHO) < 1d-16) then
                    !>    BICGiter%failed = .TRUE.
                    !>    exit
                    !endif 
                    !call M1solve(P,ilu_adjt,PT) ! P->x  PT->y
                    !call M2solve(PT,ilu_adjt,PH)
                    !
                    !call solve_Solver_BICG_OMP( .FALSE., V, SH )
                    call self%solve( .FALSE., V, SH )
                    !
                    RTV = C_ZERO; cnorm = C_ZERO
                    !
                    !$OMP PARALLEL PRIVATE(self%Ithread)
                    self%Ithread = OMP_get_thread_num()
                    !
                    !$OMP DO REDUCTION(+:RTV) PRIVATE(j)
                    do i = 1, Mii%nRow
                        !
                        V = self%VomegaMuSig2(i) * SH
                        !
                        do j = Mii%row(i), Mii%row(i+1) - 1 
                            V = V + Mii%val(j) * SH !SH(Mii%col(j)) ????
                        enddo
                        !
                        RTV = RTV + RT%conjugate() * V
                        !
                    enddo
                    !
                    !$OMP enddo
                    if( self%Ithread == 0 ) ALPHA = RHO / RTV
                    !
                    !$OMP BARRIER
                    !$OMP DO
                    do i = 1, Mii%nRow
                        xhalf = x + ALPHA * SH ! the first half of iteration
                    enddo
                    !$OMP enddo
                    !
                    !$OMP DO REDUCTION(+:cnorm) PRIVATE(j,AX1)
                    do i = 1, Mii%nRow
                        !
                        AX1 = self%VomegaMuSig2(i) * xhalf
                        !
                        do j = Mii%row(i), Mii%row(i+1) - 1 
                            AX1 = AX1+Mii%val(j) * xhalf !xhalf( Mii%col(j) ) ????
                        enddo
                        !
                        AX1 = b - AX1
                        cnorm = cnorm + AX1 * AX1%conjugate()
                        S = R - ALPHA * V !residual for the 0.5 x
                        !
                    enddo
                    !$OMP enddo
                    !$OMP END PARALLEL
                    !
                    if( RTV .EQ. 0.0 .OR. ALPHA .EQ. 0.0 ) then
                        !
                        exit
                        !
                    endif
                    !
                    cnorm = CDSQRT( cnorm ); rnorm = DREAL( cnorm )
                    !
                    self%relErr(iter) = rnorm / bnorm
                    !write( *, * ) "iter # ",iter," xhalf residual: ", self%relErr(iter)
                    !
                    if( rnorm .LT. btol ) then
                        !
                        x = xhalf
                        self%n_iter = iter
                        self%converged = .TRUE.
                        exit
                        !
                    endif
                    !
                    if( rnorm .LT. rnormin ) then
                        !
                        rnormin = rnorm
                        xmin = xhalf
                        imin = iter
                        !
                    endif
                    !
                    !> ========================== M1,M2solve ==========================
                    !> call M1solve(S,ilu_adjt,ST) = call PC_Lsolve(x,adjt,y) 
                    !>                                 call LTsolve_Cmplx(L,x,y)
                    !>                           subroutine LTsolve_Cmplx(L,b,x)
                    !> map variables:                    S->x->b, ST->y->x
                    !> matrix = Lmat => L (lower triang.)
                    !call LTsolve_Cmplx(L,x,y) S
                    !!$OMP PARALLEL PRIVATE(self%Ithread)
                    !>         self%Ithread =OMP_get_thread_num()
                    !>         if(self%Ithread==0)write(*,*)Lmat%nRow,"Lmat%nRow"
                    !!$OMP DO SCHEDULE(STATIC)
                    !>         do i = 1,Lmat%nRow
                    !>            if(iter==1)then
                    !>               if(self%Ithread==0)write(*,*)i,"THREAD0"
                    !>               if(self%Ithread==1)write(*,*)i,"THREAD1"
                    !>            endif
                    !>         enddo
                    !call solve_Solver_BICG_OMP( .FALSE., S, SH )
                    call self%solve( .FALSE., S, SH )
                    !
                    TT = C_ZERO; OMEGA = C_ZERO; cnorm = C_ZERO
                    !
                    !$OMP PARALLEL PRIVATE(self%Ithread)
                    self%Ithread = OMP_get_thread_num()
                    !
                    !$OMP DO REDUCTION(+:TT,OMEGA) PRIVATE(j)
                    do i = 1, Mii%nRow
                        !
                        T = self%VomegaMuSig2(i) * SH
                        !
                        do j = Mii%row(i), Mii%row(i+1) - 1 
                            T = T + Mii%val(j) * SH !SH(Mii%col(j)) ?????
                        enddo
                        !
                        TT = TT + T * T%conjugate() ! complex dot-product
                        OMEGA = OMEGA + T%conjugate() * S ! complex dot-product
                        !
                    enddo
                    !$OMP enddo
                    !
                    if( self%Ithread == 0 ) OMEGA = OMEGA / TT
                    !
                    !if(iter==1)print*,"OMEGA",self%Ithread,OMEGA
                    !$OMP BARRIER
                    !
                    !$OMP DO
                    do i = 1, Mii%nRow
                        x = xhalf + OMEGA * SH  ! the second half of iteration
                    enddo
                    !$OMP enddo
                    !
                    !$OMP DO REDUCTION(+:cnorm) PRIVATE(j,AX1)
                    do i = 1, Mii%nRow
                        !
                        AX1 = self%VomegaMuSig2(i) * x
                        !
                        do j = Mii%row(i), Mii%row(i+1) - 1 
                            AX1 = AX1 + Mii%val(j) * x !x(Mii%col(j)) ????
                        enddo
                        !
                        AX1 = b - AX1
                        cnorm = cnorm + AX1 * AX1%conjugate()
                        R = S - OMEGA * T  !residual for the 1.0 x
                        !
                    enddo
                    !$OMP enddo
                    !
                    !$OMP END PARALLEL
                    if( TT .EQ. 0.0 .OR. OMEGA .EQ. 0.0 ) then
                        !
                        exit
                        !
                    endif
                    !
                    rnorm = CDSQRT( cnorm )
                    self%relErr(iter) = real( rnorm / bnorm, kind=prec )
                    !> write(6,*) "iter # ",iter," x residual: ", self%relErr(2*iter)
                    !
                    if( rnorm .LT. btol ) then
                        !
                        self%n_iter = iter
                        self%converged = .TRUE.
                        exit
                        !
                    endif
                    !
                    if( rnorm .LT. rnormin ) then
                        !
                        rnormin = rnorm
                        xmin = x
                        imin = iter
                        !
                    endif
                    !
                enddo
                !
                if( .NOT. self%converged ) then 
                    !> it should be noted that this is the way my matlab version works
                    !> the bicg will return the "best" (smallest residual) iteration
                    !> x = xmin;  ! comment this line 
                    self%n_iter = self%max_iters
                    !> self%relErr(self%max_iters) = self%relErr(imin)  ! and this line
                    !> to use the last iteration result instead of the "best"
                endif
                !
                9 continue
                !
                if( lmem(1) ) deallocate( R, self%VomegaMuSig2 )
                !
                if( lmem(2) ) then
                    !
                    deallocate(xhalf,xmin)
                    deallocate(RT)
                    deallocate(P)
                    deallocate(S)
                    deallocate(SH)
                    deallocate(V)
                    deallocate(T)
                    !
                endif
                !
            class default
                call errStop( "BICG_OMP_ADJ: Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine BICG_OMP_ADJ ! BICG
    !
end module Solver_BICG_OMP
!
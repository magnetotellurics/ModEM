module MatUtils
    !
    use Constants
    !
    interface CumSum
         module procedure CumSum_real
         module procedure CumSum_complex
    end interface CumSum
    !
contains
    !
    !> Given the rank 1 array 'a', returns an array
    !> of identical type and size containing the
    !> cumulative sums of 'a'.
    function CumSum_real(a) result(r)
        implicit none
        !
        real( kind=prec ), intent( in ) :: a(:)
        !
        real( kind=prec ), allocatable :: r(:)
        !
        integer :: n, j
        !
        n = size(a)
        allocate(r(n))
        !
        r(1) = a(1)
        !
        do j = 2, n
             r(j) = r(j-1) + a(j)
        enddo
        !
    end function CumSum_Real
    !
    !> Given the rank 1 array 'a', returns an array
    !> of identical type and size containing the
    !> cumulative sums of 'a'.
    function CumSum_complex(a) result(r)
        implicit none
        !
        complex( kind=prec ), intent( in ) :: a(:)
        !
        complex( kind=prec ), allocatable :: r(:)
        integer :: n, j
        !
        n = size(a)
        allocate(r(n))
        !
        r(1) = a(1)
        !
        do j = 2, n
             r(j) = r(j-1) + a(j)
        enddo
        !
    end function CumSum_Complex
    !
end module MatUtils

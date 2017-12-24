module Poly
    use Const
!     use IO
    implicit none

contains

    function poly_value(poly, x) result(r)
        real(mpc), intent(in) :: poly(:), x
        real(mpc) :: r
        integer :: i

        r = 0

        do i = 1, size(poly)
            r = r * x + poly(i)
        end do
    end function poly_value

    !-----------------------------------------------------------------------
    ! gorner (poly, x_0) -> r_poly
    ! poly: divident 
    ! x_0 : (x-x_0) is a divisor
    ! 
    ! 	implements gorner method of poly division
    !-----------------------------------------------------------------------
    function gorner_div(poly, x_0) result(r_poly)
        real(mpc), dimension(:), intent(in):: poly
        real(mpc), intent(in) :: x_0
        real(mpc), dimension((size(poly))) :: r_poly
        integer :: n, i
!🐾
        n = size(poly)
        r_poly = 0
        r_poly(1) = poly(1)
        do i = 2, n-1
            r_poly(i) = r_poly(i-1)*x_0 + poly(i)
        end do
    end function gorner_div

    
    !-----------------------------------------------------------------------
    ! bern_solve (a) -> maxroot
    ! a: коэффициенты многочлена
    ! 
    ! 	find root of polynomial using Bernoully method
    !-----------------------------------------------------------------------
    function bern_solve(a, max_iter) result(maxroot)
        real(mpc), dimension(0:), intent(in) :: a
        integer, value, optional :: max_iter
        real(mpc) :: maxroot

        real(mpc), dimension(0:size(a)-1) :: y
        real(mpc) :: x, x_1
        integer   :: i

        if (.not. present(max_iter)) then
            max_iter = -1
        else
            max_iter = abs(max_iter) ! ну а вдруг..
        end if

        y = ini_y()
        x = huge(1.0_mpc)
        x_1 = y(0)/y(1)
        i = 0

        do while (abs(x-x_1) > eps .and. i /= max_iter )
            x = x_1
            y = shift(y)
            x_1 = y(0)/y(1)
            i = i + 1
        end do
 
        maxroot = x

        contains
            !-----------------------------------------------------------------------
            ! ini_y (.a) -> arr
            ! .a -- from outer scope
            ! 
            ! 	generates inital values for Bernoully solver.
            !-----------------------------------------------------------------------
            function ini_y() result(y)
                real(mpc), dimension(0:size(a)-1) :: y
                integer :: n
                n = size(a) -1 
                y = 1
                y(0) = -dot_product(y(1:n-1), a(2:n))/a(1)  + 1.0
            end function ini_y

            function shift(y) result(y1)
                real(mpc), dimension(0:), intent(in) :: y
                real(mpc), dimension(0:size(y)-1) :: y1
                integer :: n
                n = size(y) - 1

                ! будем орудовать с перевёрнутым массивом y_i, так можно считать
                ! скалярное произведение
                y1(1:n) = y(0:n-1)
                y1(0) = -1.0/a(0) * dot_product(a(1:), y1(1:))
            end function shift

    end function bern_solve

    function poly_roots(poly) result(roots)
        real(mpc), dimension(:), intent(in) :: poly
        real(mpc) :: roots(size(poly)-1), cpoly(size(poly))
        integer :: n,i
        n = size(poly) - 1
        cpoly = poly

        do i = 1, n
            roots(i) = bern_solve(cpoly, 100)
            cpoly = gorner_div(cpoly, roots(i))
        end do
    end function poly_roots
    
    pure function poly_by_roots(a0, roots) result(poly)
        real(mpc) :: a0, roots(:), poly(0:size(roots))
        intent(in) :: a0, roots
        integer :: i, n

        n = size(roots)

        poly    = 0
        poly(0) = 1
        ! здесь потом будет нормальное умножение многочленов
        do i = 1, n
            poly(1:n) = poly(1:n) - roots(i)*poly(0:n-1) 
        end do

        poly = poly * a0
    end function poly_by_roots

end module Poly

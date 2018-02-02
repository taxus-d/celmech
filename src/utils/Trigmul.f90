module Trigmul
    use Def_prec
    implicit none

contains
    ! тут массивы с единицы, надо подумать прежде чем передать n из остальной
    ! части проги
    function trdmatmul(a, b, n) result(c)
            real(mpc), intent(in), dimension(n, 3) :: a, b
!             real(mpc), intent(in), dimension(n), optional :: weigths
            integer, intent(in) :: n
            real(mpc), dimension (n, 5) :: c
            integer :: i
            
            c(3:n  , 1) = a(3:n, 1)*b(2:n-1, 1)
            c(2:n  , 2) = a(2:n,1)*b(1:n-1, 2) + a(2:n, 2)*b(2:n,1)
            
            c(1:n  , 3) = a(:, 2)*b(:, 2)
            c(2:n  , 3) = c(2:n, 3) + a(2:n, 1)*b(1:n-1, 3)
            c(1:n-1, 3) = c(1:n-1, 3) + a(1:n-1, 3)*b(2:n, 1)
            
            c(1:n-1, 4) = a(1:n-1, 2)*b(1:n-1, 3) + a(1:n-1, 3)*b(2:n, 2)
            c(1:n-2, 5) = a(1:n-2, 3)*b(2:n-1, 3)

    end function trdmatmul

    function transptrd(a, n) result(aT)
        integer :: n
        real(mpc), intent(in), dimension(n, 3) :: a
        real(mpc), dimension(n, 3) :: aT

        aT = a
        aT(1:n-1, 3) = a(2:n, 1)
        aT(2:n, 1) = a(1:n-1, 1)

    end function transptrd

    ! умножение трёхдиагональной матрицы на вектор (линейный оператор)
    function trdlinop(a, x, n) result(y)
        integer, intent(in) :: n
        real(mpc), intent(in) :: a(n,3), x(n)
        real(mpc), dimension (n) :: y
        integer :: i
        y(1) = sum(a(1, 2:3)*x(1:2)) 
        y(n) = sum(a(n, 1:2)*x(n-1:n)) 
        forall (i = 2:n-1) y(i) = dot_product(a(i, :), x(i-1:i+1))
    end function trdlinop

end module Trigmul

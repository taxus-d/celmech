module Deriv
    use Inival
    use Utils
    use FuncIfaces
    use IO_array
    implicit none
contains
!----------------------------------------------------------
! pder( f, i, x0) -> val
! f   : function \R^n \to \R^n
! x_0 : initial point
!       computes i-th partial derivative at given point
!----------------------------------------------------------
    function pder(f, i, x_0) result(res)
        integer, intent(in) :: i 
        real(mpc), intent(in), dimension(:) :: x_0
        real(mpc), dimension(size(x_0)) :: res
        procedure (fRnRn) :: f
        integer :: n
        real(mpc), dimension(size(x_0)) :: x_1
        real(mpc) :: delta, seps

        n = size(x_0)
        if ( i > n) stop "pder: Illegal variable number"
        x_1 = x_0

        ! тут можно сказать пару слов, почему \sqrt(\eps)
        ! С одной стороны, есть погрешность ~ \delta^2 из-за разложения
        ! С другой стороны, f(x) представляется с погрешностью f(x) * \eps, так что
        !+ при делении разности на \delta погрешность будет \delta^{-1} f(x) \eps. 
        ! Если пытаться уменьшить и то и то, как раз нужен \sqrt{\eps}
        seps = sqrt(eps)
        ! тут какая-то магия с точностью. 
        ! Оно иногда ломается, но я не знаю как лучше (
        delta = seps
        x_1(i) = x_1(i) + delta 

        res = ( f(x_1) - f(x_0) ) / (delta)
        
    end function pder

    
    !-----------------------------------------------------------------------
    ! jacobimtx (f, x_0) -> J
    ! f: \R^n -> \R^n 
    ! x_0: \in \R^n
    ! 
    ! 	computes jacobi matrix of f in x_0
    !-----------------------------------------------------------------------
    function jacobimtx(f, x_0) result(J)
        real(mpc), intent(in), dimension(:) :: x_0
        procedure (fRnRn) :: f
        real(mpc), dimension(size(x_0), size(x_0)) :: J
        integer :: n, k

        J=0

        n = size(x_0)
        do k = 1, n
            J(k, :) = pder(f, k, x_0)
        end do
        J = transpose(J)

    end function jacobimtx
    
    


end module Deriv

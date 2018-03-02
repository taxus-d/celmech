module Deriv
    use Const
    use Utils
    use FuncIfaces
    use IO_array
    implicit none
    interface pder
        module procedure :: pder_vec
        module procedure :: pder_scalar
    end interface pder
contains
!----------------------------------------------------------
! pder( f, i, x0) -> val
! f   : function \R^n \to \R^n
! x_0 : initial point
!       computes i-th partial derivative at given point
!----------------------------------------------------------
    function pder_vec(f, i, x_0) result(res)
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
        if (x_1(i) < 1.0_mpc) then 
            delta = seps
        else
            delta = seps * x_1(i)
        end if
        x_1(i) = x_1(i) + delta
        delta = x_1(i) - x_0(i)! ещё магия

        res = ( f(x_1) - f(x_0) ) / (delta)
    end function pder_vec
    function pder_scalar(f, i, x_0) result(res)
        integer, intent(in) :: i 
        real(mpc), intent(in), dimension(:) :: x_0
        real(mpc) :: res
        procedure (fRnR1) :: f
        integer :: n
        real(mpc), dimension(size(x_0)) :: x_1
        real(mpc) :: delta, seps

        n = size(x_0)
        if ( i > n) stop "pder: Illegal variable number"
        x_1 = x_0

        seps = sqrt(eps)
        if (x_1(i) < 1.0_mpc) then
            delta = seps
        else
            delta = seps * x_1(i)
        end if
        x_1(i) = x_1(i) + delta
        delta = x_1(i) - x_0(i)! ещё магия
        
        res = ( f(x_1) - f(x_0) ) / (delta)
    end function pder_scalar

    
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
    
    function ngrad(f,x0) result(grad)
        real(mpc), intent(in), dimension(:) :: x0
        procedure (fRnR1) :: f
        real(mpc), dimension(size(x0)) :: grad
        real(mpc) :: f_0, shiftsize
        integer :: n, k

        n = size(x0)
        shiftsize = opt_shiftsize(x0)
        f_0 = f(x0)
        do k = 1, n
            grad(k) = (f(x0 + shiftsize*unit_v(k)) - f_0) / shiftsize
        end do
        contains
            pure function unit_v(i) result(e)
                integer, intent(in) :: i
                real(mpc), dimension(size(x0)) :: e
                e = 0
                e(i) = 1.0_mpc
            end function
    end function ngrad
    
    
    function opt_shiftsize(x) result(s)
        real(mpc), intent(in), dimension(:) :: x
        real(mpc) :: s
        if (norm2(x) < 1.0_mpc) then
            s = sqrt(eps)
        else
            s = norm2(x)*sqrt(eps)
        end if
    end function opt_shiftsize
    
    
    function hess_mtx(f, x0) result(hessm)
        procedure (fRnR1) :: f
        real(mpc), intent(in), dimension(:) :: x0 
        real(mpc), dimension(size(x0), size(x0)) :: hessm
        real(mpc) :: f_0, f_i, shiftsize
        integer :: n,i,j

        n = size(x0)
        
        shiftsize = opt_shiftsize(x0) 
        f_0 = f(x0)

        do i = 1, n
            f_i = f(x0 + shiftsize*unit_v(i))
            do j = i, n 
                hessm(i,j) = ( f(x0 + shiftsize*(unit_v(i) + unit_v(j))) - f_i - f(x0 + shiftsize*unit_v(j)) + f_0) / (shiftsize**2)
                hessm(j,i) = hessm(i,j)
            end do
        end do
        contains
            pure function unit_v(i) result(e)
                integer, intent(in) :: i
                real(mpc), dimension(size(x0)) :: e
                e = 0
                e(i) = 1.0_mpc
            end function
    end function hess_mtx
    
    function dir_deriv(f,x0,v) result(der)
        procedure (fRnRn) :: f
        real(mpc), intent(in), dimension(:) :: x0 
        real(mpc), intent(in), dimension(size(x0)) :: v
        real(mpc), dimension(size(x0)) :: der
        real(mpc) :: f_0(size(x0)), shiftsize
        integer :: n,i,j
        
        shiftsize = opt_shiftsize(x0)
        f_0 = f(x0)
        der = (f(x0 + shiftsize*v) - f_0) / shiftsize
    end function dir_deriv
end module Deriv

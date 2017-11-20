module Samples
    use Const
    implicit none
contains 
    function linear2(t, X) result(dX)
        real(mpc) :: t, X(2), dX(size(X))
        intent(in) :: t, X
        dX(1) =  -X(2)
        dX(2) =  X(1)
    end function linear2
    function sol_linear2(t) result(X)
        real(mpc), intent(in) :: t
        real(mpc), dimension(2) :: X
        X(1) = cos(t) - sin(t)
        X(2) = sin(t) + cos(t)
    end function sol_linear2
    function linear2_nhom(t, X) result(dX)
        real(mpc) :: t, X(2), dX(size(X))
        intent(in) :: t, X
        dX(1) =  3.0_mpc*X(1) + 2.0_mpc*X(2) + sin(t)
        dX(2) = -5.0_mpc*X(1) + 1.0_mpc*X(2) + exp(-t)
    end function linear2_nhom
    function sol_linear2_nhom(t) result(X)
        real(mpc), intent(in) :: t
        real(mpc), dimension(2) :: X
        X(1) = exp(2*t)*(53.0*sin(3*t)/45.0+151.0*cos(3*t)/180.0)-sin(t)/10.0+co&
            &s(t)/20.0+exp(-t)/9.0
        X(2) = exp(2*t)*((-133.0)*sin(3*t)/72.0+97.0*cos(3*t)/72.0)+(-3.0)*sin(t&
            &)/8.0-cos(t)/8.0+(-2.0)*exp(-t)/9.0
    end function sol_linear2_nhom
end module Samples

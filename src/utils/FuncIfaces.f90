module FuncIfaces
    use Const
    implicit none
    interface

        ! Sample functions for solvers
        function fRR(x) result(y)
            import mpc
            real(mpc),intent(in) :: x
            real(mpc) :: y
        end function 

        function fRnRn(x) result(y)
            import mpc
            real(mpc), intent(in), dimension(:) :: x
            real(mpc), dimension(size(x)) ::y
        end function 
        
        function fRnR1(x) result(y)
            import mpc
            real(mpc), intent(in), dimension(:) :: x
            real(mpc) ::y
        end function 

        function genericTransform(x) result (x1)
            import mpc
            real(mpc), intent(in), dimension(:) :: x
            real(mpc), dimension(size(x)) :: x1
        end function genericTransform 
    end interface 
end module FuncIfaces

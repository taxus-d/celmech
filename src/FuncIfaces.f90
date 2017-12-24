module FuncIfaces
    use Const
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
    end interface 
end module FuncIfaces

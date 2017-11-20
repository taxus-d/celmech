module FuncIfaces
    interface

        ! Sample functions for solvers
        function fRR(x) result(y)
            use Inival
            real(mpc),intent(in) :: x
            real(mpc) :: y
        end function 

        function fRnRn(x) result(y)
            use Inival
            real(mpc), intent(in), dimension(:) :: x
            real(mpc), dimension(size(x)) ::y
        end function 
        
    end interface 
end module FuncIfaces

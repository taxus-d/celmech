module Tests
! ===<>===<>===<>===<>===<>===<>=
! не самые абстрактные тесты для многочленов. 
! Пока не выглядят нужными
! ===<>=========================
    use Prec
    use Utils
    use FuncIfaces
    use Poly
    use LegendrePoly
    use IO_gaussint
    implicit none
contains 
!     function abstract_real_test(testing, [args ...] )
! а, нет, фортран не умеет (кажется) произвольное число аргументов
! а жаль

    function test_poly_root_finder(poly, solver, max_iter) result(passedp)
        real(mpc), intent(in), dimension(:) :: poly
        procedure(poly_root_finder) :: solver
        integer, value, optional :: max_iter
        real(mpc) :: residue
        logical :: passedp
        
        residue = poly_value(poly, solver(poly, max_iter))
        passedp = residue < 8*eps
        if (.not. passedp) write(stderr, *) residue, ':',  poly
    end function test_poly_root_finder

    function test_poly_solver(poly, solver) result (passedp)
        real(mpc), intent(in), dimension(:) :: poly
        procedure(poly_solver) :: solver
        real(mpc) :: residue, roots(size(poly)-1)
        logical :: passedp
        integer :: i, n
        
        n = size(poly) - 1
        residue = 0
        roots = solver(poly)
        do i = 1, n
            residue = max(residue, poly_value(poly, roots(i)))
        end do

        passedp = residue < 8*eps
        if (.not. passedp) write(stderr, *) residue, ':',  poly
    end function test_poly_solver
    
    function test_int_coefs() result(passedp)
        real(mpc) :: residue
        logical :: passedp
        integer :: i, n
        
        residue = norm2((/ 5.0_mpc/9.0_mpc, 8.0_mpc/9.0_mpc, 5.0_mpc/9.0_mpc/) - &
            int_coefs((/ -sqrt(3._mpc)/sqrt(5._mpc), 0.0_mpc, sqrt(3._mpc)/sqrt(5._mpc)/)))
        passedp = residue < 8._mpc*eps
        if (.not. passedp) write(stderr, *) residue, ':'
    end function test_int_coefs
end module Tests

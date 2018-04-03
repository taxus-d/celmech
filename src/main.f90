program main
    use Inival
    use ODEsolve
    use Poincare
    use Utils
    use IO_array
    use ShapeTransform
    implicit none

    integer   :: fd_orbi = 10, fd_orbo = 11, fd_orbd = 12, fd_x0 = 13
    real(mpc), allocatable, dimension(:) :: x0_mod
    type(RungeKuttaInt) :: itg_rk

    procedure (eq_fun), pointer :: f => poincare_section_eq
    
    open(fd_orbo, file="data/orbit-orig.dat", action="write")
    open(fd_orbi, file="data/orbit-impr.dat", action="write")
    open(fd_orbd, file="data/orbit-ideal.dat", action="write")
    open(fd_x0  , file="data/x0.dat", action="read")
    allocate(x0_mod(size(x0)))

    itg_rk = RungeKuttaInt(f,x0,t0,h)
    current%S  => sect_x_1_body
    current%dS => sect_x_1_body_deriv
    call print_solution(itg_rk, t1, fd_orbo)
    
    ! reuse data from previous runs
    x0_mod = x0
!     read(fd_x0, *) x0_mod(1:12)
    x0_mod(13) = t0
    close(fd_x0) 
    call prarr (x0_mod)
    
    call imporove_inipos(itg_rk, t0, x0_mod, t1, stdout)
    call prarr(x0_mod)
!
    call itg_rk%set_inicond(x0_mod, t0)
    call itg_rk%crewind()
    call print_solution(itg_rk, 100*Period, fd_orbi, transf=shapecoords)
    open(fd_x0, file="data/x0.dat", action="write")
    write(fd_x0, *) x0_mod
    
    call itg_rk%set_inicond(x0_ideal, t0)
    call itg_rk%crewind()
    call print_solution(itg_rk, 10*Period, fd_orbd, transf=shapecoords)
    
    open(fd_x0, file="data/x0.dat", action="write")
    close(fd_orbo); close(fd_orbi); close(fd_orbd)
    close(fd_x0)
    deallocate(x0_mod)
end program main

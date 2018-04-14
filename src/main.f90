program main
    use Inival
    use ODEsolve
    use Poincare
    use Utils
    use IO_array
    use ShapeTransform
    implicit none

    integer   :: fd_orbi = 10, fd_orbo = 11, fd_orbd = 12, fd_x0 = 13
    real(mpc), allocatable, dimension(:) :: x0_ideal, x0
    real(mpc) :: t0, t1
    type(RungeKuttaInt) :: itg_rk

    procedure (eq_fun), pointer :: f, fp !=> shape_motion_eq
    procedure (genericTransform), pointer :: ctranform
    
    open(fd_orbo, file="data/orbit-orig.dat", action="write")
    open(fd_orbi, file="data/orbit-impr.dat", action="write")
    open(fd_orbd, file="data/orbit-ideal.dat", action="write")
    open(fd_x0  , file="data/x0.dat", action="read")
    
   
    f => shape_motion_eq
    call assign_inicond(x0_ideal, x0, t0, t1)

    itg_rk = RungeKuttaInt(poincare_section_eq,x0,t0,h)
    current%S  => sect_x_1_body
    current%dS => sect_x_1_body_deriv

    call print_solution(itg_rk, t1, fd_orbo)
    
    ! reuse data from previous runs
!     read(fd_x0, *) x0_mod(1:12)
    close(fd_x0) 
    call prarr (x0)

    call imporove_inipos(itg_rk, t0, x0, t1, stdout)
    call prarr(joinarr_sbs(x0,x0_ideal))
!
    call itg_rk%set_inicond(x0, t0)
    call itg_rk%crewind()
    call print_solution(itg_rk, 200*Period, fd_orbi)

    open(fd_x0, file="data/x0.dat", action="write")
    write(fd_x0, *) x0
    
    call itg_rk%set_inicond(x0_ideal, t0)
    call itg_rk%crewind()
    call print_solution(itg_rk, 10*Period, fd_orbd)
    
    open(fd_x0, file="data/x0.dat", action="write")
    
    close(fd_orbo); close(fd_orbi); close(fd_orbd)
    close(fd_x0)

    call cleanup_inicond(x0_ideal, x0)

contains
    function poincare_section_eq(tt,X) result(dX)
        real(mpc), intent(in) :: tt
        real(mpc), dimension(:), intent(in) :: X
        real(mpc), dimension(size(X)) :: dX
        real(mpc) :: H
        
        dX(1:finalDim) = f(tt, X(1:finalDim))
        dX(finalDim+1) = 1
        H = 1
        if (weirdstep) H=dot_product(current%dS(X(1:finalDim)),dX(1:finalDim))
        dX = dX/H
    end function
    
end program main

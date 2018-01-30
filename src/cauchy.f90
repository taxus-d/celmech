
program cauchy
    use Inival
    use ODEsolve
    use Poincare
    use Utils
    use IO_array
    implicit none

    integer   :: fd_rk = 10, fd_ae = 11, fd_ai = 12, fd_x0 = 13, fd_rk1 = 14,&
        fd_track = 15
    real(mpc) :: r_rk, r_ae, r_ai
    real(mpc), allocatable, dimension(:) :: x0_mod
    type(RungeKuttaInt) :: itg_rk
    type(ExAdamsInt) :: itg_ae
    type(ImAdamsInt) :: itg_ai

    procedure (eq_fun), pointer :: f => poincare_section_eq
    
    open(fd_rk, file="data/rk.dat", action="write")
    open(fd_rk1, file="data/rk1.dat", action="write")
    open(fd_track, file="data/track.dat", action="write")
    open(fd_ae, file="data/ae.dat", action="write")
    open(fd_ai, file="data/ai.dat", action="write")
    open(fd_x0, file="test/x0.dat", action="read")
!     write(*,*) "``Variables``"
!     write(*,'(4(a16))') "h",  "t0","t1",  "ad_ord"
!     write(*,'(3(e16.8),i7)')  h,  t0,t1,  ad_ord
!     call prarr(X0)
    allocate(x0_mod(size(x0)))

    itg_rk = RungeKuttaInt(f,x0,t0,h)
    current%S => sect_x_1_body
    current%dS => sect_x_1_body_deriv
!     itg_ae = ExAdamsInt(f, x0, t0, ad_ord, h)
!     itg_ai = ImAdamsInt(f, x0, t0, ad_ord, h)
    call print_solution(itg_rk, t1, fd_rk)
!     call print_solution(itg_ae, t1, fd_ae)
!     call print_solution(itg_ai, t1, fd_ai)
!     call print_poincare_section(itg_rk, t1, fd_rk)
!     call print_poincare_section(itg_ae, t1, fd_ae)
!     call print_poincare_section(itg_ai, t1, fd_ai)
!     x0_mod = x0
    read(fd_x0, *) x0_mod(1:12)
    close(fd_x0)
    x0_mod(13) = t0
        x0_mod = x0
    call prarr (x0_mod)
    call imporove_inipos(itg_rk, t0, x0_mod, t1, stdout)
    call itg_rk%set_inicond(x0_mod, t0)
    call itg_rk%crewind()
    call print_solution(itg_rk, 100*Period, fd_rk1)
    open(fd_x0, file="test/x0.dat", action="write")
    write(fd_x0, *) x0_mod
    
    close(fd_rk)
    close(fd_rk1)
    close(fd_ae)
    close(fd_ai)
    close(fd_x0)
    close(fd_track)
    deallocate(x0_mod)
!     write(*,*) "Residues:"
!     write(*,*) "rk", sqrt(r_rk)
!     write(*,*) "ae", sqrt(r_ae)
!     write(*,*) "ai", sqrt(r_ai)
end program cauchy

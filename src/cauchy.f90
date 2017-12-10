program cauchy
    use ODEsolve
    use Poincare
    use Utils
    use IO_array
    implicit none

    integer   :: fd_rk = 10, fd_ae = 11, fd_ai = 12
    real(mpc) :: r_rk, r_ae, r_ai
    real(mpc), allocatable, dimension(:) :: res
    type(RungeKuttaInt) :: itg_rk 
    type(ExAdamsInt) :: itg_ae
    type(ImAdamsInt) :: itg_ai

    open(fd_rk, file="data/rk.dat", action="write")
    open(fd_ae, file="data/ae.dat", action="write")
    open(fd_ai, file="data/ai.dat", action="write")

    write(*,*) "``Variables``"
    write(*,'(4(a16))') "h",  "t0","t1",  "ad_ord"
    write(*,'(3(e16.8),i7)')  h,  t0,t1,  ad_ord
    call prarr(X0)
    allocate(res(size(x0)))
    r_rk = 0;r_ae = 0; r_ai = 0
    
    itg_rk = RungeKuttaInt(f,x0,t0,h)
    itg_ae = ExAdamsInt(f, x0, t0, ad_ord, h)
    itg_ai = ImAdamsInt(f, x0, t0, ad_ord, h)
    call print_solution(itg_rk, t1, fd_rk)
    call print_solution(itg_ae, t1, fd_ae)
    call print_solution(itg_ai, t1, fd_ai)
!     call print_poincare_section(itg_rk, t1, fd_rk)
!     call print_poincare_section(itg_ae, t1, fd_ae)
!     call print_poincare_section(itg_ai, t1, fd_ai)

    close(fd_rk)
    close(fd_ae)
    close(fd_ai)
    deallocate(res)
!     write(*,*) "Residues:"
!     write(*,*) "rk", sqrt(r_rk)
!     write(*,*) "ae", sqrt(r_ae)
!     write(*,*) "ai", sqrt(r_ai)
end program cauchy

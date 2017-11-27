program cauchy
    use Poly
    use Integrators
    use Utils
    use IO_array
    implicit none

    integer   :: fd_rk = 10, fd_ae = 11, fd_ai = 12
    real(mpc) :: r_rk, r_ae, r_ai
    real(mpc), allocatable, dimension(:) :: res

    open(fd_rk, file="data/rk.dat", action="write")
    open(fd_ae, file="data/ae.dat", action="write")
    open(fd_ai, file="data/ai.dat", action="write")

    write(*,*) "``Variables``"
    write(*,'(4(a16))') "h",  "t0","t1",  "ad_ord"
    write(*,'(3(e16.8),i7)')  h,  t0,t1,  ad_ord
    call prarr(X0)
    allocate(res(size(x0)))
    r_rk = 0;r_ae = 0; r_ai = 0


    res=runge_ode(t1,fd_rk)
    write(fd_rk,*) t1, res
    res=adams_ex_ode(t1,fd_ae)
    write(fd_ae,*) t1, res
    res=adams_in_ode(t1,fd_ai)
    write(fd_ai,*) t1, res

    close(fd_rk)
    close(fd_ae)
    close(fd_ai)
    deallocate(res)
!     write(*,*) "Residues:"
!     write(*,*) "rk", sqrt(r_rk)
!     write(*,*) "ae", sqrt(r_ae)
!     write(*,*) "ai", sqrt(r_ai)
end program cauchy

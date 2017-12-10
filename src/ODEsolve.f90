module ODEsolve
    use IO_array
    use Integrators
    implicit none
    contains
        subroutine print_solution(integ,t1,fd)
            class(Integrator) :: integ
            real(mpc) :: t1, cache(integ%ord, integ%dimen)
            integer   :: i, N, fd, fd_
            intent(in):: t1, fd
            optional  :: fd
            
            fd_ = stdout
            if (present(fd)) fd_ = fd
            N=int((t1-integ%t0)/integ%h)
            cache = integ%get_cache()
            
            do i = 1, integ%cache_size-1
                write(fd,*) integ%t0 + (i-1)*integ%h, cache(i,:) 
            end do
            do i = integ%cache_size-1, N-1
                write(fd,*) integ%time(), integ%val() 
                call integ%step()
            end do
        end subroutine
end module ODEsolve

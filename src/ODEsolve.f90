module ODEsolve
    use Celmech
    use IO_array
    use Integrators
    implicit none
    contains
        subroutine print_solution(integ,t1,fd, transf)
            class(Integrator) :: integ
            real(mpc), intent(in) :: t1 
            integer  , intent(in), optional :: fd
            procedure(genericTransform), optional :: transf

            real(mpc) :: cache(integ%ord, integ%dimen), v(integ%dimen), v1(integ%dimen), t
            integer :: i,N, fd_

            fd_ = stdout
            if (present(fd)) fd_ = fd
            N=int((t1-integ%t0)/integ%h)
            cache = integ%get_cache()
            
            do i = 1, integ%cache_size-1

                v = cache(i,:)
                if (present(transf)) v(1:tDim) = transf(v(1:tDim))
                write(fd,*) integ%t0 + (i-1)*integ%h, v
            end do
            do i = integ%cache_size-1, N-1
                v = integ%val()
                if (present(transf)) v(1:tDim) = transf(v(1:tDim))
                write(fd,*) integ%time(), v
!                 t = integ%time()
                call integ%step()
!                 v1 = integ%val(); v1(1:tDim) = transf(v1(1:tDim)) - v(1:tDim)
!                 write(*,*) v1(1:tDim)/(integ%time() - t)
            end do
        end subroutine
end module ODEsolve

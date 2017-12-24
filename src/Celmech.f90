module Celmech
    use Const
    use IO_array
    implicit none
    integer, parameter :: spcdim=2, Nbodies = 3
    integer, parameter :: tDim = 2*spcdim*Nbodies
    real(mpc), parameter :: m(Nbodies) = 1.0_mpc
    
contains 
        function motion_eq(tt,X) result(f)
            real(mpc), intent(in) :: tt
            real(mpc), dimension(:), intent(in) :: X
            real(mpc), dimension(size(X)) :: f
            real(mpc) :: f_m(2,Nbodies,spcdim), R(2,Nbodies,spcdim), R_abs(Nbodies)
            real(mpc) :: R_cur(Nbodies,spcdim)
            integer   :: i,j
            ! 1 -- кординаты / скорости
            ! 2 -- тело
            ! 3 -- координаты
            R = reshape(X, (/2,Nbodies,spcdim/), order = (/3,2,1/))

            do i = 1, Nbodies
                forall (j=1:Nbodies) R_cur(j,:) =  R(1,i,:) - R(1,j,:)
                R_abs = (/ (norm2(R_cur(j,:)), j=1,Nbodies) /)

                f_m(2,i,:) = 0
                do j = 1, Nbodies
                    if (j /= i) f_m(2,i,:) = f_m(2,i,:) - m(j)*R_cur(j,:)/R_abs(j)**3 
                end do
            end do
            f_m(1,:,:) = R(2,:,:)
            f(1:tDim/2) = reshape(transpose(f_m(1,:,:)), (/Nbodies*spcdim/))
            f(tDim/2+1:tDim) = reshape(transpose(f_m(2,:,:)), (/Nbodies*spcdim/))
        end function motion_eq
        function fast_motion_eq(tt,X) result(f)
            real(mpc), intent(in) :: tt
            real(mpc), dimension(:), intent(in) :: X
            real(mpc), dimension(size(X)) :: f
            real(mpc) ::  R_abs(Nbodies), R_cur(Nbodies,spcdim)
            integer   :: i,j
            ! 1 -- кординаты / скорости
            ! 2 -- тело
            ! 3 -- координаты
            f(tDim/2+1:) = 0
            do i = 1, Nbodies
                forall (j=1:Nbodies) R_cur(j,:) =  X(ci(1,i,1):ci(1,i,spcdim)) - X(ci(1,j,1):ci(1,j,spcdim))
                R_abs = (/ (norm2(R_cur(j,:)), j=1,Nbodies) /)

                do j = 1, Nbodies
                    if (j /= i) f(ci(2,i,1):ci(2,i,spcdim)) = f(ci(2,i,1):ci(2,i,spcdim)) - m(j)*R_cur(j,:)/R_abs(j)**3 
                end do
            end do
            f(1:tDim/2) = X(ci(2,1,1):)
        contains
            pure function ci(half,body,coord) result(i)
                integer, intent(in) :: half, body, coord
                integer :: i
                i =  (half-1)*tDim/2 + (body-1)*spcdim + (coord-1) + 1
            end function ci
        end function fast_motion_eq
        
        ! dir : +/- 1
        function cyclic_shift_bodies(X, dir) result(X1)
            real(mpc), dimension(tDim), intent(in) :: X
            real(mpc), dimension(tDim) :: X1
            integer :: i, ni, dir
            do i = 1, Nbodies
                ni = mod(Nbodies + i  - 1 + dir, Nbodies) + 1
                X1(ci(1,i,1):ci(1,i,spcdim)) &
                    &= X(ci(1,ni,1): ci(1,ni,spcdim)) 
                X1(ci(2,i,1):ci(2,i,spcdim)) &
                    &= X(ci(2,ni,1):ci(2,ni,spcdim)) 
            end do
        contains 
            pure function ci(half,body,coord) result(i)
                integer, intent(in) :: half, body, coord
                integer :: i
                i =  (half-1)*tDim/2 + (body-1)*spcdim + (coord-1) + 1
            end function ci
        end function cyclic_shift_bodies
end module Celmech

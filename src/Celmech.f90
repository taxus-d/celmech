module Celmech
    use Const
    use IO_array
    implicit none
    integer, parameter :: spcdim=2, Nbodies = 3
    integer, parameter :: tDim = 2*spcdim*Nbodies
    !                                                   rotations + baricentric
    integer, parameter :: shapespDim = spcdim * Nbodies - 1       - spcdim
    real(mpc), parameter :: m(Nbodies) = 1.0_mpc
    real(mpc), parameter, dimension(3) :: c = (/&
        &sqrt((m(2)*m(3))**3 /(m(2) + m(3))),&
        &sqrt((m(1)*m(3))**3 /(m(1) + m(3))),&
        &sqrt((m(2)*m(1))**3 /(m(2) + m(1))) &
    &/)
!     real(mpc), parameter, dimension(3) :: c = 1.0/sqrt(2.0)

    real(mpc), dimension(3,3) :: b  = reshape((/&
        &0.5_mpc, -sqrt(3.0_mpc)/2.0_mpc, 0.0_mpc, & !23
        &0.5_mpc, +sqrt(3.0_mpc)/2.0_mpc, 0.0_mpc,  & !31
        &-1.0_mpc, 0.0_mpc, 0.0_mpc           & !12
        /),(/ 3,3 /))
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
                do j=1, Nbodies
                    R_cur(j,:) = X(ci(1,i,1):ci(1,i,spcdim)) - X(ci(1,j,1):ci(1,j,spcdim))
                    R_abs(j) = sqrt(sum(R_cur(j,:)**2))
                end do

                do j = 1, Nbodies
                    if (j /= i) f(ci(2,i,1):ci(2,i,spcdim)) &
                        &= f(ci(2,i,1):ci(2,i,spcdim)) - m(j)*R_cur(j,:)/R_abs(j)**3 
                end do
            end do
            f(1:tDim/2) = X(ci(2,1,1):tDim)
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


        function shape_motion_eq(t,WdW) result(DWdW)
            real(mpc), dimension(:),intent(in) :: WdW
            real(mpc), dimension(size(WdW))    :: DWdW
            real(mpc),intent(in) :: t
            
            real(mpc), dimension(shapespDim) :: w, dw, d
            real(mpc) :: I, P2, dI2
            integer :: k
            
            w  = WdW(1:shapespDim)
            dw = WdW(shapespDim+1:shapespDim*2)
            
            ! больше 3 всё равно не имеет смысла пока что
            d  = sqrt((/ &
                &norm2(w) - dot_product(w,b(:,1)),&
                &norm2(w) - dot_product(w,b(:,2)),&
                &norm2(w) - dot_product(w,b(:,3))&
            &/))
            d = 2*d**3 
            I  = 2.0_mpc*norm2(w)
            P2 = norm2(dw)**2
            dI2= dot_product(w,dw)*8.0_mpc
            
            DWdW(1:shapespDim) = dw
            do k = 1, shapespDim
                DWdW(shapespDim + k) = -dot_product(c/d, (2.0_mpc*w(k) - I*b(k,:)))&
                    & - 2.0_mpc*P2*w(k)/I**2 + dw(k)*dI2/(2.0_mpc*I**2)  
            end do
!             if (norm2(DWdW) > 1e2) stop '!'
!             write(*,*) DWdW

        contains
            pure function ci(half,body,coord) result(i)
                integer, intent(in) :: half, body, coord
                integer :: i
                i =  (half-1)*shapespDim + (body-1)*spcdim + (coord-1) + 1
            end function ci
        end function shape_motion_eq
end module Celmech

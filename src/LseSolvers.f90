module LseSolvers
    use omp_lib
    use Inival
    use Utils
    use IO_array
    implicit none
    
contains
    function solve_lse(exmtx, arg) result(sols)
        real(mpc), dimension(:, :), intent(in) :: exmtx
        real(mpc), dimension(size(exmtx, 1)) :: sols
        character(len=3) :: arg
        select case (arg)
        case('-go')
            sols = gauss_solve(exmtx, .FALSE.)
        case('-jo')
            sols = jordan_solve(exmtx)
        case('-gm')
            sols = gauss_solve(exmtx, .TRUE.)
        case default
            write (*, '(/5(a/))') "Usage: lsesol [flag]", &
                "Solves system of linear equations", &
                "Flags: -go -- ordinary Gauss method", &
                "       -jo -- with Jordan elimination", &
                "       -gm -- Gauss with impoved numerical stability"
            stop "!"
        end select
    end function solve_lse

 ! для нормировки строки, одна и та же в методе Гаусса и Жордана.
     subroutine norm_line(arr)
         real(mpc), intent(inout), dimension(:) :: arr
         integer :: j
         if (abs(arr(1)) < eps) write (stderr, *) arr(1), "Division may produce incorrect result"
         arr = arr / arr(1)
     end subroutine norm_line

 ! общая функция для вычитания строки (для диагонализации)
     subroutine substract(a)
         real(mpc), intent(inout), dimension(:, :) :: a
         integer :: i, j, n, m
         n = size(a, 1); m = size(a, 2)
         !$omp parallel
         !$omp workshare
         forall (i = 2:n)
             a(i,:) = a(i,:) - a(1,:) * a(i, 1)
         end forall
         !$omp end workshare
         !$omp end parallel
     end subroutine substract

     subroutine remove_sing(a, k, permt)
         real(mpc), intent(inout), dimension(:, :) :: a
         integer, dimension(size(a, 1)), intent(inout) :: permt
         integer :: k, n, pos(2)
         intent(in) :: k
         n = size(a,1)
         pos = maxloc(abs(a(k:n, k:n))) + (/ k-1, k-1 /)
         permt(k) = pos(2) ! номер столбца (ну то есть строки,
                           ! но оно всё транспонированное хранится)
         call swap(a, k, pos(1), a_cols=.TRUE.)
         call swap(a, k, pos(2))
     end subroutine remove_sing

     function gauss_solve(a, is_mod) result(x)
         real(mpc), intent(in), dimension(:, :) :: a
         real(mpc), dimension(size(a, 1), size(a, 2)) :: a1 ! портить матрицу некультурно
         real(mpc), dimension(size(a, 1)) :: x
         real(mpc) :: tmp
         logical, intent(in), optional :: is_mod
         logical :: c_is_mod
         integer :: k, n, permt(size(a, 1))
         n = size(a, 1)
         a1 = a
         if (.not. present(is_mod)) c_is_mod = .TRUE.

         ! тут похоже нигде не распараллелить, проход последовательный
         do k = 1, n
             if (c_is_mod) call remove_sing(a1, k, permt)
             call norm_line(a1(k, k:n+1))
             call substract(a1(k:n, k:n+1))
         end do

         ! и здесь тоже..
         do k = n,1,-1
             x(k) = a1(k, n+1) - dot_product(a1(k, k+1:n), x(k+1:n))
         end do
         if (c_is_mod) then
         do k = n, 1, -1
             tmp = x(k)
             x(k) = x(permt(k))
             x(permt(k)) = tmp
         end do
         end if
     end function gauss_solve

     function jordan_solve(a) result(x)
         real(mpc), intent(in), dimension(:, :) :: a
         real(mpc), dimension(size(a, 1), size(a, 2)) :: a1
         real(mpc), dimension(size(a, 1)) :: x
         integer :: k, n
         n = size(a, 1)
         a1 = a
         do k = 1, n
             call norm_line(a1(k, k:n+1))
             call substract(a1(k:n, k:n+1))
             call reverse_cols(a1(1:k, k:n+1))
             call substract(a1(1:k, k:n+1))
             call reverse_cols(a1(1:k, k:n+1))
         end do
         x = a1(:, n+1)
     end function jordan_solve

end module LseSolvers

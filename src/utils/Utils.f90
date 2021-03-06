module Utils
    use Const
    use omp_lib
    implicit none

contains

    pure function norm2(x) result(a)
        real(mpc), dimension(:), intent(in) :: x
        real(mpc) :: a
        a = sqrt(sum(x**2))
    end function norm2

    pure function zeron(n) result(z)
        integer, intent(in) :: n
        real(mpc), dimension(n) :: z
        z = 0.0_mpc
    end function
  
    function identity(x) result(x1) 
        real(mpc), intent(in), dimension(:) :: x
        real(mpc), dimension(size(x)) :: x1
        x1 = x
    end function

    function dyad_product(a,b) result(c)
        real(mpc), intent(in), dimension(:) :: a,b
        real(mpc), dimension(size(a),size(b)) :: c
        c = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a)) 
    end function dyad_product
    
   
    pure function symplprod(a,b) result(c)
        real(mpc), intent(in), dimension(:) :: a
        real(mpc), intent(in), dimension(size(a)) :: b
        real(mpc), dimension(size(a)) :: c
        integer :: i,j,k,n
        n = size(a)
        c = 0
        do i = 1, n
            do j = 1, n
                do k = 1, n
                    c(i) = c(i) + asimeps(i,j,k)*a(j)*b(k)
                end do
            end do
        end do
    contains
        ! вроде должно работать
        pure function asimeps(i,j,k) result(e)
            integer, intent(in) :: i,j,k
            integer :: e
            e = (k-j)*(j-i)*(k-i)
            e = e/abs(e)
        end function asimeps
    end function
    
    
    subroutine reverse_cols(a)
        real(mpc), dimension(:, :), intent(inout) :: a 
!         real(mpc), dimension(size(a,1), size(a, 2)) :: a1
        integer :: i, m
        m = size(a, 1)
        forall (i = 1:m/2) a(i, :) = a(m+1-i, :)
    end subroutine reverse_cols
    
    subroutine swap(a, i, j, a_cols)
        real(mpc), dimension(:, :) :: a
        integer :: i, j,k
        logical :: a_cols, cols
        intent(in) ::  i, j, a_cols
        optional :: a_cols
        real(mpc) :: tmp
        cols = .FALSE.
        if(present(a_cols)) cols = a_cols

        if (cols) then
            ! вроде присваивает разным кускам массива
            do k = 1, size(a,2)
                tmp = a(i,k)
                a(i,k) = a(j,k)
                a(j,k) = tmp
            end do
        else
            do k = 1, size(a,1)
                tmp = a(k,i)
                a(k,i) = a(k,j)
                a(k,j) = tmp
            end do
        end if

    end subroutine swap
    
    ! joins 2 arrays size by side
    function joinarr_sbs(a, b) result(c)
        real(mpc), intent(in), dimension(:)         :: a, b
        real(mpc), dimension(max(size(a), size(b)), 2) :: c
        integer :: i

        do i = 1, size(c,1)
            c(i,1) = a(i)
            c(i,2) = b(i)
        end do
    end function joinarr_sbs

end module Utils

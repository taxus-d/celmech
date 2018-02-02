module RandomFill
    use Const
    implicit none
    interface randomfill_arr
        module procedure randomfill_arr1
        module procedure randomfill_arr2
    end interface randomfill_arr
contains
    subroutine randomfill_arr1(arr)
        real(mpc), intent(inout) :: arr(:)
        integer, allocatable :: values(:)
        integer :: vsize

        call random_seed(size=vsize)
        allocate(values(vsize))
        call date_and_time(values=values)    

        call random_seed(put=values)
        call random_number(arr)
    end subroutine randomfill_arr1
    subroutine randomfill_arr2(arr)
        real(mpc), intent(inout) :: arr(:, :)
        integer, allocatable :: values(:)
        integer :: vsize

        call random_seed(size=vsize)
        allocate(values(vsize))
        call date_and_time(values=values)    

        call random_seed(put=values)
        call random_number(arr)
    end subroutine randomfill_arr2
end module RandomFill

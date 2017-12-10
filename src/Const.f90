module Const
    implicit none

    integer,parameter    :: mpc=8  ! Параметр разновидности типа для real
    real(mpc), parameter :: EPS = epsilon(1.0_mpc)
    
    !> error codes
    integer,parameter :: &
        &EXIT_SUCCESS = 0&
        &,EXIT_FAILURE = 1&
        &;
end module Const

module Prec
    implicit none

    integer,parameter    :: mpc=8  ! Параметр разновидности типа для real
    real(mpc), parameter :: EPS = epsilon(1.0_mpc)
end module Prec

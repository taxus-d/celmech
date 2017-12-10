module IO
    use ::  iso_fortran_env, &
        stdin => input_unit, stdout => output_unit, stderr => error_unit
    implicit none
contains
    subroutine warn(message)
        character(*) :: message
        write(stderr,*)  message 
    end subroutine warn
end module IO

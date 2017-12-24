module IO
    use ::  iso_fortran_env, &
        stdin => input_unit, stdout => output_unit, stderr => error_unit
    implicit none
contains
    subroutine warn(message)
        character(*) :: message
        write(stderr,*)  message 
    end subroutine warn

    subroutine loud_warn(message)
        character(*) :: message
        write(*,*) "=================================="
        write(*,*) 
        write(stderr,*)  message
        write(*,*) " _  _   _        ___ _"   
        write(*,*) "| \/ \ |_)/\ |\ | | /  |"
        write(*,*) "|_/\_/ | /--\| \|_|_\_ o"
        write(*,*) 
        write(*,*) "=================================="
    end subroutine

end module IO

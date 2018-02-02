module Debug
    use :: IO
    use :: IO_array 
    implicit none
contains
    subroutine step_still_alive(message)
        character(80), optional :: message
        character(80) :: message_ 
        integer :: steps_counter = 0
        
        message_ = "still alive"
        if (present(message)) message_ = message
        
        steps_counter = steps_counter + 1
        write(*,'(i4,1x,a,1x,a)') steps_counter,'--',trim(message_)
    end subroutine step_still_alive
end module Debug

!-------------- cemfDEM code: a dem simulation code ----------------------------
!      D        C enter of
!     D M       E ngineering and
!    D   M      M ultiscale modeling of
!   D     M     F luid flow    
!  EEEEEEEEM    .ir
!------------------------------------------------------------------------------
!  Copyright (C): cemf
!  website: www.cemf.ir
!------------------------------------------------------------------------------  
!  This file is part of cemfDEM code. It is a free software for simulating 
!  granular flow. You can redistribute it and/or modify it under the terms of
!  GNU General Public License version 3 or any other later versions. 
! 
!  cemfDEM code is distributed to help others in their research in the field  
!  of granular flow, but WITHOUT ANY WARRANTY; without even the implied 
!  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!------------------------------------------------------------------------------
!  file name  : g_timer.f90   
!  module name: g_timer 
! 
!  Purpose:              
!    A module to record the execution time of various parts of the code 
! 
!------------------------------------------------------------------------------

module g_timer
    
    use g_TypeDef
    implicit none
    
    private
    public:: timer
    
    type timer
        private
        integer:: numCycle  = 0
        real(8):: tot_time  = 0.0_8
        real(8):: last_time = 0.0_8
        real(8):: time_start, time_end
    contains
        procedure:: reset   => t_reset    ! resetting the timer
        procedure:: start   => t_start    ! start of execution
        procedure:: finish  => t_finish   ! end of the execution, end of event
        procedure:: average => t_average  ! average execution time per event
        procedure:: total   => t_total    ! total execution time of all events 
        procedure:: last    => t_last     ! execution time the last event 
    end type
    
contains

subroutine t_reset( this )
    implicit none
    class(timer) this
    
    this%numCycle = 0
    this%tot_time = 0.0_RK
    this%last_time = 0.0_RK
    
end subroutine

subroutine t_start( this )
    implicit none
    class(timer) this

    call CPU_TIME( this%time_start )
    
end subroutine 

subroutine t_finish( this )

    implicit none
    class(timer) this

    call CPU_TIME( this%time_end )
    this%last_time = this%time_end - this%time_start
    this%tot_time = this%tot_time+ this%last_time
    this%numCycle  = this%numCycle + 1
    
end subroutine 

real(8) function t_average(this)
    
    implicit none
    class(timer) this
    
    if(this%numCycle == 0 )then
        t_average = 0.0_8
    else
        t_average = this%tot_time/this%numCycle
    end if
    
end function

real(8) function t_total(this)
    
    implicit none
    class(timer) this
    
    t_total = this%tot_time
    
end function

real(8) function t_last(this)
    
    implicit none
    class(timer) this
    
    t_last = this%last_time
    
end function


end module

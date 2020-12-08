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
!  file name  : g_error_handling.f90
!  module name: g_error_handling
!
!  Purpose:                                   
!   A module for error handling in the program
! 
!------------------------------------------------------------------------------                                                                  
 
module g_error_handling

    
    use g_Prtcl_DefaultValues
    implicit none
    
    integer,parameter:: ErrT_NoError    = 0
    integer,parameter:: ErrT_Abort      = 1
    integer,parameter:: ErrT_Pass       = 2
    
    character(512):: g_Ch_LastReportedError
    integer          g_ET_LastReportedError
    
    character(11),parameter,dimension(0:2) :: chErrorT = (/"No error   " , "Fatal error" , "Warning    "/)
    
    type error_handling
        integer       :: iErr    = ErrT_NoError    
        character(265):: chError = ""
    contains
    procedure:: CreateMessage       => ErrHndlg_CreateMessage
    
    end type
    
    ! generic overloading 
    interface num2str
        module procedure:: num2strI, num2strR4, num2strR8, num2strInt3
    end interface
    
    
    contains

! creating an error message based on the input arguments 
subroutine ErrHndlg_CreateMessage(this, routine , message , iErr)
    
    implicit none
    class(error_handling) this
    integer,intent(in)      :: iErr
    character(*),intent(in) :: routine , message
    
    integer l_iErr
    
    !if( present(iErr) )
    
    this%iErr = or(this%iErr , iErr )
    this%chError = trim( chErrorT(iErr) )//" occurred in "// trim(routine)// " :"// trim(message)
    

end subroutine 

! integer number to string 
character(64) function num2strI( i )

    implicit none
    !integer(IK) i
    integer(IK), intent(in) :: i
    
    
    character(64) ch
    
    write(ch,"(I20)") i
    num2strI = trim(adjustl(ch))
      
end function

! float number to string 
character(64) function num2strR4( i )

    implicit none
    !integer(IK) i
    real(4), intent(in) :: i

    character(64) ch
    
    write(ch,"(En20.3)") i
    num2strR4 = trim(adjustl(ch))
      
end function

! double number to string  
character(64) function num2strR8( i )

    implicit none
    !integer(IK) i
    real(8), intent(in) :: i

    character(64) ch
    
    write(ch,"(En20.4)") i
    num2strR8 = trim(adjustl(ch))
      
end function

! integer3 variable to string 
character(64) function num2strInt3( i )

    implicit none
    !integer(IK) i
    type(integer3), intent(in) :: i

    character(64) ch
    
    write(ch , *) i%x , "," , i%y, ",", i%z
    num2strInt3 = "["// trim(adjustl(ch)) // "]"
      
end function

!    checking the input error message (Err_type) and creating a message in command window 
! and logfile, then stopping the execution of the program if necessary 
subroutine CheckForError( Err_type, chMethod, chMessage )
    implicit none
    
    integer, intent(in):: Err_type
    character(*),intent(in):: chMethod, chMessage
    
    select case(Err_type)
    case(ErrT_NoError)
        
        ! no error occurred, the program will continue its normal execution
        g_ET_LastReportedError = Err_type
        g_ch_LastReportederror = "No error occurred in :"//trim(chMethod)//". Message: "//trim(chMessage)
        return
        
    case(ErrT_Abort)
        ! a severe error occurred and program should be aborted 
        call MainLogInfo%OutInfo( "A severe error occurred in program", 1 )
        call MainLogInfo%OutInfo( "Error occurred in :"//trim(chMethod), "Error message is:"//trim(chMessage) , 1 )
        g_ET_LastReportedError = Err_type
        g_ch_LastReportederror = "A severe error occurred in :"//trim(chMethod)//". Message: "//trim(chMessage)
        
        stop
        
    case(ErrT_Pass)
        ! a warning occurred in the program, a message will appear on the screen and 
        ! a message will be sent to log file but program continues running
        
        call MainLogInfo%OutInfo( "A warning is reported in program", 1 )
        call MainLogInfo%OutInfo( "A warning occurred in :"//trim(chMethod), "Warning message is:"//trim(chMessage) , 1 )
        g_ET_LastReportedError = Err_type
        g_ch_LastReportederror = "A warning occurred in :"//trim(chMethod)//". Message: "//trim(chMessage)
        
        return
        
    end select
    
    
end subroutine


end module




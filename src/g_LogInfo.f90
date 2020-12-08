!********************************************************************!
!*     This code is a part of the program that accompanies the book *!
!* "Coupled CFD-DEM modeling: Formulation, Implementation, and      *!
!* Applications to multiphase flows" published by Wiley.            *!
!*    For bug reports please contact: Hamid.r.norouzi@gmail.com     *!
!*                                                                  *!
!********************************************************************!
!*    file name  : g_LogInfo.f90                                    *!
!*    module name: g_LogInfo                                        *!  
!*                                                                  *!
!*    purpose:                                                      *! 
!*         a module for managing the log file                       *!
!*                                                                  *!
!*                                                                  *! 
!*                                                                  *! 
!*  Author: H.R. Norouzi           Date: 00:***:2014                *!
!*  Reviewer: H.R. Norouzi         Date: 10:Apr:2015                *!
!*                                                                  *!
!*                                                                  *!
!*                                                                  *!            
!********************************************************************!
module g_LogInfo

    use g_TypeDef
    
    implicit none
   
       
    character(2),dimension(10),parameter:: bullet = (/">>", " >" , "++", " +", "--" ," -" , "**", " *", "==" , " =" /)
    
    type LogInfo
        character(256)::  chLogFileName = "Log file.txt"
        character(256)::  chTitle       = "Log file"
        integer(IK)::   rprt_lvl_file = 2
        integer(IK)::   rprt_lvl_cmdw = 3
        integer(IK)::   unit_file     = 100 
    
    contains
    
    procedure:: set_lvl_file    => LI_set_lvl_file
    procedure:: set_lvl_cmdw    => LI_set_lvl_cmdw
    procedure:: OpenFile        => LI_OpenFile
    procedure:: CloseFile       => LI_CloseFile
    
    procedure:: LI_OutInfo
    procedure:: LI_OutInfo2
    generic:: OutInfo         => LI_OutInfo, LI_OutInfo2
    
   
    end type
    
    interface LogInfo
        procedure:: LogInfo_cunst
    end interface
    
        
contains

function LogInfo_cunst(unit , Dir_Res, RunName , file_lvl, cmdw_lvl ) result(res)
    implicit none
    integer(IK),    intent(in):: unit
    character(*),   intent(in):: Dir_Res, RunName
    integer(IK),    intent(in):: file_lvl, cmdw_lvl
    type(LogInfo)   res
    
    character(256) chFile, chTitle
    
    chFile = trim(Dir_Res)//"Log file - "//trim(RunName)//".txt"
    chTitle= "Log file for run :"// trim(RunName)
    
    call res%OpenFile( unit, chFile , chTitle )
    call res%set_lvl_file(file_lvl)
    call res%set_lvl_cmdw(cmdw_lvl)

end function

subroutine LI_set_lvl_file(this, lvl )
    implicit none
    class(LogInfo) this
    integer(IK),intent(in):: lvl
    this%rprt_lvl_file = lvl
end subroutine 

subroutine LI_set_lvl_cmdw(this, lvl )
    implicit none
    class(LogInfo) this
    integer(IK),intent(in):: lvl
    this%rprt_lvl_cmdw = lvl
end subroutine 

subroutine LI_OpenFile( this, nUnit, chFile, Title )

    implicit none
    class(LogInfo) this
    integer(IK),    intent(in)  :: nUnit
    character(*),   intent(in)  :: chFile, Title
    
    this%unit_file = nUnit
    this%chLogFileName = chFile
    this%chTitle = Title
    
    open(file = chFile, unit = nUnit )
    write(this%unit_file , * ) "**************************************************************************"
    write(this%unit_file , * ) trim(Title)
    write(this%unit_file , * ) "**************************************************************************"
    write(this%unit_file , * ) 


end subroutine 

subroutine LI_CloseFile(this)
    implicit none
    class(LogInfo) this
    
    close(this%unit_file)
end subroutine

subroutine LI_OutInfo( this, chInfo, lvl , no_bull )

    implicit none
    class(LogInfo) this
    character(*),   intent(in):: chInfo
    integer(IK),    intent(in):: lvl
    logical,optional,intent(in):: no_bull
    
    logical l_no
    character(64) ch100, ch101
    
    l_no = .false.
    if( present(no_bull) ) l_no = no_bull
    
   if( lvl > 1 ) then    
      write(ch100,"(A,I3,A)")"(",2*(lvl-1), "x , A2, x ,A)" 
      write(ch101,"(A,I3,A)")"(",2*(lvl-1), "x ,A)"
   else
	
	ch100 = "(x , A2, x ,A)"
	ch101 = "(x ,A)"

    endif
    
    if( lvl <= this%rprt_lvl_file )then
        
        if(lvl == 1 )write(this%unit_file, * )
        if(l_no)then
            write(this%unit_file, ch100 ) trim(chInfo)    
        else
            write(this%unit_file, ch101 ) bullet(lvl), trim(chInfo)    
        end if
        
        
    end if
    
    
    if( lvl <= this%rprt_lvl_cmdw )then
        
        if(lvl == 1 )write(*,*)
        if(l_no)then
            write(*, ch101 ) trim(chInfo)
        else
            write(*, ch100 ) bullet(lvl), trim(chInfo)
        end if
        
        
    end if


end subroutine 


subroutine LI_OutInfo2( this, chInfo , chInfo2, lvl , no_bull )

    implicit none
    class(LogInfo) this
    character(*),   intent(in):: chInfo, chInfo2
    integer(IK),    intent(in):: lvl
    logical,optional,intent(in):: no_bull
    
    logical l_no
    
    l_no = .false.
    if( present(no_bull) ) l_no = no_bull
        
    call this%OutInfo(chInfo, lvl, l_no );
    call this%OutInfo(chInfo2, lvl, l_no );
    
end subroutine
            
end module

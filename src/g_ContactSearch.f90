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
! file name  : g_ContactSearch.f90                        
! module name: g_ContactSearch                             
!                                                               
! purpose:                                                   
!  Managing all operations related to the particle-particle contact search
! 
!------------------------------------------------------------------------------  

module g_ContactSearch
   
    use g_Geometry
    use g_Prtcl_NBS
    use g_Prtcl_NBS_Munjiza
    use g_Prtcl_Hrchl_NBS
    use g_Prtcl_DefaultValues
    use g_Prtcl_Hrchl_Munjiza
    
    implicit none
    
    !// ContactSearch
    type :: ContactSearch
        
        integer (IK) max_nPrtcl
        integer (IK) numPrtcl
        integer (IK) CS_Method
        
        class(base_ContactList),pointer :: PP_cont_List
        integer,allocatable:: dummy
        type(NBS),          pointer:: m_NBS
        type(NBS_Munjiza),  pointer:: m_NBS_Munjiza
        type(NBS_Hrchl),    pointer:: m_NBS_Hrchl
        type(NBS_Munjiza_Hrchl),pointer :: m_NBS_Munjiza_Hrchl
        
        
    contains
    
    procedure:: InitContactSearch => CS_InitContactSearch
   
    procedure:: FindContacts      => CS_FindContacts
    
    procedure:: set_numPrtcl      => CS_set_numPrtcl
    
    procedure:: get_numContact    => CS_get_numContact
    
    final:: FinalizeContactSearch
    
    end type ContactSearch
    
    !// ContactSearch
    
  
contains

!**********************************************************************
! Initializing particle-particle contact search object
!**********************************************************************
subroutine CS_InitContactSearch(this, Method,           &
                                minDomain, maxDomain,   &
                                max_nPrtcl, numPrtcl,   &
                                box_ids, flag, bndg_box,&
                                PP_cont_List,           &
                                ratio , num_lvls )
    implicit none
    class(ContactSearch) this
    integer(IK),    intent(in)  :: Method       ! contact search method
    type(real3),    intent(in)  :: minDomain, maxDomain ! minimum and maximum points of simulation domain
    integer(IK),    intent(in)  :: max_nPrtcl, numPrtcl ! maximum and current number of particles
    integer(IK),pointer,dimension(:),intent(in) :: box_ids, flag ! ids of particles and their flags
    type(real4),pointer,dimension(:),intent(in) :: bndg_box      ! position and diamter of particles 
    class(base_ContactList),pointer, intent(in) :: PP_cont_List  ! particle-particle contact list
    real(RK),       intent(in)  :: ratio                         ! cell size ratio
    integer(IK),    intent(in)  :: num_lvls                      ! number of levels
    optional ratio, num_lvls
    
    real(RK) l_ratio
    
    if( present(ratio) ) then
        l_ratio = ratio
    else
        l_ratio = 1.0_RK
    endif
        
    select case(Method)
    
    case( CSM_NBS )
        
        if( associated(this%m_NBS) ) deallocate(this%m_NBS)
        allocate( this%m_NBS )
        
        call this%m_NBS%Init_NBS(minDomain, maxDomain , max_nPrtcl , numPrtcl, box_ids, flag, bndg_box, l_ratio   )
        call this%m_NBS%setCont_List(PP_cont_List)
        
        
         !>>> log info 
        call MainLogInfo%OutInfo("Contact search method is NBS", 3 )
        call MainLogInfo%OutInfo("Cell size is [m]: "//trim( num2str(this%m_NBS%getCellSize()) ), 4)
        call MainLogInfo%OutInfo("Number of cells considered is (x,y,z) :"//trim( num2str( this%m_NBS%getNumCell()) ), 4 )
        !<<< log info 
        
    case( CSM_NBS_Munjiza )
        
	!print*, "Here", max_nPrtcl
        if( associated(this%m_NBS_Munjiza) ) deallocate(this%m_NBS_Munjiza)
	!print*, "Here2", max_nPrtcl
        allocate( this%m_NBS_Munjiza )
        call this%m_NBS_Munjiza%Init_NBS_Munjiza(minDomain, maxDomain , max_nPrtcl , numPrtcl ,box_ids, flag, bndg_box, l_ratio)
        call this%m_NBS_Munjiza%setCont_List(PP_cont_List)
        
        
        !>>> log info
        call MainLogInfo%OutInfo("Contact search method is NBS Munjiza", 3 )
        call MainLogInfo%OutInfo("Cell size is [m]: "// &
                                 trim( num2str(this%m_NBS_Munjiza%getCellSize()) ), 4)
        call MainLogInfo%OutInfo("Number of cells considered is (x,y,z) :"// &
                                 trim( num2str( this%m_NBS_Munjiza%getNumCell()) ), 4 )
        !<<< log info
        
    case( CSM_NBS_Hrchl )
        
        if( associated(this%m_NBS_Hrchl) ) deallocate( this%m_NBS_Hrchl)
        allocate( this%m_NBS_Hrchl)
        
        !>>> log info
        call MainLogInfo%OutInfo("Contact search method is NBS Hierarchial", 3 )
        !<<< log info
        
        if( present(num_lvls) ) then
            call this%m_NBS_Hrchl%Init_NBS_Hrchl (minDomain, maxDomain, &
                                                  max_nPrtcl , numPrtcl,&
                                                  box_ids, flag, bndg_box  ,num_lvls )
        else
            call this%m_NBS_Hrchl%Init_NBS_Hrchl (minDomain, maxDomain, &
                                                  max_nPrtcl , numPrtcl,&
                                                  box_ids, flag, bndg_box )
        end if
        
        call this%m_NBS_Hrchl%setCont_List(PP_cont_List)
        
                
    case( CSM_NBS_Munjiza_Hrchl )
        
        if( associated(this%m_NBS_Munjiza_Hrchl) ) deallocate( this%m_NBS_Munjiza_Hrchl)
        allocate( this%m_NBS_Munjiza_Hrchl)
        
        !>>> log info
        call MainLogInfo%OutInfo("Contact search method is NBS Munjiza Hierarchial", 3 )
        !<<< log info
        
        if( present(num_lvls) ) then
            call this%m_NBS_Munjiza_Hrchl%Init_Munjiza_Hrchl (minDomain, maxDomain, &
                                                              max_nPrtcl , numPrtcl,&
                                                              box_ids, flag, bndg_box ,num_lvls )
        else
            call this%m_NBS_Munjiza_Hrchl%Init_Munjiza_Hrchl (minDomain, maxDomain, &
                                                              max_nPrtcl , numPrtcl,&
                                                              box_ids, flag, bndg_box )
        end if
        
        call this%m_NBS_Munjiza_Hrchl%setCont_List(PP_cont_List)
        
        
    case( CSM_DESS )
        stop "It is not in the program yet! ;)"
        
        case default
        stop "It is not in the program yet! ;)"
        
    end select
    
    this%PP_cont_List => PP_cont_List
    this%CS_Method = Method
    this%numPrtcl = numPrtcl
    this%max_nPrtcl   = max_nPrtcl
    
        
end subroutine

!**********************************************************************                                
! setting current number of particles
!**********************************************************************
subroutine CS_set_numPrtcl(this, new_numPrtcl )

    implicit none
    class(ContactSearch) this
    integer(IK),intent(in)  :: new_numPrtcl
    
    select case( this%CS_Method )
    
    case( CSM_NBS )
        call this%m_NBS%set_numPrtcl(new_numPrtcl)
        !call this%
    case( CSM_NBS_Munjiza )
        call this%m_NBS_Munjiza%set_numPrtcl(new_numPrtcl)
        
    case( CSM_NBS_Hrchl )
        call this%m_NBS_Hrchl%set_numPrtcl(new_numPrtcl)
    
    case( CSM_NBS_Munjiza_Hrchl )
        call this%m_NBS_Munjiza_Hrchl%set_numPrtcl(new_numPrtcl)
        
    end select
    
    this%numPrtcl = new_numPrtcl
    
end subroutine 

!**********************************************************************
! finding contact pairs of particles
!**********************************************************************
subroutine CS_FindContacts( this )

    implicit none
    class(ContactSearch) this
    
    select case( this%CS_Method)
        
    case( CSM_NBS )
        call this%m_NBS%ContactSearch()
        !call this%
    case( CSM_NBS_Munjiza )
        call this%m_NBS_Munjiza%ContactSearch()
        
    case( CSM_NBS_Hrchl )
        call this%m_NBS_Hrchl%ContactSearch()
    
    case( CSM_NBS_Munjiza_Hrchl )
        call this%m_NBS_Munjiza_Hrchl%ContactSearch()
        
    end select
    
    ! removing previous contacts (those which are not in contact now)
    call this%PP_Cont_List%RemvReleased()

    
    
end subroutine


function CS_get_numContact( this ) result (res)
    
    implicit none
    class(ContactSearch) this
    integer(IK),dimension(2):: res
    
    res = 0
    
    select case( this%CS_Method)
        
    case( CSM_NBS )
        
        res(1) = this%m_NBS%get_num_cntct()
        
    case( CSM_NBS_Munjiza )
        
        res(1) = this%m_NBS_Munjiza%get_num_cntct()
        
    case( CSM_NBS_Hrchl )
        
        res = this%m_NBS_Hrchl%get_num_cntct()
    
    case( CSM_NBS_Munjiza_Hrchl )
        
        res = this%m_NBS_Munjiza_Hrchl%get_num_cntct()
        
    end select

end function 


!**********************************************************************
! final
!**********************************************************************
subroutine FinalizeContactSearch(this)

    implicit none
    type(ContactSearch) this
    
    if( associated(this%m_NBS) ) deallocate(this%m_NBS)
    if( associated(this%m_NBS_Munjiza) ) deallocate(this%m_NBS_Munjiza)
    if( associated(this%m_NBS_Hrchl) ) deallocate(this%m_NBS_Hrchl)

end subroutine 

    
end module

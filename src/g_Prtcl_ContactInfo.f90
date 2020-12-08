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
!  file name  : g_Prtcl_ContactInfo.f90 
!  module name: g_Prtcl_ContactInfo 
!                                   
!  Purpose: 
!   1) providing the base class for a contact pair and a derived class to store
!      contact pair data 
!                        
!    ContactInfo1 class works in conjunction with contact force calculations.
!      For simple linear and non-linear viscoelastic force models, ContactInfo1
!      is enough. This class is used in the ContactList class. But for other 
!      contact force models in which more data should be saved for each contact
!      pair, the class should be extended more to include other data and also a
!      new class should be derived from base_ContactList class and proper methods
!      should be provided for deferred boundings (see what is done for ContactList
!      which uses ContactInfo1).
! 
!------------------------------------------------------------------------------
    
module g_prtcl_ContactInfo
    
    use g_TypeDef
    
    ! the base class for the contact pairs to store their contact data 
    type,abstract:: base_ContactInfo
        integer(IK) pari
        integer(IK) parj ! wall_prop_type
        integer(IK) id_i
        integer(IK) id_j ! wall_id
    contains
    
        ! deferred: these methods should be provided in a class which is 
        ! derived from this class (see ContactInfo1)
        procedure(base_InitContact),deferred :: InitContact
        procedure(base_getIds),     deferred :: getIds
        procedure(base_getpars),    deferred :: getpars
    end type
    
    
    ! abstract interfaces for deferred methods
    abstract interface
    
        subroutine base_InitContact (this, mem_i, mem_j, id_i, id_j )
            import base_ContactInfo, IK, RK
            class(base_ContactInfo)   this
            integer(IK),intent(in) :: mem_i, mem_j, id_i, id_j
        end subroutine 
        
        function base_getIds( this )  result (res)
            import base_ContactInfo, IK
            class(base_ContactInfo)    this
            integer(IK),dimension(2):: res
        end function
        
        function base_getpars( this )  result (res)
            import base_ContactInfo, IK
            class(base_ContactInfo)   this
            integer(IK),dimension(2):: res
        end function
        
    end interface
    
    
    ! this class works with ContactList and is enough for storing 
    ! the contact pair data of linear and non-linear viscoelastic force models
    type, extends(base_ContactInfo):: ContactInfo1
        type(real3) tang_del
    contains
        ! deferred boundings 
        procedure:: InitContact => CI_InitContact
        procedure:: getIds      => CI_getIds
        procedure:: getpars     => CI_getpars
        
    end type
    
          
    
    
contains
!********************************************************************
! Initializing the contact pair data sets the 
!********************************************************************
subroutine CI_InitContact(this, mem_i, mem_j, id_i, id_j )
    implicit none
    class(ContactInfo1)        this
    integer(IK),intent(in)  :: mem_i, mem_j, id_i, id_j
    
    this%pari = mem_i
    this%parj = mem_j
    this%id_i = id_i
    this%id_j = id_j
    this%tang_del = real3( 0.0_RK, 0.0_RK, 0.0_RK)
    
end subroutine

!********************************************************************
! returning ids of contact pair 
!********************************************************************
function CI_getIds( this )  result (res)
    implicit none
    class(ContactInfo1):: this
    integer(IK),dimension(2):: res
    
    res = (/this%id_i, this%id_j/)
    
end function

!********************************************************************
! returning memory index of contact pair
!********************************************************************
function CI_getpars( this )  result (res)
    implicit none
    class(ContactInfo1):: this
    integer(IK),dimension(2):: res
    
    res = (/this%pari, this%parj/)
    
end function
    
end module

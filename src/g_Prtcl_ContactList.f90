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
!  file name  : g_Prtcl_ContactList.f90
!  module name: g_Prtcl_ContactList
!  
!  Purpose: 
!   1) To provide the base class for contact list (a contact list stores all
!      the contact pairs)
!   2) To manage operations like inserting, updating and deleting contact
!      pairs
!   3) To provide methods which enable performing a loop over all contact pairs
!      in the list
!   4) To remove previous contact pairs from list
!                                                                    
!    The first class, which is derived from base_ContactList, is ContactList.
!      This class works in conjunction with ContactInfo1 class. Since the
!      ContactInfo1 class is used for linear and nonlinear viscoelastic force
!      models, this contact list also is applicable to these models. If you wish
!      to extend it for other contact force models, you should extend a new 
!      class from the base class and provide definitions for deferred methods in
!      the new class (see what is done for ContactList).
! 
!------------------------------------------------------------------------------


module g_Prtcl_ContactList
    
    use g_LinkedList
    use g_prtcl_ContactInfo
    use g_error_handling
    
    implicit none
    
    
    ! Base class for contact list, it stores all the contact pairs and manages
    ! operations like, inserting, updating and deleting contact pairs. 
    ! It also has some methods to perform an iteration/loop on the contact pairs. 
    
    type, abstract:: base_ContactList
        
        integer(IK)         :: numPrtcl
        integer(IK)         :: numCntcts
        integer(IK)         :: max_numCntcts
        type(LinkedList)    :: List
        integer(IK),dimension(:),allocatable    :: ContctStat
        
        !for iteration
        integer(IK):: curr_item
        integer(IK):: numProcessed
        
    contains
        
        procedure,non_overridable:: getNumCntcts    => BCL_getNumCntcts  !
        procedure,non_overridable:: setNumCntcts    => BCL_setNumCntcts  !
        procedure,non_overridable:: getMaxCntcts    => BCL_getMaxCntcts  !
        
        ! for performing iteration/loop over contact pairs 
        procedure,non_overridable:: Start_iter        => BCL_Start_Iter    !
        procedure,non_overridable:: nextItem          => BCL_nextItem      !
        procedure,non_overridable:: IsFinished        => BCL_IsFinished    !
        
        ! deferred bindings 
        procedure(base_InitContactList),deferred:: InitContactList
        procedure(base_AddContact),deferred     :: AddContact
        procedure(base_RemvReleased),deferred   :: RemvReleased
    end type
    
    
    ! abstract interfaces for base_ContactList class 
    abstract interface
        subroutine base_InitContactList(this, max_numCntcts , numPrtcl )
            import base_ContactList, IK
            class(base_ContactList)   this
            integer(IK),intent(in) :: max_numCntcts, numPrtcl
        end subroutine
        
        subroutine base_AddContact(this , Mem_i , Mem_j, id_i, id_j )
            import base_ContactList, IK
            class(base_ContactList)   this
            integer(IK),intent(in) :: Mem_i , Mem_j, id_i, id_j
        end subroutine 
        
        subroutine base_RemvReleased(this)
            import base_ContactList
            class(base_ContactList)  this
        end subroutine 
        
    end interface
    
    
    type, extends(base_ContactList):: ContactList
                
        type(ContactInfo1),pointer,dimension(:) :: ContPair
        
    contains
    
        
        procedure:: getItem         => CL_getItem         ! 
        procedure:: setItem         => CL_setItem         ! 
                
        procedure:: InitContactList => CL_InitContactList ! deferred
        procedure:: AddContact      => CL_AddContact      ! deferred
        procedure:: RemvReleased    => CL_RemvReleased    ! deferred
        
        procedure,private:: DeallocateAll   => CL_DeallocateAll !
        
        final :: CL_Final
    end type
    
        
contains

!**********************************************************************
! returning number of contact pairs
!**********************************************************************
integer(IK) function BCL_getNumCntcts(this)
    implicit none
    class(base_ContactList) this
    
    BCL_getNumCntcts = this%numCntcts
    
end function

!**********************************************************************
! setting number of contact pairs
!**********************************************************************
subroutine BCL_setNumCntcts(this , val)
    implicit none
    class(base_ContactList)   this
    integer(IK),intent(in)  :: val
    
    this%numCntcts = val
    
end subroutine

!**********************************************************************
! getting maximum allowable number of contact pairs
!**********************************************************************
integer(IK) function BCL_getMaxCntcts(this)
    implicit none
    class(base_ContactList)::this
    
    BCL_getMaxCntcts = this%Max_numCntcts
    
end function

!***********************************************************************************************
! signaling the contact list to put the iteration position at the beginning of contact list and 
! returning index of contact pair in the 1D vector which stores them.
!  returns:
!      -1 : No contact pair is available
!       0 : cannot find the first contact pair in the list
!      >0 : index of contact pair in the 1D vector which stores them
!**********************************************************************************************
integer(IK) function BCL_Start_Iter(this , lnew )
    
    implicit none
    class(base_ContactList) this
    logical,intent(out)::   lnew
    
    ! locals
    integer(IK) i , nCntct, maxCntct 
    
    ! body
    
    ! getting number of contact pairs in the list
    nCntct = this%getNumCntcts()
    
    if(nCntct <= 0 ) then
         BCL_Start_Iter = -1
         return
    end if
    
    
    this%numProcessed = 0
    maxCntct = this%getMaxCntcts()
    
    !Starting the loop over the list and finding the first one in the list
    do i=1,maxCntct
        
        ! ContactStat >0 means a contact pair exists in this location
        if( this%ContctStat(i) > 0  ) then
            
            if( this%ContctStat(i) == 1 ) then
                lnew = .false.  ! old contact pair
            else
                lnew = .true.   ! new contact pair
            end if
            
            BCL_Start_Iter = i      ! index of Contact pair vector
            this%ContctStat(i) = 0  ! means processed 
            this%numProcessed = 1   ! sets number of processed items to 1 
            this%curr_item = i      ! saves it to start the loop from the rest of list - the items before i are already processed. 
            return
       
        end if
        
    end do
    
    BCL_Start_Iter = 0
    
end function

!********************************************************************************
! returning the index of next contact pair in the contact list
! it should be invoked after invoking Start_Iter method in this class
! returns
!    -1) no next item can be found
!     0) the last item has already processed and the list is finished
!   >0 ) the index of contact pair in the contact list
!********************************************************************************
integer(IK) function BCL_nextItem( this, lnew )
    implicit none
    class(base_ContactList)    this
    logical,intent(out)     :: lnew
    
    ! locals
    integer(IK) i , maxCntct 
    
    ! body
    
    ! check if the list is finished (all items are processed?)
    if( this%IsFinished() ) then
        BCL_nextItem = 0
        return
    end if
    
    ! getting maximum number of contacts
    maxCntct = this%getMaxCntcts()
    
    do i=this%curr_item+1, maxCntct
        
        select case (this%ContctStat(i) )
        case(1)
            lnew = .false.
            BCL_nextItem = i       ! index of ContPair vector
            this%ContctStat(i) = 0 ! means processed 
            this%numProcessed = this%numProcessed + 1  
            this%curr_item = i
            return
            
        case(2) 
            lnew = .true.
            BCL_nextItem = i       ! index of ContPair vector
            this%ContctStat(i) = 0 ! means processed 
            this%numProcessed = this%numProcessed + 1  
            this%curr_item = i
            return
            
        end select
        
    end do
    
    BCL_nextItem = -1       
    
end function

!**********************************************************************
! determining if all the items have been processed. 
!**********************************************************************
logical function BCL_IsFinished(this)
    implicit none
    class(base_ContactList) this
    
    if(this%numProcessed >= this%numCntcts ) then
        BCL_IsFinished = .true.
    else
        BCL_IsFinished = .false.
    end if

end function

!**********************************************************************
! setting the contact pair data in the container based on the index of 
! container vector
!**********************************************************************
subroutine CL_setItem(this, ind , item )
    implicit none
    class(ContactList) this
    integer(IK),intent(in)        :: ind
    type(ContactInfo1),intent(in) :: item
    
    this%ContPair(ind) = item

end subroutine

!**********************************************************************
! returning the contact data based on the index of conainer vector 
!**********************************************************************
type(ContactInfo1) function CL_getItem(this, ind )
    implicit none
    class(ContactList)       this
    integer(IK),intent(in) :: ind
    
    CL_getItem = this%ContPair(ind)

end function

!**********************************************************************
! Initializing the contact list for ContactList class
!**********************************************************************
subroutine CL_InitContactList(this, max_numCntcts , numPrtcl )
    implicit none
    class(ContactList)      this
    integer(IK),intent(in):: max_numCntcts, numPrtcl
    
    
    ! body
    
    this%Max_numCntcts = max_numCntcts
    this%numPrtcl = numPrtcl
    this%numCntcts = 0
    
    ! initializing the linked list to stores id pairs 
    call this%List%InitList( numPrtcl, max_numCntcts )
    
    ! deallocating the memory 
    call this%DeallocateAll()
    
    ! allocating a vector to store the status of contact pairs
    allocate( this%ContctStat(max_numCntcts) )
    this%ContctStat = -1
    
    ! allocating a vector to store the contact data of each pair 
    allocate( this%ContPair(max_numCntcts) )
    
end subroutine

!**********************************************************************
! Adding a contact pair to the contact list 
!**********************************************************************
subroutine CL_AddContact(this , Mem_i , Mem_j, id_i, id_j )

    implicit none
    class(ContactList):: this
    integer(IK),intent(in)  :: Mem_i , Mem_j, id_i, id_j
    
    integer(IK) l_mem_i, l_mem_j, l_id_i, l_id_j
    integer(IK),dimension(2):: ins_item
    
    ! this is a convention, the lower id should be the first item in the 
    ! contact pair
    if( id_i < id_j ) then
        l_mem_i = Mem_i
        l_mem_j = Mem_j
        l_id_i = id_i
        l_id_j = id_j
    else
        l_mem_i = Mem_j
        l_mem_j = Mem_i
        l_id_i = id_j
        l_id_j = id_i
    end if
    
    
        
    ins_item = this%List%Find_Insert( l_id_i , l_id_j )
    
    ! ins_item(1) is the status of the item insertion 
    ! ins_item(2) is the container index
    select case( ins_item(1) )
    
    case (1_IK) ! this is an existing contact
        
        ! the contact is kept in its place in the container
        this%ContctStat( ins_item(2) ) = 1_IK 
        this%numCntcts = this%numCntcts + 1
        
    case(2_IK)! this is new contact
        
        this%ContctStat( ins_item(2) ) = 2_IK ! this is a flag to show, this is a new contact
        this%numCntcts = this%numCntcts + 1
        
        ! the new contact is inserted into the list and the corresponding Contact pair is initiated
        call this%ContPair(ins_item(2))%InitContact( l_mem_i, l_mem_j, l_id_i , l_id_j )
        
        
        case default ! there is an error in the item insertion
        
            if (ins_Item(1) == -1 ) then
               
                call CheckForError( ErrT_Abort , "CL_AddContact" , &
                                    "The inserted item's id is greater than allowed value:"//num2str(l_id_i) )
                 
            else
                
                call CheckForError( ErrT_Abort , "CL_AddContact" , "The container is full and there is no space for new item" )
                
            endif
                       
    end select
    
end subroutine

!**********************************************************************
! Removing all released contacts (those which are not in contact in this time step) 
! from contact list.
!**********************************************************************
subroutine CL_RemvReleased(this)

    class(ContactList) this
    
    integer(IK) i, res
    integer(IK),dimension(2):: pair_id
    
    do i=1, this%Max_numCntcts
        
        if(this%ContctStat(i) == 0 ) then
            this%ContctStat(i) = -1
            pair_id = this%ContPair(i)%getIds()
            res = this%List%DeleteItem( pair_id(1) , pair_id(2) )
        end if
        
    end do
    
end subroutine

!**********************************************************************
! the final (destructor) for the class
!**********************************************************************
subroutine CL_Final(this)
    implicit none
    type(ContactList) this
    
    call this%DeallocateAll()
    
end subroutine 

!**********************************************************************
! deallocating the previously allocated memory for this class
!**********************************************************************
subroutine CL_DeallocateAll(this)
    implicit none
    class(ContactList):: this
    
    if(allocated(this%ContctStat)) deallocate(this%ContctStat)
    if(associated(this%ContPair)) deallocate(this%ContPair)
    
end subroutine
  
    
end module

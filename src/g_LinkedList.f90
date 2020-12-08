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
!  file name  : g_LinkedList.f90 
!  module name: g_LinkedList     
!                            
!  Purpose:
!   A module for creating and managing linked-lists
! 
!------------------------------------------------------------------------------

module g_LinkedList
    
    use g_TypeDef
    implicit none  
    
    type item
        integer(IK) value
        integer(IK) next
    end type
    
    type LinkedList
        
        integer(IK):: numBucket
        integer(IK):: numEntry
        integer(IK):: nextInsert
        integer(IK),dimension(:),allocatable:: Bucket
        type(item) ,dimension(:),allocatable:: Container
    
    contains
        procedure:: InitList        => LL_InitList
        procedure:: getNumEntry     => LL_getNumEntry
        procedure:: getNumBucket    => LL_getNumBucket
        procedure:: InsertItem      => LL_InsertItem
        procedure:: FindItem        => LL_FindItem
        procedure:: Find_Insert     => LL_Find_Insert
        procedure:: DeleteItem      => LL_DeleteItem
        
        procedure:: emptyContainer => LL_emptyContainer
        
        procedure::         printAll     ! prints linked list 
        procedure,private:: pushItem            => LL_pushItem
        procedure,private:: DeAllocAll
        procedure,private:: refreshContainer
    end type
    

contains

    subroutine printAll(this, bkt_id )
        
        implicit none
        class(LinkedList) this
        integer(IK) bkt_id
        
        integer(IK) n
    
        n = this%Bucket(bkt_id) 
        
        print*, " Array Loc -    Value  "
        
        do while ( n > 0 )
        
            print '(I10 , "  -" , I8)', n ,this%Container(n)%value        
            n = this%Container(n)%next
        
        end do
        
    end subroutine
     
    ! initializing the linked list for the first use 
    subroutine LL_InitList( this, nBucket, nEntry )
        
        implicit none
        class(LinkedList) this
        integer(IK) nBucket
        integer(IK) nEntry
        
        call this%DeAllocAll()
        
        allocate( this%Bucket(nBucket) )
        this%numBucket = nBucket
        this%Bucket = 0_IK
        
        
        allocate( this%Container(nEntry) )
        this%numEntry = nEntry
        call this%refreshContainer()
        
    end subroutine
    
    ! returning number of entries in the list     
    integer function LL_getNumEntry(this)
        implicit none
        class(LinkedList) this
        
        LL_getNumEntry = this%numEntry
        
    end function
    
    ! returning number of buckets in the list 
    integer function LL_getNumBucket(this)
        implicit none
        class(LinkedList) this
        
        LL_getNumBucket = this%numBucket
        
    end function
    
    !********************************************************************************  
    ! inserting an item based on bucket id (bkt_id) and a value for that bucket 
    !********************************************************************************
    integer(IK) function LL_InsertItem(this , bkt_id , value )
        
        implicit none
        class(LinkedList) this
        integer(IK) bkt_id, value
        
        
        integer(IK) nBucket, nextI
        
        nBucket = this%getNumBucket()
        if( bkt_id > nBucket ) then
            ! The bucket id is not in the range and must return 0
            ! this is an error for the linked list
            LL_InsertItem = 0 
            return
        end if
    
        nextI = this%nextInsert 
        if( nextI == 0 ) then
            !the container is full and there is no more space 
            !for this item
            LL_InsertItem = -1
            return
        endif
        
        this%nextInsert = - this%container(nextI)%next ! puts the next item 
        this%container(nextI) = item( value , this%bucket(bkt_id) ) 
        this%bucket(bkt_id) = nextI
        
        LL_InsertItem = nextI ! returning the index of container array in which the item is inserted
        
        
    end function
    
    !*************************************************************************************
    ! pushing an item to the end of the list, no check is done 
    !*************************************************************************************
    integer(IK) function LL_pushItem(this , bkt_id , value )
        
        implicit none
        class(LinkedList) this
        integer(IK) bkt_id, value
        
        
        integer(IK) nBucket, nextI
        
            
        nextI = this%nextInsert     
        this%nextInsert = - this%container(nextI)%next ! putting the next item 
        this%container(nextI) = item( value , this%bucket(bkt_id) ) 
        this%bucket(bkt_id) = nextI
        
        LL_pushItem = nextI ! returning the index of container array in which the item is inserted
        
        
    end function
    
    !***********************************************************************************************
    !* searching for an item, if it exists, it would return (/1, index of container that contain the item/),
    !     if it is new, pushing the item into the list and return (/2, index of container that
    !     contains the new item /), otherwise (-2,-2). 
    !********************************************************************************************
    function LL_Find_insert(this, bkt_id , value ) result(res)
    
        implicit none
        
        class(LinkedList) this
        integer(IK) bkt_id, value
        integer(IK),dimension(2):: res
        
        
        integer(IK) nBucket, n
                
        nBucket = this%getNumBucket()
        if( bkt_id > nBucket ) then
            ! The bucket id is not in the range and must return (-1,-1)
            ! this is an error for the linked list
            res = (-1_IK , -1_IK)
            return
        end if
        
        n = this%Bucket(bkt_id) 
    
        do while ( n > 0 )
        
            if( this%Container(n)%value == value ) then        
               ! the item already exists in the list
               ! function returns the proper code and the index of container in which the item exists 
                res = (/1_IK , n /) 
               return
            end if
            n = this%Container(n)%next
            
        end do
        
         
        if( this%nextInsert == 0 ) then
            !the container is full and there is no more space 
            !for this item
            res = (/-2_IK , -2_IK /)
            return
        endif
        
        ! the item is new and it should be pushed into the list
        ! inserting new item in the list 
        res = (/ 2_IK , this%pushItem(bkt_id , value) /)
        
    end function
    
    !***************************************************************************
    ! finding the index of the item in the container 
    !***************************************************************************
    integer(IK) function LL_FindItem(this, bkt_id , value )
    
        implicit none
    
        class(LinkedList) this
        integer(IK) bkt_id, value
    
        integer(IK) n
    
        n = this%Bucket(bkt_id) 
    
        do while ( n > 0 )
        
            if( this%Container(n)%value == value ) then
                LL_FindItem = n
                return
            end if
            
            n = this%Container(n)%next
        
        end do
    
        LL_FindItem = -1
        
    end function
    
    !*******************************************************************
    !* deleting an item from container
    !*******************************************************************
    integer(IK) function LL_DeleteItem(this, bkt_id, value )
        
    implicit none
    class(LinkedList) this
    integer(IK) bkt_id, value
        
        integer(IK) n, prev
        
        
        n = this%Bucket(bkt_id) 
        prev = 0 
        
        do while ( n > 0 )
        
            if( this%Container(n)%value == value ) then
                
                if( prev .eq. 0 ) then
                    
                    this%Bucket(bkt_id) = this%Container(n)%next
                    this%Container(n)%next = - this%nextInsert
                    this%nextInsert = n 
                    LL_DeleteItem = n
                    return
                else
                    this%Container(prev)%next = this%Container(n)%next
                    this%Container(n)%next = - this%nextInsert
                    this%nextInsert = n
                    LL_DeleteItem = n 
                    return
                end if
                
                
            endif
            
            prev = n
            n = this%Container(n)%next
            
        end do
    
        LL_DeleteItem = 0 ! means this item is not found 
    
    end function
    
    subroutine DeAllocAll(this)
        
        implicit none
        class(LinkedList) this
        
        if(allocated(this%Bucket)) deallocate(this%Bucket)
        if(allocated(this%Container)) deallocate(this%Container)
    end subroutine
    
    
    subroutine LL_emptyContainer(this)
        implicit none
        class (LinkedList) this 
        this%Bucket = 0_IK
        call this%refreshContainer()
        
    end subroutine 
    
    subroutine refreshContainer(this)
    
        implicit none
        class(LinkedList) this
        
        integer(IK) i, nEntry
        
        nEntry = this%getNumEntry()
        do i = 1, nEntry-1
            this%Container(i)%Next = -(i+1_IK)
        end do
        
        this%Container(nEntry)%Next = 0 ! blocks end of container 
        this%nextInsert = 1_IK
    end subroutine
    
    
    
end module

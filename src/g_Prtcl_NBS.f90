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
!  File name:  g_Prtcl_NBS.f90
!  Module name: g_Prtcl_NBS
!  
!  Purpose: 
!    1) Providing a basic class for NBS contact search algorithms in this code
! 
!    2) It is extended based on the CellBased class and contains necessary
!       objects and methods for doing NBS contact search
! 
!  External literature used:
!   - Munjiza, A. (2004). The Combined Finite-Discrete Element Method. West
!     Sussex, UK., John Wiley & Sons Ltd.
! 
!   - Munjiza, A. and K. R. F. Andrews (1998). "NBS contact detection algorithm
!     for bodies of similar size. International Journal for Numerical Methods
!     in Engineering 43: 131-149.
! 
!   Referenced modules: 
!       - g_Prtcl_CellBased 
! 
!------------------------------------------------------------------------------
    
    module g_Prtcl_NBS
    
    use g_Prtcl_CellBased
    use g_error_handling
    
    implicit none
    
    private
    public:: NBS
    
    ! all adjacent cells for performing no binary search
    type(integer3),parameter,dimension(0:13):: NBS_Mask &
                                      = (/ integer3( 0, 0, 0), &
                                           integer3( 0, 0,-1), &
                                           integer3(-1, 0,-1), &
                                           integer3(-1, 0, 0), &
                                           integer3(-1, 0, 1), &
                                           integer3( 0,-1,-1), &
                                           integer3( 0,-1, 0), &
                                           integer3( 0,-1, 1), &
                                           integer3(-1,-1,-1), &
                                           integer3(-1,-1, 0), &
                                           integer3(-1,-1, 1), &                                           
                                           integer3( 1,-1,-1), &
                                           integer3( 1,-1, 0), &
                                           integer3( 1,-1, 1) /)                                          

    

    !// NBS
    type, extends (CellBased):: NBS
        private
        ! head list of boxes
        integer(IK),dimension(:,:,:),allocatable:: HeadList
        ! next list of boxes
        integer(IK),dimension(:),allocatable    :: NextList
        

    contains
       
        !initializing NBS object 
        procedure:: Init_NBS        => NBS_Init_NBS
        procedure:: Allocate_NBS    => NBS_Allocate_NBS
        
        ! performing contact search based on the NBS method
        procedure:: ContactSearch   => NBS_ContactSearch
        
        ! performing broad contact search 
        procedure:: BroadSearch     => NBS_BroadSearch
        
        ! constructing contact list (HeadList and NextList)
        procedure:: CnstrctList     => NBS_CnstrctList
        
        ! returning box index in the head(i,j,k) of the contact list
        procedure:: Head            => NBS_Head
        procedure:: Head_Inc        => NBS_Head_Inc ! adds numPrtcl when returns 
        
        ! returning index of the next box in the list 
        procedure:: Next            => NBS_Next
        
        ! 
        procedure:: Ref_Head        => NBS_Ref_Head     
        
        ! freeing memory allocated for NBS object 
        final    :: FinilizeNBS

    end type NBS
    !// NBS
    
    
    
contains

!****************************************************************!   
!************************ NBS methods ***************************!
!****************************************************************!    


!******************** NBS_Init_NBS (Init_NBS) ***********************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          Initializing NBS object                                 *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in:                                                       *!
!*          this    : polymorphic NBS object                        *!
!*          minDomain: lower limit of simulation domain             *!
!           maxDomain: higher limit of simulation domain            *!
!*          max_nPrtcl  : maximum number of particles               *!
!*          numPrtcl: number of particles                           *!
!*          ids     : array pointer of bounding box ids             *!
!*          bndg_box: array pointer of bounding box positions       *!
!*          ratio   : (optional) the ratio of cell length to        *!
!*                     maximum diamter of boxes                     *!
!*                                                                  *!  
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!
    subroutine NBS_Init_NBS(this, minDomain, maxDomain, max_nPrtcl, numPrtcl, ids, flag, Bndg_box , ratio )
        implicit none
        class(NBS)  this
        type(real3),                        intent(in) :: minDomain, maxDomain
        integer(IK),                        intent(in) :: max_nPrtcl, numPrtcl
        integer(IK),pointer,dimension(:),   intent(in) :: ids, flag
        type(real4),pointer,dimension(:),   intent(in) :: Bndg_box
        real(RK),optional,                  intent(in) :: ratio
        
        ! locals
        real(RK) l_ratio
                

        if( present(ratio) ) then
            l_ratio = ratio
        else
            l_ratio = 1.0_RK
        endif
        
        
        call this%setDomain(minDomain, maxDomain)
        
        call this%InitSimWorld( max_nPrtcl, numPrtcl, ids, flag, Bndg_Box, .true. )
        
        call this%initCellBased( l_ratio)
        
        call this%Allocate_NBS()
    
    end subroutine
    
    !******************************************************************
    ! allocating memory for NBS contact search object
    !******************************************************************
    subroutine NBS_Allocate_NBS(this )
        implicit none
        class (NBS):: this
    
        ! // locals
        integer         iErr
        integer(IK)     max_nPrtcl
        type(integer3)  numCells
        
        
        numCells = this%getNumCell()
        
        ! one cell is added to the lower and upper limits of the cells in each 
        ! direction to apply the NBS mask at borders without using if condition 
        if( allocated(this%HeadList) ) deallocate( this%HeadList )        
        allocate( this%HeadList(0:numCells%x+1,0:numCells%y+1,0:numCells%z+1) , STAT = iErr )
        
        if( iErr .ne. 0)then
            call CheckForError( ErrT_Abort, "NBS_Allocate_NBS" , &
                                "Allocation of HeadList failed with size of" // num2str(numCells) )
            return
        end if
        
        ! nullifying Head
        this%HeadList = -1
        
        
        max_nPrtcl = this%get_max_nPrtcl()
        if( allocated(this%NextList) ) deallocate( this%NextList )
        allocate( this%NextList(max_nPrtcl), STAT = iErr )
        
        if( iErr .ne. 0)then
            call CheckForError( ErrT_Abort, "NBS_Allocate_NBS" , &
                                "Allocation of NextList failed with size of" // num2str(max_nPrtcl) )
            return
        end if
        
        this%NextList = -1
        

    end subroutine

!************* NBS_ContactSearch (ContactSearch) ********************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          performing all operations required for a full contact   *!
!*          search                                                  *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in: none                                                  *!
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!
    subroutine NBS_ContactSearch(this)
        implicit none
        class(NBS) this
        
        ! first setting number of conservative contacts to zero
        call this%set_num_cntct(0_IK)
        
        ! finding the integer coordinates of all boxes
        call this%BoxIndex()
        
        ! constructing the contact list 
        call this%CnstrctList()
        
        ! performing broad search (fine search is also invoked in this procedure) 
        call this%BroadSearch()

    end subroutine 
    
!**************** NBS_BroadSearch (BroadSearch) *********************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          performing no binary operations required for a full     *!
!*          contact search                                          *!
!*          also invoking the fine search routine to determine      *!
!*          particle-particle and wall-particle contacts            *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in: none                                                  *!
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!    
    subroutine NBS_BroadSearch(this)
        implicit none
        class(NBS) this
        
        ! local variables
        integer(IK) ix, iy, iz, iix, iiy, iiz
        integer(IK) m,n, i, icell
        integer(IK) nPrtcl
        type(integer3)  ind, iind
        
        
       ! getting current number of particles 
        nPrtcl = this%get_numPrtcl() 
       
        do i = 1 , nPrtcl 
            
            ! performing contact search for particles which are in simulation domain
            if ( this%bndg_flag( this%getMemIndx(i) ) >= Pflg_inDomain ) then
                
                ind = this%getIndex(i) ! getting integer coordinates of box        
            
                n   = this%Head_Inc(ind) 
            
                if( n <= nPrtcl ) then ! checking if this is the first time that this cell is considered as target cell
                                
                    !looping over all particles in the same cell (ix, iy, iz)
                    do while (n .ne. -1)
                    
                    
                        ! checking against all boxes in the same cell but not the same box 
                        m = this%Next(n)
                        do while (m .ne. -1)
                            ! performing the fine search
                            call this%FineSearch(this%getMemIndx(n),this%getMemIndx(m) )
                        
                            call this%add_num_cntct(1_IK)
                            m = this%Next(m)
                        end do

                        ! checking against all 13 adjacent cells based on the NBS mask 
                        do icell = 1, 13
                        
                            iind = ind + NBS_Mask(icell)
                            if( iind%z < 0 ) then
                                print*, i , ind, icell
                            
                            end if
                        
                            m = this%Head(iind) 
                            if( m > nPrtcl )  m = m - nPrtcl
                        
                            ! looping over all boxes in the adjacent cells based on the NBS mask 
                            do while( m .ne. -1 )
                                ! performing the fine search
                                call this%FineSearch(this%getMemIndx(n),this%getMemIndx(m) )
                                call this%add_num_cntct(1_IK)
                                m = this%Next(m)
                            end do

                        end do

                        n = this%Next(n)
                    end do

                
                end if
            
            end if
            
        end do

    end subroutine 
    
    !******************************************************************
    ! constructing the linked list of particles 
    !******************************************************************
    subroutine NBS_CnstrctList(  this )
        implicit none
        class(NBS) :: this
        

        !// local variables 
        integer n, nPrtcl
        type(integer3) ind

        !// body
        
        ! current number of particles
        nPrtcl = this%get_numPrtcl()
        
        ! first nullifying the list
        this%HeadList = -1
        
        do n = 1, nPrtcl
            
            ! for particles which are in the simulation domain
            if( this%bndg_flag( this%getMemIndx(n) ) >= Pflg_inDomain ) then
                
                ! getting integer coordinates of box
                ind = this%getIndex(n) 
                ! pushing the box number into the head of the list
                this%NextList(n) = this%HeadList(ind%x,ind%y,ind%z)
                this%HeadList(ind%x,ind%y,ind%z) = n
                
            end if
            
        enddo
       
    end subroutine

    !******************************************************************
    ! returning the head value of a cell (ind)
    !******************************************************************
    integer(IK) function NBS_Head(this, ind )
        implicit none
        class(NBS) this
        type(integer3), intent(in) :: ind
        
        NBS_Head = this%HeadList(ind%x,ind%y,ind%z)
                 
    end function
    
    !******************************************************************
    !   returning the head value of a cell (ind) and then adds number 
    ! of particles to the head value (as a flag)
    !******************************************************************
    integer(IK) function NBS_Head_Inc(this, ind )
        
        implicit none
        class(NBS) this
        type(integer3), intent(in)  :: ind
        
        integer(IK) nPrtcl, n

        nPrtcl = this%get_numPrtcl()
        n = this%HeadList(ind%x,ind%y,ind%z)
        if( n <= nPrtcl .and. n .ne. -1 ) & 
            this%HeadList(ind%x,ind%y,ind%z) = this%HeadList(ind%x,ind%y,ind%z) + nPrtcl
        
        NBS_Head_Inc = n

    end function
    
    !******************************************************************
    ! refreshing the Head
    !******************************************************************
    subroutine NBS_Ref_Head(this)
        
        implicit none
        class(NBS) this
        integer(IK) nPrtcl

        nPrtcl = this%get_numPrtcl()
        where( this%HeadList > nPrtcl )
            this%HeadList = this%HeadList - nPrtcl
        end where

    end subroutine
    
    !******************************************************************
    ! Next item in the list
    !******************************************************************
    integer(IK) function NBS_Next(this, n )
        implicit none
        class(NBS) this
        integer(IK) n
        
        NBS_Next = this%NextList(n)

    end function

    !******************************************************************
    ! destructor/final
    !******************************************************************
    subroutine FinilizeNBS(this)

        implicit none
        type(NBS):: this

        if( allocated(this%HeadList) ) deallocate( this%HeadList )
        if( allocated(this%NextList) ) deallocate( this%NextList )       

    end subroutine

     
    
!**************** end of NBS methods *********************!
!*********************************************************!   

    
end module g_Prtcl_NBS

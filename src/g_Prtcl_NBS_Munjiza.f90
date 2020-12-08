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
!  File name:  g_Prtcl_NBS_Munjiza.f90 
!  Module name: g_Prtcl_NBS_Munjiza    
! 
!  Purpose:
!    1) Providing a basic class for NBS_Munjiza contact search algorithms in
!       this code  
! 
!    2) It is extended based on the CellBased class and contains necessary
!       objects and methods for doing this contact search.
! 
!  External literature used:  
!     - Munjiza, A. (2004). The Combined Finite-Discrete Element Method. West
!       Sussex, UK., John Wiley & Sons Ltd.
!                     
!     - Munjiza, A. and K. R. F. Andrews (1998). "NBS contact detection 
!       algorithm for bodies of similar size. International Journal for 
!       Numerical Methods in Engineering 43: 131-149.
!  
!  Referenced modules:
!     - g_Prtcl_CellBased
!     
!------------------------------------------------------------------------------

module g_Prtcl_NBS_Munjiza
    
    use g_Prtcl_CellBased
    implicit none
    
    
    type,extends(CellBased) :: NBS_Munjiza
        
       
        integer(IK),dimension(:),   allocatable :: HeadY
        integer(IK),dimension(:),   allocatable :: NextY
        integer(IK),dimension(:),   allocatable :: HeadX
        
        integer(IK),dimension(:),   allocatable :: Headx0
        integer(IK),dimension(:),   allocatable :: NextX
        integer(IK) curr_xList_ind
        
        integer(IK),dimension(:,:), allocatable :: HeadZ  ! head list for iy, current row 
        integer(IK),dimension(:,:), allocatable :: HeadZ0 ! head list for (iy-1), lower row
        integer(IK),dimension(:),   allocatable :: NextZ
        integer(IK),dimension(0:2):: curr_zList_ind, curr_zList0_ind
        
    contains
        ! allocating NBS_Munjiza object
        procedure:: Init_NBS_Munjiza    => NBS_Munjiza_Init_NBS_Munjiza
        procedure:: Allocate_Munjiza    => NBS_Munjiza_Allocate_Munjiza
        
        ! constructing YList
        procedure:: YList               => NBS_Munjiza_YList
        
        ! constructing XList
        procedure:: XList               => NBS_Munjiza_XList
        
        ! constructing ZLists
        procedure:: ZList               => NBS_Munjiza_ZList
        procedure:: ZList0              => NBS_Munjiza_ZList0
        
        ! looping over cells marked by NBS mask
        procedure:: LoopNBSMask         => NBS_Munjiza_LoopNBSMask
        
        !performing a BroadSearch 
        procedure:: BroadSearch         => NBS_Munjiza_BroadSearch
        
        ! performing contact search (includes all steps)
        procedure:: ContactSearch       => NBS_Munjiza_ContactSearch
        
        ! destructor/final
        final::     FinilizeNBS_Munjiza
    
    endtype
    
contains

!**************************************************************!
!******************* NBS_Munjiza methods **********************!
!**************************************************************!

    !******************************************************************
    ! Initializing NBS_Munjiza object
    !******************************************************************
    subroutine NBS_Munjiza_Init_NBS_Munjiza(this,minDomain, maxDomain, max_nPrtcl, numPrtcl ,ids, flag, Bndg_box, ratio )
        implicit none
        class(NBS_Munjiza)  this
        type(real3),intent(in)  :: minDomain, maxDomain
        integer(IK),intent(in)  :: max_nPrtcl, numPrtcl
        integer(IK),pointer,dimension(:),intent(in):: ids, flag
        type(real4),pointer,dimension(:),intent(in):: Bndg_box
        real(RK),optional,intent(in)               :: ratio
        
        !// locals
        
        real(RK) l_ratio
        
        !// body
        if( present(ratio) ) then
            l_ratio = ratio
        else
            l_ratio = 1.0_RK
        endif
        
        call this%setDomain(minDomain, maxDomain)
        call this%InitSimWorld( max_nPrtcl, numPrtcl, ids, flag, Bndg_Box, .true. )
        call this%initCellBased( l_ratio )
        call this%Allocate_Munjiza()
        
    end subroutine

    !******************************************************************
    ! Allocating memory for this object
    !******************************************************************
    subroutine NBS_Munjiza_Allocate_Munjiza(this)
        implicit none
        class(NBS_Munjiza)  this
        
        
        !// locals 
        integer(IK) max_nPrtcl, iErr
        type(integer3):: numCell
        
        !// body
        numCell = this%getNumCell()
        max_nPrtcl = this%get_max_nPrtcl()
        
        
        !  Linked list allocation in Y direction 
        if( allocated( this%HeadY ) ) deallocate(this%HeadY)
        allocate( this%HeadY(0:numCell%y+1), STAT = iErr )
        if( iErr .ne. 0 )then
            call CheckForError( ErrT_Abort, "NBS_Munjiza_Allocate_Munjiza" , "Allocation of HeadY failed with size " )
            return
        end if
        
        this%HeadY = -1

        if( allocated( this%NextY ) ) deallocate(this%NextY)
        allocate( this%NextY( max_nPrtcl ) , STAT = iErr )
        if( iErr .ne. 0 )then
            call CheckForError( ErrT_Abort, "NBS_Munjiza_Allocate_Munjiza" , &
                                "Allocation of NextY failed with size " // num2str(max_nPrtcl) )
            return
        end if
        
        this%NextY = -1

        ! Linked list allocation in X direction
        if( allocated(this%HeadX ) ) deallocate(this%HeadX )
        allocate(this%HeadX(0:numCell%x+1), STAT = iErr )
        if( iErr .ne. 0 )then
            call CheckForError( ErrT_Abort, "NBS_Munjiza_Allocate_Munjiza" ,&
                                "Allocation of HeadX failed with size " // num2str(numCell%x) )
            return
        end if
        
        this%HeadX = -1
        
        if( allocated(this%HeadX0 ) ) deallocate(this%HeadX0 )
        allocate(this%HeadX0(0:numCell%x+1), STAT=iErr )
        if( iErr .ne. 0 )then
            call CheckForError( ErrT_Abort, "NBS_Munjiza_Allocate_Munjiza" ,&
                                "Allocation of HeadX0 failed with size " // num2str(numCell%x) )
            return
        end if
        this%HeadX0 = -1
        

        if( allocated(this%NextX ) ) deallocate(this%NextX )
        allocate(this%NextX( max_nPrtcl ), STAT = iErr )
        
        if( iErr .ne. 0 )then
            call CheckForError( ErrT_Abort, "NBS_Munjiza_Allocate_Munjiza" , &
                                "Allocation of NextX failed with size " // num2str(max_nPrtcl) )
            return
        end if
        this%NextX = -1
	
        ! Linked list allocation in Z direction
        if( allocated(this%HeadZ) ) deallocate(this%HeadZ )
        allocate(this%HeadZ(0:2 , 0:numCell%z+1 ), STAT = iErr )
        if( iErr .ne. 0 )then
            call CheckForError( ErrT_Abort, "NBS_Munjiza_Allocate_Munjiza" , &
                                "Allocation of Headz failed with size [3 ," // trim(num2str(numCell%z))//']' )
            return
        end if
        this%HeadZ = -1
        
        
        if( allocated(this%HeadZ0) ) deallocate(this%HeadZ0 )
        allocate(this%HeadZ0(0:2 , 0:numCell%z+1 ) , STAT = iErr )
        if( iErr .ne. 0 )then
            call CheckForError( ErrT_Abort, "NBS_Munjiza_Allocate_Munjiza" , &
                                "Allocation of Headz0 failed with size [3 ," // trim(num2str(numCell%z))//']' )
            return
        end if
        this%HeadZ0 = -1


        if( allocated(this%NextZ) ) deallocate(this%NextZ )
        allocate(this%NextZ(max_nPrtcl), STAT = iErr )
        if( iErr .ne. 0 )then
            call CheckForError( ErrT_Abort, "NBS_Munjiza_Allocate_Munjiza" , &
                                "Allocation of NextZ failed with size " // num2str(max_nPrtcl) )
            return
        end if
        this%NextZ = -1
        
        this%curr_xList_ind = -1
        this%curr_zList_ind = -1
        this%curr_zList0_ind = -1       
        
    
    end subroutine
    
    !******************************************************************
    ! performing a contact search on all particles (includes all steps)
    !******************************************************************
    subroutine NBS_Munjiza_ContactSearch(this)    
        implicit none
        class(NBS_Munjiza) this
        
        !// body 
        
        call this%set_num_cntct(0_IK)
        
        call this%BoxIndex()
        
        call this%BroadSearch()
    
    end subroutine
    
    
    !******************************************************************
    ! performing a broad search on all particles 
    !******************************************************************
    subroutine NBS_Munjiza_BroadSearch( this )
        implicit none
        class(NBS_Munjiza):: this
    

        !// local variables 
        integer(IK) ix, iy, iz
        integer(IK) m, n, l, lx
        integer(IK) mId, nId
        type(integer3) numCell
        
        !// body 
        
        !** creating the main linked list, Ylist
        call this%YList()
         
        ! if no particle is in the system, the contact search is skipped 
        if( this%get_numPrtcl() == 0 ) return
        
        numCell = this%getNumCell()

        this%HeadX0(:) = -1
        
        ! looping over all rows
        do iy = 1, numCell%y
        
            call this%XList( iy )
                    
            if( this%HeadY(iy) .ne. -1  ) then
            
                this%HeadZ(0,:) = -1
                this%HeadZ0(0,:) = -1
            
                !call zlists01( 1 ) ! icol = 1
                ! Creating the zlist of the lower row and current ix (column)
                call this%Zlist0( 1 , 1 )

                do ix = 1, numCell%x
                
                    call this%ZList( ix , 1 )
                    call this%Zlist0( ix , 2 )
                
                    if( this%Headx(ix) .ne. -1 ) then
                    
                        do iz = 1, numCell%z
                    
                           call this%LoopNBSMask(iz)
                           
                        end do
                    
                    end if
                    
                    ! same row, subs
                    this%HeadZ(0,:) = this%HeadZ(1,:)
            
                    ! lower row, subs
                    this%HeadZ0(0,:) = this%HeadZ0(1,:)
                    this%HeadZ0(1,:) = this%HeadZ0(2,:)

                end do
            
            end if
        
            this%Headx0(:) = this%Headx(:)

        end do

    end subroutine

!**********************************************************************
!   finding contacts between particles in the target cell and particles in
! cells determined by NBS mask.  
!**********************************************************************
subroutine NBS_Munjiza_LoopNBSMask( this, iz)
    implicit none
    class(NBS_Munjiza) this
    integer(IK),intent(in)  :: iz
    
    !// locals
    integer(IK) m, n, l, lx
  
    
    !// body
    
    m = this%HeadZ(1,iz)
    do while( m .ne. -1 )
        !over particles in the same cell but not the same particle (to prevent self-contact)
        n = this%NextZ(m)
        do while( n .ne. -1 )
            
            
            ! performing the narrow phase search
            call this%FineSearch(this%getMemIndx(n),this%getMemIndx(m) )
            
            ! Increasing number of contacts by one 
            call this%add_num_cntct(1_IK)
                                
            n = this%NextZ(n)
        end do

        ! over particles in (ix, iy , iz-1)
        n = this%HeadZ(1,iz-1)
        do while ( n .ne. -1 )
            
            ! performing the narrow phase search
            call this%FineSearch(this%getMemIndx(n),this%getMemIndx(m) )
            
            ! Increasing number of contacts by one
            call this%add_num_cntct(1_IK)

            n = this%NextZ(n)
        end do


        ! over particles in all cells located at (ix-1) and (iy)
        do l = -1,1
                            
            n = this%HeadZ(0,iz+l)
            do while( n .ne. -1 )
                
                ! performing the narrow phase search
                call this%FineSearch(this%getMemIndx(n),this%getMemIndx(m) )
                
                ! increasing number of contacts by one
                call this%add_num_cntct(1_IK)

                n = this%NextZ(n)
            end do

        end do
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
                            
            do l=-1,1
                n = this%HeadZ0(lx,iz+l)
                do while (n .ne. -1)
                    
                    ! performing narrow phase search
                    call this%FineSearch(this%getMemIndx(n),this%getMemIndx(m) )
                    
                    ! Increasing number of contacts by one
                    call this%add_num_cntct(1_IK)

                    n = this%NextZ(n)
                end do

            end do

        end do


        m = this%NextZ(m)
    end do

    
end subroutine 
    
    !******************************************************************
    ! constructing Ylist of particles
    !******************************************************************
    subroutine NBS_Munjiza_YList(this)
        implicit none
        class(NBS_Munjiza) this
        
        
        ! local variables
	    integer(IK) nPrtcl,n, iy
	    type(integer3):: box
	
        ! current number of particles
        nPrtcl = this%get_numPrtcl()
        
        ! nullifying list Y
	    this%HeadY = -1
        
	    do n = 1, nPrtcl
            
            ! if sphere is in the simulation domain
            if( this%bndg_flag( this%getMemIndx(n) ) >= Pflg_inDomain ) then
                box = this%getIndex(n)  ! getting the integer coordinates of sphere
                iy  =    box%y          ! integer y-coordinate  
                this%NextY( n ) = this%HeadY(iy)
                this%HeadY( iy ) = n
            end if
            
        end do
        
    end subroutine
    
    !******************************************************************
    ! constructing Xlist of row iy 
    !******************************************************************
    subroutine NBS_Munjiza_XList( this , iy , lcheck)
	    implicit none
        class( NBS_Munjiza )          this
	    integer(IK),intent(in)::      iy ! row index
        logical,optional,intent(in):: lcheck  ! checking if the current list at column = iy is constructed or not
     
	    !// locals
	    integer(IK) n,  ix
	    type(integer3):: box
	    
        !// body
        
	 ! first checking if the xlist of row iy has been constructed previously 
        if( present(lcheck) ) then
            if( lcheck) then
                if (iy == this%curr_xList_ind) return
            endif
        endif

        ! nullifying the xlist of current row but keeps the previous raw
	    this%HeadX(:) = -1
        this%curr_xList_ind = iy
        
	    n = this%Heady(iy)
	
	    do while ( n .ne. -1 )
    
            box = this%getIndex(n)
		    ix = box%x
			
		    this%NextX( n ) = this%HeadX(ix)
		    this%HeadX(ix) = n	

		    n = this%Nexty(n)

	    end do
	
    end subroutine 

    !*********************************************************************
    !   Constructing the ZList of column ix, ix-1, or ix+1 depending on the 
    ! value of l
    !*********************************************************************
    subroutine NBS_Munjiza_ZList( this , ix , l , lcheck )
        implicit none
        class(NBS_Munjiza ) this
	    integer(IK),intent(in) :: ix ! col index
        integer(IK),intent(in) :: l  ! the column location with resect to ix 
        ! 2: next (right) column
        ! 1: current column
        ! 0: previous (left) column   
        logical,optional,   intent(in):: lcheck
    

        !// locals
        integer(IK) n,  iz
        type(integer3) box
	    
        !// body
        
        if( present(lcheck) ) then
            if(lcheck) then
                if (ix == this%curr_zList_ind(l) ) return
            end if
        end if
        
	    ! nullifying the zlist of current col ix, but keeps the previous and next cols
        ! 2: next (right) column
        ! 1: current column
        ! 0: previous (left) column
	    this%HeadZ(l,:) = -1
	    this%curr_ZList_ind(l) = ix
        
	    n = this%HeadX(ix+l-1) ! reading from current row iy
	    
	    do while ( n .ne. -1 )
		
            box = this%getIndex(n)
	    iz = box%z
	    this%NextZ(n) = this%HeadZ(l,iz)
            this%HeadZ(l,iz) = n

            n = this%NextX(n)
                
	    end do
        
    
    end subroutine NBS_Munjiza_ZList
    
    !*********************************************************************
    !   Constructing the ZList0 of column ix, ix-1, or ix+1 depending on the 
    ! value of l
    !*********************************************************************
    subroutine NBS_Munjiza_ZList0( this , ix , l, lcheck )
        implicit none
        class(NBS_Munjiza ):: this
	    integer(IK),intent(in):: ix ! column index
        integer(IK),intent(in):: l  ! the column location with respect to ix 
        ! 2: next (right) column
        ! 1: current column
        ! 0: previous (left) column   
        logical,optional,intent(in):: lcheck

        !// locals
        integer(IK) n,  iz
        type(integer3) box
	
	    if( present(lcheck) ) then
            if(lcheck) then
                if (ix == this%curr_zList0_ind(l) ) return
            end if
        end if
        
        ! 2: next (right) column
        ! 1: current column
        ! 0: previous (left) column
	    this%HeadZ0(l,:) = -1
	    this%curr_ZList0_ind(l) = ix
	    n = this%HeadX0(ix+l-1) ! reading from the row below iy (or iy-1)
	    if( n .eq. -1 ) return

	    do while ( n .ne. -1 )
		
            box = this%getIndex(n)
		    iz = box%z
		    this%NextZ(n) = this%HeadZ0(l,iz)
            this%HeadZ0(l,iz) = n

            n = this%NextX(n)
                
	    end do
        
    end subroutine NBS_Munjiza_ZList0
   
    !******************************************************************
    ! final/destructor
    !******************************************************************
    subroutine FinilizeNBS_Munjiza(this)
        implicit none
        type(NBS_Munjiza) this
        
        if( allocated( this%HeadY ) ) deallocate(this%HeadY)
        if( allocated( this%NextY ) ) deallocate(this%NextY)
        
        ! Linked list allocation in X direction
        if( allocated(this%HeadX ) ) deallocate(this%HeadX )
        if( allocated(this%HeadX0 ) ) deallocate(this%HeadX0 )
        if( allocated(this%NextX ) ) deallocate(this%NextX )
        
        
        ! Linked list allocation in Z direction
        if( allocated(this%HeadZ) ) deallocate(this%HeadZ )
        if( allocated(this%HeadZ0) ) deallocate(this%HeadZ0 )
        if( allocated(this%NextZ) ) deallocate(this%NextZ )
        
    end subroutine
    
end module

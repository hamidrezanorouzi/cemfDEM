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
!  File name:  g_Prtcl_Hrchl_Munjiza.f90 
!  Module name: g_Prtcl_Hrchl_Munjiza 
! 
!  Purpose: 
!   1) to Provide a basic classes for hierarchical contact search algorithms 
!      in this code (Munjiza).
!                            
!   2) It is extended based on the NBS_Munjiza class and contains necessary
!      objects and methods for doing hierarchical contact search in one level.
! 
! 
!   External literature used: 
!      - Peters, J. F., et al. (2009). "A hierarchical search algorithm for
!        discrete element method of greatly differing particle sizes."
!        Engineering Computations 26(6): 621-634. 
! 
!      - Munjiza, A. (2004). The Combined Finite-Discrete Element Method.
!        West Sussex, UK., John Wiley & Sons Ltd.
!        
!      - Munjiza, A. and K. R. F. Andrews (1998). "NBS contact detection 
!        algorithm for bodies of similar size." International Journal for
!        Numerical Methods in Engineering 43: 131-149.
! 
!  Referenced modules: 
!      - g_Prtcl_NBS_Munjiza
! 
!------------------------------------------------------------------------------

module g_Prtcl_Hrchl_Munjiza
    
    use g_TypeDef
    use g_error_handling
    use g_Prtcl_NBS_Munjiza
    
    implicit none
    
    private
    public:: lvl_Munjiza, NBS_Munjiza_Hrchl
    
    
    type, extends(NBS_Munjiza) :: lvl_Munjiza
        
        integer(IK)  lvl                ! level number
        integer(IK)::lvl_max_nPrtcl = 0 ! maximum number of particles in this level
        integer(IK)::lvl_numPrtcl = 0   ! number of current particles in this level
        integer(IK)::lvl_num_cnsv_cntct = 0 ! number of conservative contact in this level
        real(RK)     lvl_minD           ! the minimum diameter of bonnding boxes in this level
        real(RK)     lvl_maxD           ! the maximum diameter of bounding boxes in this level
        
        integer(IK),allocatable,dimension(:):: memIndx ! memory index of bounding boxes in this level
        
        integer(IK),dimension(:),   allocatable :: Headx2
        integer(IK) curr_xList2_ind
        integer(IK),dimension(:,:), allocatable :: HeadZ2
        integer(IK),dimension(0:2):: curr_zList2_ind
        
        
    contains
        ! initializing the object of this level
        procedure:: Init_lvl_Munjiza    => lvl_Munjiza_Init_lvl_Munjiza
        
        ! setting and getting the minimum and maximum diameter of bounding box in this level
        procedure:: set_lvl_minmax  => lvl_Munjiza_set_lvl_minmax
        procedure:: get_lvl_minmax  => lvl_Munjiza_get_lvl_minmax
        
        ! overriding the original function defined in SimWorld class, it returns
        ! the minimum and maximum diameter in this level
        procedure:: get_minmax      => lvl_Munjiza_get_lvl_minmax ! function overriding     
        
        ! overriding the original function defined in SimWorld class, it returns 
        ! number of bounding boxes in this level
        procedure:: get_max_nPrtcl  => lvl_Munjiza_get_max_nPrtcl      ! function overriding 
        procedure:: get_numPrtcl    => lvl_Munjiza_get_numPrtcl        ! function overriding
        
        ! overriding the original function defined in SimWorld class, it returns
        ! the position of ith boxe in this level
        procedure:: getBndgBox_i    => lvl_Munjiza_getBndgBox_i   ! function overriding 
        
        ! overriding the original function defined in SimWorld class, it returns
        ! the memory index of ith box in this level
        procedure:: getMemIndx      => lvl_Munjiza_getMemIndx     ! function overriding 
        
        ! finding the corresponding index of a cell of the current level in another level 
        procedure:: mapIndex        => lvl_Munjiza_mapIndex
        
        ! calculating number of particles which exist in this level 
        procedure:: clc_lvl_numPrtcls => lvl_Munjiza_clc_lvl_numPrtcls
        
        ! calculating memory indices of bounding boxes which exist in this level 
        procedure:: clc_lvl_memIndx => lvl_NBS_clc_lvl_memIndx
        
        ! constructing extra zlists to map particles aroung 26 cells around the target cell 
        procedure:: ZList2          => lvl_Munjiza_ZList2
        procedure:: XList2          => lvl_Munjiza_XList2
        
        ! managing swaping linked lists 
        procedure:: SwapXLists      => lvl_Munjiza_SwapXLists
        procedure:: SwapZLists      => lvl_Munjiza_SwapZLists
        
        ! performing a cross-level contact search
        procedure:: Loop_CrossMask  => lvl_Munjiza_Loop_CrossMask
        
        ! number of contacts related to the cross level search         
        procedure:: set_crs_num_cntct   => lvl_Munjiza_set_num_cntct
        procedure:: get_crs_num_cntct   => lvl_Munjiza_get_num_cntct
        procedure:: add_crs_num_cntct   => lvl_Munjiza_add_num_cntct
        
        ! number of contacts related to the level (cross level + same level)
        procedure:: get_lvl_num_cntct   => lvl_Munjiza_get_lvl_num_cntct
        procedure:: set_lvl_num_cntct   => lvl_Munjiza_set_lvl_num_cntct
        
        ! freeing memory allocated for this object
        final :: Finalize_lvl_Munjiza
    end type
    
    
    type NBS_Munjiza_Hrchl
        
        integer,allocatable:: dummy
        integer(IK)::num_lvls = 1                                 ! number of levels
        type(lvl_Munjiza),allocatable,dimension(:):: lvls_Munjiza ! level
        
    contains
        
        !initializing the hierarchical search algorithm object
        procedure   :: Init_Munjiza_Hrchl=> NBS_Munjiza_Hrchl_Init_Munjiz_Hrchl 
        
        !setting the contact list (particle-particle )
        procedure   :: setCont_List     => NBS_Munjiza_Hrchl_setCont_List
        
        ! One level contact search (intra-level and cross level with next levels)
        procedure   :: OneLevelBroadSearch    =>  NBS_Munjiza_Hrchl_OneLevelBroadSearch
        
        ! performing a full contact search 
        procedure   :: ContactSearch    => NBS_Munjiza_Hrchl_ContactSearch
        
        ! methods related to number of conservative contacts 
        procedure   :: get_num_cntct => NBS_Munjiza_Hrchl_get_num_cntct
        procedure   :: set_num_cntct => NBS_Munjiza_Hrchl_set_num_cntct
        
        
        procedure   :: set_numPrtcl  => NBS_Munjiza_Hrchl_set_numPrtcl
        
	! freeing the memory allocated for this object
        final       :: Finalize_NBS_Munjiza_Hrchl

    end type
    
    
contains

!**********************************************************************
! initializing the object of this level
!**********************************************************************
subroutine lvl_Munjiza_Init_lvl_Munjiza(this, lvl )
    implicit none
    class(lvl_Munjiza) this
    integer(IK),intent(in):: lvl
    
    !// locals
    integer(IK) numBox, iErr
    type(integer3):: numCell
    
    !// body
    
    this%lvl = lvl
    
    !calculating number of particles which belongs to this level
    call this%clc_lvl_numPrtcls()
    
    ! allocating memory to keep particle ids of this level
    if(allocated(this%memIndx))deallocate(this%memIndx)
    allocate( this%memIndx(this%lvl_max_nPrtcl) )
    call this%clc_lvl_memIndx()
    
    ! number of cells
    numCell = this%getNumCell()
    
    if( allocated(this%HeadX2 ) ) deallocate(this%HeadX2 )
    allocate(this%HeadX2(0:numCell%x+1), STAT = iErr )
    if( iErr .ne. 0 )then
        call CheckForError( ErrT_Abort, "lvl_Munjiza_Init_lvl_Munjiza", &
                            "Allocation of HeadX2 failed with size" // num2str(numCell%x) )
          
        return
    end if
    
    this%HeadX2 = -1
    
    if( allocated(this%HeadZ2) ) deallocate(this%HeadZ2 )
    allocate(this%HeadZ2(0:2 , 0:numCell%z+1 ) , STAT = iErr )
    if( iErr .ne. 0 )then
        call CheckForError( ErrT_Abort, "lvl_Munjiza_Init_lvl_Munjiza", &
                            "Allocation of Headz2 failed with size [3 ," // trim(num2str(numCell%z))//']' )
        
        return
    end if
    
    this%HeadZ2 = -1
    
    this%curr_xList2_ind = -1
    this%curr_zList2_ind = -1
    
end subroutine 


!**********************************************************************
! setting the minimum and maximum diamter of bounding boex in this level
!**********************************************************************
subroutine lvl_Munjiza_set_lvl_minmax(this, minD , maxD )
    implicit none
    class(lvl_Munjiza):: this
    real(RK),intent(in) :: minD, maxD

    this%lvl_minD = minD
    this%lvl_maxD = maxD

end subroutine 

!**********************************************************************
! returning the minimum and maximum diamter of bounding boex in this level
!**********************************************************************
function lvl_Munjiza_get_lvl_minmax(this) result (res)
    implicit none
    class(lvl_Munjiza):: this
    real(RK),dimension(2):: res

    res(1) = this%lvl_minD
    res(2) = this%lvl_maxD

end function

!**********************************************************************
! overriding the original function defined in SimWorld class, it returns 
! maximum number of bounding boxes in this level
!**********************************************************************
integer(IK) function lvl_Munjiza_get_max_nPrtcl(this)
    implicit none
    class(lvl_Munjiza):: this

    lvl_Munjiza_get_max_nPrtcl = this%lvl_max_nPrtcl

end function

!**********************************************************************
! overriding the original function defined in SimWorld class, it returns 
! number of bounding boxes in this level
!**********************************************************************
integer(IK) function lvl_Munjiza_get_numPrtcl(this)
    implicit none
    class(lvl_Munjiza):: this

    lvl_Munjiza_get_numPrtcl = this%lvl_numPrtcl

end function

!**********************************************************************
! overriding the original function defined in SimWorld class, it returns
! the position of ith boxe in this level
!**********************************************************************
type(real4) function lvl_Munjiza_getBndgBox_i(this, n)
    implicit none
    class(lvl_Munjiza):: this
    integer(IK), intent(in) :: n

    lvl_Munjiza_getBndgBox_i = this%bndg_box( this%memIndx(n) ) 

end function

!*************************************************************************
! overriding the original function defined in SimWorld class, it returns
! the memory index of ith box in this level
!*************************************************************************
integer(IK) function lvl_Munjiza_getMemIndx(this, n)
    implicit none
    class(lvl_Munjiza):: this
    integer(IK), intent(in):: n

    lvl_Munjiza_getMemIndx =  this%memIndx(n)

end function

!**********************************************************************
!   finding the corresponding index of a cell of the current level in
! another level
!**********************************************************************
type(integer3) function lvl_Munjiza_mapIndex( this, ind ,crs_lvl )
    implicit none
    class(lvl_Munjiza)  this
    type(integer3), intent(in)  :: ind
    integer(IK), intent(in)     :: crs_lvl
    
    !//locals
    integer(IK) a, lvl

    !//body
        
    a = 2**( this%lvl-crs_lvl )
    
    lvl_Munjiza_mapIndex%x = (ind%x -1)/a + 1
    lvl_Munjiza_mapIndex%y = (ind%y -1)/a + 1
    lvl_Munjiza_mapIndex%z = (ind%z -1)/a + 1

end function

!**********************************************************************
! Calculating number of particles which exist in this level. 
! Each level has a valid size range. If the size of particle falls in 
! this range, that particle is assigned to this level.
!*********************************************************************   
subroutine lvl_Munjiza_clc_lvl_numPrtcls(this)
    implicit none
    class(lvl_Munjiza):: this
    
    !// locals
    integer(IK) i 
    integer(IK) lvl_max_nPrtcl ! maximum number of particles in this level
    integer(IK) lvl_numPrtcl   ! current number of particles in this level
    integer(IK) tot_Prtcls ! maximum number of particles in the entire simulation / all levels
    integer(IK) numPrtcl  ! current number of particles in all levels
    real(RK)::  minD, maxD, D(2)
    type(real4):: bndg_box

    ! calls the original function
    tot_Prtcls = this%NBS_Munjiza%CellBased%get_max_nPrtcl()
    numPrtcl = this%NBS_Munjiza%CellBased%get_numPrtcl()
        
    ! size range of the level        
    D = this%get_lvl_minmax()
    MinD = D(1)
    MaxD = D(2)
        
    lvl_max_nPrtcl = 0 
    lvl_numPrtcl = 0
        
    do i = 1, tot_Prtcls
            
        ! invokes the original function
        bndg_box = this%NBS_Munjiza%CellBased%getBndgBox_i (i)
            
        if( bndg_box%w > minD .and. bndg_box%w <= maxD  ) then
                
            lvl_max_nPrtcl = lvl_max_nPrtcl + 1
            if ( i .le. numPrtcl ) lvl_numPrtcl = lvl_numPrtcl + 1
            
        end if
    end do

    this%lvl_max_nPrtcl = lvl_max_nPrtcl
    this%lvl_numPrtcl = lvl_numPrtcl

end subroutine

!**********************************************************************
! Calculating memory indices of bounding boxes which exist in this level. 
! Each level has a valid size range. If the size of particle falls in 
! this range, that particle is assigned to this level and its memory
! index is stored in this level.
!*********************************************************************   
subroutine lvl_NBS_clc_lvl_memIndx(this)
    implicit none
    class(lvl_Munjiza) this
    
    !// locals
    integer(IK):: i, nPar
    integer(IK):: tot_nPrtcl ! total number of particles in all levels 
    real(RK)::  minD, maxD, D(2)
    type(real4):: bndg_box

    ! // body
    
    !total number of particles in simulation
    tot_nPrtcl= this%NBS_Munjiza%CellBased%get_max_nPrtcl()
    
    ! size range of the level       
    D = this%get_lvl_minmax()
    MinD = D(1)
    MaxD = D(2)

    nPar = 0 
    do i = 1, tot_nPrtcl
            
        bndg_box = this%NBS_Munjiza%CellBased%getBndgBox_i (i)
        if( bndg_box%w > minD .and. bndg_box%w <= maxD  ) then
            nPar = nPar + 1
            this%memIndx(nPar) = i
        end if
    end do

end subroutine

!**********************************************************************
! constructing Xlist of row iy
!**********************************************************************
subroutine lvl_Munjiza_XList2( this , iy , lcheck)
	    implicit none
        class( lvl_Munjiza )     this
	    integer(IK),intent(in):: iy ! row index
        logical,optional,intent(in):: lcheck  ! checking if the current list at column = iy is constructed or not
     
    
	    ! // locals
	    integer(IK) n,  ix
	    type(integer3):: box
	
	    
	    ! first checking if the xlist of row iy has been constructed previously 
        if( present(lcheck) ) then
            if( lcheck) then
                if (iy == this%curr_xList2_ind) return
            endif
        endif
        
        
        ! nullifying the xlist of current row but keeping the previous raw
	    this%HeadX2(:) = -1
        this%curr_xList2_ind = iy
	    n = this%Heady(iy)
	
	    if( n .eq. -1 ) return

	    do while ( n .ne. -1 )
    
            box = this%getIndex(n)
		    ix = box%x
			
		    this%NextX( n ) = this%HeadX2(ix)
		    this%HeadX2(ix) = n	

		    n = this%Nexty(n)

	    end do
	    
end subroutine

!**********************************************************************
!   constructing Zlists of columns ix-1, ix, or ix+1 depending on the
! value of l in the input argument. 
!*********************************************************************
subroutine lvl_Munjiza_ZList2( this , ix , l , lcheck )
    implicit none
    class(lvl_Munjiza ):: this
	integer(IK) ix ! col index
    integer(IK) l  ! the column location with respect to ix 
    ! 2: next (right) column
    ! 1: current column
    ! 0: previous (left) column   
    logical,optional,   intent(in):: lcheck
    

    !// locals
    integer(IK) n,  iz
    type(integer3) box
	
    !//body
    
    if( present(lcheck) ) then
        if(lcheck) then
            if (ix == this%curr_zList2_ind(l) ) return
        end if
    end if
        
        
	! nullifying the zlist of current col ix, but keeping the previous and next cols
    ! 2: next (right) column
    ! 1: current column
    ! 0: previous (left) column
	this%HeadZ2(l,:) = -1
	this%curr_ZList2_ind(l) = ix
    
	n = this%HeadX2(ix+l-1) ! reads from current row iy
	if( n .eq. -1 ) return

	do while ( n .ne. -1 )
		
            box = this%getIndex(n)
            iz = box%z
	    this%NextZ(n) = this%HeadZ2(l,iz)
            this%HeadZ2(l,iz) = n

            n = this%NextX(n)
                
	end do
    
end subroutine 

!******************************************************************************
!   Performing a cross-level contact search between the current level (this)
! and CrossMunjiza level. A cross level contact search performed between 
! particles from target cell in current level and 27 neighbor cells from
! CrossMunjiza level
!******************************************************************************
subroutine lvl_Munjiza_Loop_CrossMask( this , iz, crs_iz , CrossMunjiza)
    implicit none
    class(lvl_Munjiza) this
    integer(IK),        intent(in)  :: iz, crs_iz
    class(lvl_Munjiza), intent(in)  :: CrossMunjiza
    
    !// locals
    integer(IK) m, n
    integer(IK) lx,  l
    
    !// body    
    m = this%HeadZ(1,iz)
    
    
    do while( m .ne. -1 )
        
         ! first, looping all 9 cells of cross level located in the row below of the target cell
         do lx = 0,2
             do l=-1,1
                 n = CrossMunjiza%HeadZ0(lx,crs_iz+l)
                 do while( n .ne. -1)
                     ! performing the narrow phase search
                     call this%FineSearch(this%getMemIndx(m),CrossMunjiza%getMemIndx(n) )
                     call this%add_crs_num_cntct(1_IK)
                     
                     n = CrossMunjiza%NextZ(n)
                 end do
                 
             end do
         end do
         
         ! second, looping all 9 cells of cross level located in the same row as the target cell
         do lx = 0,2
             do l=-1,1
                 n = CrossMunjiza%HeadZ(lx,crs_iz+l)
                 do while( n .ne. -1)
                     ! performing the narrow phase search
                     call this%FineSearch(this%getMemIndx(m),CrossMunjiza%getMemIndx(n) )
                     call this%add_crs_num_cntct(1_IK)
                     
                     n = CrossMunjiza%NextZ(n)
                 end do
                 
             end do
         end do
         
         ! third, looping all 9 cells of cross level located in the row above the target cell
         do lx = 0,2
             do l=-1,1
                 n = CrossMunjiza%HeadZ2(lx,crs_iz+l)
                 do while( n .ne. -1)
                     ! performing the narrow phase search
                    call this%FineSearch(this%getMemIndx(m),CrossMunjiza%getMemIndx(n) )
                    call this%add_crs_num_cntct(1_IK)
                                    
                     n = CrossMunjiza%NextZ(n)
                 end do
                 
             end do
         end do
         
         m = this%NextZ(m)
    end do
    
end subroutine 

!**********************************************************************
! managing swapping linked lists 
!**********************************************************************
subroutine lvl_Munjiza_SwapZLists( this )
    implicit none
    class(lvl_Munjiza) this
    
    this%HeadZ0(0,:) = this%HeadZ0(1,:)
    this%HeadZ0(1,:) = this%HeadZ0(2,:)
    this%curr_ZList0_ind(0) = this%curr_ZList0_ind(1)
    this%curr_ZList0_ind(1) = this%curr_ZList0_ind(2)
    
    this%HeadZ(0,:) = this%HeadZ(1,:)
    this%HeadZ(1,:) = this%HeadZ(2,:)
    this%curr_ZList_ind(0) = this%curr_ZList_ind(1)
    this%curr_ZList_ind(1) = this%curr_ZList_ind(2)
    
    this%HeadZ2(0,:) = this%HeadZ2(1,:)
    this%HeadZ2(1,:) = this%HeadZ2(2,:)
    this%curr_ZList2_ind(0) = this%curr_ZList2_ind(1)
    this%curr_ZList2_ind(1) = this%curr_ZList2_ind(2)
      
end subroutine 

subroutine lvl_Munjiza_SwapXLists( this )
    implicit none
    class(lvl_Munjiza) this
    
    this%HeadX0 = this%HeadX
    !this%curr_xList0_ind = this%curr_xList_ind
    this%HeadX = this%HeadX2
    this%curr_xList_ind = this%curr_xList2_ind
end subroutine 


!**********************************************************************
!   setting number of conservative contacts in this level
!**********************************************************************
subroutine lvl_Munjiza_set_num_cntct(this , num )
    implicit none
    class(lvl_Munjiza) this
    integer(IK) num
        
    this%lvl_num_cnsv_cntct = num
        
end subroutine 

!**********************************************************************
!   returning number of conservative contacts in cross level checks
!**********************************************************************
integer(IK) function lvl_Munjiza_get_num_cntct(this)
    implicit none
    class(lvl_Munjiza) this
        
    lvl_Munjiza_get_num_cntct = this%lvl_num_cnsv_cntct
        
end function

!**********************************************************************
!   adding num to the number of conservative contacts in this level
!**********************************************************************
subroutine lvl_Munjiza_add_num_cntct(this , num )
    implicit none
    class(lvl_Munjiza) this
    integer(IK) num
        
    this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + num
    
end subroutine

!**********************************************************************
!   returning number of conservative contacts in this level (intra-level and cross level)
!**********************************************************************
function lvl_Munjiza_get_lvl_num_cntct(this) result(res)     
    implicit none
    class(lvl_Munjiza):: this
    integer(IK):: res(2)
        
    res(1) = this%NBS_Munjiza%get_num_cntct()
    res(2) = this%get_crs_num_cntct()
    
end function

!**********************************************************************
!   setting number of conservative contacts in this level (intra-level and cross level)
!**********************************************************************
subroutine lvl_Munjiza_set_lvl_num_cntct(this, num) 
    implicit none
    class(lvl_Munjiza):: this
    integer(IK):: num(2)
        
    call this%NBS_Munjiza%set_num_cntct( num(1) )
    call this%set_crs_num_cntct( num(2) )
    
end subroutine

!**********************************************************************
! final 
!**********************************************************************
subroutine Finalize_lvl_Munjiza(this)
    implicit none
    type(lvl_Munjiza) this
    
    if(allocated(this%memIndx))deallocate(this%memIndx)
    if( allocated(this%HeadX2 ) ) deallocate(this%HeadX2 )
    if( allocated(this%HeadZ2) ) deallocate(this%HeadZ2 )
    
end subroutine




!/////////////////////////////////NBS_Munjiza_Hrchl/////////////////////////////////////////////

subroutine NBS_Munjiza_Hrchl_Init_Munjiz_Hrchl (this , minDomain, maxDomain, max_nPrtcl, numPrtcl ,ids, flag, bndg_box,num_lvls )
        implicit none
        class(NBS_Munjiza_Hrchl) this
        type(real3),                     intent(in) :: minDomain, maxDomain
        integer(IK),                     intent(in) :: max_nPrtcl, numPrtcl
        type(real4),pointer,dimension(:),intent(in) :: bndg_box
        integer(IK),pointer,dimension(:),intent(in) :: ids, flag
        integer(IK),optional,            intent(in) :: num_lvls
    
        ! // locals
        integer(IK) numLevels, lvl
        real(RK)    minD, maxD, lvl_minD, lvl_maxD
        real(RK)::  dmax_min, D(2) 
        type(CellBased):: temp_simWorld

        ! //body
        ! finding the minimum and maximum diamter of bounding boxes
        call temp_simWorld%InitSimWorld( max_nPrtcl, numPrtcl ,ids, flag, bndg_box , .true. )
        D = temp_simWorld%get_minmax();
        minD = D(1)
        maxD = D(2)
        
        ! the ratio of minimum and maximum determines number of required levels 
        dmax_min = maxD/minD
        
        
        if( .not. present( num_lvls)) then
        
            select case( int(dmax_min) )
                case(0:1)
                    numLevels = 1
                case(2:3)
                    numLevels = 2
                case(4:7)
                    numLevels = 3
                case(8:15)
                    numLevels = 4
                case(16:31)
                    numLevels = 5
                case default
                    numLevels = 6
            end select

        else
            numLevels = num_lvls            
        end if

        

        this%num_lvls = numLevels
        ! allocating objects for each levels 
        if( allocated(this%lvls_Munjiza) ) deallocate(this%lvls_Munjiza)
        allocate( this%lvls_Munjiza(this%num_lvls) )
	
	
        ! setting the maximum diameter of the first level equal to maximum diameter of bounding boxes
        lvl_maxD = maxD
        
        !>>> log info
        call MainLogInfo%OutInfo(" Number of levels is :" // trim( num2str(this%num_lvls)) ,3 )
        !<<< log info
        
        do lvl = 1 ,this%num_lvls
        
        
            call this%lvls_Munjiza(lvl)%setDomain( minDomain, maxDomain) 
            call this%lvls_Munjiza(lvl)%InitSimWorld( max_nPrtcl, numPrtcl ,ids, flag, bndg_box , .true. )
            
            ! the minimum diameter of the level is half of the maximum diameter
            lvl_minD = lvl_maxD/2.0_RK 
            
            ! modifying the minimum diameter of the last level and sets it to a very small value
            if( lvl == this%num_lvls) lvl_minD = epsilon(lvl_minD) 
            
            call this%lvls_Munjiza(lvl)%set_lvl_minmax( lvl_minD, lvl_maxD )
            call this%lvls_Munjiza(lvl)%initCellBased()
            call this%lvls_Munjiza(lvl)%Init_lvl_Munjiza( lvl  )
            call this%lvls_Munjiza(lvl)%initCellBased()
            call this%lvls_Munjiza(lvl)%Allocate_Munjiza()
            
            
            !>>> log file
            call MainLogInfo%OutInfo("Level "//trim( num2str(lvl)), 4) 
            call MainLogInfo%OutInfo("   cell size is [m] :"// &
                                    trim( num2str(this%lvls_Munjiza(lvl)%getCellSize() ) ) , 4, .true.)
            call MainLogInfo%OutInfo("   number of cells (x,y,z) :"// &
                                    trim( num2str( this%lvls_Munjiza(lvl)%getNumCell()) ) , 4, .true.)
            call MainLogInfo%OutInfo("   number of particles : "// &
                                    trim( num2str(this%lvls_Munjiza(lvl)%get_max_nPrtcl() ) ) , 4, .true.)
            !<<< log file
            
            ! halving the maximum diameter of the level to be the maximum diameter for the next 
            ! level
            lvl_maxD = lvl_maxD/2.0_RK
       
        end do

end subroutine

!**********************************************************************
! setting the contact list (particle-particle )    
!**********************************************************************
subroutine NBS_Munjiza_Hrchl_setCont_List(this , PP_CL )
    implicit none 
    class(NBS_Munjiza_Hrchl) this
    class(base_ContactList),pointer, intent(in):: PP_CL
        
    integer(IK) i
        
    do i = 1, this%num_lvls
            
        call this%lvls_Munjiza(i)%setCont_List(PP_CL)
            
    end do
                
end subroutine 

!**********************************************************************
! a full contact search for all levels
!*********************************************************************
subroutine NBS_Munjiza_Hrchl_ContactSearch(this)
    implicit none
    class(NBS_Munjiza_Hrchl) this
    
    integer(IK) lvl, numlvls
    
    numlvls = this%num_lvls
    
    !Initializing number of conservative contacts 
    call this%set_num_cntct( (/0_IK, 0_IK/) )
    
    ! grid index and Ylists of all levels
    do lvl = numlvls, 1, -1
        call this%lvls_Munjiza(lvl)%BoxIndex()
        call this%lvls_Munjiza(lvl)%YList()
    end do
    
    ! contact search of all levels
    do lvl = numlvls, 1, -1
        call this%OneLevelBroadSearch( this%lvls_Munjiza(lvl) )            
    end do
    
end subroutine 

!******************************************************************************
!One level contact search (intra-level and cross level with next levels)
!******************************************************************************
subroutine NBS_Munjiza_Hrchl_OneLevelBroadSearch( this , base_lvl )
    implicit none
    class(NBS_Munjiza_Hrchl) this
    class(lvl_Munjiza) base_lvl
    
    !// locals
    integer(IK) lvl, lvl_dif
    integer(IK) ix, iy, iz
    type(integer3) numCell, ind, crs_ind
    
    !// body
    
    
    ! if number of particles in the base level is zero, contact search for this level is skipped
    if( base_lvl%get_numPrtcl() == 0 ) return
    
    numCell = base_lvl%getNumCell()
    
    
    ! same level
    
    base_lvl%HeadX0(:) = -1
    
    ! starting the loop over all rows in base level
    do iy = 1, numCell%y
        
        !
        call base_lvl%XList(iy, .true. )
        
        ! constructing the xLists of lvl at above rows
        do lvl =base_lvl%lvl-1, 1, -1
            
            crs_ind = base_lvl%mapIndex( integer3(1,iy,1) , lvl)
            
            if( crs_ind%y == 1 ) then
                
                ! for the first row, the xLists of below and top rows should be constructed
                this%lvls_Munjiza(lvl)%HeadX0(:) = -1
                call this%lvls_Munjiza(lvl)%xList(1, .true.)
                
            end if                      
            
            call this%lvls_Munjiza(lvl)%xList2(crs_ind%y+1, .true.)
            
        end do
                        
        
        ! if row is non-empty
        if( base_lvl%HeadY(iy) .ne. -1 ) then
            
            base_lvl%HeadZ(0,:) = -1
            base_lvl%HeadZ0(0,:) = -1
            
            ! Creating the zlist of the lower row and current ix (column)
            call base_lvl%Zlist0( 1 , 1 )
            
            do ix = 1, numCell%x
               
                call base_lvl%ZList( ix , 1 )
                call base_lvl%Zlist0( ix , 2 )
                
                
                ! constructing the zLists of the lvl at the right column
                do lvl =base_lvl%lvl-1, 1, -1
                    
                    crs_ind = base_lvl%mapIndex( integer3(ix,iy,1) , lvl)
                    
                    if( crs_ind%x == 1 ) then
                        
                        ! for the first column, the zLists of current and left columns should be constructed
                        this%lvls_Munjiza(lvl)%HeadZ0(0,:) = -1
                        this%lvls_Munjiza(lvl)%HeadZ(0,:) = -1
                        this%lvls_Munjiza(lvl)%HeadZ2(0,:) = -1
                        call this%lvls_Munjiza(lvl)%zList0(crs_ind%x,1, .true. )    
                        call this%lvls_Munjiza(lvl)%zList(crs_ind%x,1, .true.  )    
                        call this%lvls_Munjiza(lvl)%zList2(crs_ind%x,1, .true. )    
                        
                    end if
                        
                    call this%lvls_Munjiza(lvl)%zList0(crs_ind%x,2, .true. )    
                    call this%lvls_Munjiza(lvl)%zList(crs_ind%x,2, .true.  )    
                    call this%lvls_Munjiza(lvl)%zList2(crs_ind%x,2, .true. )  
                    
                end do                        
                              
                
                if( base_lvl%Headx(ix) .ne. -1 ) then
                    
                    do iz = 1, numCell%z
                        ! same level NBS mask check
                        call base_lvl%LoopNBSMask(iz)
                        
                        do lvl = base_lvl%lvl-1, 1, -1
                        
                            crs_ind = base_lvl%mapIndex( integer3(ix,iy,iz) , lvl)
                            ! cross level check
                            call base_lvl%Loop_CrossMask( iz, crs_ind%z , this%lvls_Munjiza(lvl) )
                        
                        end do
                        
                    end do
                    
                end if
                
                ! same row, subs
                base_lvl%HeadZ(0,:) = base_lvl%HeadZ(1,:)
            
                ! lower row, subs
                base_lvl%HeadZ0(0,:) = base_lvl%HeadZ0(1,:)
                base_lvl%HeadZ0(1,:) = base_lvl%HeadZ0(2,:)
                
                do lvl = base_lvl%lvl-1, 1, -1
                    lvl_dif = base_lvl%lvl - lvl
                    if( mod(ix, 2**lvl_dif) == 0 )then
                        call this%lvls_Munjiza(lvl)%SwapZLists()
                    endif 
                end do
                
            end do
            
        end if
        
        base_lvl%Headx0(:) = base_lvl%Headx(:)
        do lvl = base_lvl%lvl-1, 1, -1
            lvl_dif = base_lvl%lvl-lvl
            if( mod(iy, 2**lvl_dif) == 0 )then
                call this%lvls_Munjiza(lvl)%SwapXLists()
            end if
        end do
    
    end do
    
end subroutine

!**********************************************************************
! current number of particles
!**********************************************************************
 subroutine NBS_Munjiza_Hrchl_set_numPrtcl(this, new_numPrtcl )
        implicit none
        class(NBS_Munjiza_Hrchl) this
        integer(IK) new_numPrtcl
        
        integer(IK) lvl
        
        do lvl = 1, this%num_lvls
            call this%lvls_Munjiza(lvl)%set_numPrtcl(new_numPrtcl)
            call this%lvls_Munjiza(lvl)%clc_lvl_numPrtcls()
        end do
    
 end subroutine
 
!**********************************************************************
! setting number of conservative contacts for each level
!**********************************************************************
subroutine NBS_Munjiza_Hrchl_set_num_cntct(this , num)
    implicit none
    class(NBS_Munjiza_Hrchl) this
    integer(IK),dimension(2):: num
        
    integer(IK):: lvl
                
    do lvl = 1, this%num_lvls
            
        call this%lvls_Munjiza(lvl)%set_lvl_num_cntct( num )  
                        
    end do
        
end subroutine
    
!**********************************************************************
! returning sum of number of conservative contacts of all levels     
!**********************************************************************
function NBS_Munjiza_Hrchl_get_num_cntct(this) result(res)
    implicit none
    class(NBS_Munjiza_Hrchl) this
    integer(IK):: res(2)
        
    integer(IK):: lvl
        
    res = 0
        
    do lvl = 1, this%num_lvls
            
        res = res +  this%lvls_Munjiza(lvl)%get_lvl_num_cntct()   
            
    end do
        
end function

!******************************************************************
! final/destructor
!******************************************************************
subroutine Finalize_NBS_Munjiza_Hrchl(this)
    implicit none
    type(NBS_Munjiza_Hrchl)  this

    if( allocated(this%lvls_Munjiza) ) deallocate(this%lvls_Munjiza)
        
end subroutine

end module

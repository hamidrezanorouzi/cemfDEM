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
!  File name:  g_Prtcl_Hrchl_NBS.f90 
!  Module name: g_Prtcl_Hrchl_NBS 
! 
!  Purpose:
!    1) To provide a basic classes for hierarchical contact search  algorithms
!       in this code
!                                 
!    2) It is extended based on the NBS class and contains necessary objects
!       and methods for doing hierarchical contact search in one level.
! 
!  External literature used:
!      - Peters, J. F., et al. (2009). "A hierarchical search algorithm for
!        discrete element method of greatly differing particle sizes."
!        Engineering Computations 26(6): 621-634.
! 
! 
!      - Munjiza, A. (2004). The Combined Finite-Discrete Element Method. West
!        Sussex, UK., John Wiley & Sons Ltd.
!                                           
!      - Munjiza, A. and K. R. F. Andrews (1998). "NBS contact detection
!        algorithm for bodies of similar size." International Journal for
!        Numerical Methods in Engineering 43: 131-149.
! 
!  Referenced modules:
!      - g_Prtcl_CellBased
! 
!------------------------------------------------------------------------------
  
module g_Prtcl_Hrchl_NBS

    use g_TypeDef
    use g_error_handling
    use g_Prtcl_CellBased
    use g_Prtcl_NBS
    
    implicit none
    
    
    private
    public:: NBS_Hrchl
    
    
    !// lvl_NBS
    type , extends(NBS) :: lvl_NBS
        private
        integer(IK)  lvl                ! level number
        integer(IK)::lvl_max_nPrtcl = 0 ! maximum number of particles in this level
        integer(IK)::lvl_numPrtcl = 0   ! number of current particles in this level
        integer(IK)::lvl_num_cnsv_cntct = 0 ! number of conservative contact in this level
        real(RK)     lvl_minD           ! the minimum diameter of bounding boxes in this level
        real(RK)     lvl_maxD           ! the maximum diameter of bounding boxes in this level
        
        integer(IK),allocatable,dimension(:):: memIndx ! memory index of bounding boxes in this level
        
    contains
        
        ! Initializing the level 
        procedure:: Init_lvl_NBS    => lvl_NBS_Init_lvl_NBS   
        
        ! setting and getting the minimum and maximum diameter of bounding boxes in this level
        procedure:: set_lvl_minmax  => lvl_NBS_set_lvl_minmax
        procedure:: get_lvl_minmax  => lvl_NBS_get_lvl_minmax
        
        ! overriding the original function defined in SimWorld class, it returns
        ! the minimum and maximum diameter in this level
        procedure:: get_minmax      => lvl_NBS_get_lvl_minmax ! function overriding 
        
        ! overriding the original function defined in SimWorld class, it returns 
        ! number of particles in this level
        procedure:: get_max_nPrtcl  => lvl_NBS_get_max_nPrtcl      ! function overriding 
        procedure:: get_numPrtcl    => lvl_NBS_get_numPrtcl        ! function overriding
        
        ! overriding the  original function defined in SimWorld class, it returns
        ! the position of ith box in this level
        procedure:: getBndgBox_i    => lvl_NBS_getBndgBox_i   ! function overriding 
        
        ! overriding the  original function defined in SimWorld class, it returns
        ! the memory index of ith box in this level
        procedure:: getMemIndx      => lvl_NBS_getMemIndx     ! function overriding 
        
        ! performing a cross contact search between two levels
        procedure:: BroadSearchCross=> lvl_NBS_BroadSearchCross
        
        ! finding the corresponding index of a cell of the current level in another level 
        procedure:: mapIndex        => lvl_NBS_mapIndex
        
        ! calculating number of bounding boxes which exist in this level 
        procedure:: clc_lvl_numPrtcls  => lvl_NBS_clc_lvl_numPrtcls
        
        ! calculating memory indises of bounding boxes which exist in this level 
        procedure:: clc_lvl_memIndx => lvl_NBS_clc_lvl_memIndx
        
        ! number of contacts related to the cross level search         
        procedure:: set_crs_num_cntct   => lvl_NBS_set_num_cntct
        procedure:: get_crs_num_cntct   => lvl_NBS_get_num_cntct
        procedure:: add_crs_num_cntct   => lvl_NBS_add_num_cntct
        
        ! number of contacts related to the level (cross level + same level)
        procedure:: get_lvl_num_cntct   => lvl_NBS_get_lvl_num_cntct
        procedure:: set_lvl_num_cntct   => lvl_NBS_set_lvl_num_cntct
        
        ! freeing memory allocated for this object 
        final    :: Finilize_lvl_NBS 
        
    end type lvl_NBS
    !// lvl_NBS
    
    
    ! final object that contains information of all levels in hierarchical search algorithm
    type NBS_Hrchl
        
        integer(IK) num_lvls                          ! number of levels
        type(lvl_NBS),allocatable,dimension(:):: lvls_NBS ! level

    contains
        
        !initializing the hierarchical search algorithm 
        procedure   :: Init_NBS_Hrchl=> NBS_Hrchl_Init_NBS_Hrchl 
        
        ! performing broad contact search (fine search is also invoked in this procedure)
        procedure   :: ContactSearch => NBS_Hrchl_ContactSearch
        
        !setting the contact lists (particle-particle)
        procedure   :: setCont_List  => NBS_Hrchl_setCont_List
        
        ! methods related to number of conservative contacts 
        procedure   :: get_num_cntct => NBS_Hrchl_get_num_cntct
        procedure   :: set_num_cntct => NBS_Hrchl_set_num_cntct
        
        ! setting current number of particles
        procedure   :: set_numPrtcl  => NBS_Hrchl_set_numPrtcl  
        
        ! freeing the memory allocated for this object
        final       :: Finalize_NBS_Hrchl
        
    end type NBS_Hrchl
    !// NBS_Hrchl
    
    
contains


!**************************************************************!
!********************* lvl_NBS methods ************************!
!**************************************************************!

    !******************************************************************
    ! Initializing the level
    !******************************************************************
    subroutine lvl_NBS_Init_lvl_NBS(this , lvl )
        implicit none
        class(lvl_NBS) this
        integer(IK),intent(in)  :: lvl
    
    
        this%lvl = lvl
        
        ! Calculating number of particles which belongs to this level.
        ! If diameter of particle fall between the max and min size of level,
        ! it is assigned to this level
        
        call this%clc_lvl_numPrtcls()
            
        ! allocating memIndx vector to store particle indices
        if(allocated(this%memIndx))deallocate(this%memIndx)
        allocate( this%memIndx(this%lvl_max_nPrtcl) )
        call this%clc_lvl_memIndx()
    
    end subroutine
    
    !******************************************************************
    ! setting the minimum and maximum diameter of bounding boxes in this level
    !******************************************************************
    subroutine lvl_NBS_set_lvl_minmax(this, minD , maxD )
        implicit none
        class(lvl_NBS):: this
        real(RK),intent(in) :: minD, maxD

        this%lvl_minD = minD
        this%lvl_maxD = maxD

    end subroutine 
    
    !************************************************************************
    ! returning the minimum and maximum diameter of bounding boxes in this level
    !************************************************************************
    function lvl_NBS_get_lvl_minmax(this) result (res)
    
        implicit none
        class(lvl_NBS):: this
        real(RK),dimension(2):: res

        res(1) = this%lvl_minD
        res(2) = this%lvl_maxD

    end function
    
    !************************************************************************
    !   overriding the original function defined in SimWorld class, it returns 
    ! maximum number of particles in this level
    !************************************************************************
    integer(IK) function lvl_NBS_get_max_nPrtcl(this)
        implicit none
        class(lvl_NBS):: this

        lvl_NBS_get_max_nPrtcl = this%lvl_max_nPrtcl

    end function
    
    !************************************************************************
    !   Overriding the original function defined in SimWorld class. It returns 
    ! current number of particles in this level
    !************************************************************************
    
    integer(IK) function lvl_NBS_get_numPrtcl(this)
        implicit none
        class(lvl_NBS):: this

        lvl_NBS_get_numPrtcl = this%lvl_numPrtcl

    end function
    
    !**************************************************************************
    ! Overriding the original function defined in SimWorld class. It returns
    ! the position of ith box in this level
    !**************************************************************************
    type(real4) function lvl_NBS_getBndgBox_i(this, n)
        implicit none
        class(lvl_NBS):: this
        integer(IK), intent(in) :: n

        lvl_NBS_getBndgBox_i = this%bndg_box( this%memIndx(n) ) 

    end function
    
        
    !**************************************************************************
    ! Overriding the original function defined in SimWorld class. It returns
    ! the memory index of ith box in this level
    !**************************************************************************
    integer(IK) function lvl_NBS_getMemIndx(this, n)
        implicit none
        class(lvl_NBS):: this
        integer(IK), intent(in):: n

        lvl_NBS_getMemIndx =  this%memIndx(n)

    end function
    
    !**************************************************************************
    ! calculating number of bounding boxes which exist in this level
    !**************************************************************************
    subroutine lvl_NBS_clc_lvl_numPrtcls(this)
        implicit none
        class(lvl_NBS):: this
        
        !// locals
        integer(IK) i 
        integer(IK) lvl_max_nPrtcl ! maximum number of particles in this level
        integer(IK) lvl_numPrtcl   ! current number of particles in this level
        integer(IK) tot_Prtcls ! maximum number of particles in the entire simulation / all levels
        integer(IK) numPrtcl  ! current number of particles in all levels
        real(RK)::  minD, maxD, D(2)
        type(real4):: bndg_box
        
        !// body
        
        ! calling the original function
        tot_Prtcls = this%NBS%CellBased%get_max_nPrtcl()
        numPrtcl = this%NBS%CellBased%get_numPrtcl()
        
            
        D = this%get_lvl_minmax()
        MinD = D(1)
        MaxD = D(2)
        
        lvl_max_nPrtcl = 0 
        lvl_numPrtcl = 0
        
        do i = 1, tot_Prtcls
            
            ! invoking the original function
            bndg_box = this%NBS%CellBased%getBndgBox_i (i)
            
            if( bndg_box%w > minD .and. bndg_box%w <= maxD  ) then
                
                lvl_max_nPrtcl = lvl_max_nPrtcl + 1
                if ( i .le. numPrtcl ) lvl_numPrtcl = lvl_numPrtcl + 1
            
            end if
            
        end do

        this%lvl_max_nPrtcl = lvl_max_nPrtcl
        this%lvl_numPrtcl = lvl_numPrtcl
        
    end subroutine
    
    !**************************************************************************
    ! calculating memory indices of bounding boxes which exist in this level
    !**************************************************************************
    subroutine lvl_NBS_clc_lvl_memIndx(this)
        implicit none
        class(lvl_NBS):: this
        
        !// locals
        integer(IK):: i, nPar
        integer(IK):: tot_nPrtcl ! total number of particles in all levels 
        real(RK)::  minD, maxD, D(2)
        type(real4):: bndg_box

        !// body
        
        tot_nPrtcl= this%NBS%CellBased%get_max_nPrtcl()
        
        D = this%get_lvl_minmax()
        MinD = D(1)
        MaxD = D(2)

        nPar = 0 
        do i = 1, tot_nPrtcl
            
            bndg_box = this%NBS%CellBased%getBndgBox_i (i)
            if( bndg_box%w > minD .and. bndg_box%w <= maxD  ) then
                nPar = nPar + 1
                this%memIndx(nPar) = i
            end if
        end do

    end subroutine
    
    !**************************************************************************
    !    finding the corresponding index of a cell of the current level in 
    ! another level
    !**************************************************************************
    type(integer3) function lvl_NBS_mapIndex( this, ind ,crs_lvl )
        implicit none
        class(lvl_NBS)  this
        type(integer3), intent(in)  :: ind
        integer(IK), intent(in)     :: crs_lvl
    
        !//locals
        integer(IK) a, lvl

        !//body
        
        a = 2**( this%lvl-crs_lvl )
    
        lvl_NBS_mapIndex%x = (ind%x -1)/a + 1
        lvl_NBS_mapIndex%y = (ind%y -1)/a + 1
        lvl_NBS_mapIndex%z = (ind%z -1)/a + 1

    end function
    
!********** lvl_NBS_BroadSearchCross (BroadSearchCross) *************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          to perform a broad search between objects in this level *!
!           and the object in the crossNBS (input argument)         *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in:                                                       *!
!*          this    : polymorphic lvl_NBS object                    *!
!*          crossNBS: object of type lvl_NBS, the cross level       *!
!*                    contact search is performed between current   *!
!*                    level and crossNBS                            *!
!*                                                                  *!  
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!      
    subroutine lvl_NBS_BroadSearchCross (this, crossNBS)
        implicit none
        class(lvl_NBS) this
        type (lvl_NBS), intent(in) :: crossNBS

        !// locals
        integer(IK) i, n, m 
        integer(IK) numPrtcl, crs_numPrtcl, icell
        integer(IK) il, jl, kl
        type(integer3) ind, crs_ind, crs_iind

        !// body
        
        ! current number of particles in this level
        numPrtcl = this%get_numPrtcl()
        
        ! current number of particles in cross level 
        crs_numPrtcl = crossNBS%get_numPrtcl()

        ! if no particle is in the cross-level, the contact search is skipped 
        if(crs_numPrtcl == 0 ) return
        
        
        do i = 1, numPrtcl
            
            ! for particles which are in domain
            if ( this%bndg_flag( this%getMemIndx(i) ) >= Pflg_inDomain ) then
            
                ind = this%getIndex(i)
                n = this%Head_inc(ind)

                if( n <= numPrtcl ) then
            
                    !looping over all particles in the cell (ix, iy, iz)
                    do while (n .ne. -1)

                        !finding the corresponding index of ind in the cross level
                        crs_ind = this%mapIndex(ind , crossNBS%lvl)

                        !in the central cell and adjacent cells in the cross level
                        do il = -1,1
                            do jl = -1,1
                                do kl = -1,1
                                
                                    crs_iind = crs_ind + integer3(il, jl, kl)

                                    !finding the first particle 
                                    m = crossNBS%Head(crs_iind) 
                                    if( m > crs_numPrtcl ) m = m - crs_numPrtcl
                                            
                                    ! looping over particles in the central cell and adjacent cells in the cross level
                                    do while( m .ne. -1 )
                                    
                                        ! do the narrow 
                                        call this%FineSearch(this%getMemIndx(n),crossNBS%getMemIndx(m) )
                                        call this%add_crs_num_cntct(1_IK)
                                    
                                        m = crossNBS%Next(m)
                                    end do
                            
                                end do
                            end do
                        end do
                    
                        n = this%Next(n)
                    end do

                end if
                
            end if

        end do

    end subroutine

    !******************************************************************
    ! final/destructor
    !******************************************************************
    subroutine Finilize_lvl_NBS(this)
        implicit none
        type(lvl_NBS):: this
    
        if(allocated(this%memIndx))deallocate(this%memIndx)

    end subroutine
    
    !
    subroutine lvl_NBS_set_num_cntct(this , num )
        implicit none
        class(lvl_NBS) this
        integer(IK) num
        
        this%lvl_num_cnsv_cntct = num
        
    end subroutine 
    
    !
    integer(IK) function lvl_NBS_get_num_cntct(this)
        implicit none
        class(lvl_NBS) this
        
        lvl_NBS_get_num_cntct = this%lvl_num_cnsv_cntct
        
    end function
    
    !
    subroutine lvl_NBS_add_num_cntct(this , num )
        implicit none
        class(lvl_NBS) this
        integer(IK) num
        
        this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + num
    
    end subroutine
    
    !
    function lvl_NBS_get_lvl_num_cntct(this) result(res)
        
        implicit none
        class(lvl_NBS):: this
        integer(IK):: res(2)
        
        res(1) = this%NBS%get_num_cntct()
        res(2) = this%get_crs_num_cntct()
    
    end function
    
    !
    subroutine lvl_NBS_set_lvl_num_cntct(this, num) 
        
        implicit none
        class(lvl_NBS):: this
        integer(IK):: num(2)
        
        call this%NBS%set_num_cntct( num(1) )
        call this%set_crs_num_cntct( num(2) )
    
    end subroutine
    
    
!**************** end of lvl_NBS methods **********************!
!**************************************************************!




!**************************************************************!
!******************* NBS_Hrchl methods ************************!
!**************************************************************!


!************** NBS_Hrchl_Init_NBS_Hrchl (Init_NBS_Hrchl) ***********!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          Initializing NBS_Hrch object                             *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in:                                                       *!
!*          this    : polymorphic NBS object                        *!
!*          minDomain: lower limit of simulation domain             *!
!           maxDomain: higher limit of simulation domain            *!
!*          max_nPrtcl: maximum number of bounding boxes            *!
!*          numPrtcl: number of particles                           *!
!*          ids     : array pointer of bounding box ids             *!
!*          bndg_box: array pointer of bounding box positions       *!
!*          num_lvls: (optional) specifies number of levels to be   *!
!*                    established. If not present, program calculates*!
!*                    number of levels                              *!
!*                                                                  *!  
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!

    subroutine NBS_Hrchl_Init_NBS_Hrchl (this , minDomain, maxDomain, max_nPrtcl, numPrtcl ,ids, flag, bndg_box, num_lvls )
        implicit none
        class(NBS_Hrchl) this
        type(real3),                     intent(in) :: minDomain, maxDomain
        integer(IK),                     intent(in) :: max_nPrtcl, numPrtcl
        type(real4),pointer,dimension(:),intent(in) :: bndg_box
        integer(IK),pointer,dimension(:),intent(in) :: ids, flag
        integer(IK),optional,            intent(in) :: num_lvls
    
        ! locals
        integer(IK) numLevels, lvl
        real(RK)    minD, maxD, lvl_minD, lvl_maxD
        real(RK)::  dmax_min, D(2) 
        type(CellBased):: temp_simWorld

        ! body
        ! finding the minimum and maximum diameter of bounding boxes
        call temp_simWorld%InitSimWorld( max_nPrtcl, numPrtcl ,ids, flag ,bndg_box , .true. )
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
        if( allocated(this%lvls_NBS) ) deallocate(this%lvls_NBS)
        allocate( this%lvls_NBS(this%num_lvls) )
        
        ! setting the maximum diameter of the first level equal to maximum diameter of bounding boxes
        lvl_maxD = maxD
        
        call MainLogInfo%OutInfo(" Number of levels is :" // trim( num2str(this%num_lvls)) ,3 )
        
        do lvl = 1 ,this%num_lvls
        
            call this%lvls_NBS(lvl)%setDomain( minDomain, maxDomain) 
            call this%lvls_NBS(lvl)%InitSimWorld( max_nPrtcl, numPrtcl ,ids, flag, bndg_box , .true. )
            
            ! the minimum diameter of the level is half of the maximum diameter
            lvl_minD = lvl_maxD/2.0_RK 
            
            ! modifying the minimum diameter of the last level and sets it to a very small value
            if( lvl == this%num_lvls) lvl_minD = epsilon(lvl_minD) 
            
            call this%lvls_NBS(lvl)%set_lvl_minmax( lvl_minD, lvl_maxD )
            call this%lvls_NBS(lvl)%Init_lvl_NBS( lvl )
            call this%lvls_NBS(lvl)%initCellBased()
            call this%lvls_NBS(lvl)%Allocate_NBS()
            
            
            !>>> log info
            call MainLogInfo%OutInfo("Level "//trim( num2str(lvl)), 4) 
            call MainLogInfo%OutInfo("   cell size is [m] :"// &
                                     trim( num2str(this%lvls_NBS(lvl)%getCellSize() ) ) , 4, .true.)
            call MainLogInfo%OutInfo("   number of cells (x,y,z) :"//&
                                     trim( num2str( this%lvls_NBS(lvl)%getNumCell()) ) , 4, .true.)
            call MainLogInfo%OutInfo("   number of particles : "// &
                                     trim( num2str(this%lvls_NBS(lvl)%get_max_nPrtcl() ) ) , 4, .true.)
            !<<< log info
            
            
            ! halving the maximum diameter of the level to be the maximum diameter for the next 
            ! level
            lvl_maxD = lvl_maxD/2.0_RK

        end do
        
    end subroutine
    
    !******************************************************************
    ! performing broad contact search 
    !******************************************************************
    subroutine NBS_Hrchl_ContactSearch(this)
        implicit none
        class(NBS_Hrchl) this

        !// locals
        integer(IK) lvl, crs_lvl
        integer(IK) numlvls

        !// body
        
        ! number of levels 
        numlvls = this%num_lvls
        
        ! Initializing number of conservative contacts 
        call this%set_num_cntct( (/0_IK, 0_IK/) )
        
        !   looping over all levels to find bounding box index and
        ! construct the linked list of each level.
        do lvl = 1, numlvls
        
            call this%lvls_NBS(lvl)%BoxIndex()
            call this%lvls_NBS(lvl)%CnstrctList()
            
        end do
        
        
        ! looping over all level to start contact search
        do lvl= numlvls,1,-1
        
            ! same level
            call this%lvls_NBS(lvl)%Ref_Head()
            call this%lvls_NBS(lvl)%BroadSearch()

            ! cross level check with higher levels
            do crs_lvl = lvl-1,1,-1
                
                call this%lvls_NBS(lvl)%Ref_Head()
                call this%lvls_NBS(lvl)%BroadSearchCross( this%lvls_NBS(crs_lvl) )
            
            end do

        end do

    end subroutine 
    
    
    !******************************************************************
    ! setting contact list for all levels
    !******************************************************************
    subroutine NBS_Hrchl_setCont_List(this , PP_CL )
        implicit none 
        class(NBS_Hrchl) this
        class(base_ContactList),pointer, intent(in):: PP_CL        
        
        !// loclas
        integer(IK) i
        
        do i = 1, this%num_lvls
            
            call this%lvls_NBS(i)%setCont_List(PP_CL)
            
        end do
                
    end subroutine 
    
    !******************************************************************
    ! final/destructor
    !******************************************************************
    subroutine Finalize_NBS_Hrchl(this)
        implicit none
        type(NBS_Hrchl)  this
              
        if( allocated(this%lvls_NBS) ) deallocate(this%lvls_NBS)
                
    end subroutine
    
    !******************************************************************
    ! this procedure is related to number of conservative contacts 
    !******************************************************************
    subroutine NBS_Hrchl_set_num_cntct(this , num)
        class(NBS_Hrchl) this
        integer(IK),dimension(2):: num
        
        !// locals
        integer(IK):: lvl
                
        do lvl = 1, this%num_lvls
            
            call this%lvls_NBS(lvl)%set_lvl_num_cntct( num )  
                        
        end do
        
    end subroutine
    
    !******************************************************************
    ! this procedure is related to number of conservative contacts 
    !******************************************************************
    function NBS_Hrchl_get_num_cntct(this) result(res)
        implicit none
        class(NBS_Hrchl) this
        integer(IK):: res(2)
        
        integer(IK):: lvl
        
        res = 0
        
        do lvl = 1, this%num_lvls
            
            res = res +  this%lvls_NBS(lvl)%get_lvl_num_cntct()   
            
        end do
        
    end function
    
    !******************************************************************
    ! setting current number of particles for all levels
    !******************************************************************
    subroutine NBS_Hrchl_set_numPrtcl(this, new_numPrtcl )
        implicit none
        class(NBS_Hrchl) this
        integer(IK),intent(in):: new_numPrtcl
        
        integer(IK) lvl
        
        do lvl = 1, this%num_lvls
            call this%lvls_NBS(lvl)%set_numPrtcl(new_numPrtcl)
            call this%lvls_NBS(lvl)%clc_lvl_numPrtcls()
        end do
    
    end subroutine 
    
end module

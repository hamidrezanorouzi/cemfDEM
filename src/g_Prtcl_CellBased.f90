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
!  File name:  g_Prtcl_CellBased.f90
!  Module name: g_Prtcl_CellBased
!                                
!  Purposes:                      
!   1) Providing a basic class for cell-based contact search algorithms
!                                       
!   2) It is extended the SimWorld class and contains necessary  methods for
!      doing cell-based contact searches
! 
!   Referenced modules:
!       - g_Prtcl_SimWorld
!                         
!------------------------------------------------------------------------------
    
module g_Prtcl_CellBased
    
    use g_error_handling
    use g_Prtcl_SimWorld
    
    implicit none
      
    !// CellBased
    type, extends(SimWorld):: CellBased
        private
        real(RK)    dx      ! length of cells
        real(RK)    ratio   ! ratio of cell length to maximum diameter of boxes
        integer(IK) nx      ! number of divisions in x direction
        integer(IK) ny      ! number of divisions in y direction
        integer(IK) nz      ! number of divisions in z direction
        ! integer coordinate of box
        type (integer3),allocatable,dimension(:):: box_index

    contains
        
        ! initializing CellBased object
        procedure:: InitCellBased   => Cell_InitCellBased
        ! calculating integer coordinates of all boxes
        procedure:: BoxIndex        => Cell_boxIndex
        ! returning number of boxes
        procedure:: getNumBox       => Cell_getNumBox
        ! returning number of nx, ny, nz 
        procedure:: getNumCell      => Cell_getNumCell
        ! returning cell size
        procedure:: getCellSize     => Cell_getCellSize
        ! returning integer coordinate of nth box
        procedure:: getIndex        => Cell_getIndex
        ! freeing memory allocated for this object 
        final    :: FinalizeCellBased
        
    end type CellBased
    
    !// CellBased
    
contains

!*********************************************************!
!*************** CellBased methods ***********************!
!*********************************************************!

!************* Cell_InitCellBased (InitCellBased) *******************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          Initializing CellBased object                           *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in:                                                       *!
!*          this    : polymorphic CelllBased object                 *!
!*          ratio   : (optional) the ratio of cell length to        *!
!*                     maximum diameter of boxes                    *!
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!

    subroutine Cell_InitCellBased( this, ratio )
        implicit none
        class(CellBased) this
        real(RK),optional,intent(in):: ratio

        !// locals 
        integer(IK)             nx, ny, nz, max_nPrtcl
        integer                 iErr
        real(RK)                l_cell_ratio
        real(RK)                dx, maxD
        real(RK),dimension(2):: min_max
        type(real3)          :: min_Domain, max_Domain

        
        if( present(ratio) ) then
            l_cell_ratio = ratio
        else
            l_cell_ratio = 1.000_RK
        end if

        ! getting the minimum and maximum diameter of bounding boxes
        min_max = this%get_minmax()
        maxD = min_max(2)
        dx = l_cell_ratio * maxD
        
        min_Domain = this%getMinDomain()
        max_Domain = this%getMaxDomain()

        nx = max(1 , int( (max_domain%x - min_domain%x)/dx ) + 1 )
        ny = max(1 , int( (max_domain%y - min_domain%y)/dx ) + 1 )
        nz = max(1 , int( (max_domain%z - min_domain%z)/dx ) + 1 )
        
        this%dx = dx
        this%ratio = l_cell_ratio
        this%nx = nx
        this%ny = ny
        this%nz = nz
        
        ! allocating space for maximum number of particles (all particles) 
        max_nPrtcl = this%get_max_nPrtcl()
        
        if(allocated( this%box_index ) ) deallocate(this%box_index)
        allocate( this%box_index(max_nPrtcl), Stat = iErr )
        
        if( iErr .ne. 0 ) then
            call CheckForError(ErrT_Abort , "Cell_InitCellBased" , "Allocation of box_index failed")
        end if
        
    end subroutine 
    
    !******************************************************************
    ! calculating integer coordinates of all boxes
    !******************************************************************
    subroutine Cell_boxIndex(this)
        implicit none
        class(CellBased):: this

        ! locals
        integer(IK) n
        integer(IK) nPrtcl
        type(integer3) ind
        type(real3) minDomain
        type(real4) BndgBox
        
        !// body
        
        ! getting number of current particles
        nPrtcl = this%get_numPrtcl()
        
        ! getting min point of domain
        minDomain = this%getMinDomain()
	
	
	
        ! loop over all particles 
        do n = 1,nPrtcl
            
            ! doing calculations for particles which are in simulation domain 
            if( this%bndg_flag( this%getMemIndx(n) ) >= Pflg_inDomain ) then
                BndgBox = this%getBndgBox_i(n)
                
                ind%x = (BndgBox%x - minDomain%x) / this%dx + 1
                ind%y = (BndgBox%y - minDomain%y) / this%dx + 1
                ind%z = (BndgBox%z - minDomain%z) / this%dx + 1
                
                ! checking if the particle is within system boundaries 
                if( ind%x <= 0_IK .or. ind%x >= this%nx+1 .or. &
                    ind%y <= 0_IK .or. ind%y >= this%ny+1 .or. &
                    ind%z <= 0_IK .or. ind%z >= this%nz+1 ) then
                    ! if not, deleting it
                    this%bndg_flag( this%getMemIndx(n) ) = Pflg_deleted
		                        
                else
                    ! if so, returning index
                    this%box_index(n) = ind
                        
                end if
                
            end if
            
        end do
        
    end subroutine

    !******************************************************************
    ! returning total number of bounding boxes
    !******************************************************************
    integer(IK) function Cell_getNumBox ( this )
        implicit none
        class(CellBased):: this
        
        Cell_getNumBox = this%get_totnumBox()
        
    end function
    
    !******************************************************************
    ! returning number of cells in 3 directions 
    !******************************************************************
    type(integer3) function Cell_getNumCell(this) 
        implicit none
        class(CellBased):: this
        
        Cell_getNumCell%x = this%nx
        Cell_getNumCell%y = this%ny
        Cell_getNumCell%z = this%nz
    
    end function
    
    !******************************************************************
    ! returning cell size
    !******************************************************************
    function Cell_getCellSize(this) result (dx)
        implicit none
        class(CellBased) this
        real(RK) dx
        dx = this%dx
    end function

    !******************************************************************
    ! returning the grid index of nth bounding box
    !******************************************************************
    type(integer3) function Cell_getIndex(this , n)
        implicit none
        class(CellBased):: this
        integer(IK),intent(in)  :: n

        Cell_getIndex = this%box_index( n )

    end function  
    
    !******************************************************************
    ! destructor/final method for CellBased class 
    !******************************************************************
    subroutine FinalizeCellBased(this)        
        implicit none
        type(CellBased) this
        
        if(allocated( this%box_index ) ) deallocate(this%box_index)
        
    end subroutine

!*********** end of CellBased methods ********************!
!*********************************************************!

end module g_Prtcl_CellBased

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
!  File name:  g_Prtcl_SimWorld.f90
!  Module name: g_Prtcl_SimWorld
!          
!  Purpose:
!    1) Providing a basic class for all contact search algorithms in this code
! 
!    2) It is extended based on the SimDomain class and contains necessary
!       objects for doing contact search, like contact list, bounding box ids
!       and bounding box positions objects
! 
!  External literature used: none
!                       
!  Referenced modules:
!    - g_Prtcl_SimDomain
!    - g_Prtcl_ContactList
!
!------------------------------------------------------------------------------
    
    module g_Prtcl_SimWorld

    use g_Prtcl_SimDomain
    use g_Prtcl_ContactList
    use g_Prtcl_DefaultValues
    implicit none
    
    
        
    type, extends(SimDomain) :: SimWorld
        private
        integer(IK)     max_nPrtcl   ! maximum number of particles
        integer(IK)     tot_Bndg_box !
        integer(IK)     numPrtcl     ! current number of particles
        
        !number of conservative contacts in the broad search phase
        integer(IK)::   num_Cnsv_cntct = 0
        
        ! minimum and maximum diameters of bounding boxes
        real(RK)        minDiam
        real(RK)        maxDiam
        
        ! Particle-Particle contact list object 
        class(base_ContactList),pointer ::  PP_Cont_List

        
        ! Bonding box ids and positions (x,y,z,diam)
        integer(IK),public,pointer,dimension(:)::  bndg_ids
        integer(IK),public,pointer,dimension(:)::  bndg_flag
        type(real4),public,pointer,dimension(:)::  bndg_box
        
    contains
    
        ! Initializing simWorld object
        procedure :: InitSimWorld   => simWorld_InitSimWorld
        
        ! setting Ids of bounding boxes
        procedure :: setIds         => simWorld_setIds
        
        ! setting positions of bounding boxes
        procedure :: setBndgBox     => simWorld_setBndgBox
        
        ! setting contact lists 
        procedure :: setCont_List   => simWorld_setCont_List

        ! returning the ith bounding box position
        procedure :: getBndgBox_i   => simWorld_getBndgBox_i 
        procedure,non_overridable   :: getBndgBox_i_non => simWorld_getBndgBox_i_non
        
        ! returning the ith id of bounding box
        procedure :: getBndgIds_i   => simWorld_getBndgIds_i 
        procedure,non_overridable   :: getBndgIds_i_non => simWorld_getBndgIds_i_non
        
        ! returning the memory index of ith bounding box
        procedure :: getMemIndx     => simWorld_getMemIndx   
        
        ! setting, getting and calculating the minimum and maximum diameter of bounding boxes
        procedure :: set_minmax     => simWorld_set_minmax
        procedure :: get_minmax     => simWorld_get_minmax
        procedure :: clc_minmax     => simWorld_clc_minmax
        
        ! returning total number of bounding boxes
        procedure :: get_totnumBox  => simWorld_get_totnumBox
        
        ! returning current and maximum numbers of particles 
        procedure :: get_numPrtcl   => simWorld_get_numPrtcl
        procedure :: get_max_nPrtcl => simWorld_get_max_nPrtcl
               
        ! fine search for Particle-Particle and Particle-Wall contacts 
        procedure:: FineSearch      => Prtcl_FineSearch
        
        
        ! setting, getting and adding operations for number of conservative contacts 
        procedure:: set_num_cntct   => simWorld_set_num_cntct
        procedure:: get_num_cntct   => simWorld_get_num_cntct
        procedure:: add_num_cntct   => simWorld_add_num_cntct
        
        !
        procedure:: set_numPrtcl    => simWorld_set_numPrtcl
            
    end type SimWorld
    
   
    
contains

!*********************************************************!
!***************** SimWorld methods **********************!
!*********************************************************!

!************* simWorld_InitSimWorld (InitSimWorld) *****************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          Initializing SimWorld object                            *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in:                                                       *!
!*          this      : polymorphic SimWorld object                 *!
!*          max_nPrtcl: maximum number of particles                 *!
!*          numPrtcl  : current number of particles                 *!
!*          bndg_ids: array pointer of bounding box ids             *!
!*          bndg_box: array pointer of bounding box positions       *!
!*          clc_minmax: determines whether minimum and maximum      *!
!*                      diameter of bounding boxes are calculated   *!
!*                                                                  *!
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!

    subroutine simWorld_InitSimWorld( this, max_nPrtcl,     &
                                      numPrtcl , bndg_ids , &
                                      bndg_flag, bndg_box , clc_minmax ) 
        
        class(SimWorld):: this
        integer(IK), intent(in)                     :: max_nPrtcl, numPrtcl
        integer(IK),pointer,dimension(:),intent(in) :: bndg_ids, bndg_flag
        type(real4),pointer,dimension(:),intent(in) :: bndg_box
        logical,optional,intent(in)                 :: clc_minmax
        
        !// body
        
        this%max_nPrtcl = max_nPrtcl
        this%tot_Bndg_box = max_nPrtcl 
        this%numPrtcl = numPrtcl
        this%bndg_ids => bndg_ids
        this%bndg_flag=> bndg_flag
        this%bndg_box => bndg_box
        

        ! determining the minimum and maximum length of bounding boxes 
        if( present(clc_minmax) )then
            if( clc_minmax) then
                call this%clc_minmax()
            end if
        end if

    end subroutine
                                      
                                      
!******************* simWorld_SetIds (SetIds) ***********************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          setting ids of bounding boxes                           *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in:                                                       *!
!*          this    : polymorphic SimWorld object                   *!
!*          numBox  : number of bounding boxes                      *!
!*          numPrtcl: number of particles                           *!
!*          bndg_ids: array pointer of bounding box ids             *!
!*                                                                  *!
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!                
    subroutine simWorld_SetIds( this, max_nPrtcl, numPrtcl , bndg_ids ) 
        implicit none
        class(SimWorld):: this
        integer(IK),                     intent(in) :: max_nPrtcl, numPrtcl
        integer(IK),pointer,dimension(:),intent(in) :: bndg_ids
        
        this%max_nPrtcl = max_nPrtcl
        this%tot_Bndg_box = max_nPrtcl
        this%numPrtcl = numPrtcl
        this%bndg_ids => bndg_ids
        

    end subroutine

    
!*************** simWorld_setBndgBox (setBndgBox) *******************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*          setting ids of bounding boxes                           *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in:                                                       *!
!*          this    : polymorphic SimWorld object                   *!
!*          numBox  : number of bounding boxes                      *!
!*          numPrtcl: number of particles                           *!
!*          bndg_box: array pointer of bounding box positions       *!
!*          clc_minmax: determines whether minimum and maximum      *!
!*                      diameter of bounding boxes are calculated   *!
!*                                                                  *!    
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************! 
    subroutine simWorld_setBndgBox( this, max_nPrtcl, numPrtcl , bndg_box , clc_minmax ) 
        implicit none
        class(SimWorld):: this
        integer(IK),                     intent(in) :: max_nPrtcl, numPrtcl
        type(real4),pointer,dimension(:),intent(in) :: bndg_box
        logical,optional,                intent(in) ::  clc_minmax
        
        
        !// body
        
        this%max_nPrtcl = max_nPrtcl
        this%tot_Bndg_box = max_nPrtcl 
        this%numPrtcl = numPrtcl
        this%bndg_box => bndg_box
        
        if( present(clc_minmax) )then
            if( clc_minmax) then
                call this%clc_minmax()
            end if
        end if
        
    end subroutine
    
    !******************************************************************
    ! setting contact list for the object
    !******************************************************************
    subroutine simWorld_setCont_List(this, PP_CL )
        implicit none
        class(simWorld):: this
        class(base_ContactList),pointer,intent(in):: PP_CL
        
        this%PP_Cont_List => PP_CL
    
    end subroutine 

    
    !******************************************************************
    ! returning ith bounding box (real4) 
    !******************************************************************
    type(real4) function simWorld_getBndgBox_i(this , n )
        implicit none
        class(SimWorld):: this
        integer(IK),intent(in)  :: n

        simWorld_getBndgBox_i = this%bndg_box(n)

    end function

    !******************************************************************
    ! returning ith bounding box (real4), non-overridable version
    !******************************************************************
    type(real4) function simWorld_getBndgBox_i_non(this , n )
        implicit none
        class(SimWorld):: this
        integer(IK),intent(in)  :: n

        simWorld_getBndgBox_i_non = this%bndg_box(n)

    end function
    
    !******************************************************************
    ! returning id of ith bounding box
    !******************************************************************
    integer(IK) function simWorld_getBndgIds_i(this , n )
        implicit none
        class(SimWorld):: this
        integer(IK),intent(in)  :: n

        simWorld_getBndgIds_i = this%Bndg_ids(n)

    end function
    
    !******************************************************************
    ! returning id of ith bounding box, non-overridable version
    !******************************************************************
    integer(IK) function simWorld_getBndgIds_i_non(this , n )
        implicit none
        class(SimWorld):: this
        integer(IK),intent(in)  :: n

        simWorld_getBndgIds_i_non = this%Bndg_ids(n)

    end function
    
    !******************************************************************
    ! returning memory index of nth box
    !******************************************************************
    integer(IK) function simWorld_getMemIndx(this , n )
        implicit none
        class(SimWorld):: this
        integer(IK),intent(in)  :: n

        simWorld_getMemIndx = n

    end function
    
    !******************************************************************
    ! setting min and mx diameters of bounding boxes
    !******************************************************************
    subroutine simWorld_set_minmax(this, minD, maxD)
        implicit none
        class(SimWorld):: this
        real(RK),intent(in) :: minD, maxD

        this%minDiam = minD
        this%maxDiam = maxD
        
    end subroutine
    
    !******************************************************************
    ! returning min and max diameters of bounding boxes
    !******************************************************************
    function  simWorld_get_minmax(this) result(res)
        implicit none
        class(SimWorld):: this
        real(RK),dimension(2) :: res
        
        res(1) = this%minDiam 
        res(2) = this%maxDiam 

    end function
    
    !******************************************************************
    ! calculating the minimum and maximum diameters of bounding boxes
    !******************************************************************
    subroutine simWorld_clc_minmax( this )
        implicit none
        class(SimWorld):: this

        ! determines the minimum and maximum length of bounding boxes 
        this%minDiam = minval( this%bndg_box(1:this%max_nPrtcl )%w )
        this%maxDiam = maxval( this%bndg_box(1:this%max_nPrtcl )%w )

    end subroutine 
    
    !******************************************************************
    ! returning total number of bounding boxes
    !******************************************************************
    integer(IK) function simWorld_get_totnumBox(this)
    
        implicit none
        class(SimWorld):: this
        
        simWorld_get_totnumBox = this%tot_Bndg_box
        
    end function
    
    !******************************************************************
    ! returning current number of particles
    !******************************************************************
    integer(IK) function simWorld_get_numPrtcl(this)
    
        implicit none
        class(SimWorld):: this
        
        simWorld_get_numPrtcl = this%numPrtcl
        
    end function
    
    !******************************************************************
    ! returning maximum number of particles
    !******************************************************************
    integer(IK) function  simWorld_get_max_nPrtcl(this)
    
        implicit none
        class(SimWorld):: this
        
        simWorld_get_max_nPrtcl = this%max_nPrtcl
        
    end function  
    
!***************** Prtcl_FineSearch (FineSearch) ********************!
!*                                                                  *!    
!*      Purpose:                                                    *!
!*         1) performing fine contact search between two particles  *!
!*            and determines if they are in contact                 *!
!*         2) adding the contact information to particle-particle   *!
!*            Contact list                                          *!
!*                                                                  *!
!*      Arguments:                                                  *!
!*      - in:                                                       *!
!*          this    : polymorphic SimWorld object                   *!
!*          MemIndx1: Memory index of object in the vector          *!
!*          MemIndx2: Memory index of object in the vector          *!
!*                                                                  *!
!*      - out: none                                                 *!            
!*                                                                  *!    
!*      - return: none                                              *!    
!********************************************************************!
subroutine Prtcl_FineSearch( this, MemIdx1, MemIdx2 )
    implicit none
    class(SimWorld) this
    integer(IK),intent(in)  :: MemIdx1, MemIdx2 
    
    !// locals
    integer(IK) mId, nId, numPrtcl
    real(RK) ovrlp
    type(real4) pos1, pos2
    
    !// body
    
    ! ids of bounding boxes
    mId = this%getBndgIds_i_non(MemIdx1)
    nId = this%getBndgIds_i_non(MemIdx2)   
    
    ! particle-particle contact
    pos1 = this%getBndgBox_i_non(MemIdx1)
    pos2 = this%getBndgBox_i_non(MemIdx2)
    
    ! overlap of particles
    ovrlp = pos1 .ovlp. pos2
        
    if(ovrlp >= 0.0_RK ) then
    
        call this%PP_Cont_List%AddContact( MemIdx1, MemIdx2, mId, nId )
        
    end if
    
end subroutine
    
    
    ! conservative contacts in broad phase
    subroutine simWorld_set_num_cntct(this , num )
        implicit none
        class(simWorld) this
        integer(IK) num
        
        this%num_Cnsv_cntct = num
        
    end subroutine 
    
    ! conservative contacts in broad phase
    integer(IK) function simWorld_get_num_cntct(this)
        implicit none
        class(simWorld) this
        
        simWorld_get_num_cntct = this%num_Cnsv_cntct
        
    end function
    
    ! conservative contacts in broad phase
    subroutine simWorld_add_num_cntct(this , num )
        implicit none
        class(simWorld) this
        integer(IK) num
        
        this%num_Cnsv_cntct = this%num_Cnsv_cntct + num
    
    end subroutine
    
    ! setting current/new number of particles 
    subroutine simWorld_set_numPrtcl(this, new_numPrtcl)
        implicit none
        class(simWorld) this
        integer(IK),intent(in)  :: new_numPrtcl
        
        this%numPrtcl = new_numPrtcl
        
    end subroutine 
    
!*********** end of SimWorld methods *********************!
!*********************************************************!

    
end module g_Prtcl_SimWorld

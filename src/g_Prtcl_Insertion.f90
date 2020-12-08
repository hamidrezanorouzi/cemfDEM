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
!  file name  : g_Prtcl_Insertion.f90  
!  module name: g_Prtcl_Insertion 
!                                 
!  Purpose: 
!    A module for inserting particles with a predefined rate and velocity from
!      a plane
! 
!------------------------------------------------------------------------------
    
module g_Prtcl_Insertion

    use g_PlaneWall
    
    implicit none
    
    
    integer(IK):: min_iter_itv  
    integer(IK):: min_prtcl_ins 
    
    
    type Prtcl_Insertion
        
        integer(IK) numPrtcl        ! current number of particles inserted
        integer(IK) prev_numPrtcl   ! last number of particles inserted
        integer(IK) max_nPrtcl      ! maximum number of particles to be inserted 
        integer(IK) total_steps     ! total number of steps 
        integer(IK) :: max_iter = 100.0_RK 
        real(RK)    :: insertion_vel        ! insertion velocity of particles
        type(PlaneWall) insertion_box       ! insertion plane of particles 
        
        
        integer(IK) next_numPrtcl  ! number of particles to be inserted in the next time
        integer(IK) next_iteration ! iteration number at which next particles should be inserted 
                
        
    contains
    
    procedure:: InitInsertion   => PIns_InitInsertion
    procedure:: InsertParticle  => PIns_InsertParticle
    
    procedure,private:: get_pos             => PIns_get_pos
    procedure,private:: isInContact         => PIns_isInContact
    procedure,private:: Calculate_AllParams => PIns_Calculate_AllParams
    
    end type
    
    
contains

!**********************************************************************
! initializing the particle insertion object
!**********************************************************************
subroutine PIns_InitInsertion(this, box, steps , max_nPrtcl , iVel )
    
    implicit none
    class(Prtcl_Insertion) this
    type(PlaneWall),    intent(in)  :: box
    integer(IK),        intent(in)  :: steps
    integer(IK),        intent(in)  :: max_nPrtcl
    real(RK),optional,  intent(in)  :: iVel
    
    
    integer(IK) num_ins , num_iter
    
    this%numPrtcl = 0_IK
    this%prev_numPrtcl = 0_IK
    this%max_nPrtcl = max_nPrtcl
    this%total_steps= steps
    this%insertion_box = box
    
    
    min_iter_itv    = min( 300 , steps )
    if( int(steps /min_iter_itv) .lt. 10  ) min_iter_itv = max(1,steps/10)
    min_prtcl_ins   = 1
    
    if( int(1.0_RK*steps/min_iter_itv*min_prtcl_ins) > max_nPrtcl  ) then
        min_iter_itv = max( 1, int( steps/max_nPrtcl)-1 )
    end if
    
    this%next_numPrtcl  = min( min_prtcl_ins, max_nPrtcl)
    this%next_iteration = 0     
    
    if( present(iVel) ) then
        this%insertion_vel = iVel
    else
        this%insertion_vel = 0.0_RK
    end if
    
        
end subroutine 

!**********************************************************************
! inserting particles 
!**********************************************************************

function PIns_InsertParticle( this, numPrtcl, dt ,bndg_box, flag, Lin_vel , iter ) result (change)
    
    implicit none
    class(Prtcl_Insertion) this
    integer(IK),                     intent(inout) :: numPrtcl  ! number of particles inserted so far
    real(RK),                        intent(in)    :: dt        ! time step for calculations 
    type(real4),dimension(:),pointer,intent(inout) :: bndg_box  ! position of particles
    integer(IK),dimension(:),pointer,intent(inout) :: flag      ! flag
    type(real3),dimension(:),pointer,intent(inout) :: Lin_vel   ! linear velocity of particle 
    integer(IK),                     intent(in)    :: iter      ! number of iterations required for particle insertion
    logical change
    
    integer(IK) i_par, k
    real(RK) maxD
    logical finished
    type(real3) pos
    type(real4) dpos
    ! checking if it is the iteration of particle insertion
    
    change = .false.
    if( iter /= this%next_iteration ) return
    if( this%numPrtcl == this%max_nPrtcl ) return       
    
    i_par = 0
    finished = .false.
    
    maxD = 0.0_RK
    
    do k = 1, int(this%max_iter * this%next_numPrtcl)
        
        pos = this%get_pos()
        dpos = real4( pos%x, pos%y, pos%z , bndg_box( i_par + this%numPrtcl + 1 )%w )
        
        ! checking if there is a contact between new inserted particle and previous particles        
        if( .not. this%isInContact(dpos, bndg_box ,i_par) ) then
            
            i_par = i_par + 1_IK
            bndg_box(this%numPrtcl + i_par) = dpos
            flag(this%numPrtcl + i_par)     = Pflg_inDomain
            Lin_vel(this%numPrtcl + i_par)  = this%insertion_vel* this%insertion_box%getNormal()
            maxD = max( maxD , dpos%w)
            
        end if
        
        if( i_par == this%next_numPrtcl ) then
            finished = .true.
            exit
        end if
        
        
    end do
        
    this%prev_numPrtcl = max( this%numPrtcl - 50*this%next_numPrtcl , 0  ) 
    
    numPrtcl = numPrtcl + i_par
    
    this%numPrtcl = numPrtcl
    change = .true.   
    
    
    call this%Calculate_AllParams( finished, maxD, dt )
    
    
end function


subroutine PIns_Calculate_AllParams( this , finished, maxD, dt )

    implicit none
    class(Prtcl_Insertion) this 
    logical,            intent(in)  :: finished
    real(RK),           intent(in)  :: maxD, dt
    real(RK) particle_left, iter_left
    
    integer(IK) num_iter, nxt_nPrtcl
    
    if( .not. finished ) then
        
        ! if not finished, we must increase the initial velocity of insertion
        this%insertion_vel = max( 1.2_RK * this%insertion_vel, 0.03)
        
        call MainLogInfo%OutInfo("insertion velocity is increased to :"// trim( num2str(this%insertion_vel) ) , 1 )
        
    end if
    
    
    particle_left = (this%max_nPrtcl-this%numPrtcl)
    iter_left     = this%total_steps- this%next_iteration
        
    
    nxt_nPrtcl = max( min_prtcl_ins, int(1.0_RK*particle_left /(iter_left/min_iter_itv))+1 )
    num_iter   = min_iter_itv
            
    
    this%next_iteration = this%next_iteration + num_iter
        
    if( this%numPrtcl+nxt_nPrtcl > this%max_nPrtcl ) then
        nxt_nPrtcl = this%max_nPrtcl - this%numPrtcl
    end if
    
    this%next_numPrtcl = nxt_nPrtcl    
    
    
end subroutine 



logical function PIns_isInContact(this, dpos, bndg_box, nPar ) 

    implicit none
    class(Prtcl_Insertion) this
    type(real4),intent(in)  :: dpos
    type(real4),dimension(:),pointer,intent(inout) :: bndg_box
    integer(IK),intent(in)  :: nPar
    
    integer(IK) i
    
    PIns_isInContact = .false.
    
    ! searching in the current particles
    do i = this%prev_numPrtcl+1, this%numPrtcl + nPar
        
        if( ( dpos .ovlp. bndg_box(i) ) > 0.0_RK ) then
            PIns_isInContact = .true.
            return
        endif
        
    end do    
    
   

end function

function PIns_get_pos(this) result (pos)

    implicit none
    class(Prtcl_Insertion) this
    type(real3) pos
    
    real(RK) u , v, w, x, s
   
    call RANDOM_NUMBER(u)
    call RANDOM_NUMBER(v)
    call RANDOM_NUMBER(w)
    call RANDOM_NUMBER(x)
    
    s = u + v + w + x 
    
    u = u/s
    v = v/s
    w = w/s
    x = x/s
    
    pos = u*this%insertion_box%P1 + v* this%insertion_box%P2 + w* this%insertion_box%P3 + x* this%insertion_box%P4
    
        
end function

end module

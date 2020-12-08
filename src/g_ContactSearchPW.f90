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
!  file name  : g_ContactSearchPW.f90            	
!  module name: g_ContactSearchPW                                 
!                                                                   
!  Purpose:                                                        
!    Managing all operations related to the particle-wall contact search   
! 
!------------------------------------------------------------------------------ 
module g_ContactSearchPW
    
    use g_TypeDef
    use g_Prtcl_ContactList
    use g_Geometry
    
    implicit none
    
    
    !//ContactSearchPW
    type ContactSearchPW
        
        private
        integer(IK) :: start_iter_update = 0 
        integer(IK) :: iter_to_update       = 0
        class(base_ContactList),pointer:: PW_cont_list
        class(Geometry),pointer:: m_Geometry
        integer(IK) :: max_nearPrtcl_pWall = 1000
        integer(IK) :: num_nearPrtcl_pWall = 0
        integer(IK),dimension(:,:),allocatable:: pWall_NearPrtcl_range
        integer(IK),dimension(:)  ,allocatable:: pWall_NearPrtcls    
        
    contains
    
    ! Initializing the object
    procedure:: InitContactSearch   => PW_CS_InitContactSearch
    
    ! Find particles which are near all walls
    procedure:: FindNearPrtcls      => PW_CS_FindNearPrtcls
    
    ! calculating number of iterations which should be performed
    ! untill the next update of the list of particles near the wall 
    procedure:: clc_num_iter        => PW_CS_clc_num_iter
    
    ! Performing contact search to determine particle-wall contacts 
    procedure:: FindContacts        => PW_CS_FindContacts
    
    ! returning next iteration in which update should be performed
    procedure:: Next_update         => PW_CS_Next_update
    
    !// resetting the the update iteration
    procedure:: Reset_update        => PW_CS_Reset_update
    
    ! private: finding particles which are near a wall
    procedure,private:: pWall_nearPrtcl => PW_CS_pWall_nearPrtcl
    procedure,private:: Increase_1Dsize => PW_CS_Increase_1Dsize
        
    end type ContactSearchPW
    !// ContactSearchPW
    
contains

!**********************************************************************
! Initializing the object
!**********************************************************************
subroutine PW_CS_InitContactSearch(this , PW_cont_list , m_Geometry )
    implicit none
    class(ContactSearchPW) this
    class(base_Contactlist),pointer,intent(in) :: PW_cont_list
    class(Geometry),pointer,intent(in)         :: m_Geometry
    
    ! locals
    integer(IK) max_pWall
    
    !// body
    this%PW_cont_list => PW_cont_list
    this%m_Geometry   => m_Geometry
    
    max_pWall = m_Geometry%get_max_pWall()
    
    if( allocated(this%pWall_NearPrtcl_range)) deallocate( this%pWall_NearPrtcl_range )
    allocate(this%pWall_NearPrtcl_range( 2 , max_pWall) )
    
    if( allocated(this%pWall_NearPrtcls) ) deallocate(this%pWall_NearPrtcls)
    allocate( this%pWall_NearPrtcls(this%max_nearPrtcl_pWall) )    
    
end subroutine 

!**********************************************************************
! Finding particles which are near all walls
!**********************************************************************
function PW_CS_FindNearPrtcls( this, dt , iterNumber, nPrtcl, bndg_box, flag, lin_vel, lin_acc ) result( res )
    implicit none
    class(ContactSearchPW)     this
    real(RK)   ,intent(in)  :: dt
    integer(IK),intent(in)  :: iterNumber, nPrtcl
    type(real4),dimension(:),pointer,intent(in):: bndg_box
    integer(IK),dimension(:),pointer,intent(in):: flag
    type(real3),dimension(:),pointer,intent(in):: lin_vel, lin_acc
    logical res
    
    !// locals
    integer(IK) nw , i, numIter
    real(RK)    max_d, dx
    type(mvng_PlaneWall) pWall
    
    
    !// body
    ! checks if in this iteration the neighbour list should be updated
    if( .not. (iterNumber >= this%Next_update() ) ) then
        res = .false.
        return
    end if
    
    
    max_d = maxval( bndg_box(1:nPrtcl)%w ) 
    dx = max_d * d_Wall_neighbor_ratio
    
    this%num_nearPrtcl_pWall = 0 

    ! number of plane walls 
    nw = this%m_Geometry%get_num_pWall()
    
    do i = 1, nw
        
        pWall = this%m_Geometry%get_pWall(i)
        call this%pWall_nearPrtcl( pWall , i , dx, nPrtcl, bndg_box, flag )
                
    end do
    
    call this%clc_num_iter(dx-0.5*max_d, dt, nPrtcl, flag, lin_vel, lin_acc , numIter )
    
    this%start_iter_update = iterNumber
    this%iter_to_update = numIter
    
    res = .true.
    
end function

!**********************************************************************
! calculating number of iterations which should be performed
! untill the next update of the list of particles near the wall 
!**********************************************************************
subroutine  PW_CS_clc_num_iter(this, dx , dt , numPrtcl , flag,  lin_vel , lin_acc, num_iter )
    implicit none
    class(ContactSearchPW)     this
    real(RK)    ,intent(in) :: dx, dt
    integer(IK) ,intent(in) :: numPrtcl
    integer(IK),dimension(:),pointer, intent(in) :: flag
    type(real3),dimension(:),pointer, intent(in) :: lin_vel, lin_acc
    integer(IK),intent(out) :: num_iter
    
    !// locals
    integer(IK) n, nw
    real(RK) max_vel_prtcl, max_acc_prtcl, max_vel_wall
    real(RK) max_v, max_a
    real(RK) del , t
    type(mvng_PlaneWall) wall
    
    
    !// body 
    ! calculating the maximum velocity in the system
    max_vel_prtcl = 0 
    do n = 1, numPrtcl
        if( flag(n) >= Pflg_inDomain ) max_vel_prtcl = max( max_vel_prtcl , norm( lin_vel(n) ) )
    end do
    
    !calculating the maximum acceleration in the system
    max_acc_prtcl = 0 
    do n = 1, numPrtcl
        if( flag(n) >= Pflg_inDomain ) max_acc_prtcl = max( max_acc_prtcl , norm( lin_acc(n) ) )
    end do
    
    
    nw = this%m_Geometry%get_num_pWall()
    max_vel_wall = 0 
    do n = 1,nw
        wall = this%m_Geometry%get_pWall( n )
        max_vel_wall = max( max_vel_wall, norm( wall%GetVelocity( wall%P1 ) ) )
        max_vel_wall = max( max_vel_wall, norm( wall%GetVelocity( wall%P2 ) ) )
        max_vel_wall = max( max_vel_wall, norm( wall%GetVelocity( wall%P3 ) ) )
        max_vel_wall = max( max_vel_wall, norm( wall%GetVelocity( wall%P4 ) ) )
    end do
     
    
    max_v = max_vel_prtcl + max_vel_wall
    max_a = max_acc_prtcl
    
    max_v = max(max_v,0.3_RK)
    max_a = max(max_a,9.8_RK)
    
    ! solving the equation of motion for t
    ! x-x0 = v*t + 0.5*a*t^2
    ! or  0.5*a*t^2 + v*t - dx = 0 
    
    del = max_v**2.0_RK - ( 4.0_RK * (0.5_RK*max_a)*(-dx) )
    
    ! the value of del is always greater than 0. This equation has two roots, a positive and a negative
    ! we need the positive root
    
    t = ( - max_v + sqrt(del) )/(max_a) 
    
    if( t < 0.0 ) then
        call CheckForError( ErrT_Abort, "PW_CS_clc_num_iter" , "Error in root finding of equaiton of motion" )
        return
    end if
    
    num_iter = t/dt 
    num_iter = max(num_iter,1_IK)
    num_iter = min(num_iter, d_Wall_max_update_iter)
    
end subroutine 

!**********************************************************************
! Performing contact search to determine particle-wall contacts 
!**********************************************************************
subroutine PW_CS_FindContacts( this, numPrtcl, bndg_box , bndg_ids)
        
    implicit none
    class(ContactSearchPW) this
    integer(IK),intent(in):: numPrtcl
    type(real4),pointer,dimension(:),intent(in)::bndg_box
    integer(IK),pointer,dimension(:),intent(in)::bndg_ids
    
    
    integer(IK) nw, n, i
    integer(IK):: rng(2), npar
    type(real4):: prtcl
    !
    !getting the number of plane walls
    nw = this%m_Geometry%num_pWall
    
    do i =1, nw
        
         rng = this%pWall_NearPrtcl_range(:,i)
         
        do n = rng(1), rng(2)
             
             ! getting the particle index which is near wall
       
             npar = this%pWall_NearPrtcls(n)
             prtcl = bndg_box(npar)
             
             if( this%m_Geometry%pWall(i)%isInContact(prtcl) )then
                 ! mem_indx_i, mem_indx_j, id_i, id_j
                 call this%PW_cont_list%AddContact( npar, n , bndg_ids(npar) , this%m_Geometry%pWall(i)%wall_id )
             end if
                
         end do
        
    end do
    
    ! removing the prvious contacts (those which are not in contact now) 
    call this%PW_Cont_List%RemvReleased()
     
end subroutine 

!***********************************************************************
!
!***********************************************************************
function PW_CS_Next_update(this) result (res)
    implicit none
    class(ContactSearchPW) this
    integer(IK) res
    
    res = this%Start_iter_update + this%iter_to_update 
    
end function

!***********************************************************************
!
!***********************************************************************
subroutine PW_CS_Reset_update(this)
    implicit none
    class(ContactSearchPW) this
    this%Start_iter_update = 0
    this%iter_to_update = 0
    
end subroutine


!**********************************************************************
! finding particles which are near a wall
!********************************************************************** 
subroutine PW_CS_pWall_nearPrtcl(this , Wall , pw_ind, dx , nPrtcl ,  bndg_box, flag )
    implicit none 
    class(ContactSearchPW)      this 
    type(mvng_PlaneWall),intent(in):: Wall
    integer(IK) ,intent(in) :: pw_ind
    real(RK)    ,intent(in) :: dx
    integer(IK) ,intent(in) :: nPrtcl
    type(real4),dimension(:),pointer:: bndg_box
    integer(IK),dimension(:),pointer:: flag
    
    !// locals
    integer(IK) i , f_ind, l_ind
    type(integer3) indx
    type(real3) min_point, max_point, ddx , prcl
    type(PlaneWall) ext_wall, par_wall, mirror_wall
    
    
    !// body
    
   ! ext_wall = Wall%Extend_boundaries(dx)  ! extended wall
   ! par_wall = ext_wall%Parallel_Plane(dx) ! parallel wall
    
    !min_point = min( par_wall%p4, min( par_wall%p3, &
    !            min( par_wall%p2, min( par_wall%p1, &
    !            min( ext_wall%p4, min( ext_wall%p3, &
    !            min( ext_wall%p1, ext_wall%p2) ) ) ) ) ) )
    !
    !
    !max_point = max( par_wall%p4, max( par_wall%p3, &
    !            max( par_wall%p2, max( par_wall%p1, &
    !            max( ext_wall%p4, max( ext_wall%p3, &
    !            max( ext_wall%p1, ext_wall%p2) ) ) ) ) ) )
    
    
    min_point = min( wall%p4, min( wall%p3, &
                min( wall%p2,  wall%p1) ) )
                
    
    
    max_point = max( wall%p4, max( wall%p3, &
                max( wall%p2, wall%p1) ) )
    
    min_point = min_point - (dx*real3(1.0,1.0,1.0));
    max_point = max_point + (dx*real3(1.0,1.0,1.0));
                
    !if( Wall%bothside )    then
    !    
    !    mirror_wall = ext_wall%Parallel_Plane( -dx )
    !    
    !    min_point = min( min_point, min(mirror_wall%P4 , min( mirror_wall%p3 , min( mirror_wall%p2, mirror_wall%P1) ) ) )
    !    max_point = max( max_point, max(mirror_wall%P4 , max( mirror_wall%p3 , max( mirror_wall%p2, mirror_wall%P1) ) ) )
    !    
    !end if
    
    ! calculating the dimension of box
    ddx = max_point - min_point
    
    ! using multiplicaiton instead of division here after
    ddx = real3(1.0_RK/ddx%x , 1.0_RK/ddx%y, 1.0_RK/ddx%z)
    
    ! the first index
    f_ind = this%num_nearPrtcl_pWall+1
    
    do i = 1, nPrtcl
        
        !for those particles which are in the simulaiton domain
        if( flag(i) >= Pflg_inDomain ) then
        
            prcl = bndg_box(i)
            indx = (prcl-min_point)*ddx
            
            ! checking if the particle is in the box 
            if(indx%x == 0_IK .and. indx%y == 0_IK .and. indx%z == 0_IK ) then 
            
                
                ! enough space to keep new entry? 
                if( this%num_nearPrtcl_pWall >= this%max_nearPrtcl_pWall ) then
                    ! if not, extending it 
                    call this%Increase_1Dsize(this%pWall_NearPrtcls , this%max_nearPrtcl_pWall )
                end if
                
            ! add it to the list 
            this%num_nearPrtcl_pWall = this%num_nearPrtcl_pWall + 1_IK
            this%pWall_NearPrtcls( this%num_nearPrtcl_pWall ) = i 
            
            
            end if
            
        end if
            
    end do
    
    ! recording  first and last indexes of near particles (in pWall_NearPrtcls) for this wall
    l_ind = this%num_nearPrtcl_pWall
    this%pWall_NearPrtcl_range(:,pw_ind) = (/f_ind, l_ind/)
       
end subroutine

subroutine PW_CS_Increase_1Dsize(this, arry , n )
    implicit none
    class(ContactSearchPW) this
    integer(IK),allocatable,dimension(:) :: arry
    integer(IK) n
    
    !// locals
    integer(IK) old_n
    integer(IK),allocatable,dimension(:) :: temp
    
    !// body
    
    old_n = n
    
    allocate( temp(old_n) )
    temp(1:old_n) = arry(1:old_n)
    
    deallocate(arry)
    n = (1.5*n) + 1
    allocate( arry(n) )
    
    arry(1:old_n) = temp(1:old_n)
    
    deallocate(temp) 

end subroutine 


end module

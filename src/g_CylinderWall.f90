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
!  file name  : g_CylinderWall.f90 
!  module name: g_CylinderWall     
!                                  
!  Purpose:                        
!   Defining a class for a cylinder wall and required methods for creating
!   and cutting it
!                 
!------------------------------------------------------------------------------

module g_CylinderWall
    
    
    use g_PlaneWall
    use g_error_handling
    
    implicit none
    
    
    type CylinderWall
        
        private
        integer(IK) user_id
        integer(IK) prop_type
        integer(IK) nw
        type(mvng_PlaneWall),allocatable,dimension(:):: walls
        
        type(p_line) axis
        real(RK)         rad1                   ! radius at z = 0
        real(RK)         rad2                   ! radius at z = l
        logical::        bothSide   = .false.   ! 
        logical::        mvng       = .false.   ! 
        logical::        rotating   = .false.   ! 
        type(p_line) rot_line                   ! rotation line
        real(RK)         rot_vel                ! rotation velocity rad/s
        type(real3)      trans_vel              ! translational velocity   
        
    contains
    
    procedure:: CW_CreateCylinder1  ! creating cylinder wall
    procedure:: CW_CreateCylinder2  ! creating cylinder wall 
    generic:: CreateCylinder    => CW_CreateCylinder1 , CW_CreateCylinder2
    procedure:: cut_wall        => CW_cut_wall
    procedure:: get_numWall     => CW_get_numWall
    procedure:: get_Wall        => CW_get_Wall
        
        
    end type
    
    
    !// constructors 
    interface CylinderWall
        procedure CylinderWall_const1, CylinderWall_const2
    end interface
    
    
contains

!**********************************************************************
! constructor for creating a cylinder with predefined velocities
!**********************************************************************
function CylinderWall_const1( rad1, rad2, axis, nw , user_id, prop_type, t_vel, r_vel, r_line, both ) result ( this )
    implicit none
    real(RK),        intent(in) :: rad1, rad2 ! radii of the cylinder
    type(p_line),intent(in) :: axis           ! axis of the cylinder 
    integer(IK),     intent(in) :: nw     ! number of divisions 
    integer(IK),     intent(in) :: user_id, prop_type
    type(real3),     intent(in) :: t_vel  ! translational velocity
    real(RK),        intent(in) :: r_vel  ! rotational velocity
    type(p_line),intent(in) ::     r_line ! axis of rotation
    logical,optional,intent(in) :: both   ! both side active status 
    type(CylinderWall) this
    
    logical res
    
    if( present(both) ) then
        res = this%CreateCylinder(rad1, rad2, axis, nw , user_id, prop_type, t_vel, r_vel, r_line, both)
    else
        res = this%CreateCylinder(rad1, rad2, axis, nw , user_id, prop_type, t_vel, r_vel, r_line)
    end if
    
    if( .not. res ) then
        call CheckForError(ErrT_Abort, "CylinderWall_const1", "Cannot create cylindrical wall" )    
        return
    end if
    
    
end function

!**************************************************************************
! constructor for creating a stationary cylinderical wall 
!**************************************************************************
function CylinderWall_const2( rad1, rad2, axis, nw , user_id, prop_type , both ) result ( this )
    implicit none
    real(RK),        intent(in) :: rad1, rad2 ! radii of the cylinder 
    type(p_line),    intent(in) :: axis       ! axis of the cylinder 
    integer(IK),     intent(in) :: nw         ! number of divisions 
    integer(IK),     intent(in) :: user_id, prop_type
    logical,optional,intent(in) :: both       ! both side active status 
    type(CylinderWall) this
    
    logical res
    
    if( present(both) ) then
        res = this%CreateCylinder(rad1, rad2, axis, nw , user_id, prop_type, both)
    else
        res = this%CreateCylinder(rad1, rad2, axis, nw , user_id, prop_type)
    endif
    
    if( .not. res ) then
        call CheckForError(ErrT_Abort, "CylinderWall_const1", "Cannot create cylindrical wall" )    
        return
    end if
    
end function


!*******************************************************************************
! creating a cylinderical wall with specified velocity
!*******************************************************************************
logical function CW_CreateCylinder1(this, rad1, rad2, axis, nw , user_id, prop_type, t_vel, r_vel, r_line, both ) result (res)
    implicit none
    class(CylinderWall) this
    real(RK),        intent(in) :: rad1, rad2  ! radii of the cylinder 
    type(p_line),    intent(in) :: axis        ! axis of the cylinder 
    integer(IK),     intent(in) :: nw          ! number of divisions 
    integer(IK),     intent(in) :: user_id, prop_type
    type(real3),     intent(in) :: t_vel       ! translational velocity
    real(RK),        intent(in) :: r_vel       ! rotational velocity   
    type(p_line),    intent(in) :: r_line      ! axis of rotation 
    logical,optional,intent(in) :: both        ! both side active status 
        
        
    !// locals
    logical lboth
    integer(IK) n
    real(RK)    L, teta, dteta ,x, y, z
    type(real3) p1, p2, p3, p4
    type(real3),allocatable,dimension(:):: r1_p, r2_p
    type(zAxisLine) zAxis
    !// boty
        
    if( present(both) ) then
        lboth = both
    else
        lboth = .false.
    end if
        
    this%user_id   = user_id
    this%prop_type = prop_type
    this%nw = nw
        
    this%axis = axis
    this%rad1 = rad1
    this%rad2 = rad2
    this%bothSide = lboth
        
    if( norm(t_vel) > 0.0000001_RK ) then
        this%mvng = .true.
    end if
        
    if( abs(r_vel) > 0.0000001_RK )then
        this%rotating = .true.    
    end if
        
    this%trans_vel = t_vel
    this%rot_vel   = r_vel
    this%rot_line  = r_line
        
    !// translates axis to z-axis
    call zAxis%Rotate_toZ1(axis)
        
    ! now the first point of zAxis line is (0,0,0) and the second 
    ! point is (0,0,L), where L is the length of cylinder.
        
    L = zAxis%get_length()
        
    if(allocated(r1_p))deallocate(r1_p)
    allocate( r1_p(nw+1) )
        
    if(allocated(r2_p))deallocate(r2_p)
    allocate( r2_p(nw+1) )
        
        
    dteta = 2*pi/nw
    teta = 0
    ! creating points on z-axis aligned coordinates
    do n = 1, nw+1
            
        r1_p(n) = real3( rad1*cos(teta) , rad1*sin(teta) , 0.0_RK )
        r2_p(n) = real3( rad2*cos(teta) , rad2*sin(teta) , L )
        teta = teta + dteta
            
    end do
        
    ! transferring back all points to the original axis of cylinder
        
    do n = 1, nw+1
            
        r1_p(n) = zAxis%rBack_fromZ( r1_p(n) )
        r2_p(n) = zAxis%rBack_fromZ( r2_p(n) )
            
    end do
        
    ! creating plane walls from all points 
    if( allocated(this%walls) ) deallocate(this%walls)
    allocate( this%walls(nw) )
        
    do n = 1, nw
            
        if( .not. this%walls(n)%CreateWall( r1_p(n), r2_p(n), r2_p(n+1), r1_p(n+1),&
            user_id, 0 , prop_type, t_vel, r_vel, r_line ,lboth ) ) then
                
            call CheckForError(ErrT_Abort , "CW_CreateCylinder1", "Cannot create cylindrical wall" ) 
            res = .false.
            return
                
        end if
            
    end do
        
        
        
    deallocate( r1_p, r2_p)
    res = .true.
            
end function

!*******************************************************
! creating an stationary cylinderical wall 
!*******************************************************
logical function CW_CreateCylinder2(this, rad1, rad2, axis, nw , user_id, prop_type , both ) result (res)
    implicit none
    class(CylinderWall) this
    real(RK),        intent(in) :: rad1, rad2 ! radii of the cylinder 
    type(p_line),    intent(in) :: axis       ! axis of the cylinder 
    integer(IK),     intent(in) :: nw         ! number of divisions 
    integer(IK),     intent(in) :: user_id, prop_type
    logical,optional,intent(in) :: both       ! both side active status 
        
    !// locals
        
    logical lboth 
    !// body
        
    if( present(both) ) then
        lboth = both
    else
        lboth = .false.
    end if
        
    if( .not. this%CreateCylinder(rad1, rad2, axis, nw, user_id, prop_type, zero_r3, 0.0_RK , axis ,lboth)  ) then
        res = .false.
        return
    end if
        
end function

!*******************************************************************************
! cutting the cylinder wall with respect the active side of a plane
!*******************************************************************************
function CW_cut_wall(this, plane ) result( res )
    implicit none
    class(CylinderWall) this
    class(PlaneWall),intent(in):: plane   ! the cylinder is cut with respect to the active side of this plane 
    type(CylinderWall) res
        
    integer(IK) i, new_nw, n
    integer(IK),dimension(4):: off_points
    type(real3),dimension(4):: pw1_points, pw2_points
        
                
    ! first finding number of new walls after cut (new_nw)
    new_nw = 0
    do i=1, this%nw
            
        off_points = this%walls(i)%test_inactive(plane)
            
        select case( sum(off_points) )
        case (0) ! all points of wall on the active side
            new_nw = new_nw + 1
            
        case(1) ! three points are on the active side
                
            new_nw = new_nw + 2
                
        case(2,3) ! two (or one) points are on the active side
                
            new_nw = new_nw + 1
                
        case(4) ! no point is on the active side 
            new_nw = new_nw
        end select
            
    end do
        
        
    ! allocating memory for new walls
    res%nw = new_nw 
    if( allocated( res%walls) ) deallocate(res%walls)
    allocate( res%walls(new_nw) )
        
    n = 0
    do i=1, this%nw
            
        off_points = this%walls(i)%test_inactive(plane)
            
        select case( sum(off_points) )
        case (0) ! all points of wall are on the active side
            n = n+1
            ! copying  wall into the new cylinder
            res%walls(n) = this%walls(i)
            
        case(1) ! one point is on the off-side of the wall
                
            call this%walls(i)%cut_wall_1point( plane, pw1_points, pw2_points )
            n = n + 1
            if( .not. res%walls(n)%CreateWall( pw1_points(1), pw1_points(2), &
                                                pw1_points(3), pw1_points(4), &
                                                this%user_id , 0,             &
                                                this%prop_type, this%trans_vel,&
                                                this%rot_vel, this%rot_line,   &
                                                this%bothSide)   ) then
                    
                    
                call CheckForError(ErrT_Abort , "CW_cut_wall", "Cannot create plane wall after cutting" ) 
                
            end if
                
            n = n + 1
            if( .not. res%walls(n)%CreateWall( pw2_points(1), pw2_points(2), &
                                                pw2_points(3), pw2_points(4), &
                                                this%user_id , 0,             &
                                                this%prop_type, this%trans_vel,&
                                                this%rot_vel, this%rot_line,   &
                                                this%bothSide)   ) then
                    
                    
                call CheckForError(ErrT_Abort , "CW_cut_wall", "Cannot create plane wall after cutting" ) 
                
            end if
                                                   
        case(2) ! two points are on the off-side of the wall
                
            call this%walls(i)%cut_wall_2points( plane, pw1_points )
            n = n + 1
            if( .not. res%walls(n)%CreateWall( pw1_points(1), pw1_points(2), &
                                                pw1_points(3), pw1_points(4), &
                                                this%user_id , 0,             &
                                                this%prop_type, this%trans_vel,&
                                                this%rot_vel, this%rot_line,   &
                                                this%bothSide)   ) then
                    
                    
                call CheckForError(ErrT_Abort , "CW_cut_wall", "Cannot create plane wall after cutting" ) 
                
            end if
                                                   
        case(3) ! three points are on the off-side of the wall
                
            call this%walls(i)%cut_wall_3points( plane, pw1_points )
            n = n + 1
            if( .not. res%walls(n)%CreateWall( pw1_points(1), pw1_points(2), &
                                                pw1_points(3), pw1_points(4), &
                                                this%user_id , 0,             &
                                                this%prop_type, this%trans_vel,&
                                                this%rot_vel, this%rot_line,   &
                                                this%bothSide)   ) then
                    
                    
                call CheckForError(ErrT_Abort , "CW_cut_wall", "Cannot create plane wall after cutting" ) 
                
            end if
                                                   
        case(4) ! all points are on the off-side of the wall and this wall will be deleted 
                
        end select
            
            
    end do
     
        
    res%user_id   = this%user_id
    res%prop_type = this%prop_type
    res%axis      = this%axis
    res%rad1      = this%rad1
    res%rad2      = this%rad2
    res%bothSide  = this%bothSide
    res%mvng      = this%mvng
    res%rotating  = this%rotating
    res%rot_vel   = this%rot_vel
    res%trans_vel = this%trans_vel


end function
    
        
!************************************************************************
! returning number of plane walls in the cylinder 
!************************************************************************    
integer(IK) function CW_get_numWall(this)
    implicit none
    class(CylinderWall) this
        
    CW_get_numWall = this%nw
    
end function

!*************************************************************
! returning the nth plane wall in the cylinderical wall
!*************************************************************
function CW_get_Wall(this , n ) result( res )
    implicit none
    class(CylinderWall)      this
    integer(IK),intent(in):: n
    type(mvng_PlaneWall)     res
        
    if(n <= this%nw ) then
        res = this%walls(n)
    else
        call CheckForError( ErrT_Abort, "CW_getWall", "Out of range index for wall :"//num2str(n) )     
    end if
        
end function
    
    
end module

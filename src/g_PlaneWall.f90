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
!  file name  : g_PlaneWall.f90
!  module name: g_PlaneWall
!                          
!  Purpose:                
!  Defining a class for a plane wall and all the methods required for contact
!    detection, distance, contact point and etc.
! 
! 
! 
!  The wall can move with translational and rotational motions 
!    * Translational motion: 
!       - Three components of linear velocity should be specified and wall
!         moves linearly
!     * Rotational motion:
!       - A rotation axis (line) and a rotational speed around this axis 
!         should be specified.
! 
!------------------------------------------------------------------------------
    
module g_PlaneWall

   
    use g_Line
    use g_Prtcl_DefaultValues
    use g_error_handling
    
    implicit none
    
    !// stationay plane wall 
    type PlaneWall
        
        integer(IK) user_id  ! user supplied wall id
        integer(IK) wall_id  ! program generated wall id
        integer(IK) wall_prop_Type  ! property type of wall material
        logical bothSide     ! checking if both sides are active 
        real (RK)   d        ! d in implicit equation   
        type(real3) n        ! normal vector
        type(real3) P1       ! first point   
        type(real3) P2       ! second point  
        type(real3) P3       ! third point 
        type(real3) P4       ! fourth point
        type(p_line) L1      ! line connects P1 to P2
        type(p_line) L2      ! line connects P2 to P3   
        type(p_line) L3      ! line connects P3 to P4
        type(p_line) L4      ! line connects P4 to P1
        
    contains
        procedure:: PW_CreateWall
        procedure:: PW_CreateWall2
        generic::   CreateWall_nv           => PW_CreateWall, PW_CreateWall2
        
       
        procedure:: Extend_boundaries       => PW_Extend_boundaries
        procedure:: Parallel_Plane          => PW_Parallel_Plane
        procedure:: Mirror_Plane            => PW_Mirror_Plane
        procedure:: isInContact             => PW_isInContact
        procedure:: IsInPlane               => PW_IsInPlane
        procedure:: IsOnLines               => PW_IsOnLines
        procedure:: PointFromPlane          => PW_PointFromPlane    
	    procedure:: NearestPointOnPlane     => PW_NearestPointOnPlane
        procedure:: Line_Intersect          => PW_line_intersect
        procedure:: cut_wall                => PW_cut_wall
        procedure:: cut_wall_1point         => PW_cut_wall_1point
        procedure:: cut_wall_2points        => PW_cut_wall_2points
        procedure:: cut_wall_3points        => PW_cut_wall_3points
        procedure:: test_inactive           => PW_test_inactive
        procedure:: normal_dist             => PW_normal_dist
        procedure:: getPoint                => PW_getPoint
        procedure:: getNormal               => PW_getNormal
        procedure:: set_wallID              => PW_set_wallID
        
    end type !// stationary plane wall 
    
    
    !// moving plane wall 
    type, extends(PlaneWall):: mvng_PlaneWall
        logical     :: moving    = .false. ! if it linearly moves   
        logical     :: rotating  = .false. ! if it rotates  
        type(real3) :: trans_vel = zero_r3 ! translational velocity  
        real(RK)    :: rot_vel   = 0.0_RK  ! rotational velocity
        type(p_line):: rot_line            ! rotation axis  
        
    contains ! associated methods 
        
        procedure:: CreateWall              => mvng_PW_CreateWall
        procedure:: normal_dist_vel_type    => mvng_PW_normal_dist_vel_type
        procedure:: setVelocity             => mvng_PW_setVelocity
        procedure:: GetVelocity             => mvng_PW_GetVelocity
        procedure:: MoveWall                => mvng_PW_MoveWall
        procedure:: get_pointVel            => mvng_PW_get_pointVel 
    end type !// moving plane wall 
    
    
contains

!**********************************************************************
! Creating a plane wall with 4 specified vertices
!**********************************************************************
logical function PW_CreateWall( this, p1, p2, p3, p4 , user_id , wall_id, prop_type, both )
    implicit none   
    class(PlaneWall)          this
    type(real3),intent(in) :: p1, p2, p3, p4
    integer(IK),intent(in) :: user_id, wall_id, prop_type
    logical    ,intent(in) :: both
    
    !// locals
    real(RK) val
    type(real3) ln
    
    !// body
    
    ! assignments 
	this%P1 = p1
	this%P2 = p2
	this%P3 = p3
	this%P4 = p4
	
    this%wall_id = wall_id
    this%user_id = user_id
    this%wall_prop_Type = prop_type
    this%bothSide = both
    
    !normal vector of plane
    ln = (p2-p1) .cross. (p3-p1)
    val = norm( ln ) ! norm of the vector 
	this%n = ln/val  ! unit normal vector
    
    ! d in the implicit equation: ax+by+cz+d = 0
    this%d = - (this%n .dot. p1)

	!checking whether the fourth point lies on the plane!
    if( abs( (this%n .dot. p4) + this%d ) < 0.00001_RK ) then
		PW_CreateWall = .true.
	else
		PW_CreateWall = .false.
        return
    end if
    ! parametric lines 
    call this%L1%setLine(p1,p2) ! p1 --> p2
    call this%L2%setLine(p2,p3) ! p2 --> p3
    call this%L3%setLine(p3,p4) ! p3 --> p4
    call this%L4%setLine(p4,p1) ! p4 --> p1
    
end function

logical function PW_CreateWall2( this, p1, p2, p3, p4  )
    implicit none   
    class(PlaneWall)          this
    type(real3),intent(in) :: p1, p2, p3, p4
    
    PW_CreateWall2 = this%CreateWall_nv(p1,p2,p3,p4, 0, 0, 0, .false.)
    
end function

!**********************************************************************
!extending boundaries of the wall by dx 
!**********************************************************************
function PW_Extend_boundaries(this, dx ) result(res)
    class(PlaneWall)      this
    real(RK),intent(in):: dx
    type(PlaneWall)       res
    
    !// locals
    integer(IK) stat
    real(RK) t
    type(real3) P_1, P_2, P_3, P_4
    type(real3) res_p1, res_p2, res_p3, res_p4
    type(p_line) L1, L2, L3, L4, L14, L12, L23, L34
    
    ! We use barycentric coordinates to extend quadrilateral boundaries by dx
    ! The parametric equation of a plane in the barycentric form is:
    ! P = P1 + v*(P2-P1) + w*(P3-P1) + x * (P4-P1)
    !
       
    ! first finding L1, P1-P2
    call L14%setLine(this%P1,this%P4)    
    call L23%setLine(this%P2,this%P3)
    call L12%setLine(this%P1,this%P2)
    call L34%setLine(this%P3,this%P4)
    
    ! creating Line 1, P1, P2
    t = -dx/L14%get_length()
    P_1 = L14%getpoint(t)
    
    t = -dx/L23%get_length()
    P_2 = L23%getpoint(t)
    
    call L1%setLine(P_1,P_2)
    
    
    ! creating Line 2, P2, and P3
    t = 1 + dx/L12%get_length()
    P_2 = L12%getPoint(t)
    
    t = -dx/L34%get_length()
    P_3 = L34%getPoint(t)
    call L2%setLine(P_2,P_3)
    
    ! creating Line 3, P3 and P4
    t = 1 + dx/L23%get_length()
    P_3 = L23%getPoint(t)
    
    t = 1+dx/L14%get_length()
    P_4 = L14%getPoint(t)
    
    call L3%setLine(P_3, P_4)
    
    ! creating Line 4, P4, P1
    t = 1 + dx/L34%get_length()
    P_4 = L34%getPoint(t)
    
    t = - dx/L12%get_length()
    P_1 = L12%getPoint(t)
    
    call L4%setLine(P_4,P_1)
    
    
    ! Intersection of these lines creates 4 vertices of the plane
    res_p1 = L1%find_intersect(L4, stat)
    if( stat == -1 ) then
        call CheckForError( ErrT_Abort, &
                            "PW_Extend_boundaries", &
                            "Cannot find the intersect 1:"// num2str(this%wall_id) )
        return
    end if
    res_p2 = L1%find_intersect(L2, stat)
     if( stat == -1 ) then
        call CheckForError( ErrT_Abort, &
                            "PW_Extend_boundaries", &
                            "Cannot find the intersect 2:"// num2str(this%wall_id) )
        return
    end if
    res_p3 = L2%find_intersect(L3, stat)
     if( stat == -1 ) then
        call CheckForError( ErrT_Abort, &
                            "PW_Extend_boundaries", &
                            "Cannot find the intersect 3:"// num2str(this%wall_id) )
        return
     end if
    
     if( this%L3%get_vect_n() .equal. this%L4%get_vect_n() ) then
         res_p4 = 0.5_RK*(res_p3 + res_p1)
     else
         
        res_p4 = L3%find_intersect(L4, stat)
        if( stat == -1 ) then
            call CheckForError( ErrT_Abort, "PW_Extend_boundaries", "Cannot find the intersect 4:"// num2str(this%wall_id) )
            return
        end if
        
    end if
    
     ! creating a new wall
    if( .not.(res%CreateWall_nv(res_p1, res_p2, res_p3, res_p4 , &
              this%user_id, this%wall_id, this%wall_prop_Type, this%bothSide)) ) then
        call CheckForError( ErrT_Abort, &
                            "PW_Extend_boundaries", &
                            "Cannot extend wall boundaries: "// num2str(this%wall_id) ) 
    end if
    
    
end function


!**********************************************************************
!* Creating a plane which is parallel to this plane and whose distance
!* is dx. 
!**********************************************************************
function PW_Parallel_Plane(this, dx) result (res)
    implicit none
    class(PlaneWall)      this
    real(RK),intent(in):: dx
    type(PlaneWall)       res
    
    ! locals
    real(RK) t
    type(real3) p1, p2, p3, p4
    
    !// body
    
    t = dx
    
    p1 = this%p1 + this%n * t
    p2 = this%p2 + this%n * t
    p3 = this%p3 + this%n * t
    p4 = this%p4 + this%n * t
    
    ! creating parallel plane
    if( .not.(res%CreateWall_nv(p1, p2, p3, p4 , this%user_id ,this%wall_id, this%wall_prop_Type, this%bothSide)) ) then
        call CheckForError( ErrT_Abort, "PW_Parallel_Plane", "Cannot create wall boundaries in PW_Parallel_Plane")
        return 
    end if
    

end function 


subroutine PW_Mirror_Plane(this)
    implicit none
    class(PlaneWall) this
    
    type(real3) p1, p2, p3, p4
    p1 = this%P3
    p2 = this%P2
    p3 = this%P1
    p4 = this%P4
    
    if( .not. (this%CreateWall_nv( p1, p2, p3, p4, &
                                   this%user_id, this%wall_id,         &
                                   this%wall_prop_type, this%bothSide) ) ) then
        call CheckForError( ErrT_Abort, "PW_Mirror_Plane", "Cannot find a mirror wall:"// num2str(this%wall_id) )
        return
    end if
                                   
    
end subroutine 


!**********************************************************************
! determining if an sphere (box) has a contact with this wall
!**********************************************************************
logical function PW_isInContact(this, box )
    implicit none
    class(PlaneWall)         this
    type(real4),intent(in):: box
    
    !// locals
    real(RK)    rad, dist, ovrlp
    type(real3) p
    
    !// body
    rad = 0.5_RK * box%w
    p = box ! conversion from real4 to real3
    PW_isInContact = .false.
    
    ! first checking if there is an overlap between sphere and wall
    dist = this%PointFromPlane(p)
    if( this%bothside) then
        ovrlp = rad - abs(dist)
    else
        if( dist < 0.0_RK ) return
        ovrlp = rad - dist
    end if
    
    if( ovrlp .gt. 0.0_RK) then
        ! now, checking if the contact point is located within the plane boundaries
        if( this%IsInPlane( this%NearestPointOnPlane(p) ) ) then
            PW_isInContact = .true.
            return
        end if
        !checking if the particle has a contact with edges or corner points 
        if( this%IsOnLines(box) ) then
            PW_isInContact = .true.
            return
        end if   
        
    endif
    
end function

!**********************************************************************
! checking if the point lays within the boundaries of the plane
!**********************************************************************
logical function PW_IsInPlane( this, p )
    implicit none
    class(PlaneWall)         this
	type(real3),intent(in):: p
        
    !// locals
    real(RK) p1p2, p2p3, p2p2, p1p3, p2p4, p1p4, p1p1
    type(real3) p1p, p2p, p3p, p4p
        
    !// body    
    
    PW_IsInPlane = .false.    
    p1p = this%P1-p
	p2p = this%P2-p
	p3p = this%P3-p

	! first condition u.v < 0
	! u.v = [(p1-p)x(p2-p)].[(p2-p)x(p3-p)] = (p1p.p2p)(p2p.p3p) - (p2p.p2p)(p1p.p3p) 
	p1p2 = p1p .dot. p2p
	p2p3 = p2p .dot. p3p
	p2p2 = p2p .dot. p2p
	p1p3 = p1p .dot. p3p
	if( p1p2*p2p3 - p2p2*p1p3 < 0.0_RK ) return
	
	! second condition u.w < 0
	! u.w = [(p1-p)x(p2-p)].[(p3-p)x(p4-p)] = (p1p.p3p)(p2p.p4p) - (p2p.p3p)(p1p.p4p)
	p4p = this%P4-p
	p2p4 = p2p .dot. p4p
	p1p4 = p1p .dot. p4p
	if( p1p3*p2p4 - p2p3*p1p4 < 0.0_RK ) return

	! third condition u.x < 0
	! u.x = [(p1-p)x(p2-p)].[(p4-p)x(p1-p)] = (p1p.p4p)(p2p.p1p) - (p2p.p4p)(p1p.p1p)
	p1p1 = p1p .dot. p1p
	if( p1p4*p1p2 - p2p4*p1p1 < 0.0_RK ) return 
	
    PW_IsInPlane = .true.
	return
    
 
        
end function


logical function PW_IsOnLines(this, dpos, nv, cp, dist )
    implicit none
    class(PlaneWall) this
    type(real4),intent(in)  :: dpos
    type(real3),intent(out) :: nv, cp
    real(RK),   intent(out) :: dist
    optional:: nv, cp, dist
    
    !// locals
    real(RK) ldist
    type(real3) lcp, lnv
    
    
    
    if( Line_point_check(this%L1, dpos, lnv, lcp, ldist) ) goto 100
    if( Line_point_check(this%L2, dpos, lnv, lcp, ldist) ) goto 100    
    if( Line_point_check(this%L3, dpos, lnv, lcp, ldist) ) goto 100
    if( Line_point_check(this%L4, dpos, lnv, lcp, ldist) ) goto 100    
    
    
    PW_IsOnLines = .false.
    return
    
100 continue
    
    if( present(nv) ) nv = lnv
    if( present(cp) ) cp = lcp
    if( present(dist)) dist = ldist
    PW_IsOnLines = .true.
    return
end function


logical function Line_point_check( line, dpos, nv, cp , dist )
    implicit none
    type(p_line),intent(in) :: line
    type(real4), intent(in) :: dpos
    type(real3), intent(out):: nv, cp
    real(RK),    intent(out):: dist
    
    
    real(RK) t, length , r
    type(real3) pos
    
    length = line%get_length();
    
    pos = dpos
    r = 0.5_RK * dpos%w
    t = line%get_Project(pos)
    
    Line_point_check = .false.
    
    if( t>= 0.0_RK .and. t<= 1.0_RK ) then
        cp = line%getPoint(t)
    elseif( t >= (-r/length)  .and. t <0.0_RK )then
        cp = line%getPoint(0.0_RK)
    elseif( t> 1.0  .and. t>= (1.0+r/length) ) then
        cp = line%getPoint(1.0_RK)
    else
        Line_point_check = .false.
        return
    end if
    
    dist = cp .dist. pos
    
    if( (r - dist) > 0 )then
        nv = pos .nv. cp   
        Line_point_check = .true.
        return 
        
    end if
    
    ! cp --> pos
        
    
end function

!**********************************************************************
! normal distance of a point from wall
!**********************************************************************
real(RK) function PW_PointFromPlane( this,  p)
    implicit none
    class(PlaneWall)          this
    type(real3),intent(in) :: p
    
    PW_PointFromPlane = (this%n .dot. p) + this%d
    
end function

!**********************************************************************
! Nearest point (normal projection point) of point p on the
! palne
!**********************************************************************
type(real3) function PW_NearestPointOnPlane( this,  p )
    implicit none
    class(PlaneWall)         this
    type(real3),intent(in):: p
    
    !// locals
    real(RK) t
    
    !// body	
	t = - ( (this%n .dot. p)+ this%d )
    PW_NearestPointOnPlane = t * this%n + p
		
end function

!**********************************************************************
! Finding the intersection point (q) of a line with plane. If there is no 
! intersect, the function returns .false., otherwise it returns .true.
!  NOTE: plane is considered to be infinite.
!**********************************************************************
logical function PW_line_intersect( this, line, q )
    implicit none
    class(PlaneWall)    this
    type(p_line),intent(in) :: line
    type(real3),     intent(out):: q
    
    real(RK)    t, den
    type(real3) v, A
    
    ! with line equation as S(t) = vt+A where v = B-A
    ! and plane equation as nX+d = 0, where n is normal vector of plane
    ! intersection of plane and line is 
    ! t = - (d+n.A)/(n.v)
    A = line%getPoint(0.0_RK);
    v = line%get_vect();
    den = this%n .dot. v
    
    
    if( abs(den) .lt. 0.00000001_RK ) then
        ! parallel line and plane
        PW_line_intersect = .false.
        return
    end if
    
    t = -( this%d + (this%n .dot. A))/ den
    
    if( t>= 0.0_RK .and. t<= 1.0_RK ) then
        ! there is an intersection 
        q = A + (v*t)
        PW_line_intersect = .true.
        return       
    end if
        
    PW_line_intersect = .false.
    return

end function

subroutine PW_cut_wall_3points(this, plane, pw1_points )
    class(PlaneWall)                       this
    class(PlaneWall),intent(in)         :: plane
    type(real3),dimension(4),intent(out):: pw1_points
    
    !// locals
    logical res
    type(real3) q
    type(real3) p1, p2, p3, p4, p5, p6
    integer(IK),dimension(4) :: off_points
    
    off_points = this%test_inactive( plane )
    
    if( off_points(1) == 0 ) then
        ! lines l1 and l4
        p1 = this%p1
        res = plane%line_intersect(this%L1, q)
        p2 = q
        res = plane%line_intersect(this%L4, q)
        p3 = q 
        p4 = 0.5_RK*(p2+p3)
    elseif( off_points(2) == 0 ) then
        
        ! lines l1 and l2
        res = plane%line_intersect(this%L1, q)
        p1 = q
        p2 = this%p2
        res = plane%line_intersect(this%L2, q)
        p3 = q 
        p4 = 0.5_RK*(p2+p3)
    elseif( off_points(3) == 0 ) then
        
        ! lines l2 and l3
        res = plane%line_intersect(this%L2, q)
        p1 = q
        p2 = this%p3
        res = plane%line_intersect(this%L3, q)
        p3 = q 
        p4 = 0.5_RK*(p2+p3)
        
    elseif( off_points(4) == 0 ) then
        
        ! lines l3 and l4
        res = plane%line_intersect(this%L3, q)
        p1 = q
        p2 = this%p4
        res = plane%line_intersect(this%L4, q)
        p3 = q 
        p4 = 0.5_RK*(p2+p3)
    else
        call CheckForError( ErrT_Abort, "PW_cut_wall_3points", "Cannot cut the wall" )
        return
    endif
    
    pw1_points = (/p1,p2,p3,p4/)
    
end subroutine

subroutine PW_cut_wall_1point(this, plane, pw1_points, pw2_points )
    class(PlaneWall)                       this
    class(PlaneWall),intent(in)         :: plane
    type(real3),dimension(4),intent(out):: pw1_points, pw2_points
    
    !// locals
    logical res
    type(real3) q
    type(real3) p1, p2, p3, p4, p5, p6
    integer(IK),dimension(4) :: off_points
    
    !// body
    
    off_points = this%test_inactive( plane )
    
    if( off_points(1) == 1_IK ) then
       ! lines L1 and L4 
       res = plane%line_intersect(this%L1, q)
       p1 = q
       res = plane%line_intersect(this%L4, q)
       p5 = q
       
       p2 = this%p2
       p3 = this%p3
       p4 = this%p4
       p6 = 0.5_RK * (p5+p1)
       
    elseif( off_points(2) == 1_IK )then
       
       ! lines L1 and L2 
       res = plane%line_intersect(this%L1, q)
       p2 = q
       res = plane%line_intersect(this%L2, q)
       p3 = q
       p1 = this%p1
       p4 = this%p3
       p5 = this%p4
       p6 = 0.5_RK*(p5+p1)
    
    elseif(off_points(3) == 1 ) then
       ! lines L2 and L3 
       res = plane%line_intersect(this%L2, q)
       p3 = q
       res = plane%line_intersect(this%L3, q)
       p4 = q
       p1 = this%p1
       p2 = this%p2
       p5 = this%p4
       p6 = 0.5_RK*(p5+p1)

    elseif( off_points(4) == 1) then
      
        ! lines L3 and L4 
       res = plane%line_intersect(this%L3, q)
       p4 = q
       res = plane%line_intersect(this%L4, q)
       p5 = q
       p1 = this%p1
       p2 = this%p2
       p3 = this%p3
       p6 = 0.5_RK*(p5+p1)
       
    end if
    
    pw1_points = (/p1,p2,p3,p4/)
    pw2_points = (/p1,p4,p5,p6/)
    
end subroutine 


subroutine PW_cut_wall_2points( this, plane, pw1_points ) 
    implicit none
    class(PlaneWall) this
    class(PlaneWall),intent(in)         :: plane
    
    type(real3),dimension(4),intent(out):: pw1_points
    
    type(real3) q, p1, p2, p3, p4
    p1 = this%p1
    p2 = this%p2
    p3 = this%p3
    p4 = this%p4
    
        ! L1: line that connects P1 to P2
    if( plane%line_intersect(this%L1, q) ) then
        
        if( plane%PointFromPlane(p1) <= 0.0_RK ) then
            p1 = q
        else
            p2 = q
        end if
        
    end if
    
    ! L2: line that connects P2 to P3
    if( plane%line_intersect(this%L2, q) ) then
        
        if( plane%PointFromPlane(p2) <= 0.0_RK ) then
            p2 = q
        else
            p3 = q
        end if
        
    end if
    
    ! L3: line that connects P3 to P4
    if( plane%line_intersect(this%L3, q) ) then
        
        if( plane%PointFromPlane(p3) <= 0.0_RK ) then
            p3 = q
        else
            p4 = q
        end if
        
    end if
    
    ! L4: line that connects P4 to P1
    if( plane%line_intersect(this%L4, q) ) then
        
        if( plane%PointFromPlane(p4) <= 0.0_RK ) then
            p4 = q
        else
            p1 = q
        end if
        
    end if
    
    pw1_points= (/p1,p2,p3,p4/)
    
end subroutine


!**********************************************************************
! intersecting a wall with infinite plane (plane) and keeping the part 
! of the wall which is in the positive side of plane. 
!  NOTE: plane is considered to be infinite.
!**********************************************************************
subroutine PW_cut_wall(this, plane )
    implicit none
    class(PlaneWall) this
    class(PlaneWall),intent(in) :: plane
    
    !// locals
    real(RK) dist
    type(real3) q, p1, p2, p3, p4
    p1 = this%p1
    p2 = this%p2
    p3 = this%p3
    p4 = this%p4
    
    ! L1: line that connects P1 to P2
    if( plane%line_intersect(this%L1, q) ) then
        
        if( plane%PointFromPlane(p1) < 0.0_RK ) then
            p1 = q
        else
            p2 = q
        end if
        
    end if
    
    ! L2: line that connects P2 to P3
    if( plane%line_intersect(this%L2, q) ) then
        
        if( plane%PointFromPlane(p2) < 0.0_RK ) then
            p2 = q
        else
            p3 = q
        end if
        
    end if
    
    ! L3: line that connects P3 to P4
    if( plane%line_intersect(this%L3, q) ) then
        
        if( plane%PointFromPlane(p3) < 0.0_RK ) then
            p3 = q
        else
            p4 = q
        end if
        
    end if
    
    ! L4: line that connects P4 to P1
    if( plane%line_intersect(this%L4, q) ) then
        
        if( plane%PointFromPlane(p4) < 0.0_RK ) then
            p4 = q
        else
            p1 = q
        end if
        
    end if
    
    ! now creates wall from new planes
    
    
    if( .not. this%CreateWall_nv(p1,p2,p3,p4, this%user_id, this%wall_id, this%wall_prop_Type, this%bothSide) ) then
        call CheckForError( ErrT_Abort, &
                            "PW_cut_wall" , &
                            "Cannot create wall after performing cut operation on wall, wall id is :"&
                             // num2str( this%user_id) )    
    end if 
    
    
end subroutine

!***********************************************************************************
! determining number of points of this wall which are in the off-side (inactive side)
! of plane 
! NOTE: the plane is considered to be infinite
!***********************************************************************************
function PW_test_inactive( this, plane) result (res)
    implicit none
    class(PlaneWall) this
    class(PlaneWall),intent(in):: plane
    
    integer(IK),dimension(4):: res
    !// locals
    integer(IK) numP ! number of points that are in inactive side of plane
    
    !// body
    res = 0
    if( plane%PointFromPlane(this%P1) <= 0.000001_RK )  res(1) = 1
    if( plane%PointFromPlane(this%P2) <= 0.000001_RK )  res(2) = 1
    if( plane%PointFromPlane(this%P3) <= 0.000001_RK )  res(3) = 1
    if( plane%PointFromPlane(this%P4) <= 0.000001_RK )  res(4) = 1
    
end function

!**********************************************************************
! returning normal vector of plane and normal distance between center of 
! sphere(box) and plane  
!**********************************************************************
subroutine PW_normal_dist(this, box ,nv, dist, pt , cp )
    implicit none
    class(PlaneWall) this
    type(real4),intent(in)  :: box
    type(real3),intent(out) :: nv
    real(RK),intent(out)    :: dist
    integer(IK),intent(out) :: pt
    type(real3),intent(out) :: cp
    
    !// locals
    logical res
    type(real3) p
    
    !// body
    
    p = box
    cp = this%NearestPointOnPlane(p) 
    if( this%IsInPlane( cp ) ) then
        
        dist = this%PointFromPlane(p)
        
        if( this%bothside) then
            nv = sign(1.0_RK, dist) *this%getNormal()
            dist = abs(dist)
        else
            nv = this%getNormal()
        end if 
        
    else
        
        
       res = this%IsOnLines(box , nv, cp, dist )
        
        
    end if
    
    
    pt = this%wall_prop_Type
    
end subroutine

!
! returning the specified point of the wall
type(real3) function PW_getPoint(this , n )
    implicit none
    class(PlaneWall)         this
    integer(IK),intent(in):: n
    
    select case (n)
    case(1)
        PW_getPoint = this%p1
    case(2)
        PW_getPoint = this%p2
    case(3)
        PW_getPoint = this%p3
    case(4)
        PW_getPoint = this%p4
    end select
    
end function
    

! returning normal vector of the plane
type(real3) function PW_getNormal(this)
    implicit none
    class(PlaneWall) this
    
    PW_getNormal = this%n
    
    end function
    
subroutine PW_set_wallID(this , id )
    implicit none
    class(PlaneWall) this
    integer(IK),intent(in):: id
    
    this%wall_id = id
end subroutine

!**********************************************************************
! Creating a moving plane wall with 4 specified vertices, translational
! velocity, rotational velocity and rotational axis
!**********************************************************************    
logical function mvng_PW_CreateWall( this, p1, p2, p3, p4 , user_id ,wall_id, prop_type, t_vel , r_vel, r_line, both )
    implicit none   
    class(mvng_PlaneWall) this
    type(real3),intent(in):: p1, p2, p3, p4
    integer(IK),intent(in):: user_id, wall_id, prop_type
    logical,    intent(in):: both
    type(real3),intent(in):: t_vel
    real(RK),   intent(in):: r_vel
    type(p_line),intent(in):: r_line
    
    if( .not. this%CreateWall_nv(p1, p2, p3, p4, user_id ,wall_id, prop_type, both )  ) then
        mvng_PW_CreateWall = .false.
        return
    endif
    
    call this%setVelocity(t_vel , r_vel, r_line)
    mvng_PW_CreateWall = .true.
    
end function

!**********************************************************************************    
! returning normal vector of wall, normal distance of sphere center from plane,
! velocity of wall at projection point of sphere center
!**********************************************************************************
subroutine mvng_PW_normal_dist_vel_type(this, box ,nv, dist, vel , pt )
    implicit none
    class(mvng_PlaneWall) this
    type(real4),intent(in)  :: box
    type(real3),intent(out) :: nv
    real(RK),   intent(out) :: dist
    type(real3),intent(out) :: vel
    integer(IK),intent(out) :: pt
    
    !// locals
    type(real3) p, cnt_p
    
    !// body
    
    ! normal and distance
    call this%normal_dist( box, nv, dist, pt, cnt_p )
    
    if( (.not. this%moving) .and. (.not. this%rotating) ) then
        vel = zero_r3
        return
    end if    
    
    ! getting wall velocity at contact point
    vel = this%GetVelocity( cnt_p ) 

end subroutine

subroutine mvng_PW_SetVelocity(this, t_vel , r_vel, r_line )
    implicit none
    class(mvng_PlaneWall) this
    type(real3),intent(in):: t_vel
    real(RK),   intent(in):: r_vel
    type(p_line),intent(in):: r_line
    
    if( norm(t_vel) > 0.00000001_RK ) then
        this%moving = .true.        
    end if
    
    if( abs(r_vel) > 0.00000001_RK ) then
        this%rotating = .true.    
    endif
    
    this%trans_vel = t_vel
    this%rot_vel   = r_vel
    this%rot_line = r_line
    
end subroutine

!**********************************************************************
! returning the velocity of wall at the specified point on the wall
!**********************************************************************
type(real3) function mvng_PW_GetVelocity(this, p0 )
    implicit none
    class(mvng_PlaneWall)   this
    type(real3),intent(in)::p0
    
    type(real3) w_vel
    
    
    w_vel = zero_r3
    if( this%moving ) then
        w_vel = this%trans_vel
    end if
    
    if( this%rotating) then
        w_vel = w_vel + this%get_pointVel(p0)    
    end if
    
    mvng_PW_GetVelocity = w_vel
    
end function
    
!**********************************************************************    
! moving wall in the space based on its velocity    
!**********************************************************************
logical function mvng_PW_MoveWall(this, dt) 
    implicit none
    class(mvng_PlaneWall)   this
    real(RK),intent(in)::   dt
    
    !//locals
    real(RK) d_teta
    type(real3) p1, p2, p3, p4, v1, v2, v3, v4
    
    !// body
    
    ! if wall is not moving, it returns true
    if( ( .not. this%moving ) .and. ( .not. this%rotating) ) then
        mvng_PW_MoveWall = .true.
        return
    endif
    
    ! finding new positions of four vertices of the wall
    v1 = this%getVelocity( this%p1 )
    v2 = this%getVelocity( this%p2 )
    v3 = this%getVelocity( this%p3 )
    v4 = this%getVelocity( this%p4 )
    
    p1 = v1*dt + this%p1
    p2 = v2*dt + this%p2
    p3 = v3*dt + this%p3
    p4 = v4*dt + this%p4
    
    !Creating the wall again based on the new positions of vertices 
    mvng_PW_MoveWall = this%CreateWall( p1, p2, p3, p4, this%user_id, this%wall_id, &
                                   this%wall_prop_Type , this%trans_vel, &
                                   this%rot_vel, this%rot_line , this%bothSide)
    
end function

!**********************************************************************
! returning the rotational velocity of the wall at specified point
! on the wall
!**********************************************************************
type(real3) function mvng_PW_get_pointVel( this, p0 )
    implicit none
    class(mvng_PlaneWall)   this
    type(real3),intent(in)::p0
        
    !locals        
    real(RK) d_teta
    type(real3) p1
    
    ! rotating the point from current position at "t" to the next new angle at "t+dt"
    d_teta =  this%rot_vel * d_dt_wallVel
    p1 = this%rot_line%rotate_point(p0 , d_teta)
    
    ! now, calculating the velocity at p0 by taking numerical derivation
    mvng_PW_get_pointVel = (p1 - p0)/d_dt_wallVel
    
end function
    
   
end module

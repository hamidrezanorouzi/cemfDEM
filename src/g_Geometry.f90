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
!  file name  : g_Geometry.f90                          
!  module name: g_Geometry                              
!                                                       
!  Purpose:                                             
!   Managing all the actions related to the geometry, like adding walls
! 
!------------------------------------------------------------------------------

module g_Geometry
    
    
    use g_PlaneWall
    use g_CylinderWall
    use g_stl_reader
    
    
    implicit none
    
    
    integer(IK),parameter:: WT_PlaneWall = 10
    integer(IK),parameter:: WT_Cylinder  = 11
    
    type wallType
        integer(IK) WT
        integer(IK) Wall_id
        integer(IK) indx
    end type
    
    
    
    type Geometry
        
        ! number of plane walls
        integer(IK):: num_pWall = 0   
        
        ! maximum number of plane walls
        integer(IK) :: max_num_pWall = 0
        
        ! the last wall id added
        integer(IK) :: last_wall_id  = 0
        
        ! a vector that stores all plane walls 
        type(mvng_PlaneWall),allocatable,dimension(:):: pWall
        
        !
        logical :: SeparateOutput    = .false.   ! determining if program creates separate output files
        integer(IK):: numWallOut     = 0         ! number of files for output 
        character(32),dimension(:),allocatable:: outName    ! file names 
        integer(IK),dimension(:),allocatable  :: out_userid ! wall user ids to be written to each file 
        
    contains
    
    ! Initializing the geometry object
    procedure:: Initialize          => G_Initialize
    
    ! adding a plane wall to the geometry object
    procedure:: G_add_PlaneWall_nv
    procedure:: G_add_PlaneWall
    procedure:: G_add_PlaneWall_wall
    generic:: add_PlaneWall         => G_add_PlaneWall, G_add_PlaneWall_nv, G_add_PlaneWall_wall
    
    !adding a cylinderical wall to the geometry object
    procedure:: add_Cylinder        => G_add_Cylinder
    
    ! adding plane walls which are in the stl file
    procedure:: G_add_stl_file1
    procedure:: G_add_stl_file2
    generic:: add_stl_file          => G_add_stl_file1, G_add_stl_file2
    
    ! deleting a wall by user_id as input
    procedure:: delete_wall         => G_delete_wall
    
    ! setting wall velocity by user_id
    procedure:: setWallVelocity     => G_setWallVelocity
    
    ! returning number of plane walls
    procedure:: get_num_pWall       => G_get_num_pWall
    
    ! returning maximum number of plane walls
    procedure:: get_max_pWall       => G_get_max_pWall 
    
    ! returning the ith plane wall in the geometry
    procedure:: get_pWall           => G_get_pWall
    
    ! moving all walls (if they have velocity) 
    procedure:: move_walls          => G_move_walls
    
    ! determining if a sphere is in contact with specified wall in the system
    procedure:: isInContact         => G_isInContact
    
    ! normal vector and property type of wall and distance of a box from the wall 
    procedure:: normal_dist         => G_normal_dist
    
    ! normal vector and property type of wall, distance of a box from the wall and wall velocity at contact point 
    procedure:: normal_dist_vel_type=> G_normal_dist_vel_type
    
    ! setting output file name for separate geometry parts 
    procedure:: setWallOutputName   => G_setWallOutputName
    
    ! functions to create output file with VTK format 
    procedure:: GeomToVTK           => G_GeomToVTK
    procedure:: GeomToVTK_all       => G_GeomToVTK_all
    procedure:: GeomToVTK_user_id   => G_GeomToVTK_user_id
    
    
    procedure,private:: Increase_pWall
    procedure,private:: G_get_wall_index
    procedure,private:: G_delete_wall_id
    procedure,private:: G_GetWallName
    
    end type
    
    
contains


!**********************************************************************
! Initializing the geometry object 
!**********************************************************************
subroutine G_Initialize(this)
    implicit none
    class(Geometry) this
    
        
    ! it does nothing   
    
        
end subroutine

!**********************************************************************
! adding a plane wall to the geometry object
!**********************************************************************
subroutine G_add_PlaneWall_nv(this, p1, p2, p3, p4 , user_id ,prop_type, both  )
    implicit none
    class(Geometry) this
    type(real3),intent(in) :: p1, p2, p3, p4 ! corner points 
    integer(IK),intent(in) :: user_id, prop_type
    logical,intent(in)     :: both           ! both side active status
    optional both
    
    ! locals
    integer(IK) id
    logical     lboth
    type(mvng_PlaneWall) wall
    type(wallType)       code
    
    ! body 
    
    ! checking for both side spec. 
    if( present(both) )then
        lboth = both
    else
        lboth = .false.
    endif
    
    ! Id of wall will be assigned one-by-one 
    ! they are inserted into the geometry in sequence 
    id = 0 !
    
    
    if( .not. wall%CreateWall_nv(p1,p2,p3,p4 , user_id ,id , prop_type, lboth) ) then
        ! error in creating the wall, program will abort
        call CheckForError( ErrT_Abort , "G_add_PlaneWall_nv" , "Cannot create a plane wall, wall No. "//num2str(user_id ) )
        
    end if
    
    call this%add_PlaneWall(wall)
     

end subroutine

subroutine G_add_PlaneWall(this, p1, p2, p3, p4 , user_id ,prop_type, t_vel , r_vel, r_line , both  )
    implicit none
    class(Geometry) this
    type(real3),intent(in) :: p1, p2, p3, p4   ! corner points 
    integer(IK),intent(in) :: user_id, prop_type
    logical,intent(in)     :: both             !  both side active status
    optional both
    type(real3),intent(in) :: t_vel            ! translational velocity
    real(RK),intent(in)    :: r_vel            ! rotational velocity
    type(p_line),intent(in) :: r_line          ! axis of rotation
    
    !// locals
    integer(IK) id
    logical     lboth
    type(mvng_PlaneWall) wall
    
    
    ! //body 
    
    ! checking for both side spec. 
    if( present(both) )then
        lboth = both
    else
        lboth = .false.
    endif
    
    ! Id of wall will be assigned one-by-one 
    !  they are inserted into the geometry in sequence
    id = 0 
    
    
    if( .not. wall%CreateWall(p1,p2,p3,p4 , user_id ,id , prop_type, t_vel , r_vel, r_line, lboth) ) then
        ! error in creating the wall, program will abort
        call CheckForError( ErrT_Abort , "G_add_PlaneWall" , "Cannot create a plane wall, wall No. "//num2str(user_id ) )
        
    end if
    
    call this%add_PlaneWall(wall)
        

end subroutine

!****************************************************************************
! adding a plane wall 
!****************************************************************************
subroutine G_add_PlaneWall_wall(this, wall)
    implicit none
    class(Geometry) this
    type(mvng_PlaneWall),intent(in):: wall
    
    
    ! increasing the size of vector
    ! if there is not enough space to store new added wall,
    if( this%num_pWall >= this%max_num_pWall) then
        !extend sizes
        call this%Increase_pWall()    
    end if
    
    ! adding it to the vector of plane wall
    this%num_pWall = this%num_pWall  + 1 
    this%last_wall_id = this%last_wall_id + 1
    call wall%set_wallID( this%last_wall_id  + d_Base_wall_id )
    
    this%pWall(this%num_pWall) = wall
    
    
    
end subroutine

!**********************************************************************
! adding a cylinder to the geometry object
!**********************************************************************
subroutine G_add_Cylinder( this , cylinder )
    implicit none
    class(Geometry) this
    type(CylinderWall),intent(in):: cylinder
    
    !// locals
    integer(IK) n, nw
    
    nw = cylinder%get_numWall()
    
    do n=1, nw
        call this%add_PlaneWall( cylinder%get_Wall(n) )
    end do
    
end subroutine 

!****************************************************************************
! adding a wall from an stl file with predefined velocities 
!****************************************************************************
subroutine G_add_stl_file1(this, chFile, user_id ,prop_type, t_vel , r_vel, r_line , both )
    implicit none
    class(Geometry)           this
    character(*),intent(in):: chFile
    integer(IK),intent(in) :: user_id, prop_type
    type(real3),intent(in) :: t_vel
    real(RK),intent(in)    :: r_vel
    type(p_line),intent(in) :: r_line
    logical,intent(in)     :: both
    optional both
    
    !// locals
    integer(IK) numFacets, n
    logical lboth
    character(255) ch_error
    type(real3) p4
    type(p_line) L31
    type(facet),dimension(:),allocatable:: wall_facets
    
    
    !// body
    
    lboth = .false.
    if( present(both) ) lboth = both
    
    if( .not. read_from_stl_file( chFile, wall_facets, numFacets, ch_error ) )then
        
        call CheckForError( ErrT_Abort , &
                            "G_add_stl_file1" , &
                            "error in reading from stl file: " // trim(chFile)//". Error is: " // trim(ch_error) )
        return 
        
    end if

    do n=1, numFacets
        call L31%setLine(wall_facets(n)%p3, wall_facets(n)%p1)
        p4 = L31%getPoint( 0.5_RK )
        
        call this%add_PlaneWall( wall_facets(n)%p1, &
                                 wall_facets(n)%p2, &
                                 wall_facets(n)%p3, &
                                 p4 , user_id ,     &
                                 prop_type, t_vel , r_vel, r_line , lboth )
        
    end do
    
    
end subroutine 

!****************************************************************************
! adding a wall from an stl file
!****************************************************************************
subroutine G_add_stl_file2(this, chFile, user_id ,prop_type, both )
    implicit none
    class(Geometry)           this
    character(*),intent(in):: chFile
    integer(IK),intent(in) :: user_id, prop_type
    logical,intent(in)     :: both
    optional both
    
    !// locals
    logical lboth
    
    type(p_line) L
    
    !// body
    
    lboth = .false.
    if( present(both) ) lboth = both
    
    call this%add_stl_file(chFile, user_id, prop_type, zero_r3, 0.0_RK, L, lboth )    
    
end subroutine 

!********************************************************************
! deleting wall by user_id
!********************************************************************
subroutine G_delete_wall(this, user_id)
    implicit none
    class(Geometry) this
    integer(IK),intent(in):: user_id
    
    integer(iK) i
    
    
    i = 1
    do while( i <= this%num_pWall )
        if( this%pWall(i)%user_id == user_id ) then
            ! this wall should be deleted
            call this%G_delete_wall_id( i )
        else
            ! go to the next wall 
            i = i + 1
        end if
        
    end do
    
    
end subroutine

!********************************************************************
! setting the file name for walls, it works if the separate output is active
!********************************************************************
subroutine G_setWallOutputName(this, user_id, chOutputName )
    implicit none
    class(Geometry) this
    integer(IK),dimension(:) ,intent(in)::  user_id
    character(*),dimension(:),intent(in)::  chOutputName
    
    !// locals
    integer(IK):: i, numOut
    
    ! checking for size of vectors in the input
    if( size(user_id) .ne. size(chOutputName) ) then
        call CheckForError( ErrT_Abort, "G_setWallOutputName", &
        "Wrong number of elements for user_id and chOutputName vectors" )
    end if
    
    ! getting size of vector 
    numOut = size(user_id)
    
    ! allocating memory 
    if( allocated(this%outName) ) deallocate(this%outName)
    allocate( this%outName(numOut) )
    
    if( allocated(this%out_userid) ) deallocate(this%out_userid)
    allocate( this%out_userid(numOut) )
    
    ! saving user_id and output file names 
    do i= 1, numOut
        this%outName(i)    = chOutputName(i)
        this%out_userid(i) = user_id(i)
    end do
    
    this%numWallOut = numOut
    this%SeparateOutput = .true.
    
end subroutine 


!***********************************************************************
!* geometry to VTK file
!***********************************************************************
subroutine G_GeomToVTK_all( this, file_num, DEM_Options )
    implicit none
    class(Geometry) this
    integer(IK),intent(in):: file_num
    type(DEMS_Options),intent(in):: DEM_Options

    !// locals
    integer(IK) i, j, nw, nUnitFile
    integer(IK),dimension(:),allocatable:: faces
    type(real3),dimension(:),allocatable:: points
    character(256) ch, chFile, chMsg
    
    ! file name for VTK file
    write(ch,*) DEM_Options%base_Number + file_num
    chFile = trim(DEM_Options%Res_Dir)//"WallsFor"//trim(DEM_Options%RunName)//trim(adjustl(ch))//".vtk"
    
    ! Opening VTK file
    nUnitFile = d_tec_geom_unit
    open( unit =  nUnitFile , file = chfile )
    
    ! Header for VTK file
    chMsg = "Geometry for run name: " // trim(DEM_Options%RunName)   
    call vtk_header( nUnitFile, chMsg )
    
    ! getting number for plane walls
    nw = this%get_num_pWall ()
    
    if( nw < 1 ) return
    
    ! allocating memory
    if(allocated(points)) deallocate(points)
    allocate( points(nw*4) )
    if(allocated(faces)) deallocate(faces)
    allocate(faces(nw*5))
    
    
    ! starting to add walls
    do i = 1, nw
      
        do j= 1, 4
            
            points((i-1)*4+j) = this%pWall(i)%getPoint(j)
            
        end do
        faces((i-1)*5+1) = 4
        faces((i-1)*5+2) = (i-1)*4
        faces((i-1)*5+3) = (i-1)*4+1
        faces((i-1)*5+4) = (i-1)*4+2
        faces((i-1)*5+5) = (i-1)*4+3
        
    end do
    
    
    !// starting to write wall points in the file
	write( nUnitFile, "(A)" )"DATASET POLYDATA"
    write( nUnitFile, "(A8,I6,A)" )"POINTS  ", nw*4 , "  float"
    
    do i = 1, nw*4
        write(nUnitFile, *) points(i)    
    end do
    
    write( nUnitFile, "(A8, I5, A1, I6)") "POLYGONS ", nw , " " , 5*nw 
    do i =1, nw*5
        write( nUnitFile, "(I6)") faces(i)
    end do
    
    close(nUnitFile)
    
    deallocate(points)
    deallocate(faces)
    
end subroutine


!**************************************************************************************
!* getting file name associated with the user_id, it works if separated output is active
!**************************************************************************************
function G_GetWallName(this, user_id ) result(res)
    implicit none
    class(Geometry)     this
    integer(IK),intent(in)::user_id
    character(32)       res
    
    !//locals
    integer(IK) i
    
    do i = 1, this%numWallOut
        
        if( this%out_userid(i) == user_id )  then
            
            res = this%outName(i)
            return
        
        end if
        
    end do
    
    res = ""
    
end function

!*****************************************************************************
! geometry to VTK file 
!*****************************************************************************
subroutine G_GeomToVTK(this, file_num, DEM_Options)
    implicit none
    class(Geometry) this
    integer(IK),intent(in):: file_num
    type(DEMS_Options),intent(in):: DEM_Options
    
    !// locals
    integer(IK) i
    
    if(this%SeparateOutput) then
    
        do i=1, this%numWallOut
            
            call this%GeomToVTK_user_id( file_num, DEM_Options, this%out_userid(i) )
            
        end do
        
    else
        
        call this%GeomToVTK_all( file_num, DEM_Options )
        
    end if
    
end subroutine

!*******************************************************************************
!* geometry to vtk file based on user_id, it works if separate output file is active 
!*******************************************************************************
subroutine G_GeomToVTK_user_id( this, file_num, DEM_Options, user_id )
    implicit none
    class(Geometry) this
    integer(IK),intent(in):: file_num
    type(DEMS_Options),intent(in):: DEM_Options
    integer(IK),intent(in):: user_id

    !// locals
    integer(IK) i, j, nw, nUnitFile, ii
    integer(IK),dimension(:),allocatable:: faces
    type(real3),dimension(:),allocatable:: points
    character(256) ch, chFile, chMsg
    
    ! file name for VTK file
    write(ch,*) DEM_Options%base_Number + file_num
    chFile = trim(DEM_Options%Res_Dir)//"Wall-"//trim(this%G_GetWallName(user_id))// &
             "-"//trim(DEM_Options%RunName)//trim(adjustl(ch))//".vtk"
    
    ! Opening VTK file
    nUnitFile = d_tec_geom_unit
    open( unit =  nUnitFile , file = chfile )
    
    ! Header for VTK file
    chMsg = "Geometry for run name: " // trim(DEM_Options%RunName)//"-"//trim(this%G_GetWallName(user_id))  
    call vtk_header( nUnitFile, chMsg )
    
    ! getting number of plane walls
    nw = this%get_num_pWall ()
    
    if( nw < 1 ) return
    
    ! allocating memory
    if(allocated(points)) deallocate(points)
    allocate( points(nw*4) )
    if(allocated(faces)) deallocate(faces)
    allocate(faces(nw*5))
    
    
    ! starting to add walls
    ii = 0
    do i = 1, nw
      
        if( this%pWall(i)%user_id == user_id ) then
            ii = ii + 1
            do j = 1, 4
            
                points((ii-1)*4+j) = this%pWall(i)%getPoint(j)
            
            end do
            
            faces((ii-1)*5+1) = 4
            faces((ii-1)*5+2) = (ii-1)*4
            faces((ii-1)*5+3) = (ii-1)*4+1
            faces((ii-1)*5+4) = (ii-1)*4+2
            faces((ii-1)*5+5) = (ii-1)*4+3
            
        end if
        
    end do
    
    if(ii == 0 ) return
    
    !// starting to write wall points in the file
	write( nUnitFile, "(A)" )"DATASET POLYDATA"
    write( nUnitFile, "(A8,I6,A)" )"POINTS  ", ii*4 , "  float"
    
    do i = 1, ii*4
        write(nUnitFile, *) points(i)    
    end do
    
    write( nUnitFile, "(A8, I5, A1, I6)") "POLYGONS ", ii , " " , 5*ii 
    do i =1, ii*5
        write( nUnitFile, "(I6)") faces(i)
    end do
    
    close(nUnitFile)
    
    deallocate(points)
    deallocate(faces)
    
end subroutine



!**********************************************************************
! setting wall velocity by user_id
!**********************************************************************
subroutine G_setWallVelocity( this, user_id, t_vel , r_vel, r_line )
    implicit none
    class(Geometry) this
    integer(IK),intent(in):: user_id
    type(real3),intent(in) :: t_vel
    real(RK),intent(in)    :: r_vel
    type(p_line),intent(in) :: r_line
    
    integer(IK) i
    
    do i=1, this%num_pWall
    
        if( this%pWall(i)%user_id == user_id ) then
            call this%pWall(i)%setVelocity(t_vel , r_vel, r_line)    
        end if
        
    end do
    
    
end subroutine



!************************************************************************************
! determining if a sphere is in contact with specified (wall_id) wall in the system
!************************************************************************************
logical function G_isInContact( this , box, wall_id )

    implicit none
    class(Geometry) this
    type(real4),intent(in)  :: box
    integer(IK),intent(in)  :: wall_id
    
    integer(IK) i, nw
    
    nw = this%G_get_wall_index( wall_id )
      
    
    
    if( nw .ne. 0 ) then
        
        G_isInContact = this%pWall(nw)%isInContact(box)
        
    else
        
        call CheckForError( ErrT_Abort, "G_isInContact", &
                            "Wrong wall id in the input of isInContact function. Wall id is :"//num2str(wall_id) )
        
    end if
    
end function

!*****************************************************************************
! returning the normal vector (vn) and property type (pt) of a wall 
! with wall id and distance (dist) of center of a sphere(box) with this wall.
!*****************************************************************************
subroutine G_normal_dist( this , box, wall_id, nv, dist ,pt )
    implicit none
    class(Geometry)            this
    type(real4),intent(in)  :: box
    integer(IK),intent(in)  :: wall_id
    type(real3),intent(out) :: nv
    real(RK),   intent(out) :: dist
    integer(IK),intent(out) :: pt
    
    integer(IK) nw
    type(real3) cp
    
    nw = this%G_get_wall_index( wall_id )
    
    if( nw .ne. 0 ) then 
        
        call this%pWall(nw)%normal_dist(box , nv, dist, pt, cp)
    
    else
              
        call CheckForError( ErrT_Abort, "G_normal_dist", &
                            "Wrong wall id in the input of isInContact function. Wall id is :"//num2str(wall_id) )
    
    endif
         
end subroutine

!*****************************************************************************
! returning the normal vector (vn) and property type (pt) of a wall 
! with wall id and distance (dist) of center of a sphere(box) with this wall 
! and the velocity of the wall at its contact point with sphere. 
!*****************************************************************************
subroutine G_normal_dist_vel_type( this , box, wall_id, nv, dist , vel  , pt )
    implicit none
    class(Geometry)            this
    type(real4),intent(in)  :: box
    integer(IK),intent(in)  :: wall_id
    type(real3),intent(out) :: nv
    real(RK),   intent(out) :: dist
    type(real3),intent(out) :: vel
    integer(IK),intent(out) :: pt
    
    type(wallType)  WC
    
    integer(IK) nw
    
    nw = this%G_get_wall_index( wall_id )
    
    if( nw .ne. 0 ) then 
        
        call this%pWall(nw)%normal_dist_vel_type(box , nv, dist, vel, pt)
        
    else
        
        call CheckForError( ErrT_Abort, "G_normal_dist_vel_type",&
                            "Wrong wall type in the input of isInContact function. Wall id is :"//num2str(wall_id) )
        
    end if
        
end subroutine 


!************************************************************************************
!  returning number of plane walls
!************************************************************************************
function G_get_num_pWall(this) result(nw)
    implicit none
    class (Geometry) this
    integer(IK) nw
    
    nw = this%num_pWall
    
end function

!************************************************************************************
!  returning maximum number of plane walls
!************************************************************************************
function G_get_max_pWall(this) result(nw)
    implicit none
    class (Geometry) this
    integer(IK) nw
    
    nw = this%max_num_pWall
    
end function

!**********************************************************************
! moving all walls for one time step (dt)
!**********************************************************************
subroutine G_move_walls( this, dt )
    implicit none
    class(Geometry)       this
    real(RK),intent(in):: dt
    
    !// locals
    integer(IK) n, nw
    
    !//body
    
    nw = this%num_pWall
    
    do n = 1, nw
        
       if( .not. this%pWall(n)%MoveWall( dt ) ) then
           
           call CheckForError(ErrT_Abort, "G_move_walls" , "Cannot move wall correctly. Wall No."// num2str(n) )
           
       endif
           
    end do
    
end subroutine

!************************************************************************************
!  returning ith wall
!************************************************************************************
function G_get_pWall (this, i) result (res)
    implicit none
    class(Geometry)           this
    integer(IK),intent(in) :: i
    
    !// locals
    type(mvng_PlaneWall) res
    
    if( i <= this%num_pWall) then
        res = this%pWall(i)
    else
        call CheckForError( ErrT_Abort, "G_get_pWall" , "Out of range request for wall" )    
    end if
    
end function

!************************************************************************************
! increasing the size of the vector which stores plane walls
!************************************************************************************
subroutine Increase_pWall(this)
    implicit none
    class(Geometry) this
    
    !// locals
    integer(IK)  n, i, oldn
    type(mvng_PlaneWall),allocatable,dimension(:):: temp
    
    !//body
    
    if(this%max_num_pWall == 0 ) then
        n = 100
        allocate( this%pWall(n) )
        this%max_num_pWall = n
        return
    end if
    
    
    oldn = this%max_num_pWall
    
    allocate(temp (oldn) )

    do i = 1, oldn
	temp(i) = this%pWall(i)
    enddo
    

    if(allocated(this%pWall)) deallocate(this%pWall)
    
    n = 1.5*this%max_num_pWall +1  
    n = max(n,3)
    
    allocate( this%pWall(n) )
    
    do i=1, this%max_num_pWall
        this%pWall(i) = temp(i)
    end do
    
    
    
    this%max_num_pWall = n 
    
    deallocate(temp)
    
end subroutine 


!************************************************************************
! finding the index of wall in the pWalls vector based on the wall_id
!************************************************************************
function G_get_wall_index(this, wall_id) result (res)
    implicit none
    class(Geometry) this
    integer(IK),intent(in):: wall_id
    integer(IK) res
    
    integer(IK) i
    
    do i=1, this%num_pWall
        
        if( this%pWall(i)%wall_id == wall_id ) then
            
            res = i
            return
            
        end if
        
    end do
    
    ! if no wall is found
    res = 0
    
end function

!***********************************************************
! deleting a wall by index
!***********************************************************
subroutine G_delete_wall_id(this, wall_idx )
    class(Geometry) this
    integer(IK),intent(in):: wall_idx
    
    integer(IK) i
    
    ! shifting walls to the left in the vector pWall
    do i = wall_idx, this%num_pWall-1
        this%pWall(i) = this%pWall(i+1)        
    end do
    
    ! reducing number of walls by one
    this%num_pWall = this%num_pWall-1
    
end subroutine 



end module

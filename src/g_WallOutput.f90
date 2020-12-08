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
!  File name:  g_WallOutput.f90  
!  Module name: g_WallOutput
!  
!  Purpose: 
!    A module for output file for geometry in tecplot format
! 
!------------------------------------------------------------------------------

module g_WallOutput
    
    use g_TypeDef
    use g_PlaneWall
    
    
    implicit none
    
    
    integer,parameter:: WO_nDivisions = 72
    
    type WallOutput
        
        integer:: MaxPoints = 0
        integer:: MaxFaces = 0
        integer:: nPoints = 0
        integer:: nFaces  = 0 
        type(real3),dimension(:),allocatable:: points
        integer,dimension(:,:),allocatable:: faces
        
    contains
        procedure:: Add_PlaneWall  => WO_Add_PlaneWall
       
        procedure:: TecPlot_Out    => WO_TecPlot_Out
        procedure:: reset          => WO_reset
        
        procedure,private:: add_rectangle
        procedure,private:: extent_size_point
        procedure,private:: extent_size_face
    end type
    
    
    
contains


subroutine WO_Add_PlaneWall( this ,  wall )

    implicit none
    
    class(WallOutput) this
    type(mvng_PlaneWall) wall
    
    
    type(real3) p1, p2, p3, p4
    
        
    p1 = wall%getPoint(1)
    p2 = wall%getPoint(2)
    p3 = wall%getPoint(3)
    p4 = wall%getPoint(4)
   
    call this%add_rectangle(p1, p2, p3, p4)
    
    
end subroutine


subroutine WO_TecPlot_Out(this, nUnitFile, sol_time, strnd_id)
    
    implicit none
    
    class(WallOutput) this
    integer nUnitFile, strnd_id
    real(RK) sol_time
    
    integer i
    
    write(nUnitFile, *) "VARIABLES = X, Y, Z"
    write(nUnitFile, *) "ZONE NODES =", this%nPoints, ",ELEMENTS =", &
                        this%nFaces, ",ZONETYPE = FEQUADRILATERAL , DATAPACKING=POINT, SOLUTIONTIME = ", &
                        sol_time , ", STRANDID = ", strnd_id
    
    
    do i=1,this%nPoints
        write(nUnitFile, *)this%points(i)
    end do
    
    do i= 1, this%nFaces
        write(nUNitFile, *) this%faces(i,:)    
    end do
    
    
end subroutine 

subroutine WO_reset(this)

    implicit none
    class(WallOutput) this
    
    
    this%nPoints = 0
    this%nFaces  = 0 
        
end subroutine


subroutine add_rectangle( this , p1, p2, p3, p4 )
    
    implicit none
    
    class(WallOutput) this
    type(real3) p1, p2, p3, p4
    
    
    ! checking if there is enough space for new points
    call this%extent_size_point(4_IK)
    
    ! adding all four points of wall node points
    this%points( this%nPoints +1) = P1
    this%points( this%nPoints +2) = P2
    this%points( this%nPoints +3) = P3
    this%points( this%nPoints +4) = P4
        
    
    ! checking if there is enough space for new quadrilateral faces
    call this%extent_size_face(2_IK)
    
    ! adds face (P1,P2,P3, P4)
    this%faces( this%nFaces+1 , : ) = (/this%nPoints+1, this%nPoints+2, this%nPoints+3, this%nPoints+4/)
    
    ! adding second triangle (P1,P3,P4) 
   ! this%faces( this%nFaces+2 , : ) = (/this%nPoints+1, this%nPoints+3, this%nPoints+4/)
    
    this%nPoints = this%nPoints + 4
    this%nFaces = this%nFaces + 1
    
end subroutine

subroutine extent_size_point( this, extraPoints )
    
    implicit none
    class(WallOutput) this
    integer(IK) extraPoints
    
    integer(IK) numP
    type(real3),dimension(:),allocatable:: temp
        
    if( this%nPoints + extraPoints < this%maxPoints ) return
            
    numP = max( 10 , int(this%maxPoints*1.3 + extraPoints ) )
    
    if( this%maxPoints == 0 ) then
        
        if ( allocated(this%points) ) deallocate(this%points)
        allocate( this%points(numP ) )
        
    else
        
        allocate( temp(this%maxPoints) ) 
        temp(1:this%maxPoints) = this%points(1:this%maxPoints)
    
        if ( allocated(this%points) ) deallocate(this%points)
        allocate( this%points(numP ) )
    
        this%points(1:this%maxPoints) = temp(1:this%maxPoints)
        deallocate(temp)
    
    end if
    
    this%maxPoints = numP
       
end subroutine 

subroutine extent_size_face( this, extraface )
    
    implicit none
    class(WallOutput) this
    integer(IK) extraface
    
    integer(IK) numf
    integer,dimension(:,:),allocatable:: temp
        
    if( this%nFaces + extraface < this%maxFaces ) return
            
    numf = max( 10 , int(this%maxFaces*1.3 + extraface ) )
    
    if( this%maxFaces == 0 ) then
        
        if ( allocated(this%faces) ) deallocate(this%faces)
        allocate( this%faces(numf , 4 ) )
        
    else
        
        allocate( temp(this%maxFaces , 4) ) 
        temp(1:this%maxFaces,:) = this%faces(1:this%maxFaces,:)
    
        if ( allocated(this%faces) ) deallocate(this%faces)
        allocate( this%faces(numf , 4 ) )
    
        this%Faces(1:this%maxfaces,:) = temp(1:this%maxFaces,:)
        deallocate(temp)
    
    end if
    
    this%maxFaces = numf
       
end subroutine 


end module

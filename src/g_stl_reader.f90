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
!  file name  : g_stl_reader.f90  
!  module name: g_stl_reader
! 
!  Purpose:        
!    1) Reading an stl file: This file contains triangles which covers the wall
!       surface.
! 
!  The stl file should be in ASCII format 
! 
!------------------------------------------------------------------------------

module g_stl_reader

    use g_TypeDef
    use g_Prtcl_DefaultValues
    implicit none
    
    private
    public:: read_from_stl_file, facet
    
    type facet
        type(real3) p1
        type(real3) p2
        type(real3) p3
    end type
    
    
    integer numFacet
    integer max_numFacet
    
    contains

    
!******************************************************************
! taking an input string and changes all the letters to upper case  
!******************************************************************
function upperCase( ch )
    implicit none
    character,intent(in):: ch
    character upperCase
    
    integer idx

    idx = ichar ( ch )
    if ( 97 <= idx .and. idx <= 122 ) then
        upperCase = char ( idx - 32 )
    else
        upperCase = ch
    end if

end function 


!************************************************************************
! comparing two input strings to be similar, it is not case-sensitive
!************************************************************************
function string_compare( str1, str2 )
    implicit none
    character(*),intent(in):: str1, str2
    logical string_compare
    
    integer i, len1, len2
    character ch1, ch2
    character(256) s1, s2
    
    
    s1 = adjustl( str1)
    s2 = adjustl( str2)
    
    len1 = len_trim( s1 )
    len2 = len_trim( s2 )
    
    string_compare = .false.
    
    if( len1 .ne. len2 ) return
    
    do i = 1,len1
        
        if( upperCase(s1(i:i)) .ne. upperCase(s2(i:i)) ) return
        
    end do
    
    string_compare = .true.

end function 


!***************************************************************
! reading three numbers from the input string and returning them
!***************************************************************
function read_three_numbers( chNumbers , nums ) result( res )
    implicit none
    character(*),intent(in)          :: chNumbers
    real(RK),dimension(3),intent(out):: nums
    logical res
    
    !// locals
    integer istat
    
    read( chNumbers , *, IOSTAT = istat ) nums
    
    if( istat .ne. 0 ) then
        res = .false.    
    else
        res = .true. 
    end if
    
    
end function 


!**********************************************************************************************
! reading from a file with  ASCII stl format
!    it first opens file (chFile) and reads the whole file and put all triangles in all_facets
!**********************************************************************************************
function read_from_stl_file( chFile , all_facets, num_facets, ch_error ) result (res)
    implicit none
    character(*),intent(in)     :: chFile   ! address of the stl file
    type(facet),dimension(:),allocatable,intent(inout):: all_facets ! all triangles which are read from file
    integer(IK),intent( out)   :: num_facets                        ! number of triangles read from file 
    character(*),intent(out)   :: ch_error                          ! error message
    logical res                                                     ! success status   
    
    
    ! // locals
    integer istat
    integer lvl, vertex_num, line_num
    logical finished, end_f, error_read, error_fin
    real(RK),dimension(3)   :: nums
    type(real3),dimension(3):: vert
    
    character(255) chLine, word1, word2, ch
    
    ! assuming there is 100 triangles in the file and allocates memory for it 
    max_numFacet = 100
    numFacet = 0
    call re_allocate(all_facets)
    
    
    
    ! first step: opening the file 
    open( file = chFile, unit = d_stl_reader_unit , IOSTAT = istat , STATUS = "old" )
    
    if( istat .ne. 0 ) then
        res = .false.
        ch_error = "error in opening file :" // trim(chFile)
        return
    end if
    
    lvl        = 0
    line_num   = 0
    vertex_num = 0
    ch_error = ""
    error_fin = .false. 
    end_f = .false.
    error_read = .false.
    
    
    ! second step: reading all lines of the input file
    do
       
        ! reads one line 
        read( d_stl_reader_unit , "(A)" ,IOSTAT = istat ) chLine
        line_num = line_num + 1
        
        ! checking the status
        end_f = .false.
        error_read = .false.
        if( istat .gt. 0 ) then
            error_read = .true.
            exit
        elseif ( istat .lt. 0) then
            end_f = .true.
            exit
        endif
        
        ! getting next word in the line
        call get_next_word( chLine , word1 , finished )
        
        ! if there is no word left in the line, go to the next line 
        if( finished ) cycle
        
        ! converting end solid to endsolid and so on....
        if( string_compare( word1 , "END" ) .or. string_compare( word1 , "FACET" ) .or. string_compare( word1 , "OUTER" )  ) then
        
            call get_next_word( chLine, word2, finished )
            
            if( finished ) then
                error_fin = .true.
                exit
            else
                error_fin = .false.
            end if
            
            ! connecting two words 
            word1 = trim(word1)//trim(word2)
        
        end if
        
        ! this list helps to understand the code line below
        ! the first word tells what to do.
        !
        !  SOLID - begin a new solid.
        !    Valid in state 0, moves to state 1.
        !  ENDSOLID - end current solid.
        !    Valid in state 1, moves to state 0.
        !
        !  FACETNORMAL - begin a new facet.
        !    Valid in state 0 or 1, moves to state 2.
        !  ENDFACET - end current facet.
        !    Valid in state 2, moves to state 1.
        !
        !  OUTERLOOP - begin a list of vertices.
        !    Valid in state 2, moves to state 3.
        !  ENDLOOP - end vertex list.
        !    Valid in state 3, moves to state 2.
        !
        !  VERTEX - give coordinates of next vertex.
        !    Valid in state 3 if current vertex count is 0, 1 or 2.
        !
        !  End of file -
        !    Valid in state 0 or 1.
        !
        
        ! start of a new solid
        if ( string_compare ( word1, 'SOLID' ) ) then
            
            if( lvl .ne. 0 ) then
                error_fin = .true.
                ch_error = "A new SOLID statement was encountered without processing and ENDSOLID statement" 
                exit
            end if
            
            lvl = 1       
            
        ! end of the solid
        elseif( string_compare(word1,"ENDSOLID" ) )then
            
            if ( lvl .ne. 1 ) then
                error_fin = .true. 
                ch_error = '  An END SOLID statement was encountered in a wrong place'
                exit
            end if
            lvl = 0
            
        ! beginning a new facet     
        else if( string_compare ( word1, 'FACETNORMAL' ) ) then

            if ( lvl .ne. 0 .and. lvl .ne. 1 ) then
                error_fin = .true. 
                ch_error = 'Model not in right state for FACET.'
                exit
            end if
      
            lvl = 2
            ! the normal vector of facet is neglected 
        
        ! end of the facet 
        elseif (string_compare ( word1, 'ENDFACET' ) ) then

            if ( lvl .ne. 2 ) then
                error_fin = .true. 
                ch_error = '  Model not in right state for ENDFACET.'
                exit
            end if

            lvl = 1
        
        ! beginning a list of vertices 
        elseif ( string_compare ( word1, 'OUTERLOOP' ) ) then

            if ( lvl .ne. 2 ) then
        
                error_fin = .true.
                ch_error = '  Model not in right state for OUTERLOOP.'
                exit
            end if
            

            lvl = 3
            vertex_num = 0
            
        !end of the vertex list
        else if ( string_compare ( word1, 'ENDLOOP' ) ) then

            if (lvl .ne. 3 ) then
                error_fin = .true. 
                ch_error = '  Model not in right state for ENDLOOP.'
                exit
            end if
            
            if( vertex_num .ne. 3 ) then
                error_fin = .true.
                ch_error = " Model not supplied enough vertices."
                exit
            end if
            
            !** adding a facet to the facet vector (all_facets)
            
            ! first checking if there enough space for new facet 
            if( numFacet .ge. max_numFacet ) then
                call re_allocate( all_facets )
            end if
            
            numFacet = numFacet + 1 
            all_facets( numFacet ) = facet( vert(1), vert(2), vert(3) )
            
            lvl = 2
        
        ! reading a vertex
        else if ( string_compare ( word1, 'VERTEX' ) ) then

            if ( lvl .ne. 3 ) then
                error_fin = .true.
                ch_error = '  Model not in right state for VERTEX.'
                exit
            end if
            
            if( vertex_num .ge. 3 ) then
                error_fin = .true.
                ch_error = '  More than 3 vertices specified for a face.'
                exit
            end if
            
            !// reading three numbers next to the vertex keyword
    
            if( read_three_numbers( chLine, nums ) ) then
                
                vertex_num = vertex_num + 1
                vert( vertex_num ) = real3( nums(1), nums(2), nums(3) )
                
            else
                
                error_fin = .true. 
                ch_error  = " Model supplied wrong format for three numbers in vertex" 
                exit
            endif           
        
        else
            
            error_fin = .true. 
            ch_error = " unkonwn line in the file"
            exit    
        end if
        
        
    end do ! end of the loop which reads file content line-by-line 
    
    ! checking probable errors 
    if( end_f .and. lvl .ne. 0 ) then
        
        ch_error = " End of file occured without ENDSOLID statement"
        error_fin = .true.
        
    end if
    
    if( error_read ) then
        ch_error = " Error in reading line from file"
        error_fin = .true.
        
    end if
    
    ! closing the stl file 
    close(d_stl_reader_unit)
    
    
    if( error_fin ) then
        write( ch, *) line_num    
        ch_error =  trim(ch_error)// " in line " // trim(ch)
        res = .false.
        return
    end if
    
    ! returning number of facets and success status 
    num_facets = numFacet 
    res = .true.
    
end function


!****************************************************************
! getting the first next word in the input string (str), if any.
!****************************************************************
subroutine get_next_word ( str , word , finished )
    implicit none

    character(*),intent(inout) :: str
    character(*),intent(out)   :: word
    logical,intent(out)        :: finished
    
    !// locals
    character, parameter :: TAB = char ( 9 )
    integer next, ilo
    integer lenc
    
    next = 1
    finished = .false.
    lenc = len_trim ( str )
    
    ! checking if str has any character
    if ( lenc <= 0 ) then
      finished = .true.
      word = ' '
      return
    end if
    
    !
    !  searching str for the next nonblank character
    !
     ilo = next
    
  do
    ! looping until finds a non-blank character and putting it in the ilo
    if ( str(ilo:ilo) /= ' ' .and. str(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do

 
    !  If this initial nonblank is a special character,
    !  then that's the whole word as far as we're concerned,
    !  so return immediately.
    !
    if (   str(ilo:ilo) == '"' .or. &
           str(ilo:ilo) == '(' .or. &
           str(ilo:ilo) == ')' .or. &
           str(ilo:ilo) == '{' .or. &
           str(ilo:ilo) == '}' .or. &
           str(ilo:ilo) == '[' .or. &
           str(ilo:ilo) == ']' ) then

        word = str(ilo:ilo)
        str = str(ilo+1:)
        return

    end if
    
    
    !  Now searching for the last contiguous character that is not a
    !  blank, TAB, or special character.
    
    next = ilo + 1

    do while ( next <= lenc )

        if ( str(next:next) == ' ' ) then
            exit
        else if ( str(next:next) == TAB ) then
            exit
        else if ( str(next:next) == '"' ) then
            exit
        else if ( str(next:next) == '(' ) then
            exit
        else if ( str(next:next) == ')' ) then
            exit
        else if ( str(next:next) == '{' ) then
            exit
        else if ( str(next:next) == '}' ) then
            exit
        else if ( str(next:next) == '[' ) then
            exit
        else if ( str(next:next) == ']' ) then
            exit
        end if

        next = next + 1
    end do

    if ( str(next-1:next-1) == ',' ) then
        word = str(ilo:next-2)
    else
        word = str(ilo:next-1)
    end if
    
    str = str(next:)
    
    return

end subroutine 

!******************************************************
! reallocating memory for facets
!******************************************************
subroutine re_allocate( f )
    implicit none
    type(facet),dimension(:),allocatable:: f
    
    integer(IK) oldn, i
    type(facet),dimension(:),allocatable:: temp
    
    if( numFacet .eq. 0 ) then
        if( allocated(f) ) deallocate(f)
        allocate( f(max_numFacet ) )
        return
    
    end if

    oldn = max_numFacet
    if( allocated(temp) ) deallocate(temp)
    
    allocate( temp (oldn) ) 

    do i=1, oldn
	temp(i) = f(i)
    end do
  
    deallocate(f)
    
    max_numFacet = max( 10, int(max_numFacet*1.25) )
    
    allocate( f(max_numFacet) )
    f(1:oldn ) = temp(1:oldn)
    
    deallocate( temp )
    
end subroutine

    
end module

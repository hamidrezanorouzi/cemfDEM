function User_Mark( id, flag, ptype, dpos, oldMark ) result (mark)
    use g_TypeDef
    use g_Prtcl_DefaultValues
    implicit none 
    
    integer(IK),intent(in):: id, flag, ptype, oldMark
    type(real4),intent(in):: dpos 
    integer(IK) mark
    
    
    !// locals
    
   	mark = oldMark
    	   	
    return       
    
end function

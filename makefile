#-------------- cemfDEM code: a dem simulation code ----------------------------
#      D        C enter of
#     D M       E ngineering and
#    D   M      M ultiscale modeling of
#   D     M     F luid flow    
#  EEEEEEEEM    .ir
#------------------------------------------------------------------------------
#  Copyright (C): cemf
#  website: www.cemf.ir
#------------------------------------------------------------------------------  
#  This file is part of cemfDEM code. It is a free software for simulating 
#  granular flow. You can redistribute it and/or modify it under the terms of
#  GNU General Public License version 3 or any other later versions. 
# 
#  cemfDEM code is distributed to help others in their research in the field  
#  of granular flow, but WITHOUT ANY WARRANTY; without even the implied 
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  
#  makefile to build the program 
#------------------------------------------------------------------------------

FortC = gfortran # can be replaced with any standard fortran compiler 
FFLAGS = -O3

INCLUDES =
ExeName=cemfDEM
LFLAGS =
LIBS=

#current (top most directory) of the code
BASEPATH =$(CURDIR)

# Relative paths to the directories
SRCDIR =$(BASEPATH)/src
OBJDIR =$(BASEPATH)/obj
MODDIR =$(BASEPATH)/mod
RESDIR =$(BASEPATH)/Results

Exe = $(ExeName)
# Put module files in the mod dir
FFLAGS += -J$(MODDIR)

# Search path for source code
VPATH = $(SRCDIR)


ifneq ($(MODDIR),)
  $(shell test -d $(MODDIR) || mkdir -p $(MODDIR))
endif

ifneq ($(RESDIR),)
  $(shell test -d $(RESDIR) || mkdir -p $(RESDIR))
endif

ifneq ($(OBJDIR),)
  $(shell test -d $(OBJDIR) || mkdir -p $(OBJDIR))
endif


# You can also remove the full path part
modSRC = \
g_DEMSystem.f90 \
g_Prtcl_Integration.f90 \
g_Prtcl_DefaultValues.f90 \
g_TypeDef.f90 \
g_LogInfo.f90 \
g_Prtcl_OneStepIntegration.f90 \
g_Prtcl_MultiPointIntegration.f90 \
g_error_handling.f90 \
g_Prtcl_TwoStepIntegratoin.f90 \
g_ContactSearch.f90 \
g_Geometry.f90 \
g_PlaneWall.f90 \
g_Line.f90 \
g_CylinderWall.f90 \
g_stl_reader.f90 \
g_Prtcl_NBS.f90 \
g_Prtcl_CellBased.f90 \
g_Prtcl_SimWorld.f90 \
g_Prtcl_SimDomain.f90 \
g_Prtcl_ContactList.f90 \
g_LinkedList.f90 \
g_Prtcl_ContactInfo.f90 \
g_Prtcl_NBS_Munjiza.f90 \
g_Prtcl_Hrchl_NBS.f90 \
g_prtcl_Hrchl_Munjiza.f90 \
g_ContactSearchPW.f90 \
g_Prtcl_LSD_Model.f90 \
g_Prtcl_ContactForce.f90 \
g_Prtcl_Property.f90 \
g_Prtcl_pureProp.f90 \
g_MakePrtcls.f90 \
g_RandumNum.f90 \
g_Prtcl_NonLin_Model.f90 \
g_WallOutput.f90 \
g_timer.f90 \
g_Prtcl_Insertion.f90


exeSrc  =  main.f90 ProgramDefinedGeometry.f90 User_Mark.f90

# Replace .f90 with .o and append object dir
modOBJ_O = $(patsubst %.f90,%.o,$(modSRC))
modOBJ = $(addprefix $(OBJDIR)/,$(modOBJ_O))

# Replace .f90 with .O and append object dir
exeOBJ_O = $(patsubst %.f90,%.o,$(exeSrc))
exeOBJ   = $(addprefix $(OBJDIR)/,$(exeOBJ_O))

allObjects =  $(modOBJ) $(exeOBJ)


all: main

main:$(allObjects)
	$(FortC) $(INCLUDES) -o $(Exe) $+ $(LFLAGS) $(LIBS)
	@echo  The executable $(Exe) has been built

$(OBJDIR)/%.o: %.f90 
	$(FortC) $(FFLAGS) $(INCLUDES) -c -o $@ $< 	


# Dependency tree for modules and other sources
$(OBJDIR)/main.o:$(OBJDIR)/g_DEMSystem.o $(OBJDIR)/g_stl_reader.o $(OBJDIR)/g_Prtcl_pureProp.o

$(OBJDIR)/ProgramDefinedGeometry.o:$(OBJDIR)/g_Geometry.o $(OBJDIR)/g_Line.o

$(OBJDIR)/userMark.o:$(OBJDIR)/g_TypeDef.o 
	 
$(OBJDIR)/g_DEMSystem.o: $(OBJDIR)/g_Prtcl_Integration.o $(OBJDIR)/g_ContactSearch.o $(OBJDIR)/g_ContactSearchPW.o $(OBJDIR)/g_Prtcl_LSD_Model.o $(OBJDIR)/g_Prtcl_NonLin_Model.o $(OBJDIR)/g_WallOutput.o $(OBJDIR)/g_timer.o $(OBJDIR)/g_Prtcl_Insertion.o

$(OBJDIR)/g_Prtcl_Integration.o: $(OBJDIR)/g_Prtcl_DefaultValues.o $(OBJDIR)/g_Prtcl_OneStepIntegration.o $(OBJDIR)/g_Prtcl_MultiPointIntegration.o $(OBJDIR)/g_Prtcl_TwoStepIntegratoin.o

$(OBJDIR)/g_Prtcl_DefaultValues.o:$(OBJDIR)/g_TypeDef.o $(OBJDIR)/g_LogInfo.o

$(OBJDIR)/g_LogInfo.o: $(OBJDIR)/g_TypeDef.o

$(OBJDIR)/g_Prtcl_OneStepIntegration.o: $(OBJDIR)/g_TypeDef.o

$(OBJDIR)/g_Prtcl_MultiPointIntegration.o: $(OBJDIR)/g_TypeDef.o $(OBJDIR)/g_error_handling.o $(OBJDIR)/g_Prtcl_OneStepIntegration.o

$(OBJDIR)/g_error_handling.o: $(OBJDIR)/g_Prtcl_DefaultValues.o

$(OBJDIR)/g_Prtcl_TwoStepIntegratoin.o: $(OBJDIR)/g_TypeDef.o $(OBJDIR)/g_error_handling.o $(OBJDIR)/g_Prtcl_MultiPointIntegration.o

$(OBJDIR)/g_ContactSearch.o: $(OBJDIR)/g_Geometry.o $(OBJDIR)/g_Prtcl_NBS.o $(OBJDIR)/g_Prtcl_NBS_Munjiza.o $(OBJDIR)/g_Prtcl_Hrchl_NBS.o $(OBJDIR)/g_Prtcl_DefaultValues.o $(OBJDIR)/g_prtcl_Hrchl_Munjiza.o

$(OBJDIR)/g_Geometry.o: $(OBJDIR)/g_PlaneWall.o $(OBJDIR)/g_CylinderWall.o $(OBJDIR)/g_stl_reader.o

$(OBJDIR)/g_PlaneWall.o: $(OBJDIR)/g_Line.o $(OBJDIR)/g_Prtcl_DefaultValues.o $(OBJDIR)/g_error_handling.o

$(OBJDIR)/g_Line.o: $(OBJDIR)/g_TypeDef.o

$(OBJDIR)/g_CylinderWall.o: $(OBJDIR)/g_PlaneWall.o $(OBJDIR)/g_error_handling.o

$(OBJDIR)/g_stl_reader.o: $(OBJDIR)/g_TypeDef.o $(OBJDIR)/g_Prtcl_DefaultValues.o

$(OBJDIR)/g_Prtcl_NBS.o:  $(OBJDIR)/g_Prtcl_CellBased.o $(OBJDIR)/g_error_handling.o

$(OBJDIR)/g_Prtcl_CellBased.o: $(OBJDIR)/g_error_handling.o $(OBJDIR)/g_Prtcl_SimWorld.o

$(OBJDIR)/g_Prtcl_SimWorld.o: $(OBJDIR)/g_Prtcl_SimDomain.o $(OBJDIR)/g_Prtcl_ContactList.o $(OBJDIR)/g_Prtcl_DefaultValues.o

$(OBJDIR)/g_Prtcl_SimDomain.o: $(OBJDIR)/g_TypeDef.o

$(OBJDIR)/g_Prtcl_ContactList.o: $(OBJDIR)/g_LinkedList.o $(OBJDIR)/g_Prtcl_ContactInfo.o $(OBJDIR)/g_error_handling.o

$(OBJDIR)/g_LinkedList.o: $(OBJDIR)/g_TypeDef.o

$(OBJDIR)/g_Prtcl_ContactInfo.o: $(OBJDIR)/g_TypeDef.o

$(OBJDIR)/g_Prtcl_NBS_Munjiza.o: $(OBJDIR)/g_Prtcl_CellBased.o

$(OBJDIR)/g_Prtcl_Hrchl_NBS.o: $(OBJDIR)/g_TypeDef.o $(OBJDIR)/g_error_handling.o $(OBJDIR)/g_Prtcl_CellBased.o $(OBJDIR)/g_Prtcl_NBS.o

$(OBJDIR)/g_prtcl_Hrchl_Munjiza.o: $(OBJDIR)/g_TypeDef.o $(OBJDIR)/g_error_handling.o $(OBJDIR)/g_Prtcl_NBS_Munjiza.o

$(OBJDIR)/g_ContactSearchPW.o: $(OBJDIR)/g_TypeDef.o $(OBJDIR)/g_Prtcl_ContactList.o $(OBJDIR)/g_Geometry.o

$(OBJDIR)/g_Prtcl_LSD_Model.o: $(OBJDIR)/g_Prtcl_ContactForce.o

$(OBJDIR)/g_Prtcl_ContactForce.o: $(OBJDIR)/g_Prtcl_Property.o $(OBJDIR)/g_Prtcl_ContactList.o $(OBJDIR)/g_Geometry.o

$(OBJDIR)/g_Prtcl_Property.o: $(OBJDIR)/g_Prtcl_pureProp.o $(OBJDIR)/g_MakePrtcls.o

$(OBJDIR)/g_Prtcl_pureProp.o: $(OBJDIR)/g_Prtcl_DefaultValues.o

$(OBJDIR)/g_MakePrtcls.o: $(OBJDIR)/g_Prtcl_pureProp.o $(OBJDIR)/g_error_handling.o $(OBJDIR)/g_RandumNum.o

$(OBJDIR)/g_RandumNum.o: $(OBJDIR)/g_TypeDef.o

$(OBJDIR)/g_Prtcl_NonLin_Model.o: $(OBJDIR)/g_Prtcl_LSD_Model.o

$(OBJDIR)/g_WallOutput.o: $(OBJDIR)/g_TypeDef.o $(OBJDIR)/g_PlaneWall.o

$(OBJDIR)/g_timer.o: $(OBJDIR)/g_TypeDef.o

$(OBJDIR)/g_Prtcl_Insertion.o: $(OBJDIR)/g_PlaneWall.o



clean:
	rm -r -f $(allObjects) $(MODDIR)/*
	@echo The build has been cleaned

vtkClean:
	rm -f -r *.vtk $(RESDIR)/*.vtk *.plt $(RESDIR)*.plt
	@echo The content of $(RESDIR) has been cleaned

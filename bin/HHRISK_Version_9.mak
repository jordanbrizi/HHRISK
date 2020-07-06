# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=HHRISK_Version_9 - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to HHRISK_Version_9 - Win32\
 Debug.
!ENDIF 

!IF "$(CFG)" != "HHRISK_Version_9 - Win32 Release" && "$(CFG)" !=\
 "HHRISK_Version_9 - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "HHRISK_Version_9.mak" CFG="HHRISK_Version_9 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "HHRISK_Version_9 - Win32 Release" (based on\
 "Win32 (x86) Console Application")
!MESSAGE "HHRISK_Version_9 - Win32 Debug" (based on\
 "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
F90=fl32.exe
RSC=rc.exe

!IF  "$(CFG)" == "HHRISK_Version_9 - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
OUTDIR=.
INTDIR=.

ALL : "$(OUTDIR)\HHRISK_Version_9.exe"

CLEAN : 
	-@erase ".\HHRISK_Version_9.exe"
	-@erase ".\HHRISK_Version_9.obj"

# ADD BASE F90 /Ox /c /nologo
# ADD F90 /Ox /c /nologo
F90_PROJ=/Ox /c /nologo 
# ADD BASE RSC /l 0x416 /d "NDEBUG"
# ADD RSC /l 0x416 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/HHRISK_Version_9.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/HHRISK_Version_9.pdb" /machine:I386\
 /out:"$(OUTDIR)/HHRISK_Version_9.exe" 
LINK32_OBJS= \
	"$(INTDIR)/HHRISK_Version_9.obj"

"$(OUTDIR)\HHRISK_Version_9.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "HHRISK_Version_9 - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
OUTDIR=.
INTDIR=.

ALL : "$(OUTDIR)\HHRISK_Version_9.exe"

CLEAN : 
	-@erase ".\HHRISK_Version_9.exe"
	-@erase ".\HHRISK_Version_9.obj"
	-@erase ".\HHRISK_Version_9.ilk"
	-@erase ".\HHRISK_Version_9.pdb"

# ADD BASE F90 /Zi /c /nologo
# ADD F90 /Zi /c /nologo
F90_PROJ=/Zi /c /nologo /Fd"HHRISK_Version_9.pdb" 
# ADD BASE RSC /l 0x416 /d "_DEBUG"
# ADD RSC /l 0x416 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/HHRISK_Version_9.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/HHRISK_Version_9.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)/HHRISK_Version_9.exe" 
LINK32_OBJS= \
	"$(INTDIR)/HHRISK_Version_9.obj"

"$(OUTDIR)\HHRISK_Version_9.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for.obj:
   $(F90) $(F90_PROJ) $<  

.f.obj:
   $(F90) $(F90_PROJ) $<  

.f90.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "HHRISK_Version_9 - Win32 Release"
# Name "HHRISK_Version_9 - Win32 Debug"

!IF  "$(CFG)" == "HHRISK_Version_9 - Win32 Release"

!ELSEIF  "$(CFG)" == "HHRISK_Version_9 - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\HHRISK_Version_9.f90

"$(INTDIR)\HHRISK_Version_9.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################

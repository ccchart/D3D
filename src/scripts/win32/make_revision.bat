@ECHO OFF

REM =====================================
REM Arguments
REM =====================================

REM %1: Top directory of the source tree: used to define SVN_DIR and CMD_DIR
REM %2: Module directory                : svnrevision is executed in this directory
REM %3: Version number file             : containing MAJOR, MINOR and REVISION definitions
REM %4: Input file                      : containing SVN_REVISION to be replaced, normally outputfile.svn (with a double extension)
REM %5: Output file                     : contents of input file with SVN_REVISION replaced by the actual revision string
REM %6: ConfigurationName (optional)    : Name of Visual Studio configuration, to optionally skip version file generation.


REM =====================================
REM Get all directories needed
REM =====================================

SET CURDIR=%CD%
CD %1
SET TOPDIR=%CD%
CD %CURDIR%
CD %2
SET MODDIR=%CD%

SET SVN_DIR=%TOPDIR%\third_party_lgpl\subversion\bin\win32
SET CMD_DIR=%TOPDIR%\third_party_lgpl\commandline\bin\win32
SET VN_DIR=%TOPDIR%\third_party_lgpl\version_number\bin\win32

CD %CURDIR%

REM Skip generation if 6th argument is 1. PRESENT and 2. equal to "Debug" and 3. target file already exists.
IF "%6"=="Debug" (
   IF EXIST "%5" (
      echo %0: Leaving existing file '%5' as is.
      EXIT
   ) ELSE (
      echo %0: Create missing file '%5'.
   )
) ELSE (
   echo %0: Regenerating existing file '%5'.
)

CD "%2"

IF DEFINED BUILD_NUMBER (
   echo build exists
   REM =====================================
   REM BUILD_NUMBER already known
   REM =====================================

) ELSE (

   REM =====================================
   REM Execute svnrevision
   REM =====================================

   CD "%MODDIR%"
   IF EXIST "%SVN_DIR%\svnversion.exe" (
        "%SVN_DIR%\svnversion.exe" -n "%MODDIR%" > "%MODDIR%\BUILD_NUMBER"
   ) ELSE (
        ECHO exported > "%MODDIR%\BUILD_NUMBER"
   )					         
   REM also set it as an environment variable    
   SET /p BUILD_NUMBER= < "%MODDIR%\BUILD_NUMBER"	         

) 

REM ==========================================================================
REM If the source has been obtained using a svn export command, the "exported"
REM string has been generated, but this cannot be used within *.rc files
REM Replace it using 00000 (only necessary on Windows systems)
REM ==========================================================================

IF "%BUILD_NUMBER%" == "exported" (
   SET BUILD_NUMBER=00000
)
echo %BUILD_NUMBER%

REM =====================================
REM Build substitution line
REM =====================================

SET ADDLINE=%BUILD_NUMBER%

                         
                         
REM =====================================
REM Inputfile > Substitute > Outputfile
REM =====================================

CD "%CURDIR%"

"%VN_DIR%\version_number.exe" %BUILD_NUMBER% "%3" "%4" "%5"

REM =====================================
REM Clean up
REM =====================================

del /f "%MODDIR%\BUILD_NUMBER" > del.log 2>&1
del /f del.log 



REM =====================================
REM Finished
REM =====================================

REM pause

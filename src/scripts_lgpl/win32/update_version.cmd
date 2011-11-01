@echo off
rem Program will replace %1 by the %1.svn and will replace VERSION_BUILD_NUMBER by a corresponding svn version using svnversion command
rem
rem %1 - path to the target source file
rem %2 - path to the folder to be used to check svnversion
rem %3 - Single file with version number information version_number.ini

echo Generating version number in the %1 ...

set SCRIPT_DIRECTORY=%~dp0

set SED=%SCRIPT_DIRECTORY%\..\..\third_party_open\commandline\bin\win32\sed.exe
set SVNVERSION=%SCRIPT_DIRECTORY%\..\..\third_party_open\subversion\bin\win32\svnversion.exe
set VN=%SCRIPT_DIRECTORY%\..\..\third_party_open\version_number\bin\win32\version_number.exe

IF DEFINED BUILD_NUMBER (
	set version=%BUILD_NUMBER%
) ELSE (
	rem Obtain the svn version number 
	"%SVNVERSION%" %2 | "%SED%" "s/\(.*\)/set version=\1/" > setversion.bat
 	call setversion.bat & del setversion.bat > NUL
)

rem ==========================================================================
rem If the source has been obtained using a svn export command, the "exported"
rem string has been generated, but this cannot be used within *.rc files
rem Replace it using 000000 (only necessary on Windows systems)
rem ==========================================================================

set version=%version:~0,8%
IF "%version%" == "exported" (
   set version=000000
)

if exist %1 (
	del %1
)
rem Generate version number source module using version_number.exe
"%VN%" %version% "%3" "%1.svn" "%1"

:end

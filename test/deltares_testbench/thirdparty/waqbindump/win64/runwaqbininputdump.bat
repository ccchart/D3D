echo off

echo "Script      : %~dp0runwaqbininputdump.bat"
echo "In directory: %CD%"
echo "Executing   : %~dp0waqbininputdump.exe %1 %2 >%3"
%~dp0waqbininputdump.exe %1 %2 >%3

@echo off
echo Now running FCEA2 to recompile thermo_spg.inp
echo thermo_spg | FCEA2

:: Thanks ChatGPT for this simple script. I hate batch scripting.
setlocal EnableDelayedExpansion

set "file_path=thermo_spg.out"
set "keyword=FATAL"

:: Get the last line of the file
for /f %%a in ('find /c /v "" ^< "%file_path%"') do set lines=%%a
set /a last_line=!lines!-1

for /f "usebackq skip=%last_line% delims=" %%b in ("%file_path%") do set "last_line_content=%%b"

:: Check if the last line contains the keyword
echo !last_line_content! | find /i "%keyword%" > nul
if %errorlevel% equ 0 (
    echo ERROR: There was a fatal error in compilation. Look at the end of thermo_spg.out for details
) else (
    echo Recompilation successful^^! Probably... No fatal errors found in thermo_spg.out.
)

endlocal
pause
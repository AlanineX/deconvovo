@echo off
REM Build the offline WatersConsole application on Windows x64.
REM Run from the repository root; output is waters-build\dist\WatersConsole.
python -m pip install -r portable_waters_convert\build\requirements.txt
if errorlevel 1 exit /b 1
python build\build_win.py
if errorlevel 1 exit /b 1
echo Output: waters-build\dist\WatersConsole\WatersConsole.exe

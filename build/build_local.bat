@echo off
REM Build DeconVoVo Windows executable locally.
REM Run this on a Windows machine with Python installed.
REM
REM Usage:
REM   1. Open Command Prompt in the project root
REM   2. Run: build\build_local.bat
REM   3. Output: dist\DeconVoVo\DeconVoVo.exe
REM
REM Prerequisites (one-time):
REM   pip install pyinstaller pyside6 numpy scipy pandas matplotlib plotly unidec

echo === Building DeconVoVo ===
echo.

REM Install deps if needed
pip install pyinstaller pyside6 numpy scipy pandas matplotlib plotly unidec -q

REM Install the package itself
pip install -e . -q

REM Run PyInstaller
pyinstaller --noconfirm --clean build\deconvovo.spec

if %ERRORLEVEL% EQU 0 (
    echo.
    echo === Build successful ===
    echo Output: dist\DeconVoVo\DeconVoVo.exe
    echo To distribute: zip the dist\DeconVoVo folder
) else (
    echo.
    echo === Build FAILED ===
    exit /b 1
)

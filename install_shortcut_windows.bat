@echo off
rem Double-click to create the "Spectrum Simulator" shortcut (with icon) on the Desktop and Start menu.
cd /d "%~dp0"
set "PY="
py -3 -c "import sys" >nul 2>nul && set "PY=py -3"
if not defined PY python -c "import sys" >nul 2>nul && set "PY=python"
if not defined PY (
  echo No se encontro Python. Instala Python 3.10 o superior desde https://www.python.org/downloads/
  echo   y marca "Add python.exe to PATH". Luego vuelve a ejecutar este archivo.
  echo Python was not found. Install Python 3.10+ from https://www.python.org/downloads/ and run this file again.
  start "" https://www.python.org/downloads/
  pause
  exit /b 1
)
%PY% launcher\install_shortcut.py
pause

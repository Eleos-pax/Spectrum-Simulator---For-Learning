@echo off
rem Starts Spectrum Simulator (used by the desktop shortcut).
title Spectrum Simulator
cd /d "%~dp0.."
set "PY="
rem Prefer the "py" launcher; plain "python" may be the Microsoft Store alias, so test that it really runs.
py -3 -c "import sys" >nul 2>nul && set "PY=py -3"
if not defined PY python -c "import sys" >nul 2>nul && set "PY=python"
if not defined PY (
  echo No se encontro Python. Instala Python 3.10 o superior desde https://www.python.org/downloads/
  echo   ^(marca "Add python.exe to PATH" durante la instalacion^) y vuelve a abrir el acceso directo.
  echo Python was not found. Install Python 3.10+ from https://www.python.org/downloads/
  start "" https://www.python.org/downloads/
  pause
  exit /b 1
)
%PY% "%~dp0launch.py"
if errorlevel 1 pause

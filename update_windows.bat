@echo off
rem Double-click to update Spectrum Simulator to the latest version from GitHub.
cd /d "%~dp0"
set PY=
where py >nul 2>nul && set PY=py -3
if not defined PY where python >nul 2>nul && set PY=python
if not defined PY (
  echo No se encontro Python. / Python was not found.
  pause
  exit /b 1
)
%PY% launcher\update.py
pause

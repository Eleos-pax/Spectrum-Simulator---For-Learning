@echo off
rem Starts Spectrum Simulator (used by the desktop shortcut).
cd /d "%~dp0.."
where py >nul 2>nul && (py -3 "%~dp0launch.py" & goto :eof)
where python >nul 2>nul && (python "%~dp0launch.py" & goto :eof)
echo No se encontro Python. Instala Python 3.10 o superior desde https://www.python.org/downloads/
echo   (marca "Add python.exe to PATH" durante la instalacion) y vuelve a abrir el acceso directo.
echo Python was not found. Install Python 3.10+ from https://www.python.org/downloads/
start "" https://www.python.org/downloads/
pause

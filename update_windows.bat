@echo off
rem ================================================================
rem  Actualizar Spectrum Simulator / Update Spectrum Simulator
rem  Pon este archivo en la carpeta de la app (donde esta app.py)
rem  y haz doble clic. Descarga la ultima version de GitHub y
rem  reemplaza los archivos. Conserva .venv y tus archivos.
rem ================================================================
setlocal
cd /d "%~dp0"
if not exist "%~dp0app.py" (
  echo Este archivo debe estar en la carpeta de Spectrum Simulator, junto a app.py.
  echo This file must be in the Spectrum Simulator folder, next to app.py.
  pause
  exit /b 1
)
set "TMPDIR=%TEMP%\spectrum_simulator_update"
if exist "%TMPDIR%" rmdir /s /q "%TMPDIR%"
mkdir "%TMPDIR%"

echo Descargando la ultima version... / Downloading the latest version...
powershell -NoProfile -ExecutionPolicy Bypass -Command "$ErrorActionPreference='Stop'; [Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -UseBasicParsing 'https://github.com/Eleos-pax/Spectrum-Simulator---For-Learning/archive/refs/heads/main.zip' -OutFile '%TMPDIR%\update.zip'; Expand-Archive -Force '%TMPDIR%\update.zip' '%TMPDIR%'"
if errorlevel 1 goto :fail

echo Copiando archivos... / Copying files...
robocopy "%TMPDIR%\Spectrum-Simulator---For-Learning-main" "%~dp0." /E /XD .venv .git /NFL /NDL /NJH /NJS /NP
if errorlevel 8 goto :fail

rmdir /s /q "%TMPDIR%"
echo.
echo Listo. Abre Spectrum Simulator con su icono.
echo Done. Open Spectrum Simulator from its icon.
pause
exit /b 0

:fail
echo.
echo No se pudo actualizar. Revisa tu conexion a internet y vuelve a intentarlo.
echo Update failed. Check your internet connection and try again.
pause
exit /b 1

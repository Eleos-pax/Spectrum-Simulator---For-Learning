"""Create a desktop shortcut (with icon) that starts Spectrum Simulator.

* Windows: "Spectrum Simulator.lnk" on the Desktop and in the Start menu.
* macOS:   "Spectrum Simulator.app" in ~/Applications, plus an alias on the Desktop.
* Linux:   a .desktop entry in the applications menu and on the Desktop.

Run it once:  python launcher/install_shortcut.py   (or double-click the install_shortcut_* file)
"""
from __future__ import annotations

import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NAME = "Spectrum Simulator"
LAUNCH = HERE / "launch.py"


def windows():
    q = lambda p: str(p).replace("'", "''")  # PowerShell single-quote escaping
    bat = HERE / "start_windows.bat"
    ps = f"""
$s = New-Object -ComObject WScript.Shell
$targets = @([Environment]::GetFolderPath('Desktop'), (Join-Path ([Environment]::GetFolderPath('Programs')) ''))
foreach ($dir in $targets) {{
  $lnk = $s.CreateShortcut((Join-Path $dir '{NAME}.lnk'))
  $lnk.TargetPath = '{q(bat)}'
  $lnk.WorkingDirectory = '{q(ROOT)}'
  $lnk.IconLocation = '{q(HERE / "icon.ico")},0'
  $lnk.WindowStyle = 7
  $lnk.Description = 'UV-Vis, IR, MS and NMR spectrum simulator'
  $lnk.Save()
}}
"""
    subprocess.check_call(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", ps])
    return "Escritorio y menú Inicio / Desktop and Start menu"


def macos():
    app = Path.home() / "Applications" / f"{NAME}.app"
    if app.exists():
        shutil.rmtree(app)
    (app / "Contents" / "MacOS").mkdir(parents=True)
    (app / "Contents" / "Resources").mkdir()
    shutil.copy(HERE / "icon.icns", app / "Contents" / "Resources" / "icon.icns")
    (app / "Contents" / "Info.plist").write_text(f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>{NAME}</string>
  <key>CFBundleDisplayName</key><string>{NAME}</string>
  <key>CFBundleIdentifier</key><string>org.eleospax.spectrumsimulator</string>
  <key>CFBundleExecutable</key><string>launcher</string>
  <key>CFBundleIconFile</key><string>icon</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleVersion</key><string>1.0</string>
</dict></plist>
""")
    exe = app / "Contents" / "MacOS" / "launcher"
    # open a Terminal window so the first-run installation is visible; closing it quits the app
    exe.write_text(f"""#!/bin/bash
osascript -e 'tell application "Terminal" to do script "\\"{sys.executable}\\" \\"{LAUNCH}\\""' -e 'tell application "Terminal" to activate'
""")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    desk = Path.home() / "Desktop" / NAME
    if desk.is_symlink() or desk.exists():
        desk.unlink()
    desk.symlink_to(app)
    return f"{app} + Escritorio / Desktop"


def linux():
    entry = f"""[Desktop Entry]
Type=Application
Name={NAME}
Comment=UV-Vis, IR, MS and NMR spectrum simulator
Exec="{sys.executable}" "{LAUNCH}"
Path={ROOT}
Icon={HERE / "icon.png"}
Terminal=true
Categories=Education;Science;Chemistry;
"""
    places = [Path.home() / ".local/share/applications"]
    desktop = Path(subprocess.run(["xdg-user-dir", "DESKTOP"], capture_output=True, text=True).stdout.strip() or Path.home() / "Desktop") \
        if shutil.which("xdg-user-dir") else Path.home() / "Desktop"
    if desktop.exists():
        places.append(desktop)
    for d in places:
        d.mkdir(parents=True, exist_ok=True)
        f = d / "spectrum-simulator.desktop"
        f.write_text(entry)
        f.chmod(0o755)
        if shutil.which("gio") and d == desktop:
            subprocess.run(["gio", "set", str(f), "metadata::trusted", "true"], capture_output=True)
    return ", ".join(str(p) for p in places)


def main():
    where = windows() if os.name == "nt" else macos() if sys.platform == "darwin" else linux()
    print(f"Acceso directo creado / Shortcut created: {where}")
    print("La primera vez que lo abras se instalarán las dependencias (requiere internet).")
    print("The first time you open it, the dependencies are installed (internet required).")


if __name__ == "__main__":
    main()

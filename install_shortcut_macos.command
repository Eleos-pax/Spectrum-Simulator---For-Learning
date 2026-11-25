#!/bin/bash
# Double-click to create "Spectrum Simulator.app" (with icon) in ~/Applications and on the Desktop.
cd "$(dirname "$0")"
if ! command -v python3 >/dev/null; then
  echo "No se encontró Python 3. Instálalo desde https://www.python.org/downloads/ y vuelve a abrir este archivo."
  echo "Python 3 was not found. Install it from https://www.python.org/downloads/ and open this file again."
  open https://www.python.org/downloads/
  exit 1
fi
python3 launcher/install_shortcut.py

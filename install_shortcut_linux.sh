#!/bin/bash
# Creates the "Spectrum Simulator" launcher (with icon) in the applications menu and on the Desktop.
cd "$(dirname "$0")"
command -v python3 >/dev/null || { echo "Python 3 is required: sudo apt install python3 python3-venv"; exit 1; }
python3 launcher/install_shortcut.py

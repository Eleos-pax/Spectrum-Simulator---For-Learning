#!/bin/bash
# Double-click to update Spectrum Simulator to the latest version from GitHub.
cd "$(dirname "$0")"
python3 launcher/update.py
read -n 1 -s -r -p "Pulsa una tecla para cerrar / Press any key to close"

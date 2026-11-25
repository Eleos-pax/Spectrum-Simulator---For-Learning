"""Start Spectrum Simulator and open it in the browser.

First run: creates a private virtual environment in `.venv/` next to the app and installs
`requirements.txt` into it (needs internet, a few minutes). Later runs start in seconds.
If the app is already running, it only opens the browser tab.

Usage:  python launcher/launch.py        (the desktop shortcut runs this for you)
"""
from __future__ import annotations

import hashlib
import os
import socket
import subprocess
import sys
import time
import urllib.request
import venv
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VENV = ROOT / ".venv"
REQ = ROOT / "requirements.txt"
STAMP = VENV / ".requirements.sha256"
PORT = 8765  # own port range, away from Streamlit's default 8501 used by other projects
APP_ID = "spectrum-simulator-for-learning"  # served from static/app_id.txt
WIN = os.name == "nt"


def say(msg: str):
    print(f"[Spectrum Simulator] {msg}", flush=True)


def venv_python() -> Path:
    return VENV / ("Scripts/python.exe" if WIN else "bin/python")


def ensure_environment():
    if not venv_python().exists():
        say("Primera ejecución: creando el entorno de Python… / First run: creating the Python environment…")
        venv.EnvBuilder(with_pip=True).create(VENV)
    digest = hashlib.sha256(REQ.read_bytes()).hexdigest()
    if not STAMP.exists() or STAMP.read_text().strip() != digest:
        say("Instalando dependencias (solo la primera vez, unos minutos)… / Installing dependencies (first time only)…")
        py = str(venv_python())
        subprocess.check_call([py, "-m", "pip", "install", "--upgrade", "pip"])
        subprocess.check_call([py, "-m", "pip", "install", "-r", str(REQ)])
        STAMP.write_text(digest)


def is_this_app(port: int) -> bool:
    """True only if the server on this port is Spectrum Simulator (not another Streamlit project)."""
    try:
        with urllib.request.urlopen(f"http://localhost:{port}/app/static/app_id.txt", timeout=1) as r:
            return r.read().decode().strip() == APP_ID
    except Exception:
        return False


def port_free(port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("127.0.0.1", port)) != 0


def main():
    port = None
    for p in range(PORT, PORT + 50):
        if is_this_app(p):
            say(f"La app ya está abierta en el puerto {p}. / The app is already running on port {p}.")
            webbrowser.open(f"http://localhost:{p}")
            return
        if port is None and port_free(p):
            port = p  # first free port; ports used by other programs are skipped
    ensure_environment()
    say(f"Iniciando en http://localhost:{port} … (cierra esta ventana para salir / close this window to quit)")
    proc = subprocess.Popen([str(venv_python()), "-m", "streamlit", "run", str(ROOT / "app.py"),
                             "--server.port", str(port), "--server.headless", "true",
                             "--browser.gatherUsageStats", "false"], cwd=str(ROOT))
    for _ in range(120):
        if is_this_app(port):
            webbrowser.open(f"http://localhost:{port}")
            break
        if proc.poll() is not None:
            break
        time.sleep(0.5)
    try:
        proc.wait()
    except KeyboardInterrupt:
        proc.terminate()


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError:
        say("No se pudieron instalar las dependencias. Revisa tu conexión a internet y vuelve a intentarlo.")
        say("Could not install the dependencies. Check your internet connection and try again.")
        if WIN:
            input("Pulsa Enter para cerrar / Press Enter to close…")
        sys.exit(1)

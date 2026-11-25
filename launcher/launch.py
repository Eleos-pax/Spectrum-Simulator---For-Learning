"""Start Spectrum Simulator and open it in the browser.

First run: creates a private virtual environment in `.venv/` next to the app and installs
`requirements.txt` into it (needs internet, a few minutes). Later runs start in seconds.
If the app is already running, it only opens the browser tab.

Everything printed here is also written to `launcher/launch.log`, so a failed start can be diagnosed.

Usage:  python launcher/launch.py        (the desktop shortcut runs this for you)
"""
from __future__ import annotations

import hashlib
import os
import shutil
import socket
import subprocess
import sys
import time
import traceback
import urllib.request
import venv
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VENV = ROOT / ".venv"
REQ = ROOT / "requirements.txt"
STAMP = VENV / ".requirements.sha256"
LOG = Path(__file__).resolve().parent / "launch.log"
LOCK = Path(__file__).resolve().parent / ".launch.lock"
PORT = 8765  # own port range, away from Streamlit's default 8501 used by other projects
APP_ID = "spectrum-simulator-for-learning"  # served from static/app_id.txt
WIN = os.name == "nt"


def say(msg: str):
    line = f"[Spectrum Simulator] {msg}"
    print(line, flush=True)
    try:
        with LOG.open("a", encoding="utf-8") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + line + "\n")
    except OSError:
        pass


def pause():
    if sys.stdin and sys.stdin.isatty():
        try:
            input("Pulsa Enter para cerrar / Press Enter to close… ")
        except EOFError:
            pass


def venv_python() -> Path:
    return VENV / ("Scripts/python.exe" if WIN else "bin/python")


def run_logged(cmd):
    """Run a command, showing its output and copying it into the log."""
    say("$ " + " ".join(cmd))
    with LOG.open("a", encoding="utf-8") as f:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
        for line in proc.stdout:
            print(line, end="", flush=True)
            f.write(line)
        if proc.wait() != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)


def ensure_environment():
    digest = hashlib.sha256(REQ.read_bytes()).hexdigest()
    if venv_python().exists() and STAMP.exists() and STAMP.read_text().strip() == digest:
        return
    for attempt in (1, 2):
        try:
            if not venv_python().exists():
                say("Primera ejecución: creando el entorno de Python… / First run: creating the Python environment…")
                venv.EnvBuilder(with_pip=True, clear=True).create(VENV)
            say("Instalando dependencias (solo la primera vez, unos minutos; no cierres esta ventana)…")
            say("Installing dependencies (first time only, a few minutes; do not close this window)…")
            py = str(venv_python())
            run_logged([py, "-m", "pip", "install", "--upgrade", "pip"])
            run_logged([py, "-m", "pip", "install", "-r", str(REQ)])
            STAMP.write_text(digest)
            return
        except (subprocess.CalledProcessError, OSError):
            if attempt == 2:
                raise
            say("La instalación falló; se borra el entorno y se reintenta… / Install failed; rebuilding the environment…")
            shutil.rmtree(VENV, ignore_errors=True)


def get(url: str) -> str | None:
    try:
        with urllib.request.urlopen(url, timeout=1) as r:
            return r.read().decode(errors="replace")
    except Exception:
        return None


def is_this_app(port: int) -> bool:
    """True only if the server on this port is Spectrum Simulator (not another Streamlit project)."""
    return (get(f"http://localhost:{port}/app/static/app_id.txt") or "").strip() == APP_ID


def port_free(port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("127.0.0.1", port)) != 0


def another_launcher_running() -> bool:
    """A second click while the first one is still installing would break the installation."""
    try:
        return LOCK.exists() and time.time() - LOCK.stat().st_mtime < 20 * 60
    except OSError:
        return False


def main():
    say(f"Python {sys.version.split()[0]} · {sys.executable}")
    say(f"Carpeta / Folder: {ROOT}")
    for p in range(PORT, PORT + 50):
        if is_this_app(p):
            say(f"La app ya está abierta en el puerto {p}. / The app is already running on port {p}.")
            webbrowser.open(f"http://localhost:{p}")
            return
    if another_launcher_running():
        say("Spectrum Simulator ya se está iniciando en otra ventana; espera a que termine.")
        say("Spectrum Simulator is already starting in another window; please wait for it.")
        time.sleep(6)
        return
    LOCK.write_text(str(os.getpid()))
    try:
        ensure_environment()
        port = next(p for p in range(PORT, PORT + 50) if port_free(p))
        say(f"Iniciando en http://localhost:{port} … (cierra esta ventana para salir / close this window to quit)")
        proc = subprocess.Popen([str(venv_python()), "-m", "streamlit", "run", str(ROOT / "app.py"),
                                 "--server.port", str(port), "--server.headless", "true",
                                 "--browser.gatherUsageStats", "false"], cwd=str(ROOT))
        opened = False
        for i in range(240):  # up to 2 minutes
            if proc.poll() is not None:
                break
            healthy = get(f"http://localhost:{port}/_stcore/health") == "ok"
            # open as soon as the app is identified; if the identity file cannot be read, open anyway after 10 s
            if is_this_app(port) or (healthy and i > 20):
                webbrowser.open(f"http://localhost:{port}")
                say(f"Abierta en el navegador: http://localhost:{port} / Opened in the browser.")
                opened = True
                break
            time.sleep(0.5)
    finally:
        try:
            LOCK.unlink()
        except OSError:
            pass
    if not opened and proc.poll() is not None:
        raise RuntimeError(f"Streamlit se cerró al iniciar (código {proc.returncode}). / Streamlit exited on start.")
    if not opened:
        say(f"Abre manualmente / Open manually: http://localhost:{port}")
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
        say(f"Detalles en / Details in: {LOG}")
        pause()
        sys.exit(1)
    except Exception:
        say("Error inesperado / Unexpected error:\n" + traceback.format_exc())
        say(f"Envía este archivo para revisarlo / Send this file for review: {LOG}")
        pause()
        sys.exit(1)

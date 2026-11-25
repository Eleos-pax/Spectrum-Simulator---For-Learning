"""Update Spectrum Simulator to the latest version on GitHub.

* If the folder is a git clone, it runs `git pull`.
* Otherwise (installed from the ZIP) it downloads the latest ZIP of the `main` branch and copies the new
  files over the installed ones. The Python environment (.venv) and your own files are kept; the
  launcher reinstalls dependencies automatically if requirements.txt changed.

Run it with the update_* file for your system, or:  python launcher/update.py
"""
from __future__ import annotations

import io
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPO = "Eleos-pax/Spectrum-Simulator---For-Learning"
ZIP_URL = f"https://github.com/{REPO}/archive/refs/heads/main.zip"
KEEP = {".venv", ".git"}  # never touched


def say(msg: str):
    print(f"[Spectrum Simulator] {msg}", flush=True)


def update_git() -> bool:
    if not (ROOT / ".git").exists() or not shutil.which("git"):
        return False
    say("Actualizando con git… / Updating with git…")
    subprocess.check_call(["git", "-C", str(ROOT), "pull", "--ff-only"])
    return True


def update_zip():
    say("Descargando la última versión… / Downloading the latest version…")
    with urllib.request.urlopen(ZIP_URL, timeout=60) as r:
        data = r.read()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp).resolve()
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for name in z.namelist():  # refuse paths that would leave the folder
                target = (tmp / name).resolve()
                if not str(target).startswith(str(tmp)):
                    raise RuntimeError(f"unsafe path in archive: {name}")
            z.extractall(tmp)
        src = next(p for p in tmp.iterdir() if p.is_dir())
        n = 0
        for f in src.rglob("*"):
            rel = f.relative_to(src)
            if rel.parts[0] in KEEP or f.is_dir():
                continue
            dest = ROOT / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dest)
            n += 1
    say(f"{n} archivos actualizados. / {n} files updated.")


def main():
    try:
        if not update_git():
            update_zip()
        say("¡Listo! Abre Spectrum Simulator con su icono (cierra antes la ventana si estaba abierta).")
        say("Done! Open Spectrum Simulator from its icon (close its window first if it was running).")
    except Exception as e:  # noqa: BLE001 - show any problem to the user
        say(f"No se pudo actualizar / Update failed: {e}")
        say("Revisa tu conexión a internet y vuelve a intentarlo. / Check your internet connection and try again.")
        sys.exit(1)


if __name__ == "__main__":
    main()

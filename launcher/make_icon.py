"""Generate the app icon (PNG, ICO for Windows, ICNS for macOS). Requires Pillow."""
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter

OUT = Path(__file__).resolve().parent
N = 1024  # draw large, downscale for smooth edges


def make():
    img = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    # rounded-square background with a vertical gradient (deep blue -> indigo)
    grad = np.zeros((N, N, 4), dtype=np.uint8)
    t = np.linspace(0, 1, N)[:, None]
    grad[..., 0] = (15 + 35 * t).astype(np.uint8)
    grad[..., 1] = (40 + 10 * t).astype(np.uint8)
    grad[..., 2] = (110 + 60 * t).astype(np.uint8)
    grad[..., 3] = 255
    bg = Image.fromarray(grad, "RGBA")
    mask = Image.new("L", (N, N), 0)
    ImageDraw.Draw(mask).rounded_rectangle((40, 40, N - 40, N - 40), radius=210, fill=255)
    img.paste(bg, (0, 0), mask)
    d = ImageDraw.Draw(img)

    # visible-spectrum strip along the bottom
    x0, x1, y0, y1 = 150, N - 150, 760, 820
    for x in range(x0, x1):
        lam = 380 + (x - x0) / (x1 - x0) * 370
        if lam < 440: r, g, b = (440 - lam) / 60, 0, 1
        elif lam < 490: r, g, b = 0, (lam - 440) / 50, 1
        elif lam < 510: r, g, b = 0, 1, (510 - lam) / 20
        elif lam < 580: r, g, b = (lam - 510) / 70, 1, 0
        elif lam < 645: r, g, b = 1, (645 - lam) / 65, 0
        else: r, g, b = 1, 0, 0
        d.line([(x, y0), (x, y1)], fill=(int(255 * r), int(255 * g), int(255 * b), 255))

    # spectrum trace: baseline with sharp peaks (like an IR/NMR/MS plot)
    xs = np.linspace(150, N - 150, 900)
    u = (xs - 150) / (N - 300)
    peaks = [(0.18, 0.025, 0.55), (0.36, 0.018, 0.85), (0.52, 0.03, 0.40), (0.70, 0.014, 1.0), (0.84, 0.022, 0.6)]
    y = sum(h / (1 + ((u - c) / w) ** 2) for c, w, h in peaks)
    ys = 690 - 480 * y / y.max()
    pts = list(zip(xs, ys))
    glow = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    ImageDraw.Draw(glow).line(pts, fill=(120, 200, 255, 200), width=46, joint="curve")
    img.alpha_composite(glow.filter(ImageFilter.GaussianBlur(18)))
    d.line(pts, fill=(255, 255, 255, 255), width=22, joint="curve")
    # axis
    d.line([(150, 700), (N - 150, 700)], fill=(255, 255, 255, 140), width=8)

    big = img.resize((512, 512), Image.LANCZOS)
    big.save(OUT / "icon.png")
    img.save(OUT / "icon.ico", sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    img.save(OUT / "icon.icns")
    # small icon for the browser tab
    img.resize((64, 64), Image.LANCZOS).save(OUT / "favicon.png")


if __name__ == "__main__":
    make()

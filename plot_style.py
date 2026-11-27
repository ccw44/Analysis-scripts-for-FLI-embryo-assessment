"""Text setup shared by the figure scripts: one concrete font family.

Requesting a generic family makes matplotlib fall back per character, and each
fallback starts a new font run, which splits a word into several text objects in
the exported PDF.  Setting one concrete family avoids that.

Arial is tried first if it is installed; otherwise the script falls back through
metrically compatible faces.  The family actually used is printed by every
figure script.
"""
import os

from matplotlib import font_manager as fm
import matplotlib.pyplot as plt

CANDIDATES = ["Arial", "Liberation Sans", "Nimbus Sans", "Helvetica", "DejaVu Sans"]

FONT_SEARCH = [os.path.expanduser("~/.fonts")]


def register():
    """Make any locally available Arial/Liberation faces visible to matplotlib."""
    for d in FONT_SEARCH:
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn.lower().endswith((".ttf", ".otf")):
                try:
                    fm.fontManager.addfont(os.path.join(d, fn))
                except Exception:
                    pass


def family():
    register()
    names = {f.name for f in fm.fontManager.ttflist}
    for c in CANDIDATES:
        if c in names:
            return c
    return "DejaVu Sans"


def apply(base_size=None):
    """Set a single-family, fallback-free text configuration. Returns the family used."""
    fam = family()
    upd = {
        "font.family": fam,          # concrete family: no generic list, no per-glyph fallback
        "font.sans-serif": [fam],
        "pdf.fonttype": 42,          # TrueType embedded -> editable text in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "pdf.use14corefonts": False,
    }
    if base_size is not None:
        upd["font.size"] = base_size
    plt.rcParams.update(upd)
    return fam


def path():
    """Absolute path of the regular face, for ov.style(font_path=...)."""
    register()
    fam = family()
    try:
        return fm.findfont(fm.FontProperties(family=fam), fallback_to_default=False)
    except Exception:
        return None

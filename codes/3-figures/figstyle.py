#!/usr/bin/env python3
"""
Shared figure style for the ApoE reverse-QTY manuscript.

Palette: native = deep forest green, rQTY = deep gold.  Chosen to match the
green ribbon used for structures in our earlier QTY papers, and validated for
colour-vision deficiency: the pair separates by dE 18.7 (protanopia) and 27.0
(normal vision) in OKLab, and by a factor 3.0 in greyscale luminance, so the
figures survive black-and-white printing.

Identity is never carried by colour alone -- every multi-series panel also
carries a legend, and the substitution types in Figure 1a differ in marker
shape as well as hue.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))

def _pick(*cands):
    """First candidate that exists; the last is the documented default.
    Lets the scripts run unchanged from the project folder or from the
    published repository, and any of them can still be overridden by the
    matching environment variable."""
    for c in cands:
        if c and os.path.isdir(c):
            return c
    return cands[-1]

ROOT     = os.environ.get('APOE_ROOT') or os.path.normpath(os.path.join(HERE, '..', '..'))
ANALYSIS = _pick(os.environ.get('ANALYSIS_DIR'),
                 ANALYSIS,
                 os.path.join(ROOT, 'data', 'analysis_json'))
AF3DIR   = _pick(os.environ.get('AF3_DIR'),
                 os.path.join(ROOT, 'data', 'AF3', 'AlphaFold3'),
                 os.path.join(ROOT, 'data', 'AF3'))
STRUCTD  = _pick(os.environ.get('STRUCT_DIR'),
                 os.path.join(ROOT, 'data', 'Structures'))
MDDIR    = _pick(os.environ.get('MD_DIR'),
                 os.path.join(ROOT, 'data', 'MD'))

# --- series -------------------------------------------------------------
NAT = '#1F6B33'      # native            deep forest green
RQ  = '#C8941B'      # rQTY variant      deep gold

# --- Figure 1a substitution types (shape also encodes identity) ---------
SUB_COLOR = {'Q': '#C8941B',   # Q->L  gold     circle
             'T': '#2E5C9A',   # T->V  slate    square
             'Y': '#8C3A6B'}   # Y->F  plum     triangle
SUB_MARK  = {'Q': 'o', 'T': 's', 'Y': '^'}
LDLR = '#6B4B9A'     # LDLR recognition block

# --- neutrals -----------------------------------------------------------
HELIX_FILL = '#DCDBD6'
HELIX_EDGE = '#A9A8A3'
BAND       = '#EFEEEA'   # shaded region behind the N-terminal arm
RULE       = '#8F8E8A'   # zero lines, axis rules
MUTED      = '#52514E'   # secondary text
INK        = '#25292E'   # primary text

SEQ_GREEN = 'Greens'     # sequential ramp for the Supplementary S1 matrices

RC = {
    'font.family': 'DejaVu Sans', 'font.size': 8,
    'axes.linewidth': 0.8, 'axes.edgecolor': INK,
    'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': INK, 'ytick.color': INK,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'axes.spines.top': False, 'axes.spines.right': False,
    'savefig.dpi': 400, 'figure.dpi': 120,
    'legend.frameon': False,
}


def style():
    plt.rcParams.update(RC)


def panel(ax, tag, dx=-0.235, dy=1.20):
    ax.text(dx, dy, tag, transform=ax.transAxes, fontsize=10,
            fontweight='bold', va='top', color=INK)


HEL = [('H0', 11, 19), ('H1', 25, 41), ('H1b', 45, 52),
       ('H2', 55, 81), ('H3', 87, 124), ('H4', 131, 167)]
ARM = (11, 24)
SUBS = [16, 17, 18, 36, 41, 46, 48, 55, 57, 58, 67, 74, 81, 89, 98, 101,
        117, 118, 123, 156, 162, 163]
NOSUB = [21, 24, 42, 83, 128, 130]


def outdir():
    """Where the figures go: $FIGOUT if set, else ../new beside this script."""
    d = os.environ.get('FIGOUT') or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', 'new')
    os.makedirs(d, exist_ok=True)
    return os.path.normpath(d)


def save(fig, name):
    """Write one figure in every requested format, into outdir().

    FIGFORMATS   comma-separated extensions, default 'png'
    FIGDPI       raster resolution, default 400

    A vector PDF is resolution-independent, which satisfies a line-art dpi
    requirement by construction.  EPS is deliberately not offered: several
    panels use partial transparency and the PostScript backend flattens
    partially transparent artists to opaque, which would change the figures.
    """
    fmts = [f.strip().lower().lstrip('.') for f in
            os.environ.get('FIGFORMATS', 'png').split(',') if f.strip()]
    dpi = float(os.environ.get('FIGDPI', 400))
    out = []
    for f in fmts:
        path = os.path.join(outdir(), f'{name}.{f}')
        kw = dict(bbox_inches='tight', facecolor='white', dpi=dpi)
        if f in ('tif', 'tiff'):
            kw['pil_kwargs'] = {'compression': 'tiff_lzw'}
        fig.savefig(path, **kw)
        if f in ('tif', 'tiff'):
            _flatten_tiff(path, dpi)
        out.append(path)
    return out


def _flatten_tiff(path, dpi):
    """Composite the alpha channel onto white and store RGB.

    Matplotlib writes RGBA TIFF.  Editorial systems and print production
    expect RGB or CMYK without an alpha channel, and a stray alpha channel
    is a common cause of a figure being bounced or rendered with a black
    background, so it is removed here rather than left to chance.
    """
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    im = Image.open(path)
    if im.mode != 'RGBA':
        return
    bg = Image.new('RGB', im.size, (255, 255, 255))
    bg.paste(im, mask=im.split()[3])
    im.close()
    bg.save(path, dpi=(dpi, dpi), compression='tiff_lzw')

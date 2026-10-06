"""Shared figure system: 18 cm (7.09 in) page width, Liberation Sans (Arial metrics), palette from the dataviz reference."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

INK = '#0b0b0b'; INK2 = '#52514e'; MUTED = '#898781'; GRID = '#e1e0d9'; AXIS = '#c3c2b7'; SURF = '#fcfcfb'; MID = '#f0efec'
BLUE = '#2a78d6'; ORANGE = '#eb6834'; AQUA = '#1baf7a'; YELLOW = '#eda100'; MAGENTA = '#e87ba4'; GREEN = '#008300'; VIOLET = '#4a3aa7'; RED = '#e34948'
CAT = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]
SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
DIV = LinearSegmentedColormap.from_list('div', ['#184f95', '#3987e5', '#9ec5f4', MID, '#f6b79b', ORANGE, '#a8401a'])
SEQCMAP = LinearSegmentedColormap.from_list('seq', ['#f4f8fe'] + SEQ)
LDCMAP = LinearSegmentedColormap.from_list('ld', ['#d9d7cf', '#9ec5f4', '#1baf7a', '#eda100', '#e34948'])
ANC = {'EU': BLUE, 'SA': ORANGE, 'EA': AQUA, 'AF': VIOLET, 'AA': MAGENTA, 'HS': YELLOW}
ANCLAB = {'EU': 'European', 'SA': 'South Asian', 'EA': 'East Asian', 'AF': 'African', 'AA': 'African American', 'HS': 'Hispanic/Latino'}
W = 7.09  # inches, full page width

plt.rcParams.update({
    'font.family': 'Liberation Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5, 'axes.edgecolor': AXIS, 'axes.linewidth': 0.6, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.labelcolor': INK2,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.dpi': 400, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'lines.linewidth': 1.2,
    'legend.frameon': False, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.titlecolor': INK,
})


def tag(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', color=INK, va='bottom', ha='left')


def clean(ax, left=True, bottom=True):
    ax.spines['left'].set_visible(left); ax.spines['bottom'].set_visible(bottom)


PRINT_W_IN, PRINT_H_IN = 7.087, 9.45   # 180 mm x 240 mm: full-width print area used for the proof


def proof(fig, name):
    """Print-size proof: scale of the tight-cropped figure when placed at 180 mm width (and at most 240 mm height), and the
    effective size in points of every visible text element at that scale. Appends to figures/proof_report.csv."""
    import csv, os
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    bb = fig.get_tightbbox(r); w_in, h_in = bb.width, bb.height
    s = min(PRINT_W_IN / w_in, PRINT_H_IN / h_in, 1.0)
    sizes = []
    for t in fig.findobj(lambda a: hasattr(a, 'get_fontsize') and hasattr(a, 'get_text')):
        try:
            if t.get_visible() and t.get_text().strip(): sizes.append((t.get_fontsize() * s, t.get_text().strip()[:40]))
        except Exception: pass
    sizes.sort()
    row = dict(figure=name, tight_width_in=round(w_in, 2), tight_height_in=round(h_in, 2), scale_at_180mm=round(s, 3),
               min_effective_pt=round(sizes[0][0], 2) if sizes else None, n_text=len(sizes),
               n_below_6pt=sum(1 for z, _ in sizes if z < 6 - 1e-6), smallest_example=sizes[0][1] if sizes else '')
    path = 'figures/proof_report.csv'; rows = []
    if os.path.exists(path): rows = [x for x in csv.DictReader(open(path)) if x['figure'] != name]
    rows.append({k: str(v) for k, v in row.items()})
    with open(path, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(row)); wr.writeheader(); wr.writerows(rows)
    return row


MIN_PT = 6.2   # minimum text size (pt) at the 180 mm print width for main-text figures


def enforce_min_font(fig, minpt=MIN_PT):
    """Raise every visible text element (including tick labels) below minpt to minpt and let figure-level notes wrap within the figure width."""
    fig.canvas.draw()
    for ax in fig.axes:
        for axis, which in ((ax.xaxis, 'x'), (ax.yaxis, 'y')):
            labs = [l for l in axis.get_ticklabels() if l.get_text()]
            if labs and min(l.get_fontsize() for l in labs) < minpt: ax.tick_params(axis=which, labelsize=minpt)
    for t in fig.findobj(lambda a: hasattr(a, 'get_fontsize') and hasattr(a, 'get_text')):
        try:
            if t.get_text().strip() and t.get_fontsize() < minpt: t.set_fontsize(minpt)
        except Exception: pass
    for t in fig.texts:
        if len(t.get_text()) > 90: t.set_wrap(True)


MIN_PT_FIG = {'Fig2_mendelian_randomization': 7.2}   # v10: larger minimum for the dense panels (Figure 2c, 3d)
def save(fig, name):
    if name.startswith('Fig') and not name.startswith('FigS'): enforce_min_font(fig, MIN_PT_FIG.get(name, MIN_PT))
    proof(fig, name)
    fig.savefig(f'figures/{name}.pdf', bbox_inches='tight'); fig.savefig(f'figures/{name}.png', bbox_inches='tight', dpi=400)
    fig.savefig(f'figures/{name}.svg', bbox_inches='tight')

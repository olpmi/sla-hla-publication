"""Geometry shared by the structural panels of Figures 5 and 6.

These two functions were previously reached by ``ast.parse``-ing
``figures_structural_final.py`` and ``exec``-ing two of its function nodes into a
sandbox dictionary, because importing that module pulls in PyMOL. They are
reproduced here verbatim as ordinary importable functions, so the figure build
has no AST extraction, no ``exec``, and no PyMOL dependency.

The full rendering module (which *does* drive PyMOL) lives in
``slahla_pub.structural`` and is only used by the optional structural workflow.
"""
from __future__ import annotations

import math

import numpy as np

#: Pixel size of the rendered PyMOL panels. The archived rasters are 1600x1200.
PANEL_PX = (1600, 1200)


def project(view: list[float], xyz, size=PANEL_PX) -> tuple[float, float]:
    """Model coordinate -> image pixel, orthoscopic. Calibrated to 0.14 px."""
    W, H = size
    R = np.array(view[0:9]).reshape(3, 3).T
    t = np.array(view[9:12]); origin = np.array(view[12:15]); fov = view[17]
    cam = R @ (np.asarray(xyz, float) - origin) + t
    scale = (H / 2) / (-t[2] * math.tan(math.radians(fov) / 2))
    return W / 2 + scale * cam[0], H / 2 - scale * cam[1]


def place_labels(anchors: list[tuple[float, float]], mask: np.ndarray,
                 size=PANEL_PX, margin=30, gap=58, x_margin=140):
    """Label positions outside the molecule, one per anchor, without overlaps.

    Each label is pushed radially outward from the panel's ink centroid until it
    clears the silhouette, then nudged along the perimeter until it clears its
    neighbours. Deterministic: anchors are processed in angle order.
    """
    W, H = size
    ys, xs = np.nonzero(mask)
    cx, cy = xs.mean(), ys.mean()
    order = sorted(range(len(anchors)),
                   key=lambda i: math.atan2(anchors[i][1] - cy, anchors[i][0] - cx))
    placed: list[tuple[float, float]] = [(0.0, 0.0)] * len(anchors)
    taken: list[tuple[float, float]] = []
    for i in order:
        ax, ay = anchors[i]
        theta = math.atan2(ay - cy, ax - cx)
        if not math.isfinite(theta):
            theta = 0.0
        for dtheta in [0.0] + [s * d for d in np.arange(0.02, 1.2, 0.02) for s in (1, -1)]:
            th = theta + dtheta
            r = math.hypot(ax - cx, ay - cy)
            for extra in np.arange(0, max(W, H), 6):
                px = cx + (r + extra) * math.cos(th)
                py = cy + (r + extra) * math.sin(th)
                if not (x_margin < px < W - x_margin and margin < py < H - margin):
                    break
                iy, ix = int(round(py)), int(round(px))
                win = mask[max(0, iy - 14):iy + 14, max(0, ix - 30):ix + 30]
                if win.size and win.any():
                    continue
                if any(abs(px - qx) < gap * 1.95 and abs(py - qy) < gap * 0.63
                       for qx, qy in taken):
                    continue
                placed[i] = (px, py); taken.append((px, py)); break
            if placed[i] != (0.0, 0.0):
                break
        if placed[i] == (0.0, 0.0):                       # last resort: on the rim
            px = min(max(ax, x_margin), W - x_margin)
            py = min(max(ay, margin), H - margin)
            placed[i] = (px, py); taken.append((px, py))
    return placed


_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def format_permutation_p(p: float, n_at_least_as_extreme: int, n_shifts: int) -> str:
    """Render a circular-shift p-value the way the publication renders it.

    The estimator is ``(extreme + 1) / (shifts + 1)``, so when no shift is at
    least as extreme as the observed statistic the value is a *floor*, not a
    point estimate, and is shown as an upper bound at one significant figure.
    Otherwise two significant figures are shown.

    This replaces the two hard-coded strings that the delivery script carried,
    which would have gone silently stale if the statistics were ever revised.
    """
    if n_at_least_as_extreme == 0:
        bound = 1.0 / (n_shifts + 1)
        exp = math.floor(math.log10(bound))
        mant = bound / 10 ** exp
        mant = math.ceil(mant * 10) / 10          # round the bound upward
        if mant >= 10:
            mant, exp = mant / 10, exp + 1
        mant_s = f"{mant:g}"
        return f"p < {mant_s} × 10{str(exp).translate(_SUP)}"
    exp = math.floor(math.log10(p))
    mant = p / 10 ** exp
    return f"p = {mant:.1f} × 10{str(exp).translate(_SUP)}"


def format_permutation_p_table(p: float, n_at_least_as_extreme: int, n_shifts: int) -> str:
    """The plain rendering the published statistics table carries.

    The table and the figure annotate the same p-value differently: the table
    uses four significant figures (``p = 0.007099``) and an exponential bound
    (``p < 1e-04``), the figure uses scientific notation with superscripts
    (``p = 6.0 x 10^-4``). Both are reproduced rather than harmonised, so the
    regenerated outputs match the approved materials exactly.
    """
    if n_at_least_as_extreme == 0:
        return f"p < {1.0 / (n_shifts + 1):.0e}"
    return f"p = {p:.4g}"

"""Final publication figures: Figures 1-6 and S1-S10.

Derived from ``code/build_figures_final.py`` in the author-approved manuscript
package, which is the authoritative source for appearance. Every measurement,
colour, size, panel letter and layout constant is preserved verbatim, including
the deliberate omission of numeric cell labels from the large locus-confusion
matrices in S2C and S6C. No plotted quantity differs.

Four choices are not obvious from the code alone:

- Generated linkage/ordering/edge/label-audit CSVs and
  ``panel_letter_validation.json`` go to ``outputs/``, never back into
  ``source_data/``: the archived inputs are read-only and ``make verify``
  asserts they are byte-unchanged after a run.
- FigureS10's two significance labels are derived from
  ``eplet_statistics_v6.csv`` via ``figure_helpers.format_permutation_p``. The
  p-value belongs to ``profile_difference`` (the equal-position contrast), not
  to ``difference``; see ``source_data/EPLET_STATISTIC_DEFINITIONS.md``.
- Structural panel titles come from ``structural_panel_manifest_v6.csv``, not
  the superseded ``..._original_measurements.csv``. The two are asserted to
  agree on the allele fields actually consumed before the newer one is used.
- ``project()`` and ``place_labels()`` are imported from ``figure_helpers``
  rather than exec'd out of a PyMOL-importing module, so drawing a figure never
  requires PyMOL.
"""
from __future__ import annotations

import io
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from PIL import Image

from . import eplet_mask, runlog
from .figure_helpers import format_permutation_p, place_labels, project
from .paths import AUDIT, FIGURES, PANELS, SOURCE_DATA, TABLES as TABLES_DIR, ensure_outputs

D = SOURCE_DATA          # read-only archived inputs
F = FIGURES              # regenerated artwork
WORK = PANELS            # archived PyMOL rasters + render_report.json

C = {'A': '#D55E00', 'B': '#0072B2', 'C': '#009E73', 'DRB1': '#882255', 'DQB1': '#E69F00',
     'DRB3': '#CC79A7', 'DRB4': '#56B4E9', 'DRB5': '#999933', 'DPB1': '#666666'}
CLASS = ['#386A9A', '#A8682A']
EP = '#C43B3B'
OTHER = '#657786'
HLA = '#B23A94'
SLA = '#409F54'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.labelsize': 9,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.linewidth': .6,
    'xtick.major.width': .6, 'ytick.major.width': .6,
})


def letter(fig, ax, t, extras=()):
    if not hasattr(fig, '_panel_specs'):
        fig._panel_specs = []
    fig._panel_specs.append((ax, t, tuple(extras)))


def save(fig, name, tight=True):
    """Place panel letters above-and-left of all panel content, then write the figure.

    ``tight=False`` writes the full canvas. The five single-panel supplementary
    figures (S1, S3, S4, S5, S7) were delivered that way -- their approved
    rasters are exactly ``figsize x 350 dpi`` -- whereas every multi-panel figure
    is cropped to its ink. Matching this is what makes the regenerated artwork
    agree with the approved artwork rather than being a tighter crop of it.
    """
    from matplotlib.transforms import Bbox
    ensure_outputs()
    fig.canvas.draw(); renderer = fig.canvas.get_renderer(); specs = []
    for ax, t, extras in getattr(fig, '_panel_specs', []):
        artists = [ax] + list(getattr(ax, '_panel_extras', [])) + list(extras)
        boxes = [a.get_tightbbox(renderer) for a in artists]
        boxes = [b for b in boxes if b is not None]
        # Retain the complete subplot slot as well as visible text and axes.
        boxes.append(ax.get_position(original=True).transformed(fig.transFigure))
        b = Bbox.union(boxes).transformed(fig.transFigure.inverted())
        specs.append((ax, t, b))
    letters = []
    for ax, t, b in specs:
        pos = ax.get_position(original=True)
        samecol = [q for a, _, q in specs if abs(a.get_position(original=True).x0 - pos.x0) < .02]
        samerow = [q for a, _, q in specs if abs(a.get_position(original=True).y1 - pos.y1) < .02]
        x = min(q.x0 for q in samecol) - 4 / (72 * fig.get_figwidth())
        y = max(q.y1 for q in samerow) + 4 / (72 * fig.get_figheight())
        letters.append((t, fig.text(x, y, t, fontsize=12, weight='normal', ha='right', va='bottom'), b))
    fig.canvas.draw(); renderer = fig.canvas.get_renderer()
    audit = []
    for t, artist, b in letters:
        lb = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        assert lb.x1 < b.x0 and lb.y0 > b.y1, (name, t, lb, b)
        audit.append({'figure': name, 'panel': t, 'weight': artist.get_fontweight(),
                      'above_all_panel_content': True, 'left_of_all_panel_content': True})
    if audit:
        path = AUDIT / 'panel_letter_validation.json'
        previous = json.loads(path.read_text()) if path.exists() else []
        previous = [r for r in previous if r['figure'] != name]
        path.write_text(json.dumps(previous + audit, indent=2))
    crop = {'bbox_inches': 'tight', 'pad_inches': .065} if tight else {}
    for ext in ['png', 'pdf', 'svg']:
        b = io.BytesIO()
        fig.savefig(b, format=ext, dpi=350, facecolor='white', **crop)
        (F / (name + '.' + ext)).write_bytes(b.getvalue())
    plt.close(fig)


def mat_cluster(name):
    m = pd.read_csv(D / (name + '.csv'), index_col=0)
    arr = m.to_numpy()
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        # Locus pairs with no DIAMOND alignment are all-NaN by design and are
        # drawn as grey cells; nanmean's empty-slice warning is expected here.
        warnings.simplefilter('ignore', RuntimeWarning)
        sym = np.nanmean(np.stack([arr, arr.T]), axis=0)
    dist = np.where(np.isfinite(sym), 100 - sym, 100); np.fill_diagonal(dist, 0)
    z = linkage(squareform(dist, checks=True), method='average')
    return m, z, leaves_list(z)


def clustered(fig, m, z, order, rect):
    x, y, w, h = rect; left = .085; top = .065
    ax = fig.add_axes([x + left, y, w - left - .07, h - top])
    dt = fig.add_axes([x + left, y + h - top + .006, w - left - .07, top - .01])
    dl = fig.add_axes([x, y, left - .006, h - top])
    dendrogram(z, ax=dt, no_labels=True, color_threshold=0, above_threshold_color='#555555'); dt.axis('off')
    dendrogram(z, ax=dl, no_labels=True, orientation='left', color_threshold=0, above_threshold_color='#555555')
    dl.invert_yaxis(); dl.axis('off')
    cmap = plt.get_cmap('viridis').with_extremes(bad='#D9D9D9')
    im = ax.imshow(m.values[np.ix_(order, order)], vmin=0, vmax=100, cmap=cmap, aspect='auto', interpolation='none')
    labels = [str(m.index[i]) if str(m.index[i]).startswith('SLA') else 'HLA-' + str(m.index[i]) for i in order]
    ax.set(xticks=range(len(order)), yticks=range(len(order)), xticklabels=labels, yticklabels=labels)
    ax.tick_params(axis='x', rotation=90, length=0, pad=3)
    ax.tick_params(axis='y', length=0, pad=3, labelleft=False, labelright=True)
    ax.set_xticks(np.arange(-.5, len(order), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(order), 1), minor=True)
    ax.grid(which='minor', color='white', lw=.3); ax.tick_params(which='minor', length=0)
    cax = fig.add_axes([x + left, y + h + .01, w - left - .07, .012])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label('Sequence identity (%)', labelpad=3)
    cb.ax.xaxis.set_ticks_position('top'); cb.ax.xaxis.set_label_position('top')
    cb.set_ticks([0, 25, 50, 75, 100]); cb.ax.tick_params(labelsize=8)
    ax._panel_extras = [dt, dl, cax]
    return ax


def identity():
    """Figure 1 (A-C) and Figure S1."""
    fig = plt.figure(figsize=(7.2, 7.6))
    m, z, o = mat_cluster('figure1A_locus_identity_matrix')
    ax = clustered(fig, m, z, o, [.045, .45, .765, .415]); letter(fig, ax, 'A')
    vals = pd.read_csv(D / 'figure1BC_plotted_values.csv')
    for j, (panel, locs) in enumerate([('1B', ['A', 'B', 'C']),
                                       ('1C', ['DRB1', 'DQB1', 'DRB3', 'DRB4', 'DRB5', 'DPB1'])]):
        ax = fig.add_axes([.11 + j * .49, .095, .38, .20])
        p = vals[vals.panel == panel]; sl = list(p.query_locus.unique()); width = .78 / len(locs)
        for k, l in enumerate(locs):
            q = p[p.subject_locus == l].set_index('query_locus').reindex(sl)
            ax.bar(np.arange(len(sl)) + (k - (len(locs) - 1) / 2) * width, q['Maximum Identity'],
                   yerr=q.error_bar_halfwidth, width=width, color=C[l],
                   error_kw={'lw': .5, 'capsize': 1}, label=l)
        ax.set(xticks=range(len(sl)), xticklabels=sl, ylim=(0, 100),
               ylabel='Mean best-match identity (%)' if j == 0 else None)
        ax.tick_params(axis='x', rotation=50 if j == 0 else 0)
        ax.legend(frameon=False, ncol=3, bbox_to_anchor=(.5, 1.015), loc='lower center',
                  columnspacing=.8, handlelength=1.2, handletextpad=.35)
        letter(fig, ax, 'BC'[j])
    save(fig, 'Figure1')

    # FigureS1: single panel, no panel letter. The delivery script built and
    # discarded this; the save is restored here.
    fig = plt.figure(figsize=(5.8, 4.9))
    m, z, o = mat_cluster('figureS1_accessory_identity_matrix')
    clustered(fig, m, z, o, [.035, .21, .77, .55])
    save(fig, 'FigureS1', tight=False)

    for name in ['figure1A_locus_identity_matrix', 'figureS1_accessory_identity_matrix']:
        m, z, o = mat_cluster(name)
        pd.DataFrame(z, columns=['left', 'right', 'distance', 'n_loci']).to_csv(
            AUDIT / (name + '_average_linkage_v7.csv'), index=False)
        pd.DataFrame({'order': range(len(o)), 'locus': m.index[o]}).to_csv(
            AUDIT / (name + '_display_order_v7.csv'), index=False)


def classifier():
    """Figure 2: classifier performance and stability of the SLA projection."""
    pred = pd.read_csv(D / 'canonical_ensemble_predictions.csv')
    val = pd.read_csv(D / 'canonical_validation_across_seeds.csv')
    loo = pd.read_csv(D / 'leave_one_seed_out_by_locus.csv')
    fig, axs = plt.subplots(3, 2, figsize=(7.2, 7.3))
    fig.subplots_adjust(left=.105, right=.985, bottom=.085, top=.955, hspace=.65, wspace=.38)
    ax = axs[0, 0]
    for i, c in enumerate(['class1', 'class2']):
        v = val[val.analysis == c].val_accuracy
        ax.errorbar(i, v.mean(), yerr=v.std(), fmt='s', color=CLASS[i], capsize=4, ms=6)
        ax.scatter(i + np.linspace(-.14, .14, 5), v, color=CLASS[i], s=18)
    ax.set(xticks=[0, 1], xticklabels=['Class I', 'Class II'], ylabel='Held-out accuracy',
           ylim=(.86, 1.005), xlim=(-.5, 1.5))
    ax = axs[0, 1]
    for i, c in enumerate(['class1', 'class2']):
        d = pred[pred.analysis == c]; counts = d.seeds_voting_ensemble_top_locus.value_counts()
        ax.bar(np.arange(1, 6) + (i - .5) * .34, [counts.get(k, 0) / len(d) * 100 for k in range(1, 6)],
               width=.32, color=CLASS[i], label=['Class I', 'Class II'][i])
    ax.set(xticks=range(1, 6), ylabel='SLA proteins (%)',
           xlabel='Models agreeing with\nensemble-selected locus', ylim=(0, 65))
    ax.legend(frameon=False, ncol=2, loc='upper left', columnspacing=.7, handlelength=1.1)
    for j, (c, ls) in enumerate([('class1', ['A', 'B', 'C']),
                                 ('class2', ['DRB1', 'DQB1', 'DRB3', 'DRB4', 'DRB5', 'DPB1'])]):
        ax = axs[1, j]; d = loo[(loo.analysis == c) & (loo.scope == 'ALL')]
        for k, l in enumerate(ls):
            full = pred[pred.analysis == c]['mass_' + l].mean()
            dd = d[d.ensemble_id.str.startswith('drop')]['mean_mass_' + l]
            ax.plot([k, k], [dd.min(), dd.max()], color=C[l], lw=1.5)
            ax.scatter(k + np.linspace(-.12, .12, len(dd)), dd, color=C[l], s=15, alpha=.5)
            ax.scatter(k, full, color=C[l], s=35, marker='D', edgecolor='white', linewidth=.5, zorder=4)
        ax.set(xticks=range(len(ls)), xticklabels=ls, ylim=(0, .82), ylabel='Mean locus probability')
        ax.tick_params(axis='x', rotation=40 if j else 0)
        handles = [Line2D([], [], marker='D', color='#444444', lw=0, label='Five models', ms=5),
                   Line2D([], [], marker='o', color='#999999', lw=0, label='Omit one model', ms=4)]
        ax.legend(handles=handles, frameon=False, loc='upper right', handletextpad=.3)
    rng = np.random.default_rng(0)
    for j, c in enumerate(['class1', 'class2']):
        ax = axs[2, j]; d = pred[pred.analysis == c]
        for k in range(1, 6):
            ys = d[d.seeds_voting_ensemble_top_locus == k].ensemble_top_locus_probability.to_numpy()
            if len(ys):
                ax.scatter(k + rng.uniform(-.18, .18, len(ys)), ys, s=10, color=CLASS[j],
                           alpha=.4, edgecolors='none')
                q = np.quantile(ys, [.25, .5, .75])
                ax.plot([k, k], q[[0, 2]], lw=3, color=CLASS[j])
                ax.plot([k - .09, k + .09], [q[1]] * 2, lw=1.6, color='black')
            ax.text(k, 1.025, str(len(ys)), ha='center', va='bottom', fontsize=8)
        ax.text(.52, 1.025, 'n', ha='left', va='bottom', fontsize=8)
        ax.set(xticks=range(1, 6), ylim=(.25, 1.07), xlim=(.45, 5.5),
               ylabel='Highest mean probability', xlabel='Models agreeing with\nensemble-selected locus')
    for ax, t in zip(axs.ravel(), 'ABCDEF'):
        letter(fig, ax, t)
    save(fig, 'Figure2')


def relationships():
    """Figure 3: ensemble probability distributions, bubbles and networks."""
    d = pd.read_csv(D / 'figure3_unthresholded_relationships.csv')
    fig = plt.figure(figsize=(7.2, 7.4))
    axes = [fig.add_axes(r) for r in [[.13, .67, .84, .215], [.13, .43, .84, .12],
                                      [.13, .285, .84, .065], [.13, .085, .84, .10]]]
    audit = []
    for ci, cls in enumerate(['class1', 'class2']):
        q = d[d.analysis == cls]
        sl = sorted(q.sla_locus.unique(),
                    key=lambda s: (int(s.split('-')[1]) if s.split('-')[1].isdigit()
                                   else {'DRB1': 0, 'DQB1': 1}[s.split('-')[1]]))
        hl = list(q.hla_subgroup.drop_duplicates())
        hm = q.set_index('hla_subgroup').hla_locus.to_dict()
        m = q.pivot(index='sla_locus', columns='hla_subgroup',
                    values='mean_ensemble_probability').reindex(index=sl, columns=hl)
        n = q.groupby('sla_locus').n_sla_sequences.first()
        ax = axes[2 * ci]
        # Areas are proportional to mean probability with the same scale in both classes.
        for i, l in enumerate(sl):
            ax.scatter(range(len(hl)), [i] * len(hl), s=m.loc[l].to_numpy() * 230,
                       color=[C[hm[h]] for h in hl], edgecolors='none', alpha=.85)
        ax.set(xlim=(-.6, len(hl) - .4), ylim=(len(sl) - .45, -.55), xticks=range(len(hl)),
               xticklabels=[h.replace('_Other', '\nOther') if cls == 'class2' else h for h in hl],
               yticks=range(len(sl)), yticklabels=[f'{s} ({n[s]})' for s in sl])
        ax.tick_params(axis='x', rotation=90, length=0, pad=3, labelsize=8)
        ax.tick_params(axis='y', length=0, labelsize=8)
        ax.grid(axis='y', color='#eeeeee', lw=.5); ax.set_axisbelow(True)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax = axes[2 * ci + 1]
        sx = {s: v for s, v in zip(sl, np.linspace(0, len(hl) - 1, len(sl)))}
        hx = {h: i for i, h in enumerate(hl)}
        for _, r in q.sort_values('mean_ensemble_probability').iterrows():
            p = float(r.mean_ensemble_probability)
            if p > 0:
                # No numerical cutoff: a minimum stroke ensures even very small weights render.
                ax.plot([sx[r.sla_locus], hx[r.hla_subgroup]], [1, 0], color=C[r.hla_locus],
                        lw=.16 + 2.5 * np.sqrt(p), alpha=.19 + .60 * np.sqrt(p),
                        solid_capstyle='round', zorder=1)
                audit.append({'analysis': cls, 'sla_locus': r.sla_locus,
                              'hla_subgroup': r.hla_subgroup, 'mean_probability': p, 'rendered': True})
        ax.scatter(list(sx.values()), [1] * len(sl), s=15, color='#555555', zorder=3)
        ax.scatter(range(len(hl)), [0] * len(hl), s=12, color=[C[hm[h]] for h in hl], zorder=3)
        for s, x in sx.items():
            ax.text(x, 1.09, s, ha='center', va='bottom', fontsize=7.8)
        for h, x in hx.items():
            ax.text(x, -.06, h.replace('_Other', '\nOther') if cls == 'class2' else h,
                    ha='center', va='top', rotation=90, fontsize=8)
        ax.set(xlim=(-.6, len(hl) - .4), ylim=(-.02, 1.1)); ax.axis('off')
    for ax, t in zip(axes, 'ABCD'):
        letter(fig, ax, t)
    h = [Line2D([], [], marker='o', color=C[k], lw=0, label=k, ms=4)
         for k in ['A', 'B', 'C', 'DRB1', 'DQB1', 'DRB3', 'DRB4', 'DRB5', 'DPB1']]
    fig.legend(handles=h, loc='upper center', bbox_to_anchor=(.54, .999), ncol=9,
               frameon=False, columnspacing=.75, handletextpad=.25)
    handles = [plt.scatter([], [], s=p * 230, color='#888888') for p in [.05, .2, .5]]
    fig.text(.28, .944, 'Mean probability:', fontsize=8, ha='right')
    fig.legend(handles=handles, labels=['0.05', '0.20', '0.50'], loc='center left',
               bbox_to_anchor=(.29, .951), ncol=3, frameon=False, columnspacing=1.1, handletextpad=.5)
    save(fig, 'Figure3')
    pd.DataFrame(audit).to_csv(AUDIT / 'Figure3_rendered_edges_v7.csv', index=False)
    pd.DataFrame([{'analysis': c, 'n_relationships': len(q),
                   'n_rendered_edges': int((q.mean_ensemble_probability > 0).sum()),
                   'edge_rule': 'all mean probabilities > 0',
                   'edge_weight': 'mean ensemble probability per SLA sequence',
                   'retained_probability_mass_fraction': 1.0}
                  for c, q in d.groupby('analysis')]).to_csv(
        AUDIT / 'Figure3_display_rule_v7.csv', index=False)


def eplets():
    """Figure 4: global HLA/SLA mismatch at HLA eplet-associated positions."""
    recomputed = TABLES_DIR / 'Figure4_plotted_values.csv'
    a = pd.read_csv(recomputed if recomputed.exists() else D / 'Figure4_plotted_values.csv')
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.45), gridspec_kw={'width_ratios': [1, 1.18]})
    fig.subplots_adjust(left=.085, right=.99, bottom=.29, top=.91, wspace=.61)
    ax = axs[0]
    for i, c in enumerate(['class1', 'class2']):
        r = a[(a.analysis == c) & (a.reading == 'agree_gap_skip')].iloc[0]
        for j, k in enumerate(['mismatch_eplet_pct', 'mismatch_other_pct']):
            x = i + (j - .5) * .34
            ax.bar(x, r[k], width=.32, color=[EP, OTHER][j])
            ax.text(x, r[k] + 1.3, f'{r[k]:.1f}%', ha='center', fontsize=8)
    ax.set(xticks=[0, 1], xticklabels=['Class I', 'Class II'], ylim=(0, 50), ylabel='Residue mismatch (%)')
    ax = axs[1]
    labels = []
    for i, c in enumerate(['class1', 'class2']):
        for j, reading in enumerate(['agree_gap_skip', 'agree_gap_mismatch']):
            r = a[(a.analysis == c) & (a.reading == reading)].iloc[0]
            y = 3 - i * 2 - j; col = CLASS[i]
            ax.errorbar(r.excess_mismatch_pp, y,
                        xerr=[[r.excess_mismatch_pp - r.excess_ci_low_pp],
                              [r.excess_ci_high_pp - r.excess_mismatch_pp]],
                        fmt='o' if j == 0 else 's', color=col, capsize=3, ms=5)
            labels.append(f'{"I" if i == 0 else "II"}: ' + ('covered' if j == 0 else 'gaps included'))
    ax.axvline(0, color='#999999', lw=.6, ls='--')
    ax.set(yticks=[3, 2, 1, 0], yticklabels=labels, ylim=(-.5, 3.5), xlim=(-12, 24),
           xlabel='Excess mismatch\n(percentage points)')
    ax.tick_params(axis='y', labelsize=8)
    for ax, t in zip(axs, 'AB'):
        letter(fig, ax, t)
    fig.legend(handles=[Patch(color=EP, label='HLA eplet-associated positions'),
                        Patch(color=OTHER, label='Other aligned positions')],
               loc='lower left', bbox_to_anchor=(.07, .004), ncol=2, frameon=False, columnspacing=1.6)
    save(fig, 'Figure4')


def _panel_manifest():
    """The corrected manifest, after checking it against the superseded one.

    The delivery script read ``structural_panel_manifest_original_measurements.csv``
    (whose ``final_figure`` column still says "Figure 4"/"Figure 5"). Only
    ``hla_allele`` and ``sla_allele`` are consumed. We use the corrected v6
    manifest instead, and assert those two fields agree so the substitution
    cannot change a panel title.
    """
    new = pd.read_csv(D / 'structural_panel_manifest_v6.csv')
    old = pd.read_csv(D / 'structural_panel_manifest_original_measurements.csv')
    keys = ['panel_name', 'hla_allele', 'sla_allele']
    a = new[keys].drop_duplicates().sort_values('panel_name').reset_index(drop=True)
    b = old[keys].drop_duplicates().sort_values('panel_name').reset_index(drop=True)
    if not a.equals(b):
        raise AssertionError(
            'structural manifests disagree on the allele fields used for panel titles:\n'
            f'{a.compare(b) if a.shape == b.shape else (a, b)}')
    return new


def structures():
    """Figures 5 and 6: representative class I and class II structural comparisons."""
    rep = json.loads((WORK / 'render_report.json').read_text())
    lab = pd.read_csv(D / 'structural_labels_v6.csv')
    man = _panel_manifest()
    audit = []
    for cls, num, ncols in [('class1', 5, 2), ('class2', 6, 1)]:
        fig, axes = plt.subplots(3, ncols, figsize=(7.2, 7.5 if ncols == 2 else 7.4), squeeze=False)
        fig.subplots_adjust(left=.035, right=.99, top=.955, bottom=.018, wspace=.06, hspace=.14)
        for ax, name in zip(axes.ravel(), [n for n in rep if n.startswith(cls)]):
            rr = rep[name]; m = man[man.panel_name == name].iloc[0]
            full = np.array(Image.open(WORK / (name + '.png')).convert('RGB')) / 255.
            fh, fw = full.shape[:2]
            yy, xx = np.where(full.min(axis=2) < .97)
            y0 = max(0, int(yy.min()) - 175); y1 = min(fh, int(yy.max()) + 176)
            x0 = max(0, int(xx.min()) - 175); x1 = min(fw, int(xx.max()) + 176)
            img = full[y0:y1, x0:x1]; h, w = img.shape[:2]
            ax.imshow(img); ax.set(xlim=(0, w), ylim=(h, 0)); ax.axis('off')
            rows = lab[(lab.panel == name) & (lab.system == 'HLA')].sort_values('published_position')
            anchors = [tuple(np.array(project(rr['view'], rr['points'][f'{name}:{int(r.published_position)}']))
                             - [x0, y0]) for _, r in rows.iterrows()]
            placed = place_labels(anchors, img.min(axis=2) < .97, (w, h),
                                  margin=75, gap=110, x_margin=185)
            for (_, r), (x, y), (lx, ly) in zip(rows.iterrows(), anchors, placed):
                special = r.selection_reason == 'contact example'
                col = '#986610' if special else '#777777'
                ax.plot([x, lx], [y, ly], color=col, lw=.5); ax.plot(x, y, 'o', ms=1.5, color=col)
                pieces = [(r.sla_aa, SLA), ('/', '#555555'), (r.hla_aa, HLA),
                          (str(int(r.published_position)), col if special else '#303030')]
                txt = HPacker(children=[TextArea(t, textprops={'color': c, 'fontsize': 8.5,
                                                              'fontweight': 'normal'})
                                        for t, c in pieces], align='center', pad=0, sep=0)
                ab = AnnotationBbox(txt, (lx, ly), frameon=True,
                                    bboxprops={'fc': 'white', 'ec': 'none',
                                               'boxstyle': 'square,pad=0.1'},
                                    box_alignment=(.5, .5))
                ax.add_artist(ab)
                audit.append({'figure': num, 'panel': name, 'label': r.label,
                              'first_residue_species': 'SLA', 'first_residue_color': SLA,
                              'second_residue_species': 'HLA', 'second_residue_color': HLA,
                              'canonical_position': int(r.published_position)})
            # Pair names and measured values are identifiers; panel descriptions are in captions.
            bb = ax.get_position(original=True)
            hname = fig.text(bb.x0 + .10 * bb.width, bb.y1 + .007, 'HLA-' + m.hla_allele,
                             fontsize=8.5, color=HLA, va='bottom')
            sname = fig.text(bb.x0 + .56 * bb.width, bb.y1 + .007, m.sla_allele,
                             fontsize=8.5, color=SLA, va='bottom')
            rmsd = fig.text(bb.x1 - .01, bb.y0 + .007, f'{rr["ce_rmsd"]:.1f} Å',
                            fontsize=8.5, ha='right')
            letter(fig, ax, name[-1], extras=(hname, sname, rmsd))
        save(fig, f'Figure{num}')
    pd.DataFrame(audit).to_csv(AUDIT / 'Figure5_Figure6_colored_label_audit_v7.csv', index=False)


def confusion(ax, m, fig, small=False):
    cmap = plt.get_cmap('Blues').with_extremes(bad='#D9D9D9')
    im = ax.imshow(m.values, vmin=0, vmax=1, cmap=cmap, interpolation='none', aspect='equal')
    ax.set(xticks=range(len(m.columns)), xticklabels=m.columns, yticks=range(len(m)),
           yticklabels=m.index, xlabel='Predicted label', ylabel='True label')
    ax.tick_params(axis='x', rotation=90, length=0); ax.tick_params(axis='y', length=0)
    if small:
        for i in range(len(m)):
            for j in range(len(m.columns)):
                v = m.iloc[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f'{v:.3f}', ha='center', va='center', fontsize=8,
                            color='white' if v > .55 else '#222222')
    cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.035)
    cb.set_label('Within-row proportion'); ax.set_anchor('C'); ax._panel_extras = [cb.ax]


def qc():
    """Figures S2 and S6 (training + locus confusion) and S3, S4, S5, S7 (subgroup confusion).

    ``small=False`` throughout: numeric cell labels are deliberately omitted from
    the large confusion matrices, which is the approved presentation.
    """
    seedcolors = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00']
    for cls, num in [('class1', 2), ('class2', 6)]:
        h = pd.read_csv(D / f'canonical_qc_{cls}_training_histories.csv')
        m = pd.read_csv(D / f'canonical_qc_{cls}_locus_confusion_rownormalised_pooled.csv', index_col=0)
        fig = plt.figure(figsize=(7.2, 3.0))
        gs = fig.add_gridspec(1, 3, left=.085, right=.935, bottom=.22, top=.79, wspace=.60,
                              width_ratios=[1, 1, 1.1])
        axes = [fig.add_subplot(gs[0, j]) for j in range(3)]
        for j, k in enumerate(['loss', 'accuracy']):
            ax = axes[j]
            for seed, color in zip(range(42, 47), seedcolors):
                g = h[h.seed == seed]
                ax.plot(g.epoch + 1, g['train_' + k], c=color, lw=1)
                ax.plot(g.epoch + 1, g['val_' + k], c=color, lw=1, ls='--')
            ax.set(xlabel='Epoch', ylabel='Cross-entropy loss' if j == 0 else 'Accuracy',
                   xticks=sorted(h.epoch.unique() + 1)[::(1 if cls == 'class1' else 2)])
            if j:
                ax.set_ylim(0, 1.03)
        confusion(axes[2], m, fig, small=False)
        for ax, t in zip(axes, 'ABC'):
            letter(fig, ax, t)
        handles = ([Line2D([], [], color=c, label=str(s), lw=1.6)
                    for s, c in zip(range(42, 47), seedcolors)]
                   + [Line2D([], [], color='#444444', label='Training', lw=1),
                      Line2D([], [], color='#444444', ls='--', label='Held-out', lw=1)])
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.53, .995), frameon=False,
                   ncol=7, handlelength=1.7, columnspacing=1)
        save(fig, f'FigureS{num}')

    # Single-panel subgroup confusion matrices. The delivery script built and
    # discarded these four; the saves are restored here.
    for num, stem in [(3, 'figureS3_HLA-A'), (4, 'figureS4_HLA-B'),
                      (5, 'figureS5_HLA-C'), (7, 'figureS7_classII')]:
        m = pd.read_csv(D / f'{stem}_confusion_canonical_rownormalised_pooled.csv', index_col=0)
        fig, ax = plt.subplots(figsize=(6.5, 5.8))
        fig.subplots_adjust(left=.16, right=.88, bottom=.18, top=.97)
        confusion(ax, m, fig)
        save(fig, f'FigureS{num}', tight=False)


def profiles():
    """Figures S8, S9 (agreement by comparator locus) and S10 (eplet-associated positions).

    All three mark eplet-associated positions, so all three need the mask that
    ``slahla_pub.eplet_mask`` supplies locally.
    """
    mask = eplet_mask.load()
    d = pd.read_csv(D / 'figureS8_S13_UPDATED_plotted_values.csv')
    d['is_eplet_associated'] = [p in mask[a] for a, p in zip(d.analysis, d.position)]
    for cls, num, ls in [('class1', 8, ['A', 'B', 'C']), ('class2', 9, ['DRB1', 'DQB1'])]:
        fig, axes = plt.subplots(len(ls), 1, figsize=(7.2, 2.05 * len(ls)), sharex=True)
        fig.subplots_adjust(left=.105, right=.975, bottom=.105, top=.87, hspace=.32)
        for ax, l, t in zip(np.atleast_1d(axes), ls, 'ABC'):
            q = d[(d.analysis == cls) & (d.hla_comparator_locus == l)].sort_values('position')
            x = q.position
            ax.plot(x, q.agreement_gaps_excluded.rolling(10, center=True).mean(), c=C[l], lw=1.2)
            ax.plot(x, q.agreement_gaps_as_mismatch.rolling(10, center=True).mean(), c=C[l], lw=1, ls='--')
            ax.plot(x, q.coverage_share_of_pairs, c='#888888', lw=.65, alpha=.7)
            e = q[q.is_eplet_associated]
            ax.scatter(e.position, e.agreement_gaps_excluded, s=9, c=EP, edgecolors='none', zorder=3)
            ax.set(ylim=(0, 1.035), ylabel='Residue agreement')
            letter(fig, ax, t)
        axes[-1].set_xlabel('Mature-protein position')
        handles = [Line2D([], [], c='#444444', label='Gaps excluded'),
                   Line2D([], [], c='#444444', ls='--', label='Gaps = mismatch'),
                   Line2D([], [], c='#888888', lw=.8, label='Coverage'),
                   Line2D([], [], c=EP, marker='o', lw=0, ms=3, label='Eplet-associated')]
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.53, .985), frameon=False,
                   ncol=2, columnspacing=2)
        save(fig, f'FigureS{num}')

    # FigureS10. The significance labels come from the archived permutation
    # statistics, not from literals. circular_shift_p tests profile_difference
    # (the equal-position contrast), not difference.
    stats = pd.read_csv(D / 'eplet_statistics_v6.csv')
    plabel = {}
    for cls in ('class1', 'class2'):
        r = stats[(stats.analysis == cls) & (stats.scope == 'ALL')
                  & (stats.reading == 'agree_gap_skip')].iloc[0]
        plabel[cls] = format_permutation_p(float(r.circular_shift_p),
                                           int(r.n_shifts_at_least_as_extreme),
                                           int(r.n_shifts))
    d = pd.read_csv(D / 'figureS_identity_eplet_agreement_corrected_p_plotted_values.csv')
    d['is_eplet_associated'] = [p in mask[a] for a, p in zip(d.analysis, d.position)]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 4.9), gridspec_kw={'width_ratios': [2, 1]})
    fig.subplots_adjust(left=.10, right=.98, bottom=.13, top=.94, wspace=.36, hspace=.40)
    for i, cls in enumerate(['class1', 'class2']):
        q = d[d.analysis == cls]; ax = axs[i, 0]
        ax.plot(q.position, q.smoothed, c=CLASS[i], lw=1)
        e = q[q.is_eplet_associated]
        ax.scatter(e.position, e.agreement, s=9, color=EP, edgecolors='none')
        ax.set(ylim=(0, 1.04), ylabel='Residue agreement',
               xlabel='Mature-protein position' if i == 1 else None)
        ax = axs[i, 1]
        groups = [q[q.is_eplet_associated].agreement.dropna(),
                  q[~q.is_eplet_associated].agreement.dropna()]
        bp = ax.boxplot(groups, patch_artist=True, showfliers=False, widths=.55,
                        medianprops={'color': '#222222', 'linewidth': 1})
        for b, c in zip(bp['boxes'], [EP, OTHER]):
            b.set_facecolor(c); b.set_alpha(.45)
        ax.set(xticks=[1, 2], xticklabels=['Eplet', 'Other'], ylim=(0, 1.04))
        ax.text(.5, .06, plabel[cls], transform=ax.transAxes, ha='center', fontsize=8)
    for ax, t in zip(axs.ravel(), 'ABCD'):
        letter(fig, ax, t)
    save(fig, 'FigureS10')


TASKS = {'identity': identity, 'classifier': classifier, 'relationships': relationships,
         'eplets': eplets, 'structures': structures, 'qc': qc, 'profiles': profiles}


def main(argv=None):
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    ensure_outputs()
    for k, f in TASKS.items():
        if not argv or k in argv:
            try:
                f()
            except eplet_mask.MissingExternalInput as exc:
                runlog.record(exc.what, exc.blocks, exc.how)
                print(f"{k}: SKIPPED - {exc.what}", flush=True)
                continue
            print(k, flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

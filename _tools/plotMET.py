"""
plotMET.py  --  Plot MET_pt and MET_phi for QCD MC and data (2024).

Reads file paths from:
  datasets_2024_data.txt   (one xrootd path per line, lines starting with # ignored)
  datasets_2024_QCD.txt    (same format)

Samples ~1% of events at random per file, no selection applied.

Usage:
  python plotMET.py [--data datasets_2024_data.txt] [--qcd datasets_2024_QCD.txt]
                   [--frac 0.01] [--seed 42] [--outdir plots/]
"""

import argparse
import os
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_paths(filename):
    paths = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                paths.append(line)
    return paths


def load_branches(paths, branches, frac, tree_name, seed):
    """
    Read `branches` from all ROOT files in `paths`, keeping ~frac of events
    chosen via Bernoulli sampling (independent per event).  Returns a dict
    {branch: np.ndarray}.
    """
    rng = np.random.default_rng(seed)
    out = {b: [] for b in branches}
    n_files_ok = 0

    for path in paths:
        try:
            with uproot.open(path) as f:
                tree = f[tree_name]
                for batch in tree.iterate(branches, step_size="50 MB", library="np"):
                    n = len(batch[branches[0]])
                    mask = rng.random(n) < frac
                    for b in branches:
                        out[b].append(batch[b][mask])
            n_files_ok += 1
        except Exception as e:
            print(f"  [skip] {path}  →  {e}")

    print(f"  Loaded {n_files_ok}/{len(paths)} files")
    return {b: np.concatenate(out[b]) if out[b] else np.array([]) for b in branches}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_variable(ax, data_arr, qcd_arr, var, xlabel, bins):
    # normalise to unit area
    kw = dict(bins=bins, histtype='step', density=True, linewidth=1.5)
    ax.hist(data_arr, label=f'Data  (N={len(data_arr):,})', color='black', **kw)
    ax.hist(qcd_arr,  label=f'QCD MC  (N={len(qcd_arr):,})', color='royalblue', **kw)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Normalised events / bin', fontsize=11)
    ax.set_title(var, fontsize=12)
    ax.legend(fontsize=10)
    ax.set_yscale('log')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data',   default='datasets_2024_data.txt')
    parser.add_argument('--qcd',    default='datasets_2024_QCD.txt')
    parser.add_argument('--frac',   type=float, default=0.01,
                        help='Fraction of events to keep per file (default 1%%)')
    parser.add_argument('--seed',   type=int,   default=42)
    parser.add_argument('--tree',   default='Events')
    parser.add_argument('--outdir', default='plots')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    branches = ['MET_pt', 'MET_phi']

    print(f"\nLoading DATA from {args.data} ...")
    data_paths = read_paths(args.data)
    print(f"  {len(data_paths)} paths found")
    data = load_branches(data_paths, branches, args.frac, args.tree, args.seed)

    print(f"\nLoading QCD MC from {args.qcd} ...")
    qcd_paths = read_paths(args.qcd)
    print(f"  {len(qcd_paths)} paths found")
    qcd = load_branches(qcd_paths, branches, args.frac, args.tree, args.seed + 1)

    print(f"\nData events kept : {len(data['MET_pt']):,}")
    print(f"QCD  events kept : {len(qcd['MET_pt']):,}")

    # MET_pt: log-spaced bins from 0 to 500 GeV
    bins_pt  = np.linspace(0, 500, 101)
    bins_phi = np.linspace(-np.pi, np.pi, 63)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('MET distributions — 2024 (inclusive, ~1% statistics)', fontsize=13)

    plot_variable(axes[0], data['MET_pt'],  qcd['MET_pt'],
                  'MET_pt',  r'MET $p_T$ [GeV]', bins_pt)
    plot_variable(axes[1], data['MET_phi'], qcd['MET_phi'],
                  'MET_phi', r'MET $\phi$ [rad]', bins_phi)

    plt.tight_layout()
    outpath = os.path.join(args.outdir, 'MET_distributions.png')
    fig.savefig(outpath, dpi=150)
    print(f"\nSaved → {outpath}")
    plt.close()


if __name__ == '__main__':
    main()

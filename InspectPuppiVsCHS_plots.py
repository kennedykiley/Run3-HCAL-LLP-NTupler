#!/usr/bin/env python3
# Quick FWLite check of CHS vs PUPPI jets + simple plots
# Usage:
#   python inspectPuppiVsCHS_remote_plots.py files.txt --maxEv 200
# this is off of the raw EDM file!

import sys, os, argparse
from math import sqrt
from DataFormats.FWLite import Events, Handle
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Inspect CHS vs PUPPI jets remotely and plot.")
parser.add_argument("filelist", help="Text file listing EDM ROOT files (one per line)")
parser.add_argument("--maxEv", type=int, default=50, help="Max events to process")
args = parser.parse_args()

# --- read file list ---
with open(args.filelist) as f:
    files = [l.strip() for l in f if l.strip() and not l.startswith("#")]
if not files:
    sys.exit("File list empty or invalid.")

print(f"Loaded {len(files)} files; processing up to {args.maxEv} events.")
events = Events(files)

# --- handles ---
h_chs = Handle("std::vector<reco::PFJet>")
h_puppi = Handle("std::vector<reco::PFJet>")
label_chs = ("ak4PFJetsCHS", "", "RECO")
label_puppi = ("ak4PFJetsPuppi", "", "RECO")


def dr2(j1, j2):
    dphi = abs(j1.phi() - j2.phi())
    if dphi > 3.1416: dphi = 2*3.1416 - dphi
    deta = j1.eta() - j2.eta()
    return deta*deta + dphi*dphi

# --- storage for plots ---
nCHS, nPUPPI, ptRatio, ptCHS_all, ptPUPPI_all = [], [], [], [], []

for i, event in enumerate(events):
    if i >= args.maxEv: break
    event.getByLabel(label_chs, h_chs)
    event.getByLabel(label_puppi, h_puppi)
    chs = [j for j in h_chs.product() if j.pt()>20 and abs(j.eta())<2]
    puppi = [j for j in h_puppi.product() if j.pt()>20 and abs(j.eta())<2]
    nCHS.append(len(chs)); nPUPPI.append(len(puppi))

    for cj in chs:
        best, dmin = None, 9
        for pj in puppi:
            d = dr2(cj,pj)
            if d<dmin: dmin, best = d, pj
        if best and dmin<0.3**2:
            ptRatio.append(best.pt()/cj.pt())
    for j in chs: ptCHS_all.append(j.pt())
    for j in puppi: ptPUPPI_all.append(j.pt())

print(f"Processed {i+1} events. Making plots:")

# --- Plot 1: jet multiplicity per event ---
plt.figure()
bins = range(0, max(max(nCHS or [0]), max(nPUPPI or [0]))+2)
plt.hist(nCHS, bins=bins, histtype='step', label='CHS jets', linewidth=2)
plt.hist(nPUPPI, bins=bins, histtype='step', label='PUPPI jets', linewidth=2)
plt.xlabel('N jets (pT>20 GeV, |eta|<2)')
plt.ylabel('Events')
plt.legend()
plt.title('Jet multiplicity per event')

# --- Plot 2: pT ratio for matched jets ---
plt.figure()
plt.hist(ptRatio, bins=40, range=(0,2), histtype='stepfilled', alpha=0.7)
plt.xlabel('pT(PUPPI) / pT(CHS)')
plt.ylabel('Matched jets')
plt.title('pT ratio (matched jets)')

# --- Plot 3: compute statistics ---
ptRatio = np.array(ptRatio)
ptRatio_clean = ptRatio[(ptRatio > 0.01) & (ptRatio < 5)]  # remove any pathological values

mean_ratio = np.mean(ptRatio_clean)
std_ratio = np.std(ptRatio_clean)
median_ratio = np.median(ptRatio_clean)
n_matched = len(ptRatio_clean)

print(f"\n--- Matched jet pT ratio statistics ---")
print(f"  Number of matched jets: {n_matched}")
print(f"  Mean (PUPPI/CHS): {mean_ratio:.3f}")
print(f"  Median: {median_ratio:.3f}")
print(f"  RMS spread: {std_ratio:.3f}")
print(f"  Suggested systematic: ±{abs(1 - mean_ratio)*100:.1f}% if offset is non-negligible\n")

# --- histogram plot ---
plt.figure(figsize=(7,5))
plt.hist(ptRatio_clean, bins=40, range=(0,2), histtype='stepfilled', alpha=0.7, color='steelblue')
plt.xlabel('pT(PUPPI) / pT(CHS)')
plt.ylabel('Matched jets')
plt.title('pT ratio (matched jets)')

# reference lines
plt.axvline(1.0, color='gray', linestyle='--', linewidth=1)
plt.axvline(mean_ratio, color='red', linestyle='--', linewidth=1.5, label=f'Mean = {mean_ratio:.3f}')
plt.axvline(mean_ratio + std_ratio, color='orange', linestyle=':', linewidth=1)
plt.axvline(mean_ratio - std_ratio, color='orange', linestyle=':', linewidth=1)

plt.legend()
plt.grid(True, linestyle=':')
plt.tight_layout()

# --- Plot 4: pT spectra ---
plt.figure()
plt.hist(ptCHS_all, bins=60, range=(0,200), histtype='step', label='CHS', linewidth=2)
plt.hist(ptPUPPI_all, bins=60, range=(0,200), histtype='step', label='PUPPI', linewidth=2)
plt.xlabel('Jet pT [GeV]')
plt.ylabel('Jets')
plt.legend()
plt.title('Jet pT spectra')

# --- Plot 5: PUPPI/CHS pT spectrum ratio (for jets > 40 GeV) ---
pt_min, pt_max, nbins = 20, 200, 15
bins = np.linspace(pt_min, pt_max, nbins+1)

# histograms (normalized counts per bin)
chs_counts, _ = np.histogram(ptCHS_all, bins=bins)
puppi_counts, _ = np.histogram(ptPUPPI_all, bins=bins)

# compute ratio safely
ratio = np.divide(puppi_counts, chs_counts, out=np.zeros_like(puppi_counts, dtype=float), where=chs_counts>0)
bin_centers = 0.5 * (bins[1:] + bins[:-1])

plt.figure()
plt.plot(bin_centers, ratio, drawstyle='steps-mid', linewidth=2)
plt.axhline(1.0, color='gray', linestyle='--', linewidth=1)
plt.xlabel('Jet pT [GeV]')
plt.ylabel('PUPPI / CHS')
plt.title('Jet yield ratio (jets with pT > 20 GeV)')
plt.ylim(0, 2)
plt.grid(True, linestyle=':')

plt.tight_layout()
for i, num in enumerate(plt.get_fignums()):
    fig = plt.figure(num)
    fig.savefig(f"PuppiVsCHS_plot_{i+1}.png")
print("Saved all plots as PuppiVsCHS_plot_*.png")
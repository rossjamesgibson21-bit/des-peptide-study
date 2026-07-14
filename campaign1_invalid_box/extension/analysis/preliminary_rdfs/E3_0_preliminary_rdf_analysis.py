#!/usr/bin/env python3
"""
E3.0 — Preliminary RDF Analysis (Existing 10 ns Trajectories)
==============================================================

Computes atom-pair radial distribution functions between DES component
atoms and motif backbone oxygen atoms for all 9 systems. Identifies
first RDF minima for pair-specific coordination cutoffs.

Objective: Determine whether the uniform 0.35 nm cutoff used in Phase 3
is appropriate for all atom pairs, or whether pair-specific cutoffs
(particularly for bulkier DES components) would change coordination numbers.

Author: Ross Gibson
Date: June 2026
"""

import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from itertools import product as iterproduct
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ============================================================
# Configuration
# ============================================================

BASE_DIR = os.path.expanduser("~/des-peptide-study")
os.chdir(BASE_DIR)

# System definitions: (prmtop, pdb, dcd)
SYSTEMS = {
    "GGE_water":     ("systems/GGE_water/GGE_water.prmtop",
                      "systems/GGE_water/GGE_water_equilibrated.pdb",
                      "systems/GGE_water/GGE_water_prod.dcd"),
    "GGE_reline":    ("systems/GGE_reline/GGE_reline.prmtop",
                      "systems/GGE_reline/GGE_reline_equilibrated.pdb",
                      "systems/GGE_reline/GGE_reline_prod.dcd"),
    "GGE_glyceline": ("systems/GGE_glyceline/GGE_glyceline.prmtop",
                      "systems/GGE_glyceline/GGE_glyceline_equilibrated.pdb",
                      "systems/GGE_glyceline/GGE_glyceline_prod.dcd"),
    "CME_water":     ("systems/CME_water/CME_water.prmtop",
                      "systems/CME_water/CME_water_equilibrated.pdb",
                      "systems/CME_water/CME_water_prod.dcd"),
    "CME_reline":    ("systems/CME_reline/CME_reline.prmtop",
                      "systems/CME_reline/CME_reline_equilibrated.pdb",
                      "systems/CME_reline/CME_reline_prod.dcd"),
    "CME_glyceline": ("systems/CME_glyceline/CME_glyceline.prmtop",
                      "systems/CME_glyceline/CME_glyceline_equilibrated.pdb",
                      "systems/CME_glyceline/CME_glyceline_prod.dcd"),
    "YIY_water":     ("systems/YIY_water/YIY_water.prmtop",
                      "systems/YIY_water/YIY_water_equilibrated.pdb",
                      "systems/YIY_water/YIY_water_prod.dcd"),
    "YIY_reline":    ("systems/YIY_reline/YIY_reline.prmtop",
                      "systems/YIY_reline/YIY_reline_equilibrated.pdb",
                      "systems/YIY_reline/YIY_reline_prod.dcd"),
    "YIY_glyceline": ("systems/YIY_glyceline/YIY_glyceline.prmtop",
                      "systems/YIY_glyceline/YIY_glyceline_equilibrated.pdb",
                      "systems/YIY_glyceline/YIY_glyceline_prod.dcd"),
}

# Motif residue indices in topology (0-indexed, including ACE cap as residue 0)
# GGE construct: ACE-E-E-E-G-G-E-I-V-F-NME → motif GGE at residues 4,5,6
# CME construct: ACE-L-Y-Q-C-M-E-F-V-R-NME → motif CME at residues 4,5,6
# YIY construct: ACE-N-P-Y-I-Y-K-NME       → motif YIY at residues 3,4,5
MOTIF_RESIDUES = {
    "GGE": [4, 5, 6],
    "CME": [4, 5, 6],
    "YIY": [3, 4, 5],
}

# DES component atom selections
# Format: (residue_name, atom_names, label)
PARTNER_SELECTIONS = {
    "water": [("HOH", ["O"], "Water O")],
    "reline": [
        ("HOH", ["O"],       "Water O"),
        ("CHO", ["N"],       "Choline N"),
        ("CLA", ["Cl"],      "Chloride"),
        ("URE", ["N", "N1"], "Urea N"),
        ("URE", ["O"],       "Urea O"),
    ],
    "glyceline": [
        ("HOH", ["O"],            "Water O"),
        ("CHO", ["N"],            "Choline N"),
        ("CLA", ["Cl"],           "Chloride"),
        ("GOL", ["O", "O1", "O2"], "Glycerol O"),
    ],
}

# RDF parameters
RDF_BIN_WIDTH = 0.005   # nm — fine bins for accurate minimum detection
RDF_RANGE = (0.15, 1.0) # nm — range of interest
EQ_CUTOFF_NS = 1.0      # ns — equilibration cutoff (consistent with Phase 3)
SAVE_INTERVAL_PS = 1.0   # ps per frame

# Phase 3 default cutoff for comparison
DEFAULT_CUTOFF = 0.35    # nm

# ============================================================
# Helper functions
# ============================================================

def get_motif_backbone_O(topology, motif_residue_indices):
    """Get atom indices for backbone carbonyl O atoms in motif residues."""
    indices = []
    for res_idx in motif_residue_indices:
        res = topology.residue(res_idx)
        for atom in res.atoms:
            if atom.name == "O" and atom.element.symbol == "O":
                indices.append(atom.index)
                break  # one backbone O per residue
    return np.array(indices)


def get_partner_atoms(topology, residue_name, atom_names):
    """Get atom indices for specified atoms in all residues of a given name."""
    indices = []
    for res in topology.residues:
        if res.name == residue_name:
            for atom in res.atoms:
                if atom.name in atom_names:
                    indices.append(atom.index)
    return np.array(indices)


def make_pairs(group_a, group_b):
    """Generate all (i, j) pairs between two atom index arrays."""
    pairs = np.array(list(iterproduct(group_a, group_b)))
    return pairs


def find_first_minimum(r, g_r, search_start=0.2, search_end=0.6):
    """
    Find the first minimum in g(r) after the first peak.
    
    Strategy: find the first peak (max) after search_start, then find
    the first minimum (local min) after that peak.
    """
    # Restrict to search range
    mask = (r >= search_start) & (r <= search_end)
    r_sub = r[mask]
    g_sub = g_r[mask]
    
    if len(g_sub) < 5:
        return np.nan, np.nan
    
    # Smooth slightly to avoid noise-induced false minima
    kernel_size = max(3, int(0.01 / RDF_BIN_WIDTH))  # ~0.01 nm smoothing
    if kernel_size % 2 == 0:
        kernel_size += 1
    kernel = np.ones(kernel_size) / kernel_size
    g_smooth = np.convolve(g_sub, kernel, mode='same')
    
    # Find first peak
    peak_idx = None
    for i in range(1, len(g_smooth) - 1):
        if g_smooth[i] > g_smooth[i-1] and g_smooth[i] > g_smooth[i+1] and g_smooth[i] > 1.0:
            peak_idx = i
            break
    
    if peak_idx is None:
        return np.nan, np.nan
    
    # Find first minimum after peak
    for i in range(peak_idx + 1, len(g_smooth) - 1):
        if g_smooth[i] < g_smooth[i-1] and g_smooth[i] < g_smooth[i+1]:
            return r_sub[i], g_smooth[i]
    
    return np.nan, np.nan


def get_solvent_type(system_name):
    """Extract solvent type from system name."""
    if "reline" in system_name:
        return "reline"
    elif "glyceline" in system_name:
        return "glyceline"
    else:
        return "water"


def get_peptide_name(system_name):
    """Extract peptide name from system name."""
    return system_name.split("_")[0]


# ============================================================
# Main analysis
# ============================================================

print("=" * 70)
print("E3.0 — Preliminary RDF Analysis")
print("Computing DES component ↔ motif backbone O radial distribution functions")
print("=" * 70)

# Storage for results
all_results = []   # list of dicts for summary table
all_rdfs = {}      # (system, partner_label) → (r, g_r)

# Output directory
out_dir = os.path.join(BASE_DIR, "extension", "analysis", "preliminary_rdfs")
os.makedirs(out_dir, exist_ok=True)

for sys_name, (prmtop_path, pdb_path, dcd_path) in SYSTEMS.items():
    peptide = get_peptide_name(sys_name)
    solvent = get_solvent_type(sys_name)
    
    print(f"\n{'─'*70}")
    print(f"Loading {sys_name}...")
    
    # Load trajectory with equilibration cutoff
    traj = md.load(dcd_path, top=prmtop_path)
    eq_frames = int(EQ_CUTOFF_NS * 1000 / SAVE_INTERVAL_PS)  # 1000 frames for 1 ns
    traj = traj[eq_frames:]
    print(f"  Frames after {EQ_CUTOFF_NS} ns cutoff: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")
    
    # Get motif backbone O atoms
    motif_res = MOTIF_RESIDUES[peptide]
    backbone_O = get_motif_backbone_O(traj.topology, motif_res)
    print(f"  Motif backbone O atoms: {len(backbone_O)} (residues {motif_res})")
    
    # Get partner selections for this solvent type
    partners = PARTNER_SELECTIONS[solvent]
    
    for res_name, atom_names, label in partners:
        partner_idx = get_partner_atoms(traj.topology, res_name, atom_names)
        
        if len(partner_idx) == 0:
            print(f"  {label}: no atoms found — skipping")
            continue
        
        pairs = make_pairs(backbone_O, partner_idx)
        print(f"  {label}: {len(partner_idx)} atoms, {len(pairs)} pairs ... ", end="", flush=True)
        
        # Compute RDF
        n_bins = int((RDF_RANGE[1] - RDF_RANGE[0]) / RDF_BIN_WIDTH)
        r, g_r = md.compute_rdf(
            traj,
            pairs=pairs,
            r_range=RDF_RANGE,
            n_bins=n_bins,
            periodic=True,
        )
        
        # Find first minimum
        r_min, g_min = find_first_minimum(r, g_r)
        
        delta_from_default = (r_min - DEFAULT_CUTOFF) * 10 if not np.isnan(r_min) else np.nan  # in Angstrom
        
        print(f"first min = {r_min:.3f} nm (g = {g_min:.2f})" if not np.isnan(r_min) 
              else "no clear minimum found")
        
        # Store results
        all_rdfs[(sys_name, label)] = (r, g_r)
        all_results.append({
            "system": sys_name,
            "peptide": peptide,
            "solvent": solvent,
            "partner": label,
            "n_partner_atoms": len(partner_idx),
            "n_pairs": len(pairs),
            "first_min_nm": r_min,
            "first_min_g": g_min,
            "default_cutoff_nm": DEFAULT_CUTOFF,
            "delta_from_default_A": delta_from_default,
        })

# ============================================================
# Summary table
# ============================================================

print(f"\n{'='*70}")
print("SUMMARY: First RDF minima vs default 0.35 nm cutoff")
print(f"{'='*70}")
print(f"{'System':<20} {'Partner':<15} {'1st min (nm)':>12} {'g(r_min)':>10} {'Δ from 0.35':>12}")
print(f"{'─'*70}")

for res in all_results:
    r_min = res["first_min_nm"]
    g_min = res["first_min_g"]
    delta = res["delta_from_default_A"]
    
    r_str = f"{r_min:.3f}" if not np.isnan(r_min) else "N/A"
    g_str = f"{g_min:.2f}" if not np.isnan(g_min) else "N/A"
    d_str = f"{delta:+.2f} Å" if not np.isnan(delta) else "N/A"
    
    print(f"{res['system']:<20} {res['partner']:<15} {r_str:>12} {g_str:>10} {d_str:>12}")

# Flag pairs where cutoff differs substantially
print(f"\n{'─'*70}")
print("Pairs where first minimum differs from 0.35 nm by > 0.5 Å:")
flagged = [r for r in all_results if not np.isnan(r["delta_from_default_A"]) 
           and abs(r["delta_from_default_A"]) > 0.5]
if flagged:
    for r in flagged:
        print(f"  {r['system']:20s} {r['partner']:15s} → {r['first_min_nm']:.3f} nm "
              f"(Δ = {r['delta_from_default_A']:+.2f} Å)")
else:
    print("  None — 0.35 nm appears adequate for all pairs at 10 ns resolution")

# ============================================================
# Plotting — one figure per peptide, panels for each solvent
# ============================================================

fig, axes = plt.subplots(3, 3, figsize=(18, 14))
fig.suptitle("E3.0 — Preliminary RDFs: DES Components → Motif Backbone O\n"
             "(10 ns trajectories, 1 ns equilibration cutoff)", fontsize=14, fontweight="bold")

peptide_names = ["GGE", "CME", "YIY"]
solvent_names = ["water", "reline", "glyceline"]
solvent_colors = {
    "Water O":     "#1f77b4",
    "Choline N":   "#ff7f0e",
    "Chloride":    "#2ca02c",
    "Urea N":      "#d62728",
    "Urea O":      "#9467bd",
    "Glycerol O":  "#8c564b",
}

for row, peptide in enumerate(peptide_names):
    for col, solvent in enumerate(solvent_names):
        ax = axes[row, col]
        sys_name = f"{peptide}_{solvent}"
        
        partners = PARTNER_SELECTIONS[solvent]
        has_data = False
        
        for _, _, label in partners:
            key = (sys_name, label)
            if key in all_rdfs:
                r, g_r = all_rdfs[key]
                color = solvent_colors.get(label, "gray")
                ax.plot(r, g_r, color=color, linewidth=1.2, label=label, alpha=0.9)
                
                # Mark first minimum
                result = [res for res in all_results 
                         if res["system"] == sys_name and res["partner"] == label]
                if result and not np.isnan(result[0]["first_min_nm"]):
                    r_min = result[0]["first_min_nm"]
                    g_min = result[0]["first_min_g"]
                    ax.axvline(r_min, color=color, linestyle="--", alpha=0.4, linewidth=0.8)
                    ax.plot(r_min, g_min, "v", color=color, markersize=6, alpha=0.7)
                
                has_data = True
        
        # Reference line at default cutoff
        ax.axvline(DEFAULT_CUTOFF, color="black", linestyle=":", alpha=0.5, linewidth=1.0)
        ax.text(DEFAULT_CUTOFF + 0.005, ax.get_ylim()[1] * 0.9 if has_data else 2.0,
                "0.35 nm", fontsize=7, alpha=0.5)
        
        ax.axhline(1.0, color="gray", linestyle="-", alpha=0.2, linewidth=0.5)
        ax.set_xlim(RDF_RANGE)
        ax.set_title(f"{peptide} — {solvent}", fontsize=11, fontweight="bold")
        
        if row == 2:
            ax.set_xlabel("r (nm)")
        if col == 0:
            ax.set_ylabel("g(r)")
        
        if has_data:
            ax.legend(fontsize=7, loc="upper right")

plt.tight_layout(rect=[0, 0, 1, 0.95])
fig_path = os.path.join(out_dir, "preliminary_rdfs_all_systems.png")
plt.savefig(fig_path, dpi=200, bbox_inches="tight")
print(f"\nFigure saved: {fig_path}")
plt.close()

# ============================================================
# Save results to CSV
# ============================================================

import csv
csv_path = os.path.join(out_dir, "preliminary_rdf_cutoffs.csv")
fieldnames = ["system", "peptide", "solvent", "partner", "n_partner_atoms", "n_pairs",
              "first_min_nm", "first_min_g", "default_cutoff_nm", "delta_from_default_A"]
with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(all_results)
print(f"Results saved: {csv_path}")

print(f"\n{'='*70}")
print("E3.0 COMPLETE")
print(f"{'='*70}")

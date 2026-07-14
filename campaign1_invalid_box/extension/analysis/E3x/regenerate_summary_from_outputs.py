#!/usr/bin/env python3
"""Regenerate {SYSTEM}_summary.md (+ .csv) from the persisted E3.x section outputs,
WITHOUT re-running the trajectory analysis. All summary text is produced by the shared
e3x_summary module - the same code the notebook's E3.4 cell uses - so a regenerated
summary cannot differ from what a re-run would produce.

Explicit inputs (under DES_PROJECT_DIR, default ~/des-peptide-study), all on disk from
the committed runs:
  analysis/E3x/{SYSTEM}/{SYSTEM}_summary.csv                 tier, character, SASA/Rg verdicts
  analysis/E3x/{SYSTEM}/coordination/backbone_coordination.csv
  analysis/E3x/{SYSTEM}/coordination/sidechain_coordination.csv
  analysis/E3x/{SYSTEM}/hbonds/hbond_lifetimes.csv           (must be 1 ps)
  analysis/E3x/{SYSTEM}/hbonds/per_start_hbond_summary.csv   per-start counts + start labels
  analysis/E3x/{SYSTEM}/{SYSTEM}_basin_sasa.csv              tier D only (if present)
  analysis/E3x/{SYSTEM}/{SYSTEM}_trap_sasa.csv               tier C only (if present)
  analysis/convergence/convergence_pooled_{SYSTEM}.csv       carried SASA/Rg numbers

Output: {SYSTEM}_summary.md and {SYSTEM}_summary.csv, rewritten in place.
Metadata (tier, character, verdicts) is read from the committed {SYSTEM}_summary.csv and
passed straight back, so the header/verdict lines round-trip unchanged; only interpretive
prose in the body can differ from what is committed.

  python regenerate_summary_from_outputs.py GGE_reline CME_reline ...
  python regenerate_summary_from_outputs.py ALL
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e3x_summary

SAVE_INTERVAL_PS = 1.0
PROJECT_DIR = os.environ.get('DES_PROJECT_DIR', os.path.expanduser('~/des-peptide-study'))
EXT_DIR = os.path.join(PROJECT_DIR, 'extension')

SYSTEMS = ['GGE_reline', 'GGE_glyceline', 'GGE_water',
           'CME_water', 'CME_reline', 'CME_glyceline',
           'YIY_water', 'YIY_reline', 'YIY_glyceline']

COMPONENT_PARTNERS = {'reline': ['choline_N', 'chloride', 'urea_N', 'urea_O'],
                      'glyceline': ['choline_N', 'chloride', 'glycerol_O'], 'water': []}


def regenerate(SYSTEM):
    PEP, SOLVENT = SYSTEM.split('_', 1)
    OUT_DIR = os.path.join(EXT_DIR, 'analysis', 'E3x', SYSTEM)
    CONV_DIR = os.path.join(EXT_DIR, 'analysis', 'convergence')
    PARTNERS = ['water_O'] + COMPONENT_PARTNERS[SOLVENT]

    # metadata from the committed summary.csv (tier/character/verdicts round-trip unchanged)
    meta = pd.read_csv(os.path.join(OUT_DIR, f'{SYSTEM}_summary.csv')).iloc[0].to_dict()
    TIER = str(meta['tier'])
    CARRIED_VERDICT = dict(SASA=str(meta['sasa_verdict']), Rg=str(meta['rg_verdict']),
                           character=str(meta['character']))

    # H-bond fields (guard: must be the 1 ps result)
    hb = pd.read_csv(os.path.join(OUT_DIR, 'hbonds', 'hbond_lifetimes.csv')).iloc[0].to_dict()
    assert float(hb.get('resolution_ps', -1)) == SAVE_INTERVAL_PS, (
        f"{SYSTEM}: hbond_lifetimes.csv resolution_ps={hb.get('resolution_ps')} != {SAVE_INTERVAL_PS} ps "
        "- H-bond outputs are not the 1 ps result; re-run E3.3 before regenerating.")

    bb = pd.read_csv(os.path.join(OUT_DIR, 'coordination', 'backbone_coordination.csv'))
    sc = pd.read_csv(os.path.join(OUT_DIR, 'coordination', 'sidechain_coordination.csv'))
    for d in (bb, sc):
        d['structured_shell'] = d['structured_shell'].astype(str).str.lower().isin(['true', '1'])

    conv_pool = pd.read_csv(os.path.join(CONV_DIR, f'convergence_pooled_{SYSTEM}.csv'))
    def carried(obs):
        r = conv_pool[conv_pool['observable'].str.upper() == obs.upper()]
        return (float(r.iloc[0]['pooled_mean']), float(r.iloc[0]['ci95'])) if len(r) else (np.nan, np.nan)
    sasa_m, sasa_ci = carried('SASA')
    rg_m, rg_ci = carried('Rg')

    ps_path = os.path.join(OUT_DIR, 'hbonds', 'per_start_hbond_summary.csv')
    if os.path.exists(ps_path):
        ps_df = pd.read_csv(ps_path)
        per_start_hb = ps_df.to_dict('records')
        START_LABELS = list(ps_df['start'])
    else:
        per_start_hb = None
        START_LABELS = str(meta['starts']).split(';')

    def _load(name):
        p = os.path.join(OUT_DIR, f'{SYSTEM}_{name}.csv')
        return pd.read_csv(p).iloc[0].to_dict() if os.path.exists(p) else None

    def _normalise_basin(d):
        """Accept either basin schema on disk and return canonical K-state fields.
        Old binary schema (GGE_water) uses boundary/high_*/low_*; current schema uses
        n_states/states. Numbers are preserved; only the state labelling is unified."""
        if d is None or 'n_states' in d:
            return d
        hs, ls = d.get('high_starts', ''), d.get('low_starts', '')
        states = (f"HIGH:{float(d['high_mean']):.4f}+/-{float(d['high_ci95']):.4f} ({hs}) | "
                  f"LOW:{float(d['low_mean']):.4f}+/-{float(d['low_ci95']):.4f} ({ls})")
        return dict(n_states=2, states=states, per_start=d['per_start'],
                    pooled_mean=d['pooled_mean'], pooled_ci95=d['pooled_ci95'],
                    sidecar_sasa=d['sidecar_sasa'], sidecar_delta=d['sidecar_delta'])

    basin_sasa = _normalise_basin(_load('basin_sasa')) if TIER == 'D' else None
    trap_sasa = _load('trap_sasa') if TIER == 'C' else None

    text = e3x_summary.build_and_write_summary(
        OUT_DIR=OUT_DIR, SYSTEM=SYSTEM, TIER=TIER, PEP=PEP, SOLVENT=SOLVENT,
        START_LABELS=START_LABELS, CARRIED_VERDICT=CARRIED_VERDICT,
        sasa_m=sasa_m, sasa_ci=sasa_ci, rg_m=rg_m, rg_ci=rg_ci, hb=hb,
        bb_coord=bb, sc_coord=sc, per_start_hb=per_start_hb, PARTNERS=PARTNERS,
        basin_sasa=basin_sasa, trap_sasa=trap_sasa, write=True)
    print(f'[{SYSTEM}] tier {TIER} - wrote {SYSTEM}_summary.md (+ .csv)')
    return text


if __name__ == '__main__':
    args = sys.argv[1:] or ['GGE_reline']
    systems = SYSTEMS if args == ['ALL'] else args
    for s in systems:
        if s not in SYSTEMS:
            print(f'skip unknown system: {s}'); continue
        regenerate(s)

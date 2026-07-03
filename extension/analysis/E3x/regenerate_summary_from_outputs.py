#!/usr/bin/env python3
"""Regenerate {SYSTEM}_summary.md (and .csv) from the persisted E3.x section outputs,
WITHOUT re-running the trajectory analysis (~2 h). This is the decoupled-summary path:
E3.4 consumes section outputs on disk, so the readable summary is regenerable in seconds.

Reads:
  analysis/E3x/{SYSTEM}/hbonds/hbond_lifetimes.csv         -> H-bond fields (must be 1 ps)
  analysis/E3x/{SYSTEM}/hbonds/per_start_hbond_summary.csv -> start labels
  analysis/E3x/{SYSTEM}/coordination/backbone_coordination.csv
  analysis/E3x/{SYSTEM}/coordination/sidechain_coordination.csv
  analysis/convergence/convergence_pooled_{SYSTEM}.csv     -> carried SASA/Rg
  analysis/convergence/convergence_loo_{SYSTEM}.csv         -> SASA/Rg verdict

Writes:
  analysis/E3x/{SYSTEM}/{SYSTEM}_summary.md    (the readable record)
  analysis/E3x/{SYSTEM}/{SYSTEM}_summary.csv    (machine table; idempotent)

Only the three per-system knobs below change across the nine systems. Paths honour the
DES_PROJECT_DIR environment variable, defaulting to ~/des-peptide-study.
"""
import os
import numpy as np
import pandas as pd

# ── per-system knobs (only these three change across systems) ─────────────────
SYSTEM    = 'GGE_reline'
TIER      = 'A'
CHARACTER = 'clean / ergodic'     # E2.x qualitative convergence label (not on disk)

# ── paths ─────────────────────────────────────────────────────────────────────
PEP, SOLVENT = SYSTEM.split('_', 1)
PROJECT_DIR  = os.environ.get('DES_PROJECT_DIR', os.path.expanduser('~/des-peptide-study'))
EXT_DIR      = os.path.join(PROJECT_DIR, 'extension')
OUT_DIR      = os.path.join(EXT_DIR, 'analysis', 'E3x', SYSTEM)
CONV_DIR     = os.path.join(EXT_DIR, 'analysis', 'convergence')
SAVE_INTERVAL_PS = 1.0

# ── read section outputs ──────────────────────────────────────────────────────
hb = pd.read_csv(os.path.join(OUT_DIR, 'hbonds', 'hbond_lifetimes.csv')).iloc[0].to_dict()
bb_coord = pd.read_csv(os.path.join(OUT_DIR, 'coordination', 'backbone_coordination.csv'))
sc_coord = pd.read_csv(os.path.join(OUT_DIR, 'coordination', 'sidechain_coordination.csv'))
for d in (bb_coord, sc_coord):
    d['structured_shell'] = d['structured_shell'].astype(str).str.lower().isin(['true', '1'])

conv = {'pooled': pd.read_csv(os.path.join(CONV_DIR, f'convergence_pooled_{SYSTEM}.csv'))}

ps_path = os.path.join(OUT_DIR, 'hbonds', 'per_start_hbond_summary.csv')
if os.path.exists(ps_path):
    START_LABELS = list(pd.read_csv(ps_path)['start'])
else:
    START_LABELS = [f'start{i + 1}' for i in range(int(conv['pooled']['n_starts'].iloc[0]))]

# SASA/Rg verdict derived from the leave-one-out sidecar (robust unless a start dominates)
loo_path = os.path.join(CONV_DIR, f'convergence_loo_{SYSTEM}.csv')
def _verdict(obs):
    if not os.path.exists(loo_path):
        return 'robust'
    loo = pd.read_csv(loo_path)
    rows = loo[loo['observable'].str.upper() == obs.upper()]
    if rows.empty:
        return 'robust'
    dom = rows['shift_gt_ci'].astype(str).str.lower().isin(['true', '1']).any()
    return 'dominated' if bool(dom) else 'robust'
CARRIED_VERDICT = dict(SASA=_verdict('SASA'), Rg=_verdict('Rg'), character=CHARACTER)

def carried(obs):
    row = conv['pooled'][conv['pooled']['observable'].str.upper() == obs.upper()]
    return (float(row.iloc[0]['pooled_mean']), float(row.iloc[0]['ci95'])) if len(row) else (np.nan, np.nan)
sasa_m, sasa_ci = carried('SASA')
rg_m, rg_ci     = carried('Rg')

# guard: the H-bond outputs must be the 1 ps result, matching the E3.4 provenance check
assert float(hb.get('resolution_ps', -1)) == SAVE_INTERVAL_PS, (
    f"hbond_lifetimes.csv resolution_ps={hb.get('resolution_ps')} != {SAVE_INTERVAL_PS} ps; "
    "the H-bond outputs are not the full-resolution result — re-run E3.3 before regenerating.")

# ── machine summary (idempotent) ──────────────────────────────────────────────
summary = dict(
    system=SYSTEM, tier=TIER, n_starts=len(START_LABELS), starts=';'.join(map(str, START_LABELS)),
    character=CARRIED_VERDICT['character'],
    sasa_pooled_nm2=sasa_m, sasa_ci95=sasa_ci, sasa_verdict=CARRIED_VERDICT['SASA'],
    rg_pooled_nm=rg_m, rg_ci95=rg_ci, rg_verdict=CARRIED_VERDICT['Rg'],
    hbond_events=hb['n_events'], hbond_median_ps=hb['median_ps'],
    hbond_resolution_ps=hb['resolution_ps'], hbond_break_tol_frames=hb['break_tol_frames'],
    hbond_n_backbone_bonds=hb['n_distinct_bonds'],
    n_structured_backbone_pairs=int(bb_coord['structured_shell'].sum()) if len(bb_coord) else 0,
    n_structured_sidechain_pairs=int(sc_coord['structured_shell'].sum()) if len(sc_coord) else 0,
)
pd.DataFrame([summary]).to_csv(os.path.join(OUT_DIR, f'{SYSTEM}_summary.csv'), index=False)

# ── readable summary, organised by the three E3 questions (identical to E3.4) ──
IS_SIGNIFICANCE_TESTCASE = (SYSTEM == 'GGE_glyceline')
TESTS_INDIRECT_DIRECTLY  = (SOLVENT == 'glyceline')
IS_ANOMALY_SYSTEM        = (PEP == 'YIY')
BOTH_ROBUST = (CARRIED_VERDICT['SASA'] == 'robust' and CARRIED_VERDICT['Rg'] == 'robust')

def _tbl(df, cols):
    return df[cols].to_string(index=False) if len(df) else '(none)'

closed, deferred, L = [], [], []
L += [f'# E3.x system summary - {SYSTEM}', '',
      f'Tier {TIER} | {len(START_LABELS)}-start ({", ".join(map(str, START_LABELS))}) | {CARRIED_VERDICT["character"]}',
      f'Peptide {PEP} | solvent {SOLVENT}', '']

L += ['## E3.1 - Replicate-aware SASA & backbone coordination',
      f'- SASA (carried E2.2 pooled): {sasa_m} +/- {sasa_ci} nm2  [{CARRIED_VERDICT["SASA"]}]',
      f'- Rg   (carried E2.2 pooled): {rg_m} +/- {rg_ci} nm  [{CARRIED_VERDICT["Rg"]}]',
      f'- Convergence character: {CARRIED_VERDICT["character"]}',
      '- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):',
      '```', _tbl(bb_coord, ['partner', 'cutoff_nm', 'pooled_coord', 'ci95']), '```']
if IS_SIGNIFICANCE_TESTCASE:
    L += ['- CLOSES: named E3.1 significance test case (Phase-3 corrected ES 1.18); the verdict '
          'against corrected CIs is the load-bearing E3.1 result for the campaign.']
    closed.append('E3.1 significance (named test case)')
else:
    L += [f'- Replicate significance for THIS system: '
          f'{"achieved (both observables robust)" if BOTH_ROBUST else "see convergence character"}. '
          'The named E3.1 test case (GGE_glyceline, corrected ES 1.18) is a separate system; deferred there.']
    if BOTH_ROBUST:
        closed.append('E3.1 convergence for this system (robust)')
    deferred.append('E3.1 headline significance (GGE_glyceline ES 1.18)')
L += ['']

L += ['## E3.2 - Side-chain coordination (indirect-restructuring test)',
      '- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):',
      '```', _tbl(sc_coord, ['reference', 'partner', 'cutoff_nm', 'pooled_coord', 'ci95']), '```']
if TESTS_INDIRECT_DIRECTLY:
    L += ['- CLOSES (in part): glyceline system - glycerol-sidechain coordination is the direct test '
          'of the indirect hypothesis (substantial glycerol contact challenges purely-indirect; '
          'absence supports it), read against the water baseline.']
    closed.append('E3.2 glycerol-sidechain contact (this solvent)')
else:
    _base = ('urea direct-contact baseline (urea coordinates backbone AND carboxylate here)'
             if SOLVENT == 'reline' else 'aqueous baseline that rules out a peptide-intrinsic explanation')
    L += [f'- CONTRIBUTES (baseline): {SOLVENT} carries no glycerol, so it cannot test the '
          f'glycerol-sidechain question directly. It supplies the {_base}. The decisive glycerol '
          'comparison is deferred to the glyceline systems + water baseline.']
    deferred.append('E3.2 indirect-mechanism adjudication (glyceline + water)')
L += ['- Observation: motif vs flanking carboxylate coordinate near-uniformly per O (uniform acidic '
      'behaviour rather than local context) - a hypothesis the other systems test.', '']

L += ['## E3.3 - Backbone hydrogen bonds (anomaly test)',
      f'- {hb["n_distinct_bonds"]} backbone i->i-3 turns; {float(hb["resolution_ps"]):.0f} ps resolution, '
      f'{hb["n_events"]} continuous episodes, median {hb["median_ps"]} ps (break_tol={hb["break_tol_frames"]}).']
if IS_ANOMALY_SYSTEM:
    L += ['- CLOSES (in part): YIY is the anomaly system (Phase-3: 1433 reline events vs 0 in '
          'water/glyceline). Per-start consistency of the anomaly is the load-bearing E3.3 result.']
    closed.append('E3.3 YIY H-bond anomaly (per-start consistency)')
else:
    L += ['- CONTRIBUTES: the H-bond anomaly is YIY reline; this system validated the H-bond machinery '
          'and established the 1 ps resolution requirement (stride=10 aliases ~87% of episodes). '
          'Anomaly adjudication deferred to the YIY systems.']
    deferred.append('E3.3 anomaly adjudication (YIY systems)')
L += ['- Note: the continuous-episode median is largely blind to slow exchange; the persistence '
      'definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.', '']

L += ['## Carried caveats',
      '- Chloride: no structured shell (backbone or side-chain) - reportable null.']
_chol = bb_coord[(bb_coord['reference'] == 'backbone_O') & (bb_coord['partner'] == 'choline_N')] \
        if len(bb_coord) else pd.DataFrame()
if len(_chol) and bool(_chol['structured_shell'].any()):
    L += ['- Choline-N -> backbone-O registers a weak feature absent in the E3.0 preliminary and '
          'carrying the highest cutoff-sensitivity - treated as a shoulder, not an established shell.']
L += ['- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.', '']

L += ['## Question ledger', 'CLOSED by this system:']
L += [f'  - {c}' for c in closed] if closed else ['  - (none - validation/contribution system)']
L += ['CONTRIBUTES toward (deferred to cross-solvent stage):']
L += [f'  - {d}' for d in deferred] if deferred else ['  - (none)']

summary_text = '\n'.join(L)
with open(os.path.join(OUT_DIR, f'{SYSTEM}_summary.md'), 'w') as f:
    f.write(summary_text + '\n')
print(summary_text)
print(f'\n[wrote {SYSTEM}_summary.md and {SYSTEM}_summary.csv to {OUT_DIR}]')

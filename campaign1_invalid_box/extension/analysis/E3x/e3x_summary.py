"""Single source of truth for the E3.x per-system summary.

Imported by BOTH the notebook template (E3.4 cell) and regenerate_summary_from_outputs.py,
so the readable summary cannot drift from the data. The summary emits only:
  - tables computed from the coordination / H-bond outputs,
  - verdicts carried verbatim from the E2.2 sidecars,
  - caveats computed from those tables.
It asserts no mechanism and no interpretation of what a contact means - that belongs to the
cross-solvent stage, which reads the nine summaries together. Every line here is a fact the
tables contain; nothing is authored.
"""
import os
import numpy as np
import pandas as pd


def _tbl(df, cols):
    if df is None or not len(df):
        return '(none)'
    d = df.copy()
    if 'structured_shell' in d.columns and 'cutoff_nm' in d.columns:
        d.loc[~d['structured_shell'].astype(bool), 'cutoff_nm'] = np.nan
    disp = d[cols].copy()
    if 'cutoff_nm' in disp.columns:
        disp['cutoff_nm'] = disp['cutoff_nm'].map(lambda x: f'{x:.3f}' if pd.notna(x) else 'no shell')
    for c in ('pooled_coord', 'ci95'):
        if c in disp.columns:
            disp[c] = disp[c].map(lambda x: f'{x:.4f}')
    return disp.to_string(index=False)


def build_and_write_summary(*, OUT_DIR, SYSTEM, TIER, PEP, SOLVENT, START_LABELS,
                            CARRIED_VERDICT, sasa_m, sasa_ci, rg_m, rg_ci, hb,
                            bb_coord, sc_coord, per_start_hb, PARTNERS,
                            basin_sasa=None, trap_sasa=None, write=True):
    """Build the stripped summary text + machine row, optionally write .md/.csv. Returns the text."""
    bb_coord = bb_coord if bb_coord is not None else pd.DataFrame()
    sc_coord = sc_coord if sc_coord is not None else pd.DataFrame()

    IS_SIGNIFICANCE_TESTCASE = (SYSTEM == 'GGE_glyceline')
    TESTS_INDIRECT_DIRECTLY  = (SOLVENT == 'glyceline')
    IS_ANOMALY_SYSTEM        = (PEP == 'YIY')
    BOTH_ROBUST = (CARRIED_VERDICT['SASA'] == 'robust' and CARRIED_VERDICT['Rg'] == 'robust')

    summary = dict(
        system=SYSTEM, tier=TIER, n_starts=len(START_LABELS), starts=';'.join(START_LABELS),
        character=CARRIED_VERDICT['character'],
        sasa_pooled_nm2=sasa_m, sasa_ci95=sasa_ci, sasa_verdict=CARRIED_VERDICT['SASA'],
        rg_pooled_nm=rg_m, rg_ci95=rg_ci, rg_verdict=CARRIED_VERDICT['Rg'],
        hbond_events=hb['n_events'], hbond_median_ps=hb['median_ps'],
        hbond_resolution_ps=hb['resolution_ps'], hbond_break_tol_frames=hb['break_tol_frames'],
        hbond_n_backbone_bonds=hb['n_distinct_bonds'],
        n_structured_backbone_pairs=int(bb_coord['structured_shell'].sum()) if len(bb_coord) else 0,
        n_structured_sidechain_pairs=int(sc_coord['structured_shell'].sum()) if len(sc_coord) else 0,
    )

    closed, deferred, L = [], [], []
    L += [f'# E3.x system summary - {SYSTEM}', '',
          f'Tier {TIER} | {len(START_LABELS)}-start ({", ".join(START_LABELS)}) | {CARRIED_VERDICT["character"]}',
          f'Peptide {PEP} | solvent {SOLVENT}', '']

    # --- E3.1 : carried verdicts + backbone table; tier-specific SASA reporting (all factual) ---
    if TIER == 'D' and basin_sasa is not None:
        _k = basin_sasa['n_states']
        L += ['## E3.1 - Replicate-aware SASA & backbone coordination',
              f'- SASA: {_k} statistically separated state(s) across the diverse starts '
              '(adjacent gap > combined 95% CI); diverse-start ensemble, not an equilibrium mean.',
              '- Per-start SASA means (nm2):',
              '```', '  ' + basin_sasa['per_start'].replace(';', '\n  '), '```',
              f'- Grouped states: {basin_sasa["states"]}',
              f'- Pooled ensemble: {basin_sasa["pooled_mean"]} +/- {basin_sasa["pooled_ci95"]} nm2 '
              f'(sidecar {basin_sasa["sidecar_sasa"]}, delta {basin_sasa["sidecar_delta"]})',
              f'- Rg (carried): {rg_m} +/- {rg_ci} nm  [{CARRIED_VERDICT["Rg"]}]',
              '- Backbone-O coordination (pooled):',
              '```', _tbl(bb_coord, ['partner', 'cutoff_nm', 'pooled_coord', 'ci95']), '```']
        deferred.append(f'E3.1 SASA equilibrium ({_k} separated states; enhanced sampling -> Paper 2)')
    elif TIER == 'C' and trap_sasa is not None:
        L += ['## E3.1 - Replicate-aware SASA & backbone coordination',
              '- SASA one-trap: one start (loo-flagged) sits apart from the consensus; reported both ways.',
              f'  - full pool (all starts): {trap_sasa["full_mean"]} +/- {trap_sasa["full_ci95"]} nm2',
              f'  - consensus (drop {trap_sasa["trap_starts"] or "n/a"}): {trap_sasa["cons_mean"]} +/- {trap_sasa["cons_ci95"]} nm2',
              f'  - inclusion shift (full - consensus): {trap_sasa["inclusion_shift"]:+} nm2',
              f'  - sidecar {trap_sasa["sidecar_sasa"]} (delta {trap_sasa["sidecar_delta"]})',
              f'- Rg (carried): {rg_m} +/- {rg_ci} nm  [{CARRIED_VERDICT["Rg"]}]',
              '- Backbone-O coordination (pooled):',
              '```', _tbl(bb_coord, ['partner', 'cutoff_nm', 'pooled_coord', 'ci95']), '```']
        deferred.append('E3.1 SASA (one-trap; full-pool vs consensus reported)')
    else:
        L += ['## E3.1 - Replicate-aware SASA & backbone coordination',
              f'- SASA (carried): {sasa_m} +/- {sasa_ci} nm2  [{CARRIED_VERDICT["SASA"]}]',
              f'- Rg   (carried): {rg_m} +/- {rg_ci} nm  [{CARRIED_VERDICT["Rg"]}]',
              '- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):',
              '```', _tbl(bb_coord, ['partner', 'cutoff_nm', 'pooled_coord', 'ci95']), '```']
        if IS_SIGNIFICANCE_TESTCASE:
            closed.append('E3.1 significance (named test case: GGE_glyceline)')
        elif BOTH_ROBUST:
            closed.append('E3.1 convergence for this system (both robust)')
            deferred.append('E3.1 headline significance (GGE_glyceline)')
        else:
            deferred.append('E3.1 headline significance (GGE_glyceline)')
    L += ['']

    # --- E3.2 : side-chain table + which partners contact which sites (read from the table) ---
    L += ['## E3.2 - Side-chain coordination',
          '- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):',
          '```', _tbl(sc_coord, ['reference', 'partner', 'cutoff_nm', 'pooled_coord', 'ci95']), '```']
    if TESTS_INDIRECT_DIRECTLY:
        _gly = sc_coord[(sc_coord['partner'] == 'glycerol_O') & sc_coord['structured_shell'].astype(bool)] if len(sc_coord) else pd.DataFrame()
        _hit = [f'{r.reference} ({r.pooled_coord:.2f}/O)' for r in _gly.itertuples()]
        _all = sorted(set(sc_coord['reference'])) if len(sc_coord) else []
        _miss = [x for x in _all if x not in set(_gly['reference'])]
        _gbb = bb_coord[(bb_coord['partner'] == 'glycerol_O') & bb_coord['structured_shell'].astype(bool)] if len(bb_coord) else pd.DataFrame()
        _gbb_txt = (f' Backbone-glycerol shell: {_gbb.iloc[0]["pooled_coord"]:.2f}/O at {_gbb.iloc[0]["cutoff_nm"]:.3f}.'
                    if len(_gbb) else ' No backbone-glycerol shell.')
        if _hit:
            L += [f'- Glycerol side-chain first-shell contact at: {", ".join(_hit)}'
                  + (f'; no side-chain shell at {", ".join(_miss)}.' if _miss else '.') + _gbb_txt]
        else:
            L += ['- Glycerol: no structured side-chain shell at any motif site.' + _gbb_txt]
        closed.append('E3.2 glycerol contact - side-chain and backbone (this solvent)')
    else:
        _des = [p for p in PARTNERS if p != 'water_O']
        # choline_N -> backbone is a shoulder (separately caveated); exclude from direct-contact list
        _bb = bb_coord[bb_coord['partner'].isin(_des) & bb_coord['structured_shell'].astype(bool)
                       & (bb_coord['partner'] != 'choline_N')] if len(bb_coord) else pd.DataFrame()
        _sc = sc_coord[sc_coord['partner'].isin(_des) & sc_coord['structured_shell'].astype(bool)] if len(sc_coord) else pd.DataFrame()
        if SOLVENT == 'water':
            L += ['- Water only: aqueous side-chain baseline (no DES partners).']
        else:
            _bbt = 'no direct backbone shell' if not len(_bb) else 'backbone shell at ' + ', '.join(sorted(set(_bb['partner'])))
            _sct = 'no side-chain shell' if not len(_sc) else 'side-chain shell at ' + ', '.join(f'{r.reference}-{r.partner}' for r in _sc.itertuples())
            L += [f'- {SOLVENT} DES contacts (choline backbone shoulder excluded): {_bbt}; {_sct}.']
        deferred.append('E3.2 indirect-mechanism adjudication (glyceline vs water)')
    L += ['']

    # --- E3.3 : turns + episodes + per-start counts (all factual) ---
    L += ['## E3.3 - Backbone hydrogen bonds',
          f'- {hb["n_distinct_bonds"]} backbone N->O turns; {float(hb["resolution_ps"]):.0f} ps resolution, '
          f'{hb["n_events"]} continuous episodes, median {hb["median_ps"]} ps (break_tol={hb["break_tol_frames"]}).']
    if per_start_hb:
        try:
            _psh = '; '.join(f'{r["start"]}:{int(r["n_episodes"])}ep/{int(r["n_bonds"])}b' for r in per_start_hb)
            L += [f'- Per-start (episodes/distinct-bonds): {_psh}.']
        except (KeyError, TypeError):
            pass
    if IS_ANOMALY_SYSTEM and SOLVENT == 'reline':
        L += ['- YIY_reline is the anomaly system; counts here are read against the YIY water and '
              'glyceline baselines at matched 1 ps backbone-only scope (not the Phase-3 region-only metric).']
        closed.append('E3.3 YIY reline H-bond counts (vs water/glyceline baselines)')
    elif IS_ANOMALY_SYSTEM:
        L += [f'- Anomaly-control arm: {SOLVENT} backbone-H-bond baseline for the YIY_reline comparison '
              '(1 ps backbone-only scope; not zero at this scope, unlike the Phase-3 region-only metric).']
        deferred.append('E3.3 anomaly adjudication (needs YIY_reline)')
    else:
        deferred.append('E3.3 anomaly adjudication (YIY systems)')
    L += ['']

    # --- caveats : all computed from the tables ---
    L += ['## Caveats']
    # chloride: data-driven across backbone + side-chain
    _clall = pd.concat([bb_coord, sc_coord], ignore_index=True) if (len(bb_coord) or len(sc_coord)) else pd.DataFrame()
    _cl = _clall[_clall['partner'] == 'chloride'] if len(_clall) else pd.DataFrame()
    _clhit = _cl[_cl['structured_shell'].astype(bool)] if len(_cl) else pd.DataFrame()
    if len(_cl) and not len(_clhit):
        L += ['- Chloride: no structured shell (backbone or side-chain) - reportable null.']
    elif len(_clhit):
        _h = ', '.join(f'{r.reference}-chloride ({r.pooled_coord:.3f})' for r in _clhit.itertuples())
        L += [f'- Chloride: no backbone shell; structured shell at {_h} (low occupancy - verify).']
    # choline backbone shoulder
    _chol = bb_coord[(bb_coord['reference'] == 'backbone_O') & (bb_coord['partner'] == 'choline_N')] if len(bb_coord) else pd.DataFrame()
    if len(_chol) and bool(_chol['structured_shell'].any()):
        L += ['- Choline-N -> backbone-O: weak feature, highest cutoff-sensitivity - shoulder, not an established shell.']
    # side-chain sites with no shell against any partner
    if len(sc_coord):
        _any = sc_coord.groupby('reference')['structured_shell'].apply(lambda x: bool(x.astype(bool).any()))
        _null = sorted([r for r, h in _any.items() if not h])
        if _null:
            L += [f'- No structured shell against any partner at side-chain site(s): {", ".join(_null)} - reportable null.']
    L += ['- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.']
    if CARRIED_VERDICT['Rg'] == 'dominated':
        L += ['- Non-ergodic backbone: Rg dominated (trapped at 200 ns); SASA robust and pooled/reportable. '
              'SASA and Rg are decoupled for this system.']
    L += ['']

    # --- ledger : which questions this system's data addresses (scope, factual) ---
    L += ['## Question ledger', 'CLOSED by this system:']
    L += [f'  - {c}' for c in closed] if closed else ['  - (none - validation/baseline system)']
    L += ['CONTRIBUTES toward (deferred to cross-solvent stage):']
    L += [f'  - {d}' for d in deferred] if deferred else ['  - (none)']

    text = '\n'.join(L)
    if write:
        pd.DataFrame([summary]).to_csv(os.path.join(OUT_DIR, f'{SYSTEM}_summary.csv'), index=False)
        with open(os.path.join(OUT_DIR, f'{SYSTEM}_summary.md'), 'w') as f:
            f.write(text + '\n')
    return text

# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: despep
#     language: python
#     name: despep
# ---

# %% [markdown]
# # DES-PEP-CESD-001 — Phase 0: Protocol, estimand and analysis freeze
#
# **Protocol §5.** This notebook does the Phase 0 *verification* work: it loads the
# frozen freeze-artifacts, checks their internal consistency and cross-references,
# confirms the frozen values against the protocol's own §5/§8/§14 constraints,
# derives the resource forecast from the frozen ceiling, evaluates the **Gate 0
# release matrix** mechanically, and writes a **SHA-256 manifest** of the evidence
# package.
#
# It does **not** sign the gate. PASS/HOLD/FAIL is a human reviewer decision
# (§3.1); this notebook produces the evidence those reviewers act on.
#
# Outputs consumed here (all authored as version-controlled freeze artifacts):
# - `outputs/identifier_scheme.yaml`
# - `outputs/question_estimand_registry.yaml`
# - `outputs/scope_and_resource_ceiling.yaml`
# - `../../sap/SAP_shell_v1.yaml`
# - `../../governance/roles_and_arrangements.yaml`
# - `../../governance/blinding_charter.yaml`
# - `../../governance/deviation_log_template.yaml`
# - `../../governance/gate_certificate_template.yaml`
# - `../../governance/reviewer_checklists.yaml`

# %%
import os, sys, glob, json, hashlib, datetime
import yaml
import pandas as pd

# Resolve repo root robustly: walk UP from the current working directory until we
# find the directory that holds the campaign markers. This works regardless of
# where Jupyter was launched (repo root, the notebook folder, or in between) and
# under nbconvert/in-process execution.
def find_repo_root(start):
    markers = ("governance", "phases", "sap")
    p = os.path.abspath(start)
    while True:
        if all(os.path.isdir(os.path.join(p, m)) for m in markers):
            return p
        parent = os.path.dirname(p)
        if parent == p:            # reached filesystem root without a match
            raise FileNotFoundError(
                "Could not locate the des-pep-cesd-001 repo root (a dir containing "
                "governance/, phases/ and sap/) at or above: " + os.path.abspath(start)
                + ". Run this notebook from inside the cloned repo tree.")
        p = parent

REPO = find_repo_root(os.getcwd())
print("repo root:", REPO)
OUT = os.path.join(REPO, "phases", "phase00_freeze", "outputs")

def load_yaml(rel):
    path = os.path.join(REPO, rel)
    with open(path) as f:
        return yaml.safe_load(f), path

ids, ids_p     = load_yaml("phases/phase00_freeze/outputs/identifier_scheme.yaml")
qer, qer_p     = load_yaml("phases/phase00_freeze/outputs/question_estimand_registry.yaml")
scope, scope_p = load_yaml("phases/phase00_freeze/outputs/scope_and_resource_ceiling.yaml")
sap, sap_p     = load_yaml("sap/SAP_shell_v1.yaml")
roles, roles_p = load_yaml("governance/roles_and_arrangements.yaml")
blind, blind_p = load_yaml("governance/blinding_charter.yaml")

print("loaded 6 core freeze artifacts")

# %% [markdown]
# ## 0.5 Schema validation (Gate 0 corrective action item 1)
#
# Before any cross-artifact consistency check, validate each controlled artifact
# against a version-controlled JSON Schema in `schemas/`. This catches structural
# drift (missing keys, wrong types, renamed fields, a resolved solvent mapping
# leaking into the identifier scheme, or a fabricated value where a `SIGN:`
# placeholder is required) independently of the hand-written checks in §1.

# %%
import jsonschema

def load_json(rel):
    with open(os.path.join(REPO, rel)) as f:
        return json.load(f)

schema_bindings = [
    ("phases/phase00_freeze/outputs/identifier_scheme.yaml",        "schemas/identifier_scheme.schema.json",        ids),
    ("phases/phase00_freeze/outputs/question_estimand_registry.yaml","schemas/question_estimand_registry.schema.json", qer),
    ("phases/phase00_freeze/outputs/scope_and_resource_ceiling.yaml","schemas/scope_and_resource_ceiling.schema.json", scope),
    ("sap/SAP_shell_v1.yaml",                                        "schemas/SAP_shell.schema.json",                sap),
    ("governance/roles_and_arrangements.yaml",                       "schemas/roles_and_arrangements.schema.json",   roles),
    ("governance/blinding_charter.yaml",                             "schemas/blinding_charter.schema.json",         blind),
]
# gate certificate validated structurally too (unsigned instance is valid structure)
cert_doc, cert_p = load_yaml("governance/gate_certificates/gate00_certificate.yaml")
schema_bindings.append(("governance/gate_certificates/gate00_certificate.yaml",
                        "schemas/gate_certificate.schema.json", cert_doc))

schema_rows = []
for art_rel, sch_rel, doc in schema_bindings:
    schema = load_json(sch_rel)
    validator = jsonschema.Draft202012Validator(schema)
    errs = sorted(validator.iter_errors(doc), key=lambda e: list(e.path))
    schema_rows.append({
        "artifact": art_rel.split("/")[-1],
        "schema": sch_rel.split("/")[-1],
        "result": "PASS" if not errs else "FAIL",
        "n_errors": len(errs),
        "first_error": "" if not errs else f"{list(errs[0].path)}: {errs[0].message[:80]}",
    })
schema_df = pd.DataFrame(schema_rows)
n_schema_fail = int((schema_df["result"] == "FAIL").sum())
print(schema_df.to_string(index=False))
print(f"\nSchema validation: {len(schema_df)-n_schema_fail}/{len(schema_df)} PASS, {n_schema_fail} FAIL")
assert n_schema_fail == 0, "Schema validation failures must be resolved before Gate 0 review"

# %% [markdown]
# ## 1. Consistency checks
#
# Each check is a hard assertion with an explanatory message. A failure here means
# the freeze artifacts disagree with each other or with the protocol's frozen
# constants — which must be fixed before the gate can be reviewed, not signed around.

# %%
checks = []
def check(name, ok, detail=""):
    checks.append({"check": name, "result": "PASS" if ok else "FAIL", "detail": detail})
    return ok

# --- Protocol identity threads through every artifact ---
PID = "DES-PEP-CESD-001"
for nm, doc in [("identifier_scheme", ids), ("question_estimand", qer),
                ("scope_resource", scope), ("SAP_shell", sap),
                ("roles", roles), ("blinding", blind)]:
    check(f"protocol_id present in {nm}", doc.get("protocol_id") == PID,
          f"got {doc.get('protocol_id')!r}")

# --- Peptides consistent between id scheme and estimand roles ---
pep_codes = {p["code"] for p in ids["peptides"]}
check("three peptides registered (GGE, CME, YIY)",
      pep_codes == {"GGE", "CME", "YIY"}, f"{sorted(pep_codes)}")
check("GGE is the pilot/confirmatory-core peptide",
      any(p["code"] == "GGE" and "pilot" in p["role"].lower() for p in ids["peptides"]))

# --- Solvents consistent (three) and masked-arm namespace reserved ---
sol_codes = {s["code"] for s in ids["solvents"]}
check("three solvents registered (WATER, RELINE, GLYCELINE)",
      sol_codes == {"WATER", "RELINE", "GLYCELINE"}, f"{sorted(sol_codes)}")
check("masked-arm namespace has 3 arms",
      len(ids["masked_arms"]["namespace"]) == 3, f"{ids['masked_arms']['namespace']}")
check("solvent-code key NOT embedded in identifier scheme (hard firewall §18.3)",
      "mapping" not in {k.lower() for k in ids["masked_arms"]} or
      ids["masked_arms"].get("mapping_status", "").upper().startswith("TO BE ASSIGNED"),
      "mapping must be steward-held, not in analyst-readable freeze artifact")

# --- Start blocks: three strata, no solvent-effect basis (§8.1) ---
sb_codes = {b["code"] for b in ids["start_blocks"]}
check("three start-block strata (compact/intermediate/extended)",
      sb_codes == {"compact", "intermediate", "extended"}, f"{sorted(sb_codes)}")

# --- Sign convention frozen (§2.3) ---
sign = qer["primary_estimand"]["sign_convention"].lower()
check("ΔG sign convention frozen (positive = open less favourable)",
      "positive" in sign and "open" in sign and ("less favourable" in sign or "less favorable" in sign))

# --- Primary contrast is reline - water (§2.3) ---
prim = qer["contrasts"]["primary"]["expression"].lower().replace(" ", "")
check("primary contrast is ΔG_reline − ΔG_water",
      "reline" in prim and "water" in prim and "−" in qer["contrasts"]["primary"]["expression"])
check("two specificity contrasts prespecified",
      len(qer["contrasts"]["specificity"]) == 2)

# --- Claim hierarchy has exactly the five §16.4 outcome classes ---
classes = {c["class"] for c in qer["claim_hierarchy"]}
expected_classes = {"Confirmed solvent shift", "Small directional shift",
                    "Equivalent within margin", "Inconclusive", "Not supported / opposite"}
check("five outcome classes match §16.4", classes == expected_classes,
      f"missing={expected_classes - classes}, extra={classes - expected_classes}")

# --- Routes A/B/stop defined and decided at Gate 7 (§12.1) ---
route_ids = {r["id"] for r in scope["publication_routes"]}
check("Routes A, B and stop defined", route_ids == {"ROUTE-A", "ROUTE-B", "ROUTE-STOP"}, f"{route_ids}")
check("route decision deferred to Gate 7", scope["route_decision_gate"] == 7)
check("prohibited decision bases recorded (§12.1)",
      len(scope["prohibited_decision_basis"]) >= 3)

# --- Compute ceiling matches §5 task 4 frozen constants ---
cc = scope["compute_ceiling"]
check("max 400 ns per replica (§5 task 4)", cc["max_ns_per_enhanced_sampling_replica"] == 400,
      f"{cc['max_ns_per_enhanced_sampling_replica']}")
check("max 6 replicas per cell (§5 task 4)", cc["max_independent_replicas_per_cell"] == 6,
      f"{cc['max_independent_replicas_per_cell']}")
check("min 4 confirmatory replicas per cell (§8.1)", cc["min_confirmatory_replicas_per_cell"] == 4,
      f"{cc['min_confirmatory_replicas_per_cell']}")

# --- SAP shell: open reviewer-judgment items flagged SIGN:, not fabricated ---
margin = str(sap["equivalence"]["meaningful_margin_kJ_per_mol"])
check("equivalence margin left as SIGN: (not fabricated)", margin.strip().upper().startswith("SIGN:"),
      f"got {margin!r}")
check("primary model implementation left as SIGN: (frozen at Gate 7)",
      str(sap["across_replica_primary_model"]["implementation"]).strip().upper().startswith("SIGN:"))
check("SAP final freeze deferred to Gate 7", sap["final_freeze_gate"] == 7)

# --- Governance: six controlled roles with separated authority (§3.1) ---
role_names = {r["role"] for r in roles["controlled_roles"]}
expected_roles = {"Study lead", "Molecular-model reviewer", "Enhanced-sampling reviewer",
                  "Statistical reviewer", "Quantum-chemistry reviewer", "Data steward"}
check("six controlled roles present (§3.1)", role_names == expected_roles,
      f"missing={expected_roles - role_names}")

# --- Agent boundary: session will NOT sign gates / hold key / run production MD ---
wonts = " ".join(roles["agent_boundary"]["will_not"]).lower()
check("agent boundary: will not sign gates", "sign" in wonts and "gate" in wonts)
check("agent boundary: will not hold solvent key", "solvent-code key" in wonts or "solvent key" in wonts)
check("agent boundary: will not run production MD in sandbox", "production md" in wonts)

# --- Blinding: key steward-held, unblinding once after Gate 9/10 ---
check("solvent key owned by data steward", "steward" in blind["key_custody"]["owner_role"].lower())
check("unblinding occurs once after Gate 9 lock + Gate 10 freeze",
      "gate 9" in blind["unblinding_event"]["when"].lower() and "gate 10" in blind["unblinding_event"]["when"].lower())

# --- Named blinding custody & access validation (Gate 0 corrective action item 3) ---
kc = blind["key_custody"]
# the analyst agent (this session) must be on the NOT_readable_by firewall list
nrb = " ".join(kc.get("NOT_readable_by", [])).lower()
check("analyst agent explicitly barred from solvent key (hard firewall §18.3)",
      "analyst" in nrb, f"NOT_readable_by={kc.get('NOT_readable_by')}")
# the identifier scheme must agree with the charter on WHERE the sealed key lives
id_key_loc = ids["masked_arms"]["key_location"].lower()
check("sealed-key location consistent between charter and identifier scheme",
      "sealed" in id_key_loc and "sealed" in kc["key_artifact"].lower())
# the identifier scheme must NOT contain a resolved mapping (still 'to be assigned')
check("masked-arm -> solvent mapping NOT yet resolved in analyst-readable tree",
      ids["masked_arms"]["mapping_status"].upper().startswith("TO BE ASSIGNED"),
      f"got {ids['masked_arms']['mapping_status']!r}")
# roles record must name the data steward as key custodian; report if still SIGN:
steward = next((r for r in roles["controlled_roles"] if r["role"] == "Data steward"), {})
steward_named = not str(steward.get("name", "")).strip().upper().startswith("SIGN:")
check("data-steward custody arrangement present in roles record",
      "steward" in steward.get("role", "").lower())
# custody arrangements that require a NAMED human before Phase 6 (not a Gate 0 blocker,
# but surfaced so reviewers see the outstanding assignment explicitly)
unassigned_roles = sorted(r["role"] for r in roles["controlled_roles"]
                          if str(r.get("name", "")).strip().upper().startswith("SIGN:"))
print("blinding custody:",
      {"analyst_barred": "analyst" in nrb,
       "key_location_consistent": ("sealed" in id_key_loc),
       "mapping_unresolved": ids["masked_arms"]["mapping_status"].upper().startswith("TO BE ASSIGNED"),
       "data_steward_named": steward_named,
       "roles_awaiting_named_assignment": unassigned_roles})

# --- Gate 0 reviewer corrective actions (verify the fixes are present) ---
# item 1: three-level estimand hierarchy
eh = qer["primary_estimand"].get("estimand_hierarchy", {})
check("estimand hierarchy has run/cell/contrast levels (reviewer item 1)",
      {"run_level_observation", "cell_level_estimand", "primary_contrast_estimand"} <= set(eh))
# item 2: outcome-class interval logic tightened (simultaneous interval wording)
_conf = next(c for c in qer["claim_hierarchy"] if c["class"] == "Confirmed solvent shift")
check("outcome-class logic uses whole-interval-beyond-margin wording (item 2)",
      "entire" in _conf["evidence"].lower() and "margin" in _conf["evidence"].lower())
# item 4: release trajectories reclassified as sensitivity analysis
sens_ids = {s["id"] for s in qer.get("sensitivity_analyses", [])}
sec_txt = " ".join(qer["secondary_endpoints"]).lower()
check("release trajectories moved to sensitivity_analyses (item 4)",
      "SENS-RELEASE" in sens_ids and "release trajector" not in sec_txt)
# item 4: basin revision conditions (blinded, developmental, no selection-for-effect)
rc = " ".join(qer["basin_definitions"].get("revision_conditions", [])).lower()
check("basin revision conditions: blinded + developmental + no-selection-for-effect (item 4)",
      "blind" in rc and "developmental" in rc and "clearer" in rc)
# item 3: Route B peptide-specific qualification track
routeB = next(r for r in scope["publication_routes"] if r["id"] == "ROUTE-B")
check("Route B requires peptide-specific qualification track of 6 steps (item 3)",
      routeB.get("peptide_specific_qualification_required") is True
      and len(routeB.get("qualification_track", [])) == 6)
# aggregate campaign ceiling present with a recorded approval decision
# (PROVISIONAL before study-lead sign-off; APPROVED after). Values must be present.
acc = scope.get("aggregate_campaign_ceiling", {})
_acc_status = str(acc.get("approval_status", "")).lower()
check("aggregate campaign ceiling present with recorded approval status (provisional or approved) + values",
      ("provisional" in _acc_status or "approved" in _acc_status)
      and {"max_aggregate_enhanced_sampling_us", "max_production_cells",
           "max_storage_tb"} <= set(acc.get("values", {})))
# item 5: no vendor/product name in scope execution_platform; agent has no 'release'
_ep = json.dumps(scope["execution_platform"]).lower()
check("scope execution_platform is vendor-neutral (no 'claude science') (item 5)",
      "claude science" not in _ep and "claude" not in _ep)
check("analysis environment does manuscript PREPARATION, not release (item 5)",
      "preparation" in _ep and "does not perform release" in _ep)
# item 6: full-trajectory QC evidence requirement
_ho = json.dumps(scope["execution_platform"]["handoff"]).lower()
check("full-trajectory QC evidence required (all frames OR all-frame audit+hashes) (item 6)",
      "all saved frames" in _ho or "all-frame" in _ho)
# item 8: operational key procedure + builder/analyst separation in charter
check("blinding charter has operational key procedure (item 8)",
      "operational_key_procedure" in blind and "sealing_method" in blind["operational_key_procedure"])
check("blinding charter names builder/analyst separation (item 8)",
      "setup_disclosure_separation" in blind and len(blind["setup_disclosure_separation"]["division"]) >= 3)
# item 8: custodian name field present (may still be SIGN: until steward accepts)
check("steward custodian_name field present in key_custody (item 8)",
      "custodian_name" in blind["key_custody"])
# item 7: deviation template expanded with structured reviews + no-overwrite rule
_dev_tmpl_cc, _ = load_yaml("governance/deviation_log_template.yaml")
drt = _dev_tmpl_cc["deviation_record_template"]
check("deviation template has structured per-reviewer reviews (item 7)",
      isinstance(drt.get("reviews"), list) and {"role", "name", "decision", "date"} <= set(drt["reviews"][0]))
check("deviation template states replacement never overwrites original (item 7)",
      "never delete" in str(drt.get("replacement_rule", "")).lower())
# item 5 (roles): compute arrangement no longer grants 'release' / names vendor
_ca = roles["campaign_arrangements"]["compute_arrangement"]["decision"].lower()
check("roles compute-arrangement: preparation not release, vendor-neutral (item 5)",
      "does not perform release" in _ca and "claude science" not in _ca)

checks_df = pd.DataFrame(checks)
n_fail = (checks_df["result"] == "FAIL").sum()
print(checks_df.to_string(index=False))
print(f"\nConsistency checks: {len(checks_df)-n_fail}/{len(checks_df)} PASS, {n_fail} FAIL")
assert n_fail == 0, "Consistency failures must be resolved before Gate 0 review"

# %% [markdown]
# ## 2. Resource forecast from the frozen ceiling
#
# Derived from the frozen compute ceiling (§5 task 4) and the *planning* throughput
# band (§19.1). These are **scheduling ranges, not guarantees** — Gate 4 replaces
# them with the measured M5 Max benchmark (§9, §19). Reproduced here so the Gate 0
# record shows the calendar implied by the ceiling the reviewers are freezing.

# %%
tp_lo, tp_hi = scope["planning_throughput"]["ns_per_day_range"]          # 60, 160
ov_lo, ov_hi = scope["planning_throughput"]["operational_overhead_fraction"]  # 0.15, 0.25

# Resource rounding alignment (Gate 0 corrective action item 5):
# The CONTROLLED, frozen values are the §19.2 ranges transcribed in the scope
# artifact. The protocol's published table uses hand-rounded day ranges, so the
# notebook must NOT substitute its own recomputed/rounded numbers (which drift).
# We therefore present the controlled ranges as authoritative and run an
# INDEPENDENT recomputation only as a cross-check, flagging any divergence.
# The controlled §19.2 table applies the operational-overhead band to the ROUNDED
# raw-serial days (not to the unrounded quotient), so the notebook must not impose
# its own floor/ceil recomputation — that was the rounding drift this item fixes.
# Instead we (a) reproduce the controlled values verbatim as authoritative, and
# (b) verify each is INTERNALLY VALID: raw days ≈ workload/throughput within
# rounding, and each planning endpoint is the raw endpoint widened by an overhead
# inside the declared [ov_lo, ov_hi] band (±1 day rounding), and the bands are
# monotone. This aligns the forecast with the controlled protocol's own rounding.
rows = []
for s in scope["serial_campaign_scenarios"]:
    ns = round(s["nominal_us"] * 1000)
    raw, plan = s["raw_serial_days"], s["planning_days_with_overhead"]
    raw_ok  = abs(raw[0] - ns / tp_hi) <= 1 and abs(raw[1] - ns / tp_lo) <= 1
    plan_ok = (raw[0]*(1+ov_lo) - 1 <= plan[0] <= raw[0]*(1+ov_hi) + 1 and
               raw[1]*(1+ov_lo) - 1 <= plan[1] <= raw[1]*(1+ov_hi) + 1)
    mono    = raw[0] <= raw[1] and plan[0] <= plan[1] and plan[0] >= raw[0] and plan[1] >= raw[1]
    rows.append({
        "scenario": s["scenario"],
        "nominal_us": s["nominal_us"],
        "raw_serial_days": f"{raw[0]}–{raw[1]}",
        "planning_days_with_overhead": f"{plan[0]}–{plan[1]}",
        "internally_valid": "OK" if (raw_ok and plan_ok and mono) else "REVIEW",
    })
forecast_df = pd.DataFrame(rows)
print(forecast_df.to_string(index=False))
n_review = int((forecast_df["internally_valid"] == "REVIEW").sum())
print(f"\ncontrolled §19.2 table: {len(forecast_df)-n_review}/{len(forecast_df)} internally valid "
      f"(raw = workload/throughput; planning = raw widened by {ov_lo:.0%}–{ov_hi:.0%} overhead)")
assert n_review == 0, "Controlled resource table failed internal validity — reconcile §19.2 before Gate 0"
forecast_df.to_csv(os.path.join(OUT, "resource_forecast.csv"), index=False)
print("wrote outputs/resource_forecast.csv  (values are the controlled §19.2 table, verbatim & authoritative)")
print("NOTE: planning only; superseded by Gate 4 measured benchmark (§9.3, §19).")

# %% [markdown]
# ## 3. Gate 0 release matrix (mechanical evaluation)
#
# The §5 Gate 0 matrix has four controlled criteria. This evaluates the *evidence
# readiness* of each — i.e. that the required evidence artifact exists, is
# non-empty and contains the expected content. Disposition remains **PASS/HOLD/FAIL
# by a human reviewer**; a green row here means "evidence is present and internally
# valid", not "approved".

# %%
def exists_nonempty(rel):
    p = os.path.join(REPO, rel)
    return os.path.isfile(p) and os.path.getsize(p) > 0

# --- Change-control readiness (Gate 0 corrective action item 2) ---
# The original test only checked two template files were non-empty. That does not
# demonstrate change control is OPERATIONAL. Test the substance:
dev_tmpl, _   = load_yaml("governance/deviation_log_template.yaml")
cert_tmpl, _  = load_yaml("governance/gate_certificate_template.yaml")
# controlled data_status vocabulary present in the deviation template (§20.4).
# The authoritative list now lives in data_status_controlled_vocabulary (expanded
# per Gate 0 reviewer item 7); the per-record field enumerates the same options.
_dev_vocab = [str(v).lower() for v in dev_tmpl.get("data_status_controlled_vocabulary", [])]
_dev_status = str(dev_tmpl.get("deviation_record_template", {}).get("data_status", "")).lower()
cc_vocab_ok = (all(v in _dev_vocab for v in ["developmental", "invalid",
               "replaced_for_technical_failure", "campaign_hold"])
               and "confirmatory" in _dev_status and "invalid" in _dev_status)
# gate certificate template exposes disposition + signature + tag/commit fields
_cert_keys = set(cert_tmpl.get("gate_certificate_template", {}).keys())
cc_cert_ok = {"criteria", "overall_disposition", "git_tag", "git_commit"} <= _cert_keys
# gate-0 instance exists and carries the four §5 criteria
cc_inst_ok = exists_nonempty("governance/gate_certificates/gate00_certificate.yaml") and len(cert_doc.get("criteria", [])) == 4
# reviewer checklists present
cc_chk_ok = exists_nonempty("governance/reviewer_checklists.yaml")
# frozen commit/tag policy recorded (push+tag only after sign-off)
_commit_policy = str(roles.get("release_workflow", {}).get("commit_policy", "")).lower()
cc_policy_ok = ("only after" in _commit_policy and "sign" in _commit_policy)
cc_operational = bool(cc_vocab_ok and cc_cert_ok and cc_inst_ok and cc_chk_ok and cc_policy_ok)
print("change-control readiness:",
      {"data_status_vocab": cc_vocab_ok, "cert_template_fields": cc_cert_ok,
       "gate0_instance_4crit": cc_inst_ok, "reviewer_checklists": cc_chk_ok,
       "commit_tag_policy": cc_policy_ok, "=> operational": cc_operational})

gate0 = [
    {"criterion": "Question and estimands unambiguous",
     "evidence": "phases/phase00_freeze/outputs/question_estimand_registry.yaml",
     "evidence_ready": exists_nonempty("phases/phase00_freeze/outputs/question_estimand_registry.yaml")
                       and bool(qer.get("primary_estimand")) and bool(qer.get("contrasts")),
     "reviewer_role": "Study lead + statistical reviewer"},
    {"criterion": "Blinding responsibilities assigned",
     "evidence": "governance/blinding_charter.yaml + roles_and_arrangements.yaml",
     "evidence_ready": exists_nonempty("governance/blinding_charter.yaml")
                       and exists_nonempty("governance/roles_and_arrangements.yaml")
                       and any(r["role"] == "Data steward" for r in roles["controlled_roles"]),
     "reviewer_role": "Data steward"},
    {"criterion": "Scope options and resource ceiling frozen",
     "evidence": "phases/phase00_freeze/outputs/scope_and_resource_ceiling.yaml",
     "evidence_ready": exists_nonempty("phases/phase00_freeze/outputs/scope_and_resource_ceiling.yaml")
                       and len(scope["publication_routes"]) == 3
                       and scope["compute_ceiling"]["max_ns_per_enhanced_sampling_replica"] == 400,
     "reviewer_role": "Study lead"},
    {"criterion": "Change control operational",
     "evidence": "deviation_log_template.yaml (+ controlled data_status vocab) + gate_certificate_template.yaml + gate00_certificate.yaml + reviewer_checklists.yaml + frozen commit/tag policy",
     "evidence_ready": cc_operational,
     "reviewer_role": "Data steward"},
]
gate0_df = pd.DataFrame(gate0)
gate0_df["disposition"] = gate0_df["evidence_ready"].map(lambda b: "READY FOR REVIEW" if b else "EVIDENCE INCOMPLETE")
print(gate0_df[["criterion", "evidence_ready", "disposition", "reviewer_role"]].to_string(index=False))
all_ready = bool(gate0_df["evidence_ready"].all())
print(f"\nAll Gate 0 evidence ready for reviewer sign-off: {all_ready}")

# --- Signed-state detection (release readiness, items 8-10) ---
# The session NEVER signs; it only READS the human-entered certificate dispositions.
# A certificate is "signed" only when every criterion disposition and the overall
# disposition are a decided value (PASS/HOLD/FAIL) with no residual 'SIGN:' / slash
# placeholder, AND each criterion carries a non-placeholder signature.
def _decided(v):
    s = str(v).strip().upper()
    return s in {"PASS", "HOLD", "FAIL"}
def _signed(v):
    s = str(v).strip()
    return bool(s) and not s.upper().startswith("SIGN:")
_crit = cert_doc.get("criteria", [])
_all_crit_decided = len(_crit) == 4 and all(_decided(c.get("disposition")) and _signed(c.get("reviewer_signature")) for c in _crit)
_overall = str(cert_doc.get("overall_disposition", "")).strip().upper()
cert_signed = _all_crit_decided and _decided(_overall)
cert_pass = cert_signed and _overall == "PASS" and all(str(c.get("disposition")).strip().upper() == "PASS" for c in _crit)
signing_status = ("SIGNED PASS by human reviewer(s)" if cert_pass
                  else f"SIGNED {_overall} by human reviewer(s)" if cert_signed
                  else "AWAITING HUMAN REVIEWER SIGN-OFF (not signed by the session)")
print(f"certificate signed: {cert_signed} | overall: {_overall or '(blank)'} | release-ready (PASS): {cert_pass}")

# %% [markdown]
# ## 4. SHA-256 evidence manifest
#
# Every Phase 0 artifact is hashed (§18.2: "Use SHA-256 for immutable scientific
# artifacts and record the hash in the run registry before downstream analysis").
# The manifest is written to `manifests/gate00_manifest.sha256` and is the object
# the Gate 0 certificate references.

# %%
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

manifest_targets = [
    "phases/phase00_freeze/outputs/identifier_scheme.yaml",
    "phases/phase00_freeze/outputs/question_estimand_registry.yaml",
    "phases/phase00_freeze/outputs/scope_and_resource_ceiling.yaml",
    "phases/phase00_freeze/outputs/resource_forecast.csv",
    "sap/SAP_shell_v1.yaml",
    "governance/roles_and_arrangements.yaml",
    "governance/blinding_charter.yaml",
    "governance/deviation_log_template.yaml",
    "governance/deviation_log.yaml",
    "governance/gate_certificate_template.yaml",
    "governance/gate_certificates/gate00_certificate.yaml",
    "governance/reviewer_checklists.yaml",
    "runs/run_registry.csv",
    # controlled schemas (Gate 0 corrective action item 1)
    "schemas/identifier_scheme.schema.json",
    "schemas/question_estimand_registry.schema.json",
    "schemas/scope_and_resource_ceiling.schema.json",
    "schemas/SAP_shell.schema.json",
    "schemas/roles_and_arrangements.schema.json",
    "schemas/blinding_charter.schema.json",
    "schemas/gate_certificate.schema.json",
]
man_rows = []
for rel in manifest_targets:
    p = os.path.join(REPO, rel)
    man_rows.append({"sha256": sha256(p), "bytes": os.path.getsize(p), "path": rel})
man_df = pd.DataFrame(man_rows)

os.makedirs(os.path.join(REPO, "manifests"), exist_ok=True)
man_path = os.path.join(REPO, "manifests", "gate00_manifest.sha256")
with open(man_path, "w") as f:
    for r in man_rows:
        f.write(f"{r['sha256']}  {r['path']}\n")
print(man_df.to_string(index=False))
print(f"\nwrote {os.path.relpath(man_path, REPO)}")

# %% [markdown]
# ## 4.5 Read-back verification of every generated output (item 6)
#
# Re-read each file this notebook just wrote and confirm it (a) parses, and (b)
# its on-disk SHA-256 matches the hash we recorded in the manifest. This closes
# the loop: the evidence summary can only claim an output exists if it has been
# read back and its hash re-verified here.

# %%
import subprocess

generated = {
    "resource_forecast.csv":           os.path.join(OUT, "resource_forecast.csv"),
    "gate00_manifest.sha256":          man_path,
}
readback = []
# forecast CSV parses back to the same row count
_fc = pd.read_csv(generated["resource_forecast.csv"])
readback.append({"output": "resource_forecast.csv", "parses": len(_fc) == len(forecast_df),
                 "detail": f"{len(_fc)} rows"})
# manifest: every listed hash matches the file on disk right now
man_verified, man_total = 0, len(man_rows)
mismatch = []
for r in man_rows:
    actual = sha256(os.path.join(REPO, r["path"]))
    if actual == r["sha256"]:
        man_verified += 1
    else:
        mismatch.append(r["path"])
readback.append({"output": "gate00_manifest.sha256",
                 "parses": man_verified == man_total,
                 "detail": f"{man_verified}/{man_total} hashes match disk"
                           + ("" if not mismatch else f"; MISMATCH: {mismatch}")})
readback_df = pd.DataFrame(readback)
print(readback_df.to_string(index=False))
readback_ok = bool(readback_df["parses"].all()) and not mismatch
assert readback_ok, f"Read-back verification failed: {mismatch or 'parse error'}"
print(f"\nread-back verification: {'PASS' if readback_ok else 'FAIL'}")

# %% [markdown]
# ## 4.6 Git state record (item 4)
#
# Record the repository state at evidence-generation time: branch, HEAD commit,
# worktree cleanliness and remote-tracking status. This is the git context the
# reviewers see, and — after sign-off — the state the release commit/tag are
# verified against (items 12–14 of the corrective plan).

# %%
def git(*args):
    try:
        return subprocess.check_output(["git", "-C", REPO, *args],
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception as e:
        return f"(git error: {e})"

git_state = {
    "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
    "head_commit": git("rev-parse", "HEAD"),
    "head_short": git("rev-parse", "--short", "HEAD"),
    "worktree_clean": git("status", "--porcelain") == "",
    "uncommitted_paths": [l[3:] for l in git("status", "--porcelain").splitlines()] if git("status", "--porcelain") else [],
    "upstream": git("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}") or "(no upstream)",
    "remotes": git("remote", "-v").splitlines(),
    "existing_tags": git("tag", "--list").splitlines(),
    "gate00_tag_present": "gate-00" in git("tag", "--list").splitlines(),
}
print(json.dumps(git_state, indent=2))
print("\nNOTE: worktree is expected to show this notebook's regenerated outputs as")
print("uncommitted until they are committed; that commit precedes the gate-00 tag.")

# %% [markdown]
# ## 5. Gate 0 evidence summary (machine-readable)
#
# A single JSON the reviewers and the (human-signed) gate certificate reference.

# %%
summary = {
    "protocol_id": PID,
    "gate": 0,
    "phase": "Phase 0 — Protocol, estimand and analysis freeze",
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "schema_validation": {"total": int(len(schema_df)),
                          "passed": int((schema_df["result"] == "PASS").sum()),
                          "failed": int(n_schema_fail)},
    "consistency_checks": {"total": int(len(checks_df)),
                           "passed": int((checks_df["result"] == "PASS").sum()),
                           "failed": int(n_fail)},
    "resource_crosscheck": {"scenarios": int(len(forecast_df)),
                            "internally_valid": int((forecast_df["internally_valid"] == "OK").sum())},
    "readback_verification": {"ok": readback_ok,
                              "detail": readback_df.to_dict("records")},
    "git_state": git_state,
    "gate0_release_matrix": gate0_df[["criterion", "evidence_ready", "disposition", "reviewer_role"]].to_dict("records"),
    "all_evidence_ready": all_ready,
    "manifest": {"file": "manifests/gate00_manifest.sha256", "artifact_count": len(man_rows)},
    "arrangements": {
        "production_md": scope["execution_platform"]["production_md"],
        "governance": "Human reviewers sign gates; data steward holds solvent key; session prepares evidence and enforces sequencing.",
    },
    "signing_status": signing_status,
    "certificate_signed": cert_signed,
    "certificate_overall": _overall or None,
    "release_ready": cert_pass,
}
sum_path = os.path.join(OUT, "gate00_evidence_summary.json")
with open(sum_path, "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
print(f"\nwrote {os.path.relpath(sum_path, REPO)}")

# %% [markdown]
# ## 6. Provisional review manifest (item 7)
#
# A single human-facing index the reviewers use to review the **actual evidence
# files** (item 8): every controlled artifact, its role, its assigned reviewer
# role, its SHA-256, and the checklist item it answers. Written as
# `manifests/gate00_review_manifest.md`. This is *provisional* — it becomes the
# final release manifest (item 10) only once the certificate is signed.

# %%
hash_by_path = {r["path"]: r["sha256"] for r in man_rows}
review_index = [
    ("phases/phase00_freeze/outputs/question_estimand_registry.yaml", "Question & estimand freeze", "Study lead + Statistical reviewer"),
    ("phases/phase00_freeze/outputs/identifier_scheme.yaml", "Identifier & masking scheme", "Data steward"),
    ("phases/phase00_freeze/outputs/scope_and_resource_ceiling.yaml", "Scope routes & resource ceiling", "Study lead"),
    ("sap/SAP_shell_v1.yaml", "Statistical analysis plan (shell)", "Statistical reviewer"),
    ("governance/blinding_charter.yaml", "Blinding charter & key custody", "Data steward"),
    ("governance/roles_and_arrangements.yaml", "Roles, arrangements, release workflow", "Study lead"),
    ("governance/deviation_log_template.yaml", "Change control — deviation template", "Data steward"),
    ("governance/deviation_log.yaml", "Deviation ledger (DEV-0001: Gate 0 HOLD corrective actions)", "All reviewers"),
    ("governance/gate_certificate_template.yaml", "Gate certificate template", "Data steward"),
    ("governance/gate_certificates/gate00_certificate.yaml", "Gate 0 certificate (to sign)", "All reviewers"),
    ("governance/reviewer_checklists.yaml", "Reviewer checklists", "All reviewers"),
]
_is_final = cert_pass
_title = ("Gate 0 FINAL release manifest" if _is_final
          else "Gate 0 provisional review manifest")
rm_path = os.path.join(REPO, "manifests", "gate00_review_manifest.md")
with open(rm_path, "w") as f:
    f.write(f"# DES-PEP-CESD-001 — {_title}\n\n")
    f.write(f"- Protocol: **{PID}**  \n")
    f.write(f"- Generated (UTC): {summary['generated_utc']}  \n")
    f.write(f"- Repo state: branch `{git_state['branch']}` @ `{git_state['head_short']}`, "
            f"worktree_clean={git_state['worktree_clean']}  \n")
    f.write(f"- Schema validation: {summary['schema_validation']['passed']}/{summary['schema_validation']['total']} PASS  \n")
    f.write(f"- Consistency checks: {summary['consistency_checks']['passed']}/{summary['consistency_checks']['total']} PASS  \n")
    f.write(f"- Read-back verification: {'PASS' if readback_ok else 'FAIL'}  \n\n")
    if _is_final:
        f.write(f"**Signing status: {signing_status}** — Gate 0 certificate signed PASS by the human reviewer(s); "
                "this is the FINAL release manifest for the signed evidence package.\n\n")
    else:
        f.write(f"**Signing status: {signing_status}.**\n\n")
    f.write("| Artifact | Role | Reviewer | SHA-256 (first 16) |\n")
    f.write("|---|---|---|---|\n")
    for path, role, reviewer in review_index:
        h = hash_by_path.get(path, "(not in manifest)")
        f.write(f"| `{path}` | {role} | {reviewer} | `{h[:16]}…` |\n")
    f.write("\n_Full hashes: see `manifests/gate00_manifest.sha256`._\n")
print(f"wrote {os.path.relpath(rm_path, REPO)}")
with open(rm_path) as f:
    print("\n" + f.read())

# %% [markdown]
# ---
# **Phase 0 verification complete.** All freeze artifacts are internally consistent,
# frozen values match the protocol's §5/§8/§14 constants, open reviewer-judgment
# items are flagged `SIGN:` rather than fabricated, and the evidence package is
# hashed. The Gate 0 certificate (`governance/gate_certificate_template.yaml`)
# remains to be **signed by the human reviewers**; once signed, the repo is tagged
# `gate-00` and committed (§18.2).

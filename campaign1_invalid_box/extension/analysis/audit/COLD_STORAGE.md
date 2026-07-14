# Trajectory cold storage — 2026-07-14

All 49 DCDs (245 GB) moved off local disk. Verified sha256, 49/49 match.

Location:   /Volumes/M3_Backup/des-peptide-study-cold/20260714
Manifest:   extension/analysis/audit/trajectory_manifest_sha256.txt
Topologies: RETAINED on local disk. A DCD is unreadable without them.

## STATUS: INVALID FOR THEIR SCIENTIFIC PURPOSE — RETAINED AS EVIDENCE

BOX-01 (30 A cubic box) is violated in all 35 production runs. Closest
solute-image approach 1.77 A; 35/35 runs breach the 10 A minimum-image
criterion. The box suppresses extended conformations — i.e. it biases the
SASA ensemble against the very observable the study measures.

See assumptions.yaml (BOX-01, SAMPLING-01, RDF-01, PHASE4-01) and
extension/analysis/audit/minimum_image_audit.csv.

DO NOT DELETE. These are the substrate of the methods post-mortem.

## THREE DISTINCT SETS — DO NOT CONFLATE

| Path                              | N  | Construct              | Notes |
|-----------------------------------|----|------------------------|-------|
| systems/                          |  9 | **43-mer parent** (L2-GGE, MNPVDHPHGGGEGRAPIGRKKPATPWGYPALSKCFFFYYLNIQ) | Phase-1/3 pilot. 3x3 systems. |
| extension/trajectories/           |  5 | 9-mer constructs       | E1 rep1. Only 5 of 9 systems. Superseded intermediate. |
| extension/trajectories_extended/  | 35 | 9-mer constructs       | Production. 4 starts x 9 systems MINUS GGE_reline/open (never run) = 35. |

The pilot ran a DIFFERENT MOLECULE from the extension campaign. This is the
cause of the 5.8-8.9x SASA baseline discrepancy previously logged as an
unresolved "motif-atom-selection definition flag". It is not an indexing bug.
Resolved as SASA-01.

Constructs (extension): Ace-EEEGGEIVF-NMe (11 units), Ace-LYQCMEFVR-NMe (11),
Ace-NPYIYK-NMe (8).

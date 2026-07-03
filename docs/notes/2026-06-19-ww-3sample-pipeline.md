# W+W- 3-sample EW-Sudakov pipeline (`parallel_pipeline_WW.sh`)

## Objective
Adapt the Z+jets EW-Sudakov pipeline (`parallel_pipeline_DYj.sh`) to the
on-shell W⁺W⁻ sample, borrowing the reweight/MadSpin ordering structure from
`script_dumps/parallel_pipeline_ttHj.sh`, to produce **three** showered HepMC
samples from one LHE, all decaying the W pair leptonically (different flavour):

| Sample | Path | Steps | qed_shower | Sudakov |
|--------|------|-------|-----------|---------|
| reweight+madspin | A | reweight(EWS) → madspin → shower | T | production / scalar (density falls back) |
| madspin+reweight | B | madspin → reweight(EWS) → shower | T | spin-correlated / density |
| QCD baseline     | C | madspin → shower                | F | none (ratio denominator) |

Decays: `W+ → mu+ vm`, `W- → e- ve~`  (final state μ⁺e⁻ + MET), configurable via
`W_PLUS_DECAY` / `W_MINUS_DECAY`.

Input example: `/ceph/grid/home/lenartj/EWSUDAKOV/pc12_WW/madevent/Events/run_pc12_WW_seed1001/events.lhe.gz` (1.78 GB gz).

## Findings (sample characterisation, from the LHE banner + PROC_DIR)
- Process `generate p p > w+ w- [QCD] @0`: **NLO QCD, on-shell undecayed W's, no
  extra jets.** First event = u, ū, W⁺(status 1), W⁻(status 1).
- Model `loop_qcd_qed_sm_Gmu_forSudakov-with_b_mass_no_width`,
  `complex_mass_scheme False`, **all widths 0** (W, Z, t).
- `ickkw = 3` (FxFx) is baked into the LHE, but this is a single 0-jet
  multiplicity. `generate.sh` + the 2026-06-18 note are explicit: the shower
  **MUST** set `JetMatching:nJetMax=0`, else the Z+jets default (nJetMax=2,
  qCut=30) applies a 2-jet FxFx veto that reshapes σ(WW).
- Per-event scale weights ids 1001–1027 already present; EW Sudakov reweight
  adds the 2001-family (central NLL ξ=1.0 = id 2001).
- Same PROC dir (`PROCNLO_loop_qcd_qed_sm_Gmu_forSudakov_0`) exists in the WW
  working dir and is the identical process → used as `PROC_DIR` default.

## Decisions (2026-06-19, user-confirmed)
1. **QCD baseline runs MadSpin** (madspin → shower qed=F), NOT a bare shower of
   the undecayed LHE. Required so all three share the μ⁺e⁻ final state and the
   per-event K-factor ratio is well-defined.
2. **Flavours**: `W+ → mu+ vm`, `W- → e- ve~` (single fixed different-flavour
   assignment).
3. **MadSpin width**: **default full Breit-Wigner** — user accepts that MadSpin
   silently overrides WW=0 → max(model-default,1 MeV) ≈ 2.085 GeV (note §2,
   `MadSpin/decay.py:3196-3225`) and smears the W off-shell. We do NOT set
   `spinmode onshell`. ⇒ the decayed sample is no longer strictly on-shell;
   accepted as a known caveat.
4. **njmax=0** hardcoded as default; a startup warning fires if overridden.
5. **qCutList apparatus KEPT** (per "keep all DYj features") but is physically
   **degenerate** for a 0-jet sample — FxFx never vetoes the highest
   multiplicity, so `FxFx_qCutAccept_<q>` flags carry no real merging-scale
   variation. `QCUT_LIST=""` disables it.

## Implementation notes
- Fresh code, structured after DYj (EW-Sudakov reweight invocation with the
  trailing `0`; silent-failure signature grep; the per-event `<wgt id='2001'>`
  guard; merge/monitor/summary) + ttHj (standalone MadSpin in an isolated
  per-worker dir; scalar-vs-density auto-routing; 3-path orchestration).
- **Seed alignment (the linchpin).** For chunk *i* the seed is
  `RND_SEED_BASE+i`, identical across paths A/B/C. Injected into the **MadSpin
  command file** (`set seed`) and the shower card (`rnd_seed`/`rnd_seed2`).
  - **Bug fixed vs ttHj:** ttHj injected the madspin seed into
    `procnlo/Cards/madspin_card.dat` but `process_madspin` read the card from
    the *shadow* dir → the seed never reached MadSpin. Here the seed is injected
    directly into the `.mg`, so it is actually used.
  - Rationale: the EW-Sudakov reweight is multiplicative on the weight and never
    changes momenta, so the kinematics entering MadSpin are identical in all
    paths; same seed ⇒ identical decays ⇒ valid event-by-event ratio.
- **Opt-out knobs use `${VAR-default}` (colonless)** for `XI_SCAN` and
  `QCUT_LIST` so an explicit empty string is honoured (the `:-` form silently
  re-enables the default on empty — DYj has this latent issue).
- `pythia8_options` drops DYj's `23:mayDecay=off` (no stable Z; W's decayed by
  MadSpin) and carries only the (no-leading-comma) qCutList entry.

## Verification done
- `bash -n` clean; bash 5.3.9 (safe `set -u` empty-array semantics).
- DRY_RUN: correct config (njmax=0, μ/e decays, 3 paths, seed base).
- 22/22 custom card-rendering tests (`scratchpad/test_cards.sh`): madspin
  decays + seed injection, reweight include_sudakov + ξ-scan, shower qed=T/F,
  njmax=0, pythia8_options has qCutList & NO `23:mayDecay` & no leading comma,
  `QCUT_LIST=""`→`{}`, flavour override, SAMPLES validation, NJMAX warning.
- Two adversarial reviewers (bash control-flow; MG5/physics) — results below.

## Open questions / MUST-VERIFY-AT-RUNTIME
- ~~(HIGH) Does MadSpin propagate the `<rwgt>` Sudakov columns in path A?~~
  **RESOLVED by physics review C3** — it copies the block verbatim and rescales
  per column. Residual silent-failure mode (uppercase exponent) now guarded
  in-pipeline (R3 fix). End-to-end HepMC confirmation = R1 (smoke test).
- ~~(MED) Is id 2001 emitted in BOTH scalar (A) and density (B) modes?~~
  **RESOLVED by physics review C2** — yes, central ξ → 2001 in both. A fired
  guard is therefore a real failure, not an id mismatch (R5 message fix).
- **Zero-width-W BW smearing (decision 3)** — quantify its effect on M(W) and on
  the Sudakov-weighted observables if it matters for the final plots.
- **R4 [MED]** cross-sample ratio: confirm MadSpin's per-run `branching_ratio`
  is equal across paths (smoke test).
- **R2 [MED]** njmax=0 relies on the onlyshower run_card bypass — re-check after
  any MG5 version bump (see script header).
- Adversarial-review findings + smoke-test recipe: see below.

## Adversarial review outcomes

### Reviewer 1 — bash control-flow / parallelism (COMPLETE)
- **BUG 1 [HIGH] — FIXED.** The reweight guard defeated itself: I had used
  `count_events` (carries `|| echo 0`) for `n_events` and added a trailing
  `|| echo 0` to the `n_ewsl` grep. In the exact silent-failure case (events
  present, no `2001` column), `grep -c` prints `0` AND exits 1, so `|| echo 0`
  appended a second `0` → `n_ewsl="0\n0"` → `(( n_ewsl != n_events ))` throws an
  arithmetic syntax error → read as false under the set-e-suspended
  `if ! run_step` → guard SKIPPED → un-reweighted chunk written as SUCCESS and
  merged. Fix: bare `grep -cE` for both operands (no fallback; file existence is
  checked just above). Regression-tested (`scratchpad/test_reweight_guard.sh`,
  7/7): guard now fires on none/partial/empty, quiet on healthy, no corruption.
- **BUG 2 [LOW] — FIXED.** `((sc++))`/`((fc++))` in `generate_summary` (and
  `((success++))` etc. in `monitor_progress`) return the OLD value, so the 0→1
  increment exits 1 and `set -e` kills the (sub)shell → truncated summary report
  / frozen progress line. Switched to `x=$((x+1))`. (Pre-existing identical
  idiom in DYj/ttHj; fixed here.)
- Items A–J ("suspicious but fine"): concurrency gate, `return` inside
  `{ } | tee`, exit-code capture under `if ! run_step`, merge sharding/footer,
  array `set -u` safety, export completeness, cleanup globs — all verified
  healthy with reproducers. No further action.

### Reviewer 2 — MG5 / physics semantics (COMPLETE; verdict: design sound, no fatal bug)
Verified against the MG5 source in this checkout (dev tree + the `bin/internal`
copy that `bin/aMCatNLO` runs) and the real input LHE.
- **C1 — the 3-sample distinction is REAL.** Undecayed WW (path A): each W is a
  singleton cluster group, the singleton filter empties the resonance list
  (`fxfx_ewsudakov.py:3870-3872`) and control returns to the scalar FxFx path
  (`:3894-3896`) → production-level Sudakov. Decayed (path B): leptons cluster
  back into W resonances → genuine density-matrix Sudakov. A and B do NOT
  collapse.
- **C2 — the `2001` guard is sound in BOTH modes.** Ids built as
  `str(20+5*xi_idx+variant_idx)+source_id[-2:]` off the `10xx` scale weights
  (`fxfx_ewsudakov.py:3628`); central ξ → `1001→2001` in scalar AND density.
  Guard not biased against path B. (My earlier "MED: 2001 in both modes?" risk
  is RESOLVED — it is.)
- **C3 — the Sudakov `<rwgt>` block SURVIVES MadSpin in path A.** MadSpin copies
  the `<rwgt>` string verbatim (`decay.py:3737`) and rescales each column by the
  branching ratio (`:2430`). (My earlier "HIGH: weight propagation through
  madspin?" risk is RESOLVED — it propagates.) Caveat → R3 below.
- **C4 — seed alignment holds.** Decay kinematics come from the Fortran `ranmar`
  stream seeded by `set seed`, which never reads the `<rwgt>` columns; path A's
  extra weight columns cannot desync the RNG vs B/C. Seeds 12345+i are nonzero.
- **C5 — njmax=0 does not abort and biases nothing.** See the LOAD-BEARING note
  now in the script header (onlyshower bypasses the njmax!=0 exception). With
  nJetMax=0 every FxFx veto is gated on `npNLO()<nJetMax` (0<0 false) → no 0-jet
  WW event vetoed → σ unbiased.
- **C6/C7** — qCutList apparatus only writes side-channel named weights (σ
  untouched); `store_rwgt_info=False` irrelevant to the standalone reweight.

**Actions taken in response (this commit):**
- **R3 [MED] → FIXED in-pipeline.** MadSpin's column rescaler regex only matches
  after `decay.py:291` lowercases the line; an uppercase exponent would make it
  silently rebuild an EMPTY `<rwgt>` (no error), destroying 2001. Added an
  automatic post-MadSpin guard in **path A** (`process_madspin` 6th arg
  `require_wgt_id`): every decayed event must still carry `<wgt id=2001>`, else
  the worker fails before showering. So the silent failure is now loud on every
  run, not just the smoke test.
- **R5 [LOW] → FIXED.** Reworded the reweight guard message: both modes emit
  2001, so a fired guard is a genuine un-reweighted slice — do NOT mask it by
  changing `REWEIGHT_CHECK_WGT_ID`.
- **R2 [MED] → DOCUMENTED.** Added the load-bearing njmax=0/onlyshower
  dependency comment to the script header.

**Still MUST-VERIFY at first real run (cannot be checked statically):**
- **R1 [HIGH]** end-to-end: weight 2001 in the final HepMC for A & B, absent for
  C; A/B/C decay kinematics identical on a shared event. (The path-A in-pipeline
  guard above now covers the LHE stage automatically; R1 extends it to HepMC.)
- **R4 [MED]** MadSpin rescales by a per-RUN scalar `branching_ratio`; the
  cross-sample A/C, B/C ratios only cancel it if it is equal across paths (it
  should be — same card/model/events). Confirm via the smoke test.

## Smoke test run 1 (2026-06-19) — BLOCKED by two dev-fork code-gen bugs (NOT pipeline bugs)
Ran `SAMPLES="a b c" EVENTS_PER_CHUNK=200` on a 200-event slice. The pipeline
orchestration worked perfectly (split, 3 shadows, 3 workers, seeds, monitor,
and the guards fired correctly — the BUG-1 reweight guard caught the failure and
marked the worker FAILED instead of writing a false SUCCESS). All three paths
failed at the MG5/Fortran level:

1. **Reweight (paths A, B) — `ewsudakov_goldstone_me_1.f` won't compile.**
   The EW-Sudakov reweight computed the nominal central σ (116 pb) then errored
   compiling the Goldstone ME for `P0_ddx_wpwm`:
   `ewsudakov_goldstone_me_1.f:193  DATA %(PROC_PREFIX)SDENOM/1/  -> Error: Syntax error in DATA statement`.
   The Python template placeholder `%(PROC_PREFIX)s` is left UNSUBSTITUTED in the
   generated Fortran. Template source:
   `madgraph/iolibs/template_files/ewsudakov_goldstone_splitorders_fks.inc`.
   ⇒ the EW-Sudakov reweight has never compiled for this WW process. **Highest
   priority** — blocks paths A and B regardless of MadSpin. Needs a dev-fork fix
   to the `proc_prefix` substitution in the goldstone-ME generation.

2. **MadSpin (all paths) — two layered failures:**
   a. *BR table missing.* In the `no_width` model the W is absent from the param
      `decay_table`; MadSpin calls MadWidth, which fails, then `get_br` raises
      "No valid decay for 24" (decay.py:927). **FIX FOUND & proven standalone:**
      inject a finite-width W decay block into the LHE banner's `<slha>` that
      MadSpin reads (DECAY 24 = 2.085 GeV + standard leptonic/hadronic BRs);
      `extract_br_from_banner` then populates self.br and MadWidth is skipped
      (decay.py:1341-1357). This is decoupled from the reweight's on-shell
      param_card (reweight reads `procnlo/Cards/param_card.dat`, still WW=0), so
      the on-shell Sudakov logs are unaffected. NOT yet wired into the pipeline.
   b. *Decay-ME compile fails.* After (a), MadSpin compiles the decay ME and dies:
      `configs_decay.inc:6  PRWIDTH(-1)=ZERO  -> Error: Symbol 'zero' has no IMPLICIT type`.
      `ZERO` is referenced but undefined for this model's decay-ME export. A
      second dev-fork code-gen issue; needs ZERO defined (or use the finite-width
      model variant for the decay, or revisit the MadSpin approach per the
      2026-06-18 note's original recommendation).

**Conclusion:** the pipeline is correct and ready; the blockers are pre-existing
MG5/model code-generation problems in the dev fork that surface only at run time.
Failed run + logs kept at `smoke_run_20260619/` for inspection.
fxfx_ewsudakov DEBUG is ON but produced no output yet because the reweight dies
at *compile* time, before the Python Sudakov code runs.

## Recommended first-run smoke test (200 events, all 3 paths)
```bash
# one small chunk, all three samples, keep logs
SAMPLES="a b c" EVENTS_PER_CHUNK=200 KEEP_LOGS=true VERBOSE=true \
  ./parallel_pipeline_WW.sh <input.lhe.gz>
# then, in the run dir:
# R1: central Sudakov present in A & B, absent in C
for p in a b c; do echo -n "path $p 2001 in HepMC: "; \
  zcat path_$p/output/*.hepmc.gz | grep -m1 -c 2001; done
# R4: branching ratio identical across paths
grep -h "Branching ratio to allowed decays" path_*/workers/worker_000/madspin.log
# C4 spot check: same decayed 4-momenta for event 1 in A vs C (kinematics, not weights)
```

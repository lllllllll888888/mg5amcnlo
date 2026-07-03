# WW (4FS, on-shell) generation setup — adapting the Z+jets EW-Sudakov area

## Objective
Repurpose `mg5amcnlo_devewsudakov_on_valentin_WW/` (a fresh Jun-18 copy of the
Z+jets EW-Sudakov working area) to generate `p p > w+ w-`:
- on-shell W's, **no extra-jet multiplicities**;
- **4-flavour scheme** (massive b);
- W width = 0;
- biased toward large M(W⁺W⁻) (Sudakov logs dominate the high-mass tail);
- EW Sudakov logs applied downstream by reweighting (as on the Z leg).

`generate.sh` is currently byte-identical to the Zjj copy (`diff` = empty); it still
generates `p p > Z / Z j / Z j j [QCD]`. Nothing is WW-specific yet.

## Findings (4 code-grounded sub-audits, read-only)

### 1. 4FS / massive-b import — SOLVED, no new card needed
- Use `import model loop_qcd_qed_sm_Gmu_forSudakov-with_b_mass_no_width`.
  That ships card sets **MB = 4.5** and **WW = WZ = WT = WH = 0** in one shot
  (`models/loop_qcd_qed_sm_Gmu_forSudakov/restrict_with_b_mass_no_width.dat`).
  Suffix→file mapping verified: `models/import_ufo.py:221,223,231-232`.
- The current `generate.sh:47` imports with **no suffix** ⇒ `restrict_default.dat`
  ⇒ massless b + finite widths (i.e. silently 5FS). So the switch is real.
- **CAVEAT:** the `with_b_mass_no_width` card also moves **MT 173.3→174.3** and
  **MH 125→120** vs default. Must reconcile with the other legs or override.
- `define p`/`j` auto-drop the b when MB≠0 (`madgraph/interface/madgraph_interface.py:6001-6013`)
  — no manual edit, just regenerate the proc dir.
- **`maxjetflavor` is NOT auto-synced to MB** — must set `maxjetflavor 4` by hand
  (`generate.sh:79` currently `5`). Nothing errors if left at 5 → silent 5FS clustering.
- Do **not** carry `set mass 5 0.0` (that is the 5FS line in `generate_FO.sh:196`).
- **PHYSICS RED FLAG:** lhaid 324900 (NNPDF31_nlo_as_0118_luxqed) is a **5-flavour PDF**
  (`…/NNPDF31_nlo_as_0118_luxqed.info`: `NumFlavors: 5`, `MBottom: 4.92`, has a b-distribution).
  4FS matrix elements + 5F PDF double-counts g→bb̄ collinear logs — a genuine scheme
  inconsistency. No 4F luxqed set is installed. b-mass is triple-valued: ME 4.5 / model
  default 4.7 / PDF 4.92. **Decision required.** (For inclusive pp→WW the 4FS↔5FS difference
  is sub-percent — light qq̄→WW dominates — so the 4FS request is likely motivated by
  *consistency with the Sudakov reweighting*, not by the WW rate itself. UNVERIFIED whether
  the reweight strictly needs massive b.)

### 2. MadSpin with zero width — safe here, but a footgun
- MadSpin does **not** honour WW=0: it overwrites it with `max(model-default-width, 1 MeV)`
  ≈ 2.085 GeV, rewrites `param_card.dat`, and only logs a warning
  (`MadSpin/decay.py:3196-3225`, `1847-1861`). Intent (on-shell W) would be silently lost.
- BUT MadSpin is **not invoked anywhere** in this pipeline
  (`parallel_pipeline_DYj.sh:14` "No MadSpin — Z kept stable throughout"; it only fires if
  `Cards/madspin_card.dat` exists — only the inert `madspin_card_default.dat` is present;
  gate at `amcatnlo_run_interface.py:1354`, `madevent_interface.py:861`).
- ⇒ Keep W's stable with WW=0; do **not** add a `madspin_card.dat`. If leptonic decays are
  ever needed, let Pythia8 decay them (`wp_stable=F`/`wm_stable=F`), accepting approximate
  spin correlations. Zero width is correct & harmless for stable-boson Sudakov observables.

### 3. Biasing toward large M(WW) — mechanism found, already partly active
- NLO has **no `bias_module`** (that is LO/madevent-only: `madgraph/various/banner.py:4357`;
  `RunCardNLO` has none, `banner.py:5697+`). The NLO bias is a Fortran routine
  **`bias_weight_function(p,ipdg,bias_wgt)`** injected via `set custom_fcts`
  (mapped at `banner.py:5685`), consumed in the MINT integrand
  (`Template/NLO/SubProcesses/fks_singular.f:2243-2310`), and divided back out of the stored
  event weight (`include_inverse_bias_wgt`, `fks_singular.f:2317-2339`; `driver_mintMC.f:305-308`).
  Cross section is preserved (importance sampling).
- **`set event_norm bias` (`generate.sh:78`) is the on/off GATE + a normalization convention**
  — not itself a phase-space bias. Without a non-trivial `bias_weight_function` it does nothing
  to sampling; with any non-`bias` value the bias routine is neutralized
  (`fks_singular.f:2259-2264`). σ off the banner is tagged "do not use"; recover as Σw/N.
- **A bias is ALREADY ACTIVE in this dir's `dummy_fct.f:112-116`: `bias_wgt = H_T**2`**
  (the Template ships these lines commented out). So "the quadratic" = the existing H_T² bias.
  H_T = Σ final-state transverse masses; for WW it correlates with but is **not** M(WW). The
  clean Sudakov variable is ŝ = M(W⁺W⁻)². "Something smart" = a steeper power/exponential of
  M(WW) to flatten dσ/dM so MC stats are ~uniform per bin.
- Built-in alternative (no Fortran): generation cut `mxx_min_pdg {24:X}` +
  `mxx_only_part_antipart {24:True}` ⇒ M(W⁺W⁻)>X (`banner.py:5833-5834`, `cuts.f:713-757`,
  `setcuts.f:130-142`). **Changes σ** (drops low mass); no `mxx_max_pdg` exists.
- Bias × Sudakov reweight **factorize cleanly**: reweight is multiplicative on the stored
  (bias-compensated) weight, `w_final = w_phys × K_Sudakov`
  (`fxfx_ewsudakov.py:3750,4418`; `reweight_interface.py:1679-1696`).

### 4. `ickkw 3` on a single 0-jet WW sample — CONTRADICTS the "keep it 3" instruction
- **The EW Sudakov reweight does NOT require ickkw 3.** It routes on ickkw
  (`reweight_interface.py:1300-1344`: 0→scalar, 3→density, else→Exception) and computes the
  logs from hard-process momenta + α_s only — kernel signature `ewsudakov(pdgs, p, g)`, no
  scale/merging input. For undecayed WW the density path falls back to the scalar path anyway
  (`fxfx_ewsudakov.py:3818-3825,3870,3894-3896`). ⇒ `ickkw 0` is sufficient and correct.
- **ickkw 3 is harmful for a single-multiplicity sample:**
  (a) at generation it forces `dynamical_scale_choice=-1`, `jetalgo=jetradius=1.0`,
      `fixed_*_scale=False` (`banner.py:5880-5901`) — overrides the intended QCD scale;
  (b) downstream it sets Pythia `JetMatching:doFxFx=on` with the pipeline's hardcoded Z values
      `nJetMax=2`, `qCut=30` (`MCatNLO_MadFKS_PYTHIA8.Script:613-633`;
      `amcatnlo_run_interface.py:4750`; `parallel_pipeline_DYj.sh:56,57,71`) — an FxFx veto
      applied to a 0-jet sample reshapes the inclusive WW rate.
- The only reason to keep 3: the `FxFx_qCutAccept_<q>` named-weight bookkeeping is gated on
  `isFxFx` (`Pythia83_hep.cc:70-71`), so ickkw 0 ⇒ no qcut flags. But that machinery is
  meaningless for a 0-jet diboson sample. ⇒ recommend **ickkw 0** + drop the qcut apparatus.
  **Needs Rikkert confirmation if the stated reason for "keep 3" was the qcut weights.**

## Decisions (2026-06-18)
- **ickkw: kept at 3** per Rikkert's instruction, despite the audit showing it
  distorts a single-multiplicity WW sample. `generate.sh` carries a loud comment;
  the downstream shower MUST set Pythia `JetMatching:nJetMax=0` (pipeline still
  hardcodes 2/qCut=30). Reconcile with Rikkert (likely the qcut-weight bookkeeping).
- **Bias: M(W⁺W⁻)² quadratic.** Replaced `bias_wgt=H_T**2` in `dummy_fct.f` with the
  W-pair invariant mass squared (identify W's by |PDG|=24, sum 4-momenta, M²).
- **4FS import: `loop_qcd_qed_sm_Gmu_forSudakov-with_b_mass_no_width`** (MB=4.5,
  all widths 0); `set maxjetflavor 4`. No new restriction card needed.
- **No MadSpin; W's stable, width 0.** (Do not add a `madspin_card.dat`.)
- (still pending — user's call) PDF/4FS scheme; MT/MH shift; the FO "full EW" leg.

## Edits made
- `dummy_fct.f:99-135` — `bias_weight_function` now returns M(W⁺W⁻)² (IR-safe; guard
  returns 1 if the two W's are not both present).
- `generate.sh` — header; process block (`generate p p > w+ w- [QCD] @0`, b-mass/no-width
  model); `maxjetflavor 5→4`; `ickkw 3` kept with warning; PDF-scheme warning comment.
- `models/loop_qcd_qed_sm_Gmu_forSudakov/restrict_with_b_mass_no_width.dat` — **BUG FIX.**
  The shipped card was missing the `TADPOLE` block that this NLO-EW model requires, so
  `import model …-with_b_mass_no_width` aborted with "Invalid restriction card … Missing
  block: tadpole". Appended the exact block the default card carries:
  `Block TADPOLE / 1 1.0001e+00` (MH/MT-independent `ntadpole` tracker; nonzero ⇒ keep the
  tadpole sector). Backup: `…_no_width.dat.pretadpole_backup`. NOTE: the sibling
  `restrict_with_b_mass.dat` (finite-width) has the identical bug — fix it the same way if used.

## Verification (2026-06-18)
- `bash -n generate.sh` clean.
- `dummy_fct.f` compiles (`gfortran -Wall`); only pre-existing unused-arg warnings in
  the template `dummy_cuts`/`user_dynamical_scale` stubs, none in the bias routine.
- Numerical unit test (`/tmp/ww_bias_test/`): bias = 479100 = M² for a hand-built
  W⁺W⁻+gluon event (gluon correctly ignored); guard returns 1.0 when nW≠2. Both PASS.
- Adversarial counter-agent review: verdict FIX-FIRST, **no defect in the edits**.
  Confirmed the (E,px,py,pz) layout by tracing `fks_singular.f:2271→momenta_m(0:3,…)`,
  index 0 = energy (`fks_singular.f:2027,3910`); |PDG|=24 W-id robust (real radiation is
  gluon/quark only under [QCD], never a 3rd W); bias can't NaN (`include_inverse_bias_wgt`
  `stop 1` on 0, `fks_singular.f:2345-2348`); `w+`/`w-` exist (`particles.py:37-39`);
  Born is O(α²) EW with a real NLO QCD correction (vertex `d~ u W-`, `vertices.py:655`);
  `readlink -f dummy_fct.f` resolves to THIS dir (fixes the stale `_on_valentin` path).
  FIX-FIRST driven only by the (pre-existing) PDF mismatch + "not yet run at runtime".
- **Runtime validation (after the TADPOLE fix):**
  - `import model …-with_b_mass_no_width` → exit 0, no restriction error.
  - 4FS confirmed: final `p = j = g u c d s u~ c~ d~ s~ a` (photon added, **no b**).
  - `generate p p > w+ w- [QCD] @0` → exit 0. Born auto-set to `QED²≤4 QCD=0` (pure EW);
    9 Born subprocesses = 8 light qq̄ + `a a > w+ w-` (γγ→WW), **no bb̄**; 29 Born +
    168 real + 48 virtual diagrams. (Benign `DEBUG: invalid command for loop` = the photon
    channel has no QCD loop.)
  - NOT yet run: the `output -f` + `launch --parton` (compile + biased integration). The
    M(WW)² bias compiles into the process only at `launch`; standalone unit test already PASS.

## Open Questions / next steps
1. **Run the full `generate.sh`** — import + `generate` are validated (see Verification);
   still need `output` + `launch --parton` to confirm the process compiles and the biased
   integration completes (the M(WW)² bias is only exercised at `launch`).
6. **γγ→W⁺W⁻ is INCLUDED** (the model auto-adds the photon to `p`, so `a a > w+ w-` is the
   9th Born subprocess). Inert in the old Z setup (no γγ→Z) but ACTIVE for WW (~1-2% with
   luxqed). Decide: keep it, or override with `define p = g u c d s u~ c~ d~ s~` (no photon)
   after import to drop it. Must match whatever the comparison/closure leg assumes.

## Runtime fixes — gridpack build attempt (2026-06-18 PM)
Goal clarified by user: generate.sh is meant to build a **gridpack**. Key fact: NLO has
**no `gridpack=True`** (LO/MadEvent only; `banner.py:4300` is RunCardLO). The NLO "gridpack"
is the **integrated process dir**, reused via `generate_events --only_generation -p`
(`amcatnlo_run_interface.py:5974`). Wrapper added: `regenerate_events.sh`.

Three blockers hit when actually running generate.sh — all fixed:
1. **Import aborted** — `restrict_with_b_mass_no_width.dat` missing the TADPOLE block
   (see Edits made). Appended it; import now exits 0, 4FS confirmed (p has no b).
2. **`launch … --parton` invalid** — `-p/--parton` is NOT on the top-level mg5_aMC `launch`
   parser; it lives on the **aMCatNLO** `generate_events`/`launch` parser
   (`amcatnlo_run_interface.py:5971,6007`). Fix: part2 drives the process-local
   `PROC/bin/aMCatNLO` with `generate_events -p` (mirrors the shower/reweight pipeline).
   `set custom_fcts` is honored there — it is a RunCardNLO param (`banner.py:5823`) applied
   via `edit_dummy_fct_from_file` (`:3687`), so the M(WW)² bias compiles in.
3. **Wrong process-dir name** — `output -f` auto-derives the dir name from the MODEL string,
   so the `-with_b_mass_no_width` suffix produced `PROCNLO_…-with_b_mass_no_width_0` instead
   of the hard-coded `$PROC` (which the pipeline also hard-codes). Fix: `output ${PROC} -f`
   (explicit name). VERIFIED: part1 rebuilds `PROCNLO_…_forSudakov_0` + `bin/aMCatNLO`
   (exit 0); misnamed orphan removed.

**State:** part1 (output) verified working; PROC dir built and ready. **part2 integration
NOT yet run** — command flow is code-verified but no full integration executed yet, and it
should wait on the open PDF / MT-MH decisions (a real run bakes those in).

## Runtime fixes — deeper cascade (2026-06-18, later)
Root cause of most of this: the `_WW` dir was COPIED from `_on_valentin` without re-pointing
its toolchain. Additional blockers found by actually running, all addressed:
4. **Interactive `set` lines don't reach the run_card editor** when piped to aMCatNLO (they
   hit the config `set` at the main prompt → InvalidCmd). FIX: generate.sh part2 now
   PRE-WRITES the run_card via a `set_rc` sed helper, then runs a single bare
   `generate_events -p` (the pipeline's proven pattern). NOTE: `maxjetflavor` is absent from
   this single-multiplicity run_card and auto-derives to 4 from the massive-b model (do NOT set it).
5. **Stale `_on_valentin/HEPTools` paths** in `input/mg5_configuration.txt` (lhapdf_py3, ninja,
   pythia8_path, hepmc_path, mg5amc_py8_interface_path) + baked prefix in `lhapdf-config`/`lhapdf.pc`.
   FIX: swept 92 HEPTools text files `_on_valentin → _on_valentin_WW`.
6. **Bundled `libLHAPDF.so` will not link** under this host's binutils (`ld: failed to set
   dynamic section sizes: bad value`); only bites with `pdlabel=lhapdf` (luxqed 324900). FIX:
   point MG5 at the **system** `/usr/bin/lhapdf-config` (LHAPDF 6.5.6; test-linked OK) in
   `input/mg5_configuration.txt`. 324900 resolves via the system `pdfsets.index`; the SET
   data comes from this checkout via `LHAPDF_DATA_PATH` (the system data dir lacks luxqed).
7. **`set custom_fcts` does NOT inject the bias**: `edit_dummy_fct_from_file` (`banner.py:3401`,
   regex `:3432`) rewrote `SubProcesses/dummy_fct.f` but dropped `bias_weight_function`. FIX:
   generate.sh part2 **direct-copies `dummy_fct.f`** into `SubProcesses/` and every `P0_*/`
   (what actually compiles), with `custom_fcts` left empty so the editor is a no-op on a fresh build.

**Verified:** import + part1 (correct dir name) OK; a built-in-`nn23nlo` run produced events
(engine works); system libLHAPDF links. **NOT yet validated end-to-end:** a full luxqed run that
(a) compiles+integrates+writes events AND (b) confirms the direct-copied M(WW)^2 bias survives
`edit_dummy_fct` and is active. generate.sh fully rewritten + `bash -n` clean.
2. **PDF/4FS:** accept 4FS-ME × 5F-PDF hybrid, install a 4F luxqed set, or reconsider 4FS
   (sub-% for inclusive WW)? Confirm whether the Sudakov reweight truly needs massive b.
3. **MT/MH shift** (174.3/120, bundled in the b-mass card) — intended, or override to
   173.3/125 for cross-leg consistency? (Negligible for the [QCD] WW rate; matters for the
   reweight / full-EW leg.)
4. **"Full EW corrections" (Rikkert):** most likely the FO `[QCD QED]` closure leg for WW
   (mirror of `generate_FO.sh` Zj/Zjj) — validate Sudakov-approx vs full NLO EW. Confirm scope.
5. **ickkw 3 (kept):** if showering, set Pythia `JetMatching:nJetMax=0` in the pipeline
   (currently hardcoded 2/qCut=30 for Z+jets). Confirm Rikkert's reason for keeping 3.

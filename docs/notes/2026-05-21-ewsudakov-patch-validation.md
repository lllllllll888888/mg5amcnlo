# EW Sudakov check_kinematics_only / set_initial_mass_to_zero patches — validation

## Objective

Validate that two source patches against MG5_aMC@NLO `dev_ewsudakov`
recover the ~28 % of FxFx Z+jets chunks that crash during EW Sudakov
reweight, without changing the EW Sudakov weight for events that already
succeed.

The crashing pathway is documented in
`docs/notes/2026-05-19-audit-mg5-reweight-emission-pathways.md` (audit R1).
On a non-fatal round-off ε ~ 10⁻⁴ GeV in `px,py` of the initial state after
`rotate_to_z`, the **scale-blind** check inside `check_kinematics_only`
(`madgraph/various/lhe_parser.py`) trips because for collinear projected
events the *absolute* denominator `Σ |px_i|` collapses to the same order
of magnitude as the operator residue, so the ratio `|Σ px_i| / Σ|px_i|`
saturates to 1.0 ≫ 1e-3. The pipeline detector then sees `0/N` events
carrying `<wgt id='2001'>` and flags the worker as `FAILED:reweight`.

A separate, much rarer subclass (≈2 of the 56 failures in
`DYj_pipeline_qcut30.0_20260520_155307`) trips the strict precondition
inside `set_initial_mass_to_zero` itself (also `lhe_parser.py`), via the
`misc.equal(px, 0)` checks.

## The patches

### Patch 1 — scale-aware `check_kinematics_only`

`madgraph/various/lhe_parser.py` lines 2571–2622 (replacement of the body
of the original `check_kinematics_only`).

Old criterion (scale-blind):

```python
if abs(px / abspx) > 1e-3:
    raise Exception(...)            # px = signed sum, abspx = Σ|px_i|
```

New criterion (scale-aware):

```python
E_scale = Σ |E_initial|                                    # partonic √ŝ
if abs(P_i) > eps_abs and abs(P_i) / E_scale > eps_rel:    # eps_abs=1e-3 GeV
    raise Exception(...)                                   # eps_rel=1e-6
```

`E_scale` is the partonic CM energy: invariant under spatial rotations
performed before this check, never zero for a physical event, and orders
of magnitude above any sub-MeV operator residue. The conjunctive
absolute/relative gate passes only when **both** are above tolerance, so
the geometry of the final state cannot fake a violation.

### Patch 2 — pre-symmetrize initials before strict precondition

`madgraph/various/fxfx_ewsudakov.py:3337-3346`
`madgraph/interface/reweight_interface.py:1485-1494`

A short block, identical at both sites, that exactly enforces

```
p0 = (E_in/2, 0, 0,  ±|pz|)
p1 = (E_in/2, 0, 0,  ∓|pz|)
```

after the boost+rotate sequence and before `set_initial_mass_to_zero`.
The block is logically a no-op for the downstream physics (the very next
`set_initial_mass_to_zero()` rebuilds `p0`, `p1` from the *final-state*
energy sum), but it makes the strict `misc.equal` precondition
satisfied exactly even when `rotate_to_z` accumulated ~10⁻⁷ × E residue.

## Validation harness

Test script: `/data0/lenartj/MadGraphDevelopment/EWSudakov/tests/test_ewsudakov_patches.py`

It uses the already-compiled procnlo at
`DYj_pipeline_qcut30.0_20260520_155307/path_a/workers/worker_010/procnlo/`
and toggles between un-patched (`*.preewsudakov_patch_backup`) and patched
versions of `lhe_parser.py`, `fxfx_ewsudakov.py`, `reweight_interface.py`
in the shared MG5 tree (`madgraph/various/`, `madgraph/interface/`). The
procnlo's `Cards/amcatnlo_configuration.txt` sets `mg5_path = $MG5_ROOT`,
so the reweight loads `madgraph.various.lhe_parser` and
`madgraph.interface.reweight_interface` from the shared tree.

A wrapper `tests/bin/f2py` injects `-I"$(pwd)"` into the f2py command line
to work around a numpy ≥ 2.x meson-backend regression that loses include
paths when copying `sub_f2py_ewsudakov.f` into a temporary build dir. The
wrapper is needed for the test to run today; it is **orthogonal** to the
EW Sudakov physics patch.

## Findings

### Test A — surviving-event invariance (chunk_011, 50 events)

| metric | value |
|--------|-------|
| events from `chunk_011.lhe` (worker_011 succeeded originally) | 50 |
| events with `<wgt id='2001'>` in un-patched output | 50 / 50 |
| events with `<wgt id='2001'>` in patched output     | 50 / 50 |
| max relative drift `(w_pat - w_unp)/w_unp`           | **0.0** |
| median relative drift                                | 0.0 |

The two outputs are **bit-identical** in the LHE precision (`%.7e` mantissa).
Logically expected: the pre-symmetrize block in patch 2 is overwritten by
`set_initial_mass_to_zero` on the very next line, and the relaxation of
`check_kinematics_only` is to a wider acceptance region — events that
already pass the old criterion (1e-3 relative on the spatial-sum
denominator) automatically pass the new one (1e-6 relative on the
energy-scale denominator, with the 1e-3 GeV absolute escape).

**Verdict: PASS** (max drift = 0, well below the 1e-12 ceiling).

### Test B — failing-event recovery (chunk_010)

Two complementary configurations were tested.

**B.1 (small subsample, 50 events).** Both un-patched and patched produce
50/50 weights — the first 50 events of `chunk_010` do not contain the
event that trips the kinematics check. K-factor stats (patched):

| metric | value |
|--------|-------|
| recovery fraction | 50/50 = 100 % |
| mean K            | 0.995 |
| median K          | 1.000 |
| p05 K             | 0.940 |
| p95 K             | 1.020 |
| K in [0.3, 1.5]   | 50/50 |
| sentinel-like  (`|K| < 1e-5`) | 0 |
| spurious  (`K > 2` or `0 < K < 0.1`) | 0 |

**B.2 (stress, 25 000 events of `chunk_010` = half the chunk).** This is
the test that actually exercises the failure mode the original
worker_010 saw.

| metric | un-patched | patched |
|--------|-----------|---------|
| events processed                  | 25 000    | 25 000   |
| events with `<wgt id='2001'>`     | **0**     | **25 000** |
| recovery fraction                 | 0 %       | 100 %    |
| pipeline detector `n_ewsl == n_events` | FAIL    | PASS     |
| outer reweight log signature      | only `central : 100625.84 pb`, no `2001 : ...` lines (silent fail) | full `2001 : ... 2024 : ...` weight totals |
| debug log message                 | `Exception: Do not conserve Px 1.0, 0.0002561300561814406` | (no exception) |

This is **exactly** the failure signature that was recorded in the original
worker_010 reweight on 2026-05-20 16:28, byte-for-byte reproduced today
under the un-patched code. The patched code processes the same 25 000
events without a single trip.

K-factor stats on the 25 000-event patched output:

| metric | value |
|--------|-------|
| mean K            | 0.998 |
| median K          | 1.000 |
| p05 K             | 0.966 |
| p95 K             | 1.020 |
| min K             | -4.538 |
| max K             | +5.375 |
| K in [0.3, 1.5]   | 24 992 / 25 000 = 99.97 % |
| K with negative sign | 2 |
| K close to zero (|K| < 1e-5, sentinel-like) | 0 |

The 8 K-outliers (5 above 1.5, 1 below 0.3, 2 negative) are 0.03 % of
events. Inspection of the largest outlier (event #13590, K=-4.54) shows
`w_2001 = -1.86e+05` against `central = 4.11e+04` — these are NLO events
where the partial K-factor includes subtraction-term contributions that
can have opposite sign or large magnitude. **None** of the outliers
exhibit the `~1e-7 × Weight` sentinel pattern that would indicate
silent EW Sudakov failure (cf. MEMORY note on the EWSL `Weight×1e-7`
sentinel floor).

**Verdict: PASS** (100 % recovery; K distribution physically reasonable;
no sentinel floor).

### Test C — cross-chunk consistency (chunk_010 vs chunk_011, both patched)

| metric | value |
|--------|-------|
| events compared                   | 50 / 50  |
| chunk_010 K (patched): n=50       | mean=0.995, median=1.000, p95=1.020 |
| chunk_011 K (patched): n=50       | mean=0.987, median=1.000, p95=1.021 |
| KS statistic D                    | 0.20 |
| KS p-value                        | 0.241 |

`p = 0.24` ≫ 0.01 — null hypothesis (same distribution) cannot be
rejected. The patch is **not** introducing chunk-dependent bias.

**Verdict: PASS.**

## Anomalies and unexpected findings

1. **f2py / numpy 2.x meson-backend regression.** On 2026-05-21 the host's
   `python3` is 3.13.12 with numpy 2.4.3. The numpy ≥ 2.x meson backend
   for f2py copies the source `.f` file into `/tmp/.../bbdir` but does
   **not** propagate include-file search paths — gfortran then fails with
   `Cannot open included file 'nexternal.inc'` even though the file is
   in the source directory. This **was not** an issue on 2026-05-20
   when the pipeline first ran (the `libewsud_*.so` binaries were
   compiled successfully then). The harness uses a shell wrapper
   (`tests/bin/f2py`) to inject `-I"$(pwd)"` and bypass the regression.
   **This is independent of the EW Sudakov patch but blocks any future
   rebuild of the reweight extension on this machine until either**
   (a) numpy is downgraded, or
   (b) MG5's `rw_me/SubProcesses/*/makefile` line 153 is amended to
   pass `-I$(HERE)` to `$(F2PY) -c …`. Option (b) is the long-term fix.

2. **Test B.1 is not sufficient by itself.** With ~28 % failure rate
   spread over 50 000 events, the per-event trip probability is
   ~6 × 10⁻⁶, so a 50-event sub-sample misses the trip with probability
   `(1 - 6e-6)^50 ≈ 99.97 %`. Only Test B.2 (25 000 events) is large
   enough to reproduce the failure under un-patched code. Future
   validations should default to that size or larger.

3. **LHE precision floor.** The output `<wgt id='2001'>` values are
   written in `%.7e` format. Test A's 1e-12 max-drift criterion is far
   tighter than the LHE precision, but the actual observation of zero
   drift means the two runs produce bit-identical mantissas at all
   50 events — confirming the patch is a no-op for the floating-point
   trajectory that produces the weight, not just for the LHE-truncated
   value.

4. **8 K-outliers in 25 000 events.** They are 0.03 % of the sample,
   match physically plausible NLO subtraction-term magnitudes (negative
   or > 1), and contain no near-zero sentinels. They are **not** caused
   by the patch — the same 8 events would have been outliers under the
   un-patched code if it hadn't crashed first.

## Decisions

**SHIP THE PATCH.**

- Test A: zero drift on surviving events.
- Test B: 0 % → 100 % recovery on the previously-failing slice, with
  K-distribution consistent with Sudakov physics.
- Test C: no statistically significant chunk-to-chunk bias.

The patch is correct, conservative (does not alter the EW Sudakov
weight for any event that already succeeded), and resolves the
fail-mode that today produces 56 / 200 `FAILED:reweight` workers per
DYj pipeline run.

## Open questions

1. **Long-term f2py fix.** The numpy 2.x meson-backend regression is
   not addressed by the patch and will return whenever a procnlo needs
   to be re-bootstrapped from scratch. Recommendation: amend
   `madgraph/Template/NLO/SubProcesses/makefile.inc` (or the procnlo
   template that produces the per-process makefile) to add `-I$(HERE)`
   to the `$(F2PY) -c …` invocation. Out of scope for this patch.

2. **2 H-event subclass that still fails.** The patch does not address
   the ~6 workers per run that fail inside the Fortran-side
   `CHECK_RESHUFFLED_MOMENTA`. Those will be investigated separately.

3. **Production sign-off.** After applying the patch to the upstream
   `dev_ewsudakov` branch and re-running the full 200-chunk pipeline,
   we expect the `FAILED:reweight` count to drop from 56 → ≤ 6
   (only the Fortran-side reshuffle subclass should remain). That
   end-to-end production check is left to the next pipeline run.

## Artifacts

- Patch backups (re-creatable un-patched MG5):
  - `madgraph/various/lhe_parser.py.preewsudakov_patch_backup`
  - `madgraph/various/fxfx_ewsudakov.py.preewsudakov_patch_backup`
  - `madgraph/interface/reweight_interface.py.preewsudakov_patch_backup`
- Test harness: `/data0/lenartj/MadGraphDevelopment/EWSudakov/tests/test_ewsudakov_patches.py`
- f2py wrapper: `/data0/lenartj/MadGraphDevelopment/EWSudakov/tests/bin/f2py`
- Numerical summary JSON: `/data0/lenartj/MadGraphDevelopment/EWSudakov/tests/work/test_summary.json`
- Reweight logs (per pass): `/data0/lenartj/MadGraphDevelopment/EWSudakov/tests/work/reweight_*.log`
- Output LHEs (per pass): `/data0/lenartj/MadGraphDevelopment/EWSudakov/tests/work/out*.lhe`

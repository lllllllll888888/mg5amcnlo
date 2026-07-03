# Sudakov reweighting factor analysis + KeyError silent-skip fix

## Objective

Two things, on top of the running 10M Z+jets EW Sudakov pipeline:

1. Compute the **literal** Sudakov reweighting factor distribution per event
   (mean, std, min, max), independent of bias / scale / PDF mixing.
2. Fix a silent-skip bug where some FxFx-clustered events exited the Sudakov
   reweighting code without writing any 2XYZ weights into the output.

Source data: 70 path_a worker outputs from the prior 10M run
`DYj_pipeline_qcut30_20260509_120137` (started 12:01:37, killed 14:42 after
silent-skip fix). 5 worker files = 250k events were used for stats.

## Findings

### 1. Weight schema in the HEPMC2 output

275 weights per event, partitioned as:

- `"0"`                       — 1 weight, bias-normalized
- `"1010"`–`"1045"`           — 36 weights, scale × PDF variations (no Sudakov)
- `"2X01"`–`"2X45"`, X=0..4   — 5 × 45 = 225 Sudakov-applied weights
- `"FxFx_qCutAccept_{20,30,40}"` — 3 weights, qCut acceptance flags
- `"MUR{0.5,1,2}_MUF{0.5,1,2}"`  — 9 weights, named QCD scale variations
- `"Weight"`                   — 1 weight, nominal pre-Sudakov reference

The 5 X-prefixes correspond to the 5 NLL Sudakov variants:
`X=0 central, X=1 s_to_rij_off, X=2 LL, X=3 both_off, X=4 rij_ge_mw_off`.
The 45 sub-indices are the (PDF tag × μR × μF) combinations that
`event.parse_reweight()` returns; only 1010-1045 of those (36) are written
into the HEPMC's 1XXX block — the 1001-1009 endings are dropped because
they duplicate the named MUR×MUF set. So 2XYZ for ending 01-09 has no
matching 1XYZ in the file and only 36 endings out of 45 are useful for
ratio analysis.

### 2. Literal Sudakov reweight factor — `w[2XYZ] / w[1XYZ]`

Aggregated over **199,999 events** × **180 (variant, ending) pairs** (only
counting cases where `w[1XYZ] != 0` so the ratio is defined):

| Variant | N (event×ending) | Mean | Std | Min | Max |
|---------|----:|----:|----:|----:|----:|
| central        | 2,914,308 | 0.99879 | 0.33008 | −6.329 | +186.88 |
| s_to_rij_off   | 2,914,308 | 1.00686 | 0.65514 | −11.53 | +371.40 |
| LL             | 2,914,308 | 0.99914 | 0.33011 | −6.328 | +186.89 |
| both_off       | 2,914,308 | 1.00702 | 0.65522 | −11.53 | +371.40 |
| rij_ge_mw_off  | 2,914,308 | 0.99877 | 0.33007 | −6.329 | +186.88 |

Observations:

- **Means ≈ 1.0** for all five variants, consistent with the σ-level result
  (Sudakov correction at σ-level is ~1.0024 for central per the prior 100k
  run, ~0.5% spread across variants).
- **Stds are O(0.3–0.7) per event**, dominated by NLO H-event configurations
  where `1 + δ_h/Born` blows up when Born is small. Most events sit in a
  narrow core near 1.0; the std is heavy-tailed, not gaussian.
- **The damping cap (`abs() > 200` → set to 1`) only protects sudrat0 (central)**
  in `_compute_ewsudakov_reweight`. The other variants (s_to_rij_off,
  both_off) leak through and reach +371 in the data. This is a separate
  bug — the cap should be applied to all 5 variants symmetrically.
- The `central / LL / rij_ge_mw_off` triplet has identical std (0.33008,
  0.33011, 0.33007) → these three variants are physically near-degenerate
  for this process; the discriminating switch is `s_to_rij`.

### 3. Silent-skip pathology in the FxFx Sudakov code

Two `except KeyError` paths in `madgraph/various/fxfx_ewsudakov.py` returned
only `{"orig": event.wgt}` instead of the full 2XYZ weight set:

- Line 3564 in `_compute_ewsudakov_fxfx_reweight` (scalar FxFx path)
- Line 3890 in `_compute_density_ewsudakov_reweight` (density-matrix FxFx)

Trigger: `_prepare_fxfx_sudakov_inputs` raises KeyError when the clustered
event's PDG ordering is not registered in `sud_mod.pdg2ewsud_dict`. This
happens for clustering paths the EW Sudakov module has no amplitude for.

Effect on output: those events' 225 `2XYZ` HEPMC slots got default-filled
(often 0). When taking ratios `2XYZ/1XYZ` downstream, you get either
`0/positive = 0` or `negative/positive = negative`, contaminating the
event-level distributions with bogus zero/negative ratios.

### 4. The 2→1 topology branch was a red herring

The user's mental model was that 2→1 events skip writing weights. They don't:
lines 3556 and 3680 already call `_build_sudakov_rwgt_dict(event, [1.0]*5)`
which writes `2XYZ = 1XYZ × 1.0` (i.e. nominal pass-through). Same for
small-invariant cases at lines 3572 and 3899. Only the KeyError paths were
the actual silent-skip.

## Decisions

### Fix (applied 2026-05-09 14:42, before relaunch)

Both KeyError paths in `madgraph/various/fxfx_ewsudakov.py` now call
`self._build_sudakov_rwgt_dict(event, [1.0, 1.0, 1.0, 1.0, 1.0])` and write
nominal weights, matching the existing 2→1 / small-invariant branches:

```python
# scalar FxFx path (line 3563)
except KeyError as exc:
    _dbg(f"ERROR: {exc}, returning weight=1 pass-through (nominal)")
    return self._build_sudakov_rwgt_dict(event, [1.0, 1.0, 1.0, 1.0, 1.0])

# density-matrix FxFx path (line 3889)
except KeyError as exc:
    _dbg(f"  -> KeyError: {exc}, returning weight=1 pass-through (nominal)")
    density_logger.log_event("passthrough_keyerror", [], 0, 1.0)
    return self._build_sudakov_rwgt_dict(event, [1.0, 1.0, 1.0, 1.0, 1.0])
```

Sync points checked: only one `fxfx_ewsudakov.py` exists in the tree
(no PROCNLO snapshot, no HEPTools copy). `.pyc` cache cleared after edit.

### Pipeline restart

Killed the in-flight 10M run (PID 249255, PGID 249255, 344 descendant
processes) cleanly via SIGTERM to the process group. Removed 649 GB
working dir. Disk: 347 GB free for the relaunch.

Relaunched at 14:45:22 with identical parameters:
```
QCUT=30 EVENTS_PER_CHUNK=50000 MAX_CONCURRENT=32 \
  ./parallel_pipeline_DYj.sh \
  /data0/lenartj/MadGraphDevelopment/EWSudakov/lhe_dumps/DYj/pc16_seed10077.lhe.gz
```

New work_dir: `DYj_pipeline_qcut30_20260509_144522`. PID 103343.

## Open Questions

1. **Damping cap leak**: The `abs(sudrat) > 200 → 1.0` damping in
   `_compute_ewsudakov_reweight` (reweight_interface.py:1505) only checks
   `sudrat1` (central). `sudrat2/sudrat3/sudrat4` can reach ±371 in the
   data. Should the cap be applied per-variant?
2. **`"0"`/`"Weight"` ratio == 0.0 for every event** in the prior run.
   The bias-normalized "0" weight should be ~10^-7 of "Weight" (cross-section
   ratio: 985 / 9.8e9). Either the bias normalization is broken or the "0"
   slot isn't getting populated correctly. Investigate after the 10M finishes.
3. **9 missing 1XXX endings** (1001-1009 with tag=0): these are dropped by
   the LHE → HEPMC writer because they duplicate the named MUR×MUF set, but
   the corresponding 2X01-2X09 entries DO get written. There's no way to
   compute Sudakov ratios for those endings from the HEPMC alone. Acceptable,
   but worth documenting for downstream consumers.
4. **Heavy-tailed Sudakov factors**: per-event std 0.33–0.65 with min/max
   reaching ±300 is dominated by H-event configurations with small Born.
   Whether the σ-level result is robust against these tails (jackknife, or
   tail-trimming) is worth checking before quoting final cross-sections.

## Reproducer

```bash
# Stats script (after the 10M run produces output):
python3 /tmp/sudakov_factor_stats.py
# Edit ROOT and N_FILES at the top.
```

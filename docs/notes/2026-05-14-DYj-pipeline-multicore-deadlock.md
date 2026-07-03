# DYj reweight pipeline — MG `MultiCore` compile-phase deadlock

## Objective

Two consecutive runs of `parallel_pipeline_DYj.sh` over `lhe_dumps/DYj/pc12seed1001.lhe.gz`
failed to make any progress past the reweight stage. Yesterday's run (`DYj_pipeline_qcut30.0_20260513_141958`)
took the host down hard; today's (`DYj_pipeline_qcut30.0_20260514_100107`) left the host
alive but every worker stuck. Root-cause both failures and define a minimal, verified fix.

## Findings

### Live state (worker_013 representative, PID 13936)

- `/proc/13936/wchan` = `futex_wait_queue`, single thread, no children.
- `voluntary_ctxt_switches: 7978` — flat for tens of minutes → truly stuck, not slow.
- `grep -c libewsud /proc/13936/maps` = **0**. The EW Sudakov reweight extension is
  not even mapped — the process is upstream of the reweight execution.
- Last log line in `procnlo/Events/run_01/reweight.log`: `INFO: Compiling on 32 cores...`
  (`reweight_interface.py:2526`). The `Idle: 0, Running: 1, Completed: 0` line that I
  initially fixated on belongs to the *outer* MultiCore wrapping the whole reweight
  task — `Running: 1` there means "the reweight controller is running," not "one
  Fortran worker is computing."
- fd 4 and fd 5 point to `/dev/null`. No pipes or sockets open. Stdin (`reweight.mg`)
  fully consumed.

### Call chain to the hang

1. `common_run_interface.py:686-694` — default `run_mode = 2`, `nb_core = None`.
2. `common_run_interface.py:3724-3727` — `configure_run_mode(2)` promotes
   `nb_core = multiprocessing.cpu_count()` = 32 on this 32-core box.
3. `reweight_interface.py:2510-2533` — for each reweight, MG must compile ~402
   `P*_*` subprocess libraries:
   ```python
   nb_core = int(self.options['nb_core'])            # = 32
   compile_cluster = cluster.MultiCore(nb_core=32)
   for p_dir in p_dirs: compile_cluster.submit(make, ...)   # 402 submits
   compile_cluster.wait(self.me_dir, update_status)         # ← hangs here
   ```
   with `update_status = lambda i, r, f: (i, r, f)` — a deliberate no-op. That is
   why the log goes silent at `Compiling on 32 cores...`.
4. `cluster.py:893-925` — `MultiCore.wait()` accounting bug:
   ```python
   Idle = self.queue.qsize()
   Done = self.nb_done + self.done.qsize()
   Running = max(0, self.submitted.qsize() - Idle - Done)
   if Idle + Running <= 0 and not force_one_more_loop: break
   use_lock = self.lock.wait(300)  # threading.Event
   ```
   The daemon worker threads exit when the queue drains (`Threads: 1` post-hoc
   confirms it), but the main thread's exit predicate doesn't trip reliably under
   load — and once all daemons are gone, no one signals `self.lock` anymore.

### Outer concurrency × inner concurrency

`parallel_pipeline_DYj.sh` defaults `MAX_CONCURRENT=$NCPU` = 32 (line 44). With
the default config, each of those 32 outer workers also instantiates a 32-thread
`MultiCore` for compile. **32 × 32 = 1024 concurrent `make`s**, on a 32-core box.
This matches the observed load=207 spike during the compile burst.

Yesterday's host went down at exactly this transition — most likely OOM/IO
exhaustion under the 1024-make storm. Today the host survived because we have
slightly more headroom on this boot session, but the same `MultiCore.wait()`
livelock then made all workers wedge instead of crash.

### Yesterday vs. today

- Same driver script (mtime 2026-05-13 11:23, not touched since).
- Same input, same chunking, identical procnlo template.
- Same hang point — `Compiling on N cores...` in `reweight.log`.
- Different visible failure mode only because of host load at the moment of
  livelock (compile-storm exhaustion → kernel killing things vs. straightforward
  livelock with idle box afterwards).

### Why no make-fail log line

Grep across all 32 `reweight.log`s for "non zero status" / "ends with non zero" /
"Compilation .* failed" returned nothing. So whatever broke the accounting was
either an exit timing race (worker drained queue and exited *just* as main thread
re-entered `lock.wait()`) or a make whose stderr was eaten elsewhere. The fix
applies regardless.

## Decisions

**Fix:** edit
`PROCNLO_loop_qcd_qed_sm_Gmu_forSudakov_0/Cards/amcatnlo_configuration.txt`:

- Line 164: `# nb_core = None` → `nb_core = 1`
- Line 123: `# run_mode = 2` → `run_mode = 0`  (belt-and-suspenders)

Effect: `cluster.MultiCore(nb_core=1)` for the compile path — sequential makes,
no accounting race. `run_mode = 0` also disables the cluster path for the event
reweight phase (`common_run_interface.py:2178` shortcut). Outer parallelism
(32 workers, one MG run each, one make each) ≈ 32 concurrent serial compiles,
matching the core count exactly.

**Concurrency:** keep `MAX_CONCURRENT=16` for the first verification run
(half-bandwidth, leaves headroom).

**Process management:** the 32 zombie controllers + ~600 ancestor bash/pythons
hold ~20 GiB RAM and contribute nothing. Kill the driver tree (`kill 9050`,
sweep `parallel_pipeline_DYj` and `amcatnlo_run_interface.py` stragglers) before
restart. Pending user OK.

**Pipeline dir** `DYj_pipeline_qcut30.0_20260514_100107` is unsalvageable — no
output produced; delete and start fresh after fix.

## Open questions

- Whether `MultiCore.wait()` should be patched upstream (proper fix, but invasive)
  — punted in favor of the config workaround.
- Whether the same livelock can fire during the shower stage (Pythia8). Same
  `cluster.MultiCore` infrastructure used; the `nb_core = 1` fix should cover it
  too. Verify after first successful path_a worker.
- The mysterious 312500-event default bunch size — not the problem here, but
  still worth understanding for future tuning.

## Diagnostic agents fielded

1. **Config-tracer (Explore)** — pinned `run_mode = 2` default + `nb_core` promotion
   and the `--multicore=create` injection at `common_run_interface.py:2267`.
2. **MG MultiCore source-diver (Explore)** — identified the `cluster.py:920`
   futex point, the `submitted/done` accounting bug, and the missing SIGCHLD reaper.
3. **Pipeline-script auditor (Explore)** — read the entire driver, exposed the
   `MAX_CONCURRENT` knob, no-resume semantics, no hang detection, and the precise
   line to patch.
4. **Yesterday-vs-today forensic (Explore)** — confirmed identical pipeline
   templates, identical chunks, identical hang point. (Note: this agent
   overcounted active workers as 64 — actually 32 per `MAX_CONCURRENT` cap.)
5. **Adversarial reviewer (general-purpose)** — demolished my original
   "multiprocessing.Queue.get() with dead worker" framing by reading
   `/proc/13936/maps` (no libewsud), the actual reweight.log tail, and the
   `MultiCore` source. Reattributed the hang to `cluster.py:925`'s
   `lock.wait(300)` after worker-thread exit.
6. **Live-process forensic (general-purpose)** — could not attach gcore/py-spy
   under `ptrace_scope=1`, but via `/proc/$pid/{maps,fd,status}` reached the
   same Python frame conclusion: stuck in `MultiCore.wait`.

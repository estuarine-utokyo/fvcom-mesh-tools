# Time-efficient debugging (project rule)

Owner directive (2026-07-09, stated repeatedly): **time-to-decision
is the FIRST design criterion for every debugging step.** When
choosing how to diagnose a problem, pick the approach that reaches
the next decision in the least wall-clock time — not the most
thorough one. Thoroughness applies to fixes and verification of
results, not to waiting.

## Design rules

1. **Identify the earliest decisive signal.** Before launching
   anything, name the single observable that decides the next
   action (a log line, a timing, a stack frame) and how early it
   appears. Watch THAT; never wait for job completion when the
   verdict is available mid-run.
2. **Time-to-decision includes QUEUE WAIT, and queue wait
   dominates for large jobs on a congested cluster.** During
   debugging, design the SMALLEST job that decides the question
   (fewest cores, smallest dataset that exhibits the effect) and
   clear questions ONE BY ONE — a 4-core probe backfills in
   minutes while a 32-core production job waits hours-days.
   Never use a production-scale run as the debug vehicle unless
   it is ALREADY RUNNING; a running job's early log lines are
   free, a queued one's are not (owner correction, 2026-07-09).
3. **Instrument before re-running.** Every long-running stage gets
   elapsed-time logs; every runner arms
   `faulthandler.dump_traceback_later(900, repeat=True)` so a hang
   names its line automatically. Never run blind twice.
4. **Mini-repro first.** If the suspect code is scale-independent,
   reproduce on the login node in <30 s (small array, cProfile)
   before any batch job.
5. **Judge in ~1 minute.** "Fast enough vs too slow" is decidable
   from the first one or two per-stage timings. Monitors poll at
   30-60 s and EXIT on first evidence, not on completion.
6. **Probe sizing: small cores, generous elapse.** Elapse is free;
   a too-short limit wastes an entire cycle (a 20-min limit killed
   an 8-core probe that needed 25 min). Small cores + honest
   elapse also backfills fastest on a congested cluster.
7. **Verify every patch landed before depending on it.**
   `grep` the edited file for the new symbol after writing; two
   silent anchor-mismatch failures cost full cycles on 2026-07-09.
8. **Kill doomed jobs immediately.** A job that cannot finish
   within its limit (rate × remaining work > remaining time) is
   dead capital — `pjdel`, fix, resubmit.

## Reference outcomes (why these rules exist)

- Rossby unit bug: 2-h blind profile job → replaced by stack-dump
  short job → line-level diagnosis in minutes, slope stage
  4 h → 58 s.
- Signed-clamp dt bug: dedicated 40-min probe was unnecessary —
  the production run prints the deciding `dt_auto` at +9 min.
- DEM cache: profiling showed 97 % of sizing time was one GDAL
  read; a subset cache turned 81 s into 0.5 s for every later run.
